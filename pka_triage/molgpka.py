#!/usr/bin/env python3
"""
molgpka.py -- MolGpKa per-site pKa prediction, plain PyTorch, CPU, batched.

MolGpKa (Pan et al., J. Chem. Inf. Model. 2021, 61, 3159; github.com/Xundrug/MolGpKa,
MIT licence) is a 5-layer GCN with attention pooling that predicts the pKa of ONE
ionizable atom at a time from the 2D graph. This file re-implements its network
and featurizer so the trained weights run without torch_geometric/torch_scatter,
and fixes three things that make the upstream code unusable at library scale:

  * upstream reloads both 6 MB weight files for every molecule
    (predict_acid/predict_base call load_model each time);
  * upstream finds smarts_pattern.tsv relative to the CURRENT directory
    (ionization_group.py: os.path.abspath("")), so it only works from src/;
  * upstream runs one forward pass per site; here all sites of a molecule share
    one dense normalized adjacency and go through in one batch.

Weights and SMARTS are read from a MolGpKa checkout (MOLGPKA_DIR in config), not
copied here. `check_against_upstream()` runs the original code on the same
molecules and reports the largest difference -- run it once per install
(`python molgpka.py --check SMILES_FILE`); it needs torch_geometric.

Semantics kept from upstream, on purpose:
  * the molecule is uncharged (rdMolStandardize.Uncharger) and given explicit Hs
    before prediction, as in predict();
  * every site is predicted on that SAME neutral molecule, independently.
    Neither MolGpKa nor QupKake (as we run it: one prediction per input
    molecule, each site's conjugate made from the input state) models the
    shift a neighbouring site's ionization causes. See TRIAGE.md.
  * acid sites are reported by upstream as the index of the acidic HYDROGEN;
    here they are mapped to its heavy atom, which is what apply_changes and
    QupKake's `idx` use.
"""
from __future__ import annotations

import csv
import hashlib
import math
import os
from dataclasses import dataclass
from pathlib import Path

import torch
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

N_FEATURES = 29

# sha256 of the weight files at Xundrug/MolGpKa 4dc8352 (2024-01-11)
KNOWN_WEIGHTS = {
    "weight_acid.pth": "e0cb05874937f17a4df6943f37beff1c76d14d760787c57d553ac56390147b39",
    "weight_base.pth": "212021f43168208bd67a5c15750f9a176cf4d5906d73a85718f4b15ff8d4b95d",
}


# --------------------------------------------------------------------------
# network: same parameter names as upstream utils/net.py GCNNet, so the
# state dicts load unchanged
# --------------------------------------------------------------------------
class _GCNLayer(torch.nn.Module):
    """Upstream utils/gcn_conv.py GCNConv: X' = D^-1/2 (A+I) D^-1/2 X W + b."""

    def __init__(self, n_in: int, n_out: int):
        super().__init__()
        self.weight = torch.nn.Parameter(torch.empty(n_in, n_out))
        self.bias = torch.nn.Parameter(torch.empty(n_out))

    def forward(self, x, a_hat):            # x: (k, n, in), a_hat: (n, n)
        return a_hat @ (x @ self.weight) + self.bias


class _Gate(torch.nn.Module):
    """torch_geometric GlobalAttention(gate_nn=Linear(1024, 1)), one graph per row."""

    def __init__(self, hidden: int):
        super().__init__()
        self.gate_nn = torch.nn.Linear(hidden, 1)

    def forward(self, x):                   # (k, n, h) -> (k, h)
        g = torch.softmax(self.gate_nn(x), dim=1)
        return (g * x).sum(dim=1)


class GCNNet(torch.nn.Module):
    def __init__(self):
        super().__init__()
        dims = [N_FEATURES, 1024, 512, 256, 512, 1024]
        for i in range(5):
            setattr(self, f"conv{i + 1}", _GCNLayer(dims[i], dims[i + 1]))
            setattr(self, f"bn{i + 1}", torch.nn.BatchNorm1d(dims[i + 1]))
        self.att = _Gate(1024)
        self.fc2 = torch.nn.Linear(1024, 128)
        self.fc3 = torch.nn.Linear(128, 16)
        self.fc4 = torch.nn.Linear(16, 1)

    def forward(self, x, a_hat):            # x: (k sites, n atoms, 29)
        for i in range(1, 6):
            x = torch.relu(getattr(self, f"conv{i}")(x, a_hat))
            k, n, h = x.shape
            x = getattr(self, f"bn{i}")(x.reshape(k * n, h)).reshape(k, n, h)
        x = self.att(x)
        x = torch.relu(self.fc2(x))
        x = torch.relu(self.fc3(x))
        return self.fc4(x).squeeze(-1)


def normalized_adjacency(mol: Chem.Mol) -> torch.Tensor:
    n = mol.GetNumAtoms()
    a = torch.eye(n)
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        a[i, j] = a[j, i] = 1.0
    d = a.sum(1).pow(-0.5)
    return d[:, None] * a * d[None, :]


# --------------------------------------------------------------------------
# featurizer: upstream utils/descriptor.py get_atom_features, reproduced
# column for column -- quirks included -- because the weights were trained on it
# --------------------------------------------------------------------------
_SYMBOLS = ["C", "H", "O", "N", "S", "Cl", "F", "Br", "P", "I"]
_HYB = [Chem.rdchem.HybridizationType.SP, Chem.rdchem.HybridizationType.SP2,
        Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP3D,
        Chem.rdchem.HybridizationType.SP3D2]


def _one_hot(x, allowed):
    if x not in allowed:
        x = allowed[-1]
    return [float(x == s) for s in allowed]


def base_features(mol: Chem.Mol) -> torch.Tensor:
    """Columns 0..26, which do not depend on the site (computed once per molecule)."""
    Chem.AssignStereochemistry(mol)
    ring = mol.GetRingInfo()
    rows = []
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        o = _one_hot(atom.GetSymbol(), _SYMBOLS)
        o += [float(atom.GetDegree())]
        o += _one_hot(atom.GetHybridization(), _HYB)
        o += [float(atom.GetImplicitValence()), float(atom.GetIsAromatic())]
        o += [float(ring.IsAtomInRingOfSize(i, r)) for r in (3, 4, 5, 6, 7, 8)]
        # H-bond donor / acceptor flags: ALWAYS 0 upstream, which tests
        # `atom_idx in GetSubstructMatches(...)` -- an int against a list of
        # 1-tuples -- so the weights were trained with these columns zero.
        # Computing them "correctly" moves predictions by up to 3 pKa units.
        o += [0.0, 0.0, float(atom.GetFormalCharge())]
        rows.append(o)
    return torch.tensor(rows, dtype=torch.float32)


def site_features(base: torch.Tensor, dist: torch.Tensor, sites: list[int]) -> torch.Tensor:
    """Append (path length to site, is-site) for each site -> (k, n, 29).

    Upstream uses len(GetShortestPath(a, site)), which counts ATOMS on the path:
    bond distance + 1, and 0 for the site itself.
    """
    k, n = len(sites), base.shape[0]
    x = torch.empty(k, n, N_FEATURES)
    x[:, :, :27] = base
    for j, s in enumerate(sites):
        d = dist[s] + 1.0
        d[s] = 0.0
        x[j, :, 27] = d
        x[j, :, 28] = 0.0
        x[j, s, 28] = 1.0
    return x


# --------------------------------------------------------------------------
# ionizable sites: upstream utils/ionization_group.py, read from the checkout
# --------------------------------------------------------------------------
def _load_smarts(path: Path):
    acid, base = [], []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r = {k.strip(): v.strip() for k, v in r.items()}
            patt = Chem.MolFromSmarts(r["SMARTS"])
            idx = [int(i) for i in r["Index"].split(",")]
            (acid if r["Acid_or_base"] == "A" else base).append((patt, idx))
    return acid, base


def _unique(matches: list[list[int]]) -> list[int]:
    """Upstream unique_acid_match + flatten, order included."""
    singles = list({m[0] for m in matches if len(m) == 1})
    out = [m for m in matches if len(m) == 2] + [[j] for j in singles]
    return [j for m in out for j in m]


def find_sites(mol: Chem.Mol, acid_patterns, base_patterns):
    acid, base = [], []
    for patt, idx in acid_patterns:
        for m in mol.GetSubstructMatches(patt):
            acid.append([m[idx[0]], m[idx[1]]] if len(idx) > 1 else [m[idx[0]]])
    for patt, idx in base_patterns:
        for m in mol.GetSubstructMatches(patt):
            base.extend([[m[i]] for i in idx])
    return _unique(acid), _unique(base)


# --------------------------------------------------------------------------
# predictor
# --------------------------------------------------------------------------
@dataclass
class SitePka:
    idx: int            # heavy-atom index in the H-less, uncharged parent
    kind: str           # "acidic" (proton loss) | "basic" (proton gain)
    pka: float
    h_idx: int | None   # upstream's index (the H for acids) in the H-explicit mol


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class MolGpKa:
    def __init__(self, molgpka_dir: str | os.PathLike | None = None, threads: int = 1):
        root = Path(molgpka_dir or os.environ.get("MOLGPKA_DIR", ""))
        if not (root / "models" / "weight_acid.pth").exists():
            raise SystemExit(f"MolGpKa checkout not found at {root!s} (set MOLGPKA_DIR "
                             f"in config/hpc.env to a clone of github.com/Xundrug/MolGpKa)")
        torch.set_num_threads(threads)
        self.root = root
        self.weights_sha = {}
        self.models = {}
        for kind, fname in (("acidic", "weight_acid.pth"), ("basic", "weight_base.pth")):
            path = root / "models" / fname
            self.weights_sha[fname] = _sha256(path)
            if self.weights_sha[fname] != KNOWN_WEIGHTS[fname]:
                raise SystemExit(f"{path} is not the published MolGpKa weight file "
                                 f"(sha256 {self.weights_sha[fname]}); the port was "
                                 f"checked against Xundrug/MolGpKa 4dc8352 only")
            net = GCNNet()
            net.load_state_dict(torch.load(path, map_location="cpu", weights_only=True))
            net.eval()
            self.models[kind] = net
        self.acid_patterns, self.base_patterns = _load_smarts(
            root / "src" / "utils" / "smarts_pattern.tsv")
        self._uncharger = rdMolStandardize.Uncharger()

    def prepare(self, mol: Chem.Mol) -> Chem.Mol:
        """Upstream predict(): uncharge, round-trip through SMILES, add Hs."""
        mol = self._uncharger.uncharge(mol)
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        return AllChem.AddHs(mol)

    @torch.no_grad()
    def predict_prepared(self, molh: Chem.Mol, collapse: bool = True) -> list[SitePka]:
        acid, base = find_sites(molh, self.acid_patterns, self.base_patterns)
        if not acid and not base:
            return []
        feats = base_features(molh)
        dist = torch.tensor(Chem.GetDistanceMatrix(molh), dtype=torch.float32)
        a_hat = normalized_adjacency(molh)
        out: list[SitePka] = []
        for kind, sites in (("basic", base), ("acidic", acid)):
            if not sites:
                continue
            pkas = self.models[kind](site_features(feats, dist, sites), a_hat).tolist()
            seen: set[int] = set()
            for s, p in zip(sites, pkas):
                atom = molh.GetAtomWithIdx(s)
                heavy = atom.GetNeighbors()[0].GetIdx() if atom.GetAtomicNum() == 1 else s
                # an NH2 / SO2NH2 matches once per H; equivalent Hs give the
                # same pKa, and the state model needs one site per heavy atom
                if collapse and heavy in seen:
                    continue
                seen.add(heavy)
                out.append(SitePka(heavy, kind, float(p), s))
        return out

    def predict(self, mol: Chem.Mol) -> tuple[Chem.Mol, list[SitePka]]:
        """Return (heavy-atom parent the indices refer to, sites).

        AddHs appends H after every heavy atom, so heavy-atom indices in the
        H-explicit mol are the same as in RemoveHs(mol).
        """
        molh = self.prepare(mol)
        return Chem.RemoveHs(molh), self.predict_prepared(molh)


# --------------------------------------------------------------------------
# equivalence check against the upstream code
# --------------------------------------------------------------------------
def check_against_upstream(smiles: list[str], molgpka_dir: str) -> float:
    """Run upstream predict() and this port on the same molecules; return max |diff|."""
    import sys
    import types

    try:
        import torch_scatter  # noqa: F401
    except ImportError:      # upstream gcn_conv only needs scatter_add
        from torch_geometric.utils import scatter
        shim = types.ModuleType("torch_scatter")
        shim.scatter_add = lambda src, index, dim=0, dim_size=None: scatter(
            src, index, dim=dim, dim_size=dim_size, reduce="sum")
        sys.modules["torch_scatter"] = shim
    import torch_geometric.data as tgd
    if not hasattr(tgd, "DataLoader"):          # removed from torch_geometric.data
        from torch_geometric.loader import DataLoader
        tgd.DataLoader = DataLoader

    src = str(Path(molgpka_dir) / "src")
    here = os.getcwd()
    sys.path.insert(0, src)
    os.chdir(src)                               # upstream reads SMARTS from cwd
    try:
        import predict_pka as up
        _orig_load = torch.load
        torch.load = lambda f, map_location=None, **kw: _orig_load(
            f, map_location=map_location, weights_only=True)
        ours = MolGpKa(molgpka_dir)
        worst, n_sites = 0.0, 0
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            ub, ua = up.predict(mol)
            mine = {(s.kind, s.h_idx): s.pka
                    for s in ours.predict_prepared(ours.prepare(mol), collapse=False)}
            theirs = {("basic", k): float(v) for k, v in ub.items()}
            theirs.update({("acidic", k): float(v) for k, v in ua.items()})
            if set(mine) != set(theirs):
                raise SystemExit(f"site sets differ for {smi}:\n  upstream {sorted(theirs)}"
                                 f"\n  port     {sorted(mine)}")
            for k in theirs:
                worst = max(worst, abs(mine[k] - theirs[k]))
                n_sites += 1
        torch.load = _orig_load
    finally:
        os.chdir(here)
    print(f"{len(smiles)} molecules, {n_sites} sites, max |pKa(port) - pKa(upstream)| "
          f"= {worst:.2e}")
    return worst


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("smi", help="SMILES file (first column SMILES, optional name)")
    ap.add_argument("--molgpka-dir", default=os.environ.get("MOLGPKA_DIR"))
    ap.add_argument("--check", action="store_true",
                    help="compare against the upstream code (needs torch_geometric)")
    ap.add_argument("--tol", type=float, default=1e-3)
    a = ap.parse_args()
    smis = [line.split()[0] for line in open(a.smi) if line.strip()]
    if a.check:
        raise SystemExit(0 if check_against_upstream(smis, a.molgpka_dir) <= a.tol else 1)
    m = MolGpKa(a.molgpka_dir)
    for line in open(a.smi):
        if not line.strip():
            continue
        parts = line.split()
        mol = Chem.MolFromSmiles(parts[0])
        if mol is None:
            print(f"{parts[0]}\tUNPARSABLE")
            continue
        _, sites = m.predict(mol)
        desc = " ".join(f"{s.kind[0]}{s.idx}:{s.pka:.2f}" for s in sites) or "-"
        print(f"{parts[1] if len(parts) > 1 else parts[0]}\t{desc}")
