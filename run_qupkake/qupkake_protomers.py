#!/usr/bin/env python3
"""
qupkake_protomers.py
====================

Upstream protonation-state enumeration for a DOCK db2 build.

Pipeline
--------
    input .smi  ->  [standardize]  ->  qupkake  ->  micro-pKa sites
                ->  enumerate protomers at each pH
                ->  dedupe + label  ->  protomers.smi  ->  db2 builder

QupKake writes one SDF *record per reaction site*, with SD tags:
    idx       heavy-atom index of the site (0-based, into the mol in that record)
    pka_type  "acidic" = proton LOSS at that atom
              "basic"  = proton GAIN at that atom (pKa of the conjugate acid)
    pka       predicted micro-pKa
and the record's _Name is the molecule name.  The embedded mol carries explicit
hydrogens and is the tautomer QupKake actually used (which is NOT necessarily
the input tautomer when -t/--tautomerize is on), so all protomers are built
from that embedded mol rather than from the input SMILES.

State model
-----------
Independent-site (Henderson-Hasselbalch) approximation.  For each site the
fraction in the "changed" form at a given pH is

    acidic (proton loss):  f_deprot = 1 / (1 + 10^(pKa - pH))
    basic  (proton gain):  f_prot   = 1 / (1 + 10^(pH - pKa))

Sites with |pKa - pH| > margin are locked to their dominant form; sites inside
the margin are enumerated both ways.  A state's population is the product of
its per-site fractions over ALL sites.  This ignores site-site coupling, which
is the main systematic error for polyprotic molecules -- see README.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shlex
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import product
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.EnumerateStereoisomers import (EnumerateStereoisomers,
                                                StereoEnumerationOptions)
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

# db2 M-line writes name[-16:] (mol2db2/hierarchy.py), i.e. the LAST 16 chars.
DB2_NAME_LIMIT = 16

# db2_converter uses the ligand name unquoted as a directory name, a filename,
# and inside `subprocess.run(f"rm -r {zinc} {zinc}.smi", shell=True)`
# (db2_converter/db2_converter.py).  Anything outside this set is unsafe.
SAFE_NAME = re.compile(r"^[A-Za-z0-9_-]+$")

# Fixed-width name fields use base62 so one character covers 62 slots.
BASE62 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

# db2_converter refuses molecules with more than 2**5 stereoisomers
# (pipeline.py: "0enumerate_many").  Enumerating upstream bypasses that, but
# matching the cap keeps library growth sane.
DEFAULT_STEREO_CAP = 32


# Sites no predictor should be trusted on, dropped before any state is built
# (both paths: QupKake here, MolGpKa in pka_triage/triage.py).
#
#  amide_nh_acid  N-H of an amide, anilide, urea, carbamate, hydrazide or
#                 thioamide as an ACID. QupKake put 44 of them at 2.4-6.7 on
#                 Owen's set (2026-10-08), making an amide anion the dominant
#                 state at pH 7.4 in 22 of 50 parents; MolGpKa put some anilides
#                 at 6-7. Reference values (sources in pka_triage/TRIAGE.md):
#                 acetamide 15.1 in water (extrapolated) and 25.5 in DMSO;
#                 DMSO: N-methylacetamide 25.9, benzamide 23.3, acetanilide
#                 21.5, urea 26.9, ethyl carbamate 24.2, thioacetamide 18.5 --
#                 all well above the imides that ARE acidic at physiological
#                 pH (succinimide 14.7 in DMSO; phthalimide 8.30 in water).
#                 Kept: imide N-H (between two C=O) and N-sulfonyl N-H
#                 (acylsulfonamides, sulfonylureas; acidic, value not sourced).
#  amide_n_base   any N on C=O/C=S or S(=O)=O as a BASE. Amides protonate on
#                 O, not N, and even that is pKa -0.62 (acetamide, water);
#                 MolGpKa's patterns put anilides at 3.9-6.
# --keep-amide-sites turns both off.
IMPLAUSIBLE_SITES = {
    ("acidic", "amide_nh_acid"): Chem.MolFromSmarts(
        "[#7;!H0;$([#7]-[#6]=[#8,#16]);!$([#7](-[#6]=[#8,#16])-[#6]=[#8,#16]);"
        "!$([#7]-[#16](=[#8])=[#8])]"),
    ("basic", "amide_n_base"): Chem.MolFromSmarts(
        "[#7;$([#7]-[#6]=[#8,#16]),$([#7]-[#16](=[#8])=[#8])]"),
}


def implausible_sites(mol: Chem.Mol) -> dict[tuple[int, str], str]:
    """{(atom index, kind): reason} for sites on `mol` that are dropped."""
    out = {}
    for (kind, reason), patt in IMPLAUSIBLE_SITES.items():
        for m in mol.GetSubstructMatches(patt):
            out[(m[0], kind)] = reason
    return out


# --------------------------------------------------------------------------
# data model
# --------------------------------------------------------------------------
@dataclass
class Site:
    idx: int
    kind: str          # "acidic" | "basic"
    pka: float

    def changed_fraction(self, ph: float) -> float:
        """Fraction of molecules where this site differs from the reference."""
        if self.kind == "acidic":                       # proton loss
            return 1.0 / (1.0 + 10.0 ** (self.pka - ph))
        return 1.0 / (1.0 + 10.0 ** (ph - self.pka))    # proton gain

    @property
    def delta(self) -> int:
        return -1 if self.kind == "acidic" else +1


@dataclass
class Protomer:
    parent: str
    smiles: str
    charge: int
    ph_values: list[float] = field(default_factory=list)
    population: float = 0.0          # best population across the pH values
    changes: tuple = ()
    name: str = ""
    note: str = ""
    stereo_index: int = 0
    protomer_index: int = 0

    @property
    def n_pos(self) -> int:
        return self._count(lambda c: c > 0)

    @property
    def n_neg(self) -> int:
        return self._count(lambda c: c < 0)

    def _count(self, pred) -> int:
        m = Chem.MolFromSmiles(self.smiles)
        if m is None:
            return 0
        return sum(1 for a in m.GetAtoms() if pred(a.GetFormalCharge()))

    @property
    def is_zwitterion(self) -> bool:
        return self.charge == 0 and self.n_pos > 0 and self.n_neg > 0


# --------------------------------------------------------------------------
# input / standardization
# --------------------------------------------------------------------------
# protomers.csv, read by assign_names.py and pka_triage/merge_protomers.py.
# pka_triage/triage.py writes the same columns plus its own at the end.
PROTOMER_COLUMNS = ["protomer_name", "parent_name", "smiles", "net_charge",
                    "n_pos_atoms", "n_neg_atoms", "zwitterion", "ph_values",
                    "population_estimate", "n_changes", "protomer_index",
                    "stereo_index", "note"]


def protomer_row(p: Protomer) -> list:
    return [p.name, p.parent, p.smiles, p.charge, p.n_pos, p.n_neg,
            int(p.is_zwitterion), ";".join(f"{v:g}" for v in sorted(set(p.ph_values))),
            f"{p.population:.4f}", len(p.changes), p.protomer_index, p.stereo_index,
            p.note]


def read_smi(path: Path, delim: str | None = None) -> list[tuple[str, str]]:
    """Read a 2-column SMILES file (SMILES NAME), no header."""
    out, seen = [], Counter()
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(delim) if delim else line.split()
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"mol_{lineno}"
            name, changed = sanitize_name(name)
            if changed:
                print(f"  ! unsafe name on line {lineno}: {changed}",
                      file=sys.stderr)
            seen[name] += 1
            if seen[name] > 1:
                name = f"{name}_dup{seen[name]}"
                print(f"  ! duplicate name on line {lineno}, renamed to {name}",
                      file=sys.stderr)
            out.append((smi, name))
    return out


_UNCHARGER = rdMolStandardize.Uncharger()


def standardize(smi: str, neutralize: bool = True, largest_fragment: bool = True):
    """Clean an input SMILES.  QupKake expects neutral, single-fragment input."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    if largest_fragment:
        mol = rdMolStandardize.FragmentParent(mol)
    mol = rdMolStandardize.Cleanup(mol)
    if neutralize:
        mol = _UNCHARGER.uncharge(mol)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


# --------------------------------------------------------------------------
# QupKake invocation + parsing
# --------------------------------------------------------------------------
class QupKakeProgress:
    """Periodic one-line progress for a running QupKake, read off its files.

    QupKake's tqdm bars name only the molecule in hand, never the stage, and
    render as carriage-return noise in a SLURM log. Its working tree says more,
    in the same files audit_qupkake() reads:

      stage 1  featurize molecules  processed/{name}.pt per molecule done
                                    (xtb --opt --alpb water --lmo, + --vfukui)
      stage 2  site pairs           raw/{output} = every predicted site (the
                                    total); processed/{name}_{idx}_{type}_pair.pt
                                    per site done (molecule + conjugate
                                    re-featurized: ~4 xtb runs each)
      errors                        'Error processing' lines in logs/error_log.txt

    Prints a line every `interval` seconds and a per-stage timing summary at
    the end -- the numbers a timing preflight is run for.
    """

    def __init__(self, root: Path, output: str, n_mols: int, interval: float):
        import threading
        self.root, self.output, self.n_mols, self.interval = root, output, n_mols, interval
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self._loop, daemon=True)
        self.t0 = self._now()
        self.t_stage2 = None
        self.base1, self.base2 = self._counts()[:2]   # work reused from an earlier run

    @staticmethod
    def _now() -> float:
        # wall clock, so it compares with file mtimes: the stage boundary is
        # when QupKake WROTE the site file, not when this loop next looked
        import time
        return time.time()

    def _counts(self) -> tuple[int, int, int | None, int]:
        proc = self.root / "processed"
        mols = pairs = 0
        if proc.is_dir():
            for f in os.scandir(proc):
                if not f.name.endswith(".pt") or f.name.startswith("pre_"):
                    continue
                if f.name.endswith("_pair.pt"):
                    pairs += 1
                else:
                    mols += 1
        raw = self.root / "raw" / self.output
        sites = None
        if raw.exists():
            with open(raw, errors="replace") as fh:
                sites = sum(1 for line in fh if line.startswith("$$$$"))
        errs = 0
        log = self.root / "logs" / "error_log.txt"
        if log.exists():
            with open(log, errors="replace") as fh:
                errs = sum(1 for line in fh if "Error processing" in line)
        return mols, pairs, sites, errs

    @staticmethod
    def _fmt(sec: float) -> str:
        if sec < 60:
            return f"{sec:.1f}s"
        sec = int(round(sec))
        return f"{sec // 3600}:{sec % 3600 // 60:02d}:{sec % 60:02d}"

    def _rate(self, done: int, since: float) -> str:
        if done <= 0:
            return ""
        per = (self._now() - since) / done
        return f", {per:.1f} s/item"

    def line(self) -> str:
        mols, pairs, sites, errs = self._counts()
        now = self._now()
        if sites is None:
            stage = (f"stage 1/2 featurize molecules {mols}/{self.n_mols}"
                     f"{self._rate(mols - self.base1, self.t0)}")
        else:
            self._mark_stage2()
            stage = (f"stage 1/2 done ({mols}/{self.n_mols}) | stage 2/2 site pairs "
                     f"{pairs}/{sites}{self._rate(pairs - self.base2, self.t_stage2)}")
            done = pairs - self.base2
            if 0 < done < sites:
                eta = (now - self.t_stage2) / done * (sites - pairs)
                stage += f", ~{self._fmt(eta)} left"
        err = f" | errors logged: {errs}" if errs else ""
        return f"  [qupkake {self._fmt(now - self.t0)}] {stage}{err}"

    def _mark_stage2(self):
        if self.t_stage2 is None:
            raw = self.root / "raw" / self.output
            # a site file left by an earlier run predates t0: stage 2 starts now
            self.t_stage2 = max(self.t0, raw.stat().st_mtime) if raw.exists() else self._now()

    def _loop(self):
        while not self.stop.wait(self.interval):
            print(("\n" if sys.stderr.isatty() else "") + self.line(),
                  file=sys.stderr, flush=True)

    def __enter__(self):
        if self.interval > 0:
            self.thread.start()
        return self

    def __exit__(self, *exc):
        self.stop.set()
        if self.thread.is_alive():
            self.thread.join()
        mols, pairs, sites, errs = self._counts()
        end = self._now()
        if sites is not None:
            self._mark_stage2()
        s1_end = self.t_stage2 if self.t_stage2 is not None else end
        print(f"  [qupkake timing] stage 1 featurize: {mols}/{self.n_mols} molecules "
              f"in {self._fmt(s1_end - self.t0)}; stage 2 site pairs: {pairs}/"
              f"{sites if sites is not None else 0} in "
              f"{self._fmt(end - s1_end)}; total {self._fmt(end - self.t0)}; "
              f"errors logged {errs}", file=sys.stderr, flush=True)
        return False


def run_qupkake(csv_path: Path, root: Path, output: str, tautomerize: bool,
                nproc: int | None, exe: str = "qupkake",
                n_mols: int = 0, progress_interval: float = 60) -> Path:
    cmd = [exe, "file", str(csv_path),
           "-o", output,
           "-s", "smiles",
           "-n", "name",
           "-r", str(root)]
    if tautomerize:
        cmd.append("-t")
    if nproc:
        cmd += ["-mp", str(nproc)]
    print("  $ " + " ".join(shlex.quote(c) for c in cmd))
    with QupKakeProgress(root, output, n_mols, progress_interval):
        # xtb --vfukui takes no -P and torch sizes itself to the node, so
        # pin them; -mp/-P already sets xtb's main calculation.
        env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        res = subprocess.run(cmd, env=env)
    if res.returncode != 0:
        raise RuntimeError(f"qupkake exited with status {res.returncode}")
    out = root / "output" / output
    if not out.exists():
        raise FileNotFoundError(
            f"expected {out} but it was not created -- QupKake prints "
            "'No protonation/deprotonation sites were found' and skips the "
            "output file when nothing is ionizable in the whole input set.")
    return out


def parse_qupkake_sdf(path: Path) -> dict[str, tuple[Chem.Mol, list[Site]]]:
    """Group the per-site SDF records by molecule name."""
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    grouped: dict[str, tuple[Chem.Mol, list[Site]]] = {}
    bad = 0
    for mol in supplier:
        if mol is None:
            bad += 1
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else "unnamed"
        try:
            site = Site(idx=int(float(mol.GetProp("idx"))),
                        kind=mol.GetProp("pka_type").strip().lower(),
                        pka=float(mol.GetProp("pka")))
        except KeyError:
            continue
        if name not in grouped:
            grouped[name] = (mol, [])
        grouped[name][1].append(site)
    if bad:
        print(f"  ! {bad} SDF record(s) failed to parse", file=sys.stderr)
    return grouped


def audit_qupkake(root: Path, sdf_name: str, names: list[str],
                  grouped: dict) -> dict[str, str] | None:
    """Which molecules did QupKake FAIL on, as opposed to find nothing in?

    QupKake catches per-molecule and per-site exceptions, logs them to
    logs/error_log.txt and drops the row (mol_dataset.py, _process_chunk).
    Absent from the output SDF therefore means either "no ionizable site" or
    "crashed", and treating a crash as no-sites carries an amine through
    neutral -- a wrong protomer that looks like a clean one.  QupKake's own
    working files separate the cases:

      processed/{name}.pt       written for every molecule it featurized
      raw/{sdf_name}            every predicted site, BEFORE pKa prediction
      output/{sdf_name}         the sites that survived pKa prediction

    no .pt                      -> qupkake_failed        (molecule never processed)
    site in raw, not in output  -> qupkake_site_failed   (states would be built
                                   from an incomplete site list)

    Returns None when the working tree is not there to audit (--from-sdf on a
    bare SDF), so the caller can say the audit did not run.
    """
    processed, raw = root / "processed", root / "raw" / sdf_name
    if not processed.is_dir():
        return None
    failed: dict[str, str] = {}
    for n in names:
        if not (processed / f"{n}.pt").exists():
            failed[n] = "qupkake_failed"
    if raw.exists():
        predicted = defaultdict(set)
        for mol in Chem.SDMolSupplier(str(raw), removeHs=False, sanitize=False):
            if mol is None or not mol.HasProp("idx"):
                continue
            predicted[mol.GetProp("_Name")].add(
                (int(float(mol.GetProp("idx"))), mol.GetProp("pka_type").strip().lower()))
        for n, sites in predicted.items():
            got = {(s.idx, s.kind) for s in grouped.get(n, (None, []))[1]}
            if sites - got and n not in failed:
                failed[n] = "qupkake_site_failed"
    return failed


# --------------------------------------------------------------------------
# protomer construction
# --------------------------------------------------------------------------
def reconcile_stereo(mol: Chem.Mol, input_mol: Chem.Mol | None, mode: str):
    """Make the protomer's stereochemistry match the *input* SMILES.

    QupKake embeds every molecule in 3D, so reading its SDF back gives chiral
    tags derived from an arbitrary conformer -- RDKit will happily report
    C[C@@H](O)C for an input that specified no stereocentre at all.  Feeding
    that to a db2 build silently commits the whole screen to one enantiomer.

    Heavy-atom order is preserved through QupKake (MolFromSmiles -> AddHs ->
    embed appends H at the end), which is the same assumption that makes the
    `idx` tag meaningful, so tags can be copied index-for-index.
    """
    if mode == "keep":
        return mol
    mol.RemoveAllConformers()

    if mode == "strip" or input_mol is None:
        Chem.RemoveStereochemistry(mol)
        return mol

    heavy = [a for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    if len(heavy) != input_mol.GetNumAtoms() or any(
            a.GetAtomicNum() != input_mol.GetAtomWithIdx(i).GetAtomicNum()
            for i, a in enumerate(heavy)):
        # tautomerization or a rewritten skeleton broke the correspondence
        Chem.RemoveStereochemistry(mol)
        return mol

    Chem.RemoveStereochemistry(mol)
    for i, atom in enumerate(heavy):
        atom.SetChiralTag(input_mol.GetAtomWithIdx(i).GetChiralTag())
    for bond in input_mol.GetBonds():
        if bond.GetStereo() == Chem.BondStereo.STEREONONE:
            continue
        tgt = mol.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if tgt is None:
            continue
        tgt.SetStereoAtoms(*bond.GetStereoAtoms())
        tgt.SetStereo(bond.GetStereo())
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def apply_changes(ref: Chem.Mol, changes: Iterable[tuple[int, int]],
                  input_mol: Chem.Mol | None = None, stereo: str = "input"):
    """Return a new mol with H count / formal charge adjusted at each site.

    Mirrors QupKake's own Conjugate class: change the formal charge by +/-1 and
    the hydrogen count by the same amount on the heavy atom at `idx`.
    """
    rw = Chem.RWMol(ref)
    doomed_hs: list[int] = []

    for idx, delta in changes:
        if idx >= rw.GetNumAtoms():
            return None
        atom = rw.GetAtomWithIdx(idx)

        if delta > 0:                                   # protonate
            if not atom.GetNoImplicit():
                base = atom.GetTotalNumHs(includeNeighbors=False)
                atom.SetNoImplicit(True)
            else:
                base = atom.GetNumExplicitHs()
            atom.SetNumExplicitHs(base + 1)
            atom.SetFormalCharge(atom.GetFormalCharge() + 1)

        else:                                           # deprotonate
            h = next((nb.GetIdx() for nb in atom.GetNeighbors()
                      if nb.GetAtomicNum() == 1
                      and nb.GetDegree() == 1
                      and nb.GetIdx() not in doomed_hs), None)
            if h is not None:
                doomed_hs.append(h)
            else:
                if not atom.GetNoImplicit():
                    base = atom.GetTotalNumHs(includeNeighbors=False)
                    atom.SetNoImplicit(True)
                else:
                    base = atom.GetNumExplicitHs()
                if base <= 0:
                    return None                         # no proton to remove
                atom.SetNumExplicitHs(base - 1)
            atom.SetFormalCharge(atom.GetFormalCharge() - 1)

        atom.SetNumRadicalElectrons(0)

    for h in sorted(doomed_hs, reverse=True):
        rw.RemoveAtom(h)

    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
        mol = reconcile_stereo(mol, input_mol, stereo)
        mol = Chem.RemoveHs(mol)
        smi = Chem.MolToSmiles(mol)
        round_trip = Chem.MolFromSmiles(smi)
    except Exception:
        return None
    if round_trip is None:
        return None
    return Chem.MolToSmiles(round_trip)


def enumerate_states(sites: list[Site], ph: float, margin: float,
                     max_ambiguous: int):
    """Yield (changes, population) for every state worth building at this pH."""
    locked: list[tuple[int, int]] = []
    ambiguous: list[Site] = []
    base_pop = 1.0

    for s in sites:
        f = s.changed_fraction(ph)
        if abs(s.pka - ph) > margin:
            if f >= 0.5:
                locked.append((s.idx, s.delta))
                base_pop *= f
            else:
                base_pop *= (1.0 - f)
        else:
            ambiguous.append(s)

    # cap the combinatorics: keep the sites closest to the pH
    if len(ambiguous) > max_ambiguous:
        ambiguous.sort(key=lambda s: abs(s.pka - ph))
        overflow, ambiguous = ambiguous[max_ambiguous:], ambiguous[:max_ambiguous]
        for s in overflow:
            f = s.changed_fraction(ph)
            if f >= 0.5:
                locked.append((s.idx, s.delta))
                base_pop *= f
            else:
                base_pop *= (1.0 - f)

    for combo in product([False, True], repeat=len(ambiguous)):
        changes = list(locked)
        pop = base_pop
        for flip, s in zip(combo, ambiguous):
            f = s.changed_fraction(ph)
            if flip:
                changes.append((s.idx, s.delta))
                pop *= f
            else:
                pop *= (1.0 - f)
        yield tuple(sorted(changes)), pop


def unchanged_protomers(name: str, isomers: list[Chem.Mol],
                        ph_values: list[float]) -> list[Protomer]:
    """No ionizable site predicted: carry each isomer through unchanged."""
    return [Protomer(parent=name, smiles=Chem.MolToSmiles(iso),
                     charge=Chem.GetFormalCharge(iso), ph_values=list(ph_values),
                     population=1.0, note="no_sites_predicted", stereo_index=k)
            for k, iso in enumerate(isomers)]


def build_protomers(name: str, isomers: list[Chem.Mol], ref: Chem.Mol,
                    sites: list[Site], ph_values: list[float], margin: float,
                    min_population: float, max_states: int, max_ambiguous: int,
                    stereo: str, fallback_note: str = "built_on_qupkake_tautomer"):
    """One parent's sites -> its protomers at every pH, for every stereoisomer.

    `ref` is the molecule the site indices refer to (QupKake's SDF mol, or
    MolGpKa's re-parsed parent); `isomers` are the parent's input stereoisomers.
    Returns (protomers, states dropped by the filters, states that failed to
    build, whether the isomers themselves could be used).
    """
    # Build on the isomer itself where possible: it already carries the
    # right stereo in the site predictor's atom order, so nothing has to be
    # transferred and E/Z survives.  Fall back to `ref` (with stereo
    # reconciliation) when the predictor rewrote the molecule (QupKake -t).
    use_isomer = skeleton_matches(ref, isomers[0])
    dropped_states = build_failures = 0

    # 1. collect unique protonation STATES for the parent, judged on the
    #    first isomer so a state gets one index across every isomer
    states: dict[tuple, Protomer] = {}
    for ph in ph_values:
        found = []
        for changes, pop in enumerate_states(sites, ph, margin, max_ambiguous):
            if pop < min_population:
                dropped_states += 1
                continue
            probe = apply_changes(isomers[0] if use_isomer else ref, changes,
                                  isomers[0], stereo)
            if probe is None:
                build_failures += 1
                continue
            found.append((changes, pop, probe))
        found.sort(key=lambda x: -x[1])
        if len(found) > max_states:
            dropped_states += len(found) - max_states
            found = found[:max_states]
        for changes, pop, probe in found:
            if changes in states:
                states[changes].ph_values.append(ph)
                states[changes].population = max(states[changes].population, pop)
            else:
                states[changes] = Protomer(
                    parent=name, smiles=probe,
                    charge=Chem.GetFormalCharge(Chem.MolFromSmiles(probe)),
                    ph_values=[ph], population=pop, changes=changes)

    # 2. apply every surviving state to every stereoisomer
    emitted: list[Protomer] = []
    ordered = sorted(states.values(), key=lambda x: (-x.population, x.smiles))
    for p_idx, state in enumerate(ordered):
        seen_here: set = set()
        for k, iso in enumerate(isomers):
            out_smi = apply_changes(iso if use_isomer else ref, state.changes,
                                    iso, stereo)
            if out_smi is None or (k, out_smi) in seen_here:
                continue
            seen_here.add((k, out_smi))
            emitted.append(Protomer(
                parent=name, smiles=out_smi,
                charge=Chem.GetFormalCharge(Chem.MolFromSmiles(out_smi)),
                ph_values=list(state.ph_values), population=state.population,
                changes=state.changes, stereo_index=k, protomer_index=p_idx,
                note="" if use_isomer else fallback_note))
    return emitted, dropped_states, build_failures, use_isomer


# --------------------------------------------------------------------------
# naming
# --------------------------------------------------------------------------
def enumerate_stereoisomers(canonical_smiles: str, cap: int,
                            try_embedding: bool = False):
    """Stereoisomers of the *re-parsed canonical* SMILES, in QupKake's order.

    QupKake builds its molecule with MolFromSmiles(smiles) -> AddHs -> embed, so
    heavy atoms 0..N-1 of a mol parsed from the same string line up with the
    `idx` tags in its output.  Enumerating from that same re-parsed string means
    every isomer inherits that ordering, so a protonation change can be applied
    straight to the isomer -- no stereo transfer, and E/Z survives, which it
    does not when chiral tags are copied atom-by-atom onto a rebuilt skeleton.

    onlyUnassigned=True leaves centres the input already specified alone, so a
    partially-defined molecule keeps what it declared and only fills the gaps.
    """
    ref = Chem.MolFromSmiles(canonical_smiles)
    if ref is None:
        return [], 0
    opts = StereoEnumerationOptions(onlyUnassigned=True,
                                    tryEmbedding=try_embedding)
    seen, isomers = set(), []
    for iso in EnumerateStereoisomers(ref, options=opts):
        smi = Chem.MolToSmiles(iso)
        if smi in seen:
            continue
        seen.add(smi)
        isomers.append(iso)
        if len(isomers) > cap:
            return isomers[:cap], len(isomers)
    return isomers, len(isomers)


def _heavy_signature(m: Chem.Mol):
    """Per-heavy-atom fingerprint used to confirm two mols share atom ORDER.

    Element sequence alone is not enough: benzoic acid written as OC(=O)c1ccccc1
    and as its canonical O=C(O)c1ccccc1 have identical element sequences but
    swapped oxygens, so a site index meaning "the hydroxyl O" in one points at
    the carbonyl O in the other.  Including degree and H count catches that.
    """
    # heavy-neighbour count, not GetDegree(): the QupKake SDF carries explicit
    # H atoms in the graph while a SMILES-parsed mol does not, so raw degree
    # would differ for every molecule
    return [(a.GetAtomicNum(),
             sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1),
             a.GetTotalNumHs(includeNeighbors=True), a.GetFormalCharge())
            for a in m.GetAtoms() if a.GetAtomicNum() > 1]


def skeleton_matches(a: Chem.Mol, b: Chem.Mol) -> bool:
    """Do two mols agree atom-for-atom, so site indices mean the same thing?"""
    return _heavy_signature(a) == _heavy_signature(b)


def sanitize_name(name: str) -> tuple[str, str | None]:
    """Make a parent name safe for db2_converter.

    The name becomes a directory, a filename, and an argument to an unquoted
    `rm -r` run with shell=True, so anything outside [A-Za-z0-9_-] is a
    correctness-or-worse hazard.  Dots are stripped too: --mergeiso does
    `name.split(".")[0]`, and --checkstereo appends `.0`/`.1` for enumerated
    stereoisomers, so a dot in the parent name silently truncates the db2 name.
    """
    clean = re.sub(r"[^A-Za-z0-9_-]", "_", name)
    if clean != name:
        return clean, f"{name!r} -> {clean!r}"
    return clean, None


def format_ph(ph: float) -> str:
    return f"{ph:g}".replace(".", "p")


def base62(n: int) -> str:
    return BASE62[n] if n < 62 else "?"


def base62_serial(n: int) -> str:
    if n == 0:
        return "0"
    out = ""
    while n:
        n, r = divmod(n, 62)
        out = BASE62[r] + out
    return out


def plan_names(parent_names: list[str], style: str) -> tuple[dict[str, str], str]:
    """Choose a base stem per parent that leaves 2 characters for the suffix.

    The db2 M line stores `name[-16:]`, so a 16-character ZINC ID plus any
    suffix loses its front: ZINC000012345678_p1 truncates to 000012345678_p1,
    or worse, two parents collide.  'short' rewrites the leading ZINC to a
    single C -- C + 12 digits = 13 -- leaving room for a protomer character and
    a stereo character inside 16.

    The suffix carries NO dot.  --mergeiso sets the db2 name to
    name.split(".")[0] (pipeline.py), which would collapse every protomer and
    every stereoisomer of a molecule into a single db2 name.  With no dot,
    --mergeiso is a harmless no-op.

    If any parent will not fit, all parents switch to a compact serial, so the
    whole set stays one fixed-width scheme instead of a mix.
    """
    if style != "short":
        return {n: n for n in parent_names}, style

    stems, ok = {}, True
    for n in parent_names:
        stem = "C" + n[4:] if n.upper().startswith("ZINC") else n
        stems[n] = stem
        if len(stem) + 2 > DB2_NAME_LIMIT:
            ok = False
    if ok:
        return stems, "short"

    width = max(2, len(base62_serial(max(0, len(parent_names) - 1))))
    return ({n: "C" + base62_serial(i).rjust(width, "0")
             for i, n in enumerate(sorted(parent_names))}, "serial")


def make_name(stem: str, p_idx: int, s_idx: int, style: str,
              protomer: "Protomer", n_isomers: int) -> str:
    """Build a protomer name.  Never contains a dot (see plan_names)."""
    if style in ("short", "serial"):
        return f"{stem}{base62(p_idx)}{base62(s_idx)}"

    suffix = f"_p{p_idx + 1}"
    if n_isomers > 1:
        suffix += f"s{base62(s_idx)}"
    if style == "index":
        return f"{stem}{suffix}"
    if style == "charge":
        q = protomer.charge
        tag = "0" if q == 0 else f"{'p' if q > 0 else 'm'}{abs(q)}"
        return f"{stem}_q{tag}{suffix}"
    if style == "ph":
        phs = "-".join(format_ph(x) for x in sorted(set(protomer.ph_values)))
        return f"{stem}_pH{phs}{suffix}"
    raise ValueError(style)


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------
def build_stats(protomers: list[Protomer], grouped, inputs, ph_values,
                site_records) -> dict:
    per_parent = defaultdict(list)
    for p in protomers:
        per_parent[p.parent].append(p)

    counts = [len(v) for v in per_parent.values()] or [0]
    charges = Counter(p.charge for p in protomers)
    ph_counter = Counter()
    for p in protomers:
        for ph in p.ph_values:
            ph_counter[ph] += 1

    all_ph = {p.smiles for p in protomers if len(p.ph_values) == len(ph_values)}
    single_ph = [p for p in protomers if len(p.ph_values) == 1]
    pkas = [s.pka for _, sites in grouped.values() for s in sites]
    heavy = [Chem.MolFromSmiles(p.smiles).GetNumHeavyAtoms() for p in protomers]
    long_names = [p.name for p in protomers if len(p.name) > DB2_NAME_LIMIT]

    # names that collide once db2 truncates to the last 16 characters
    trunc = Counter(p.name[-DB2_NAME_LIMIT:] for p in protomers)
    collisions = {k: v for k, v in trunc.items() if v > 1}

    return {
        "input": {
            "molecules_read": len(inputs),
            "molecules_with_predicted_sites": len(grouped),
            "molecules_without_sites": len(inputs) - len(grouped),
            "site_records_parsed": site_records,
        },
        "sites": {
            "total": len(pkas),
            "acidic": sum(1 for _, s in grouped.values() for x in s
                          if x.kind == "acidic"),
            "basic": sum(1 for _, s in grouped.values() for x in s
                         if x.kind == "basic"),
            "pka_min": round(min(pkas), 2) if pkas else None,
            "pka_median": round(statistics.median(pkas), 2) if pkas else None,
            "pka_max": round(max(pkas), 2) if pkas else None,
            "sites_per_molecule_mean": round(len(pkas) / len(grouped), 2)
                                       if grouped else 0.0,
        },
        "protomers": {
            "unique_total": len(protomers),
            "per_ph": {str(ph): ph_counter[ph] for ph in ph_values},
            "present_at_all_ph": len(all_ph),
            "unique_to_one_ph": len(single_ph),
            "expansion_factor": round(len(protomers) / len(inputs), 2)
                                if inputs else 0.0,
        },
        "charge": {
            "cationic": sum(1 for p in protomers if p.charge > 0),
            "anionic": sum(1 for p in protomers if p.charge < 0),
            "neutral": sum(1 for p in protomers if p.charge == 0),
            "zwitterionic": sum(1 for p in protomers if p.is_zwitterion),
            "distribution": {str(k): v for k, v in sorted(charges.items())},
            "net_charge_min": min(charges) if charges else None,
            "net_charge_max": max(charges) if charges else None,
        },
        "per_molecule": {
            "min": min(counts),
            "max": max(counts),
            "mean": round(sum(counts) / len(counts), 2),
            "median": statistics.median(counts),
            "histogram": {str(k): v for k, v in sorted(Counter(counts).items())},
            "single_protomer_molecules": sum(1 for c in counts if c == 1),
        },
        "size": {
            "heavy_atoms_min": min(heavy) if heavy else None,
            "heavy_atoms_mean": round(sum(heavy) / len(heavy), 1) if heavy else None,
            "heavy_atoms_max": max(heavy) if heavy else None,
        },
        "db2_naming": {
            "name_limit": DB2_NAME_LIMIT,
            "names_over_limit": len(long_names),
            "truncation_collisions": collisions,
        },
    }


def print_stats(stats: dict) -> None:
    def section(title, d, indent="  "):
        print(f"\n{title}")
        for k, v in d.items():
            if isinstance(v, dict):
                if not v:
                    continue
                print(f"{indent}{k}:")
                for k2, v2 in v.items():
                    print(f"{indent}  {k2:<24} {v2}")
            else:
                print(f"{indent}{k:<28} {v}")

    print("\n" + "=" * 64)
    print("PROTOMER ENUMERATION SUMMARY")
    print("=" * 64)
    for title, key in [("INPUT", "input"), ("PREDICTED SITES", "sites"),
                       ("PROTOMERS", "protomers"), ("CHARGE", "charge"),
                       ("STEREO", "stereo"),
                       ("PER MOLECULE", "per_molecule"), ("SIZE", "size"),
                       ("DB2 NAMING", "db2_naming")]:
        if key in stats:
            section(title, stats[key])
    print()


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="QupKake -> protonation-state enumeration -> db2 input",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ap.add_argument("input", help="2-column SMILES file: 'SMILES name', no header")
    ap.add_argument("-o", "--outdir", default="protomers", help="output directory")
    ap.add_argument("--ph", type=float, nargs="+", default=[6.4, 7.4, 8.4],
                    help="pH values to enumerate")
    ap.add_argument("--margin", type=float, default=1.0,
                    help="enumerate both forms when |pKa - pH| is within this")
    ap.add_argument("--min-population", type=float, default=0.01,
                    help="drop states below this estimated fractional population")
    ap.add_argument("--max-states", type=int, default=8,
                    help="max protomers kept per molecule per pH (highest population)")
    ap.add_argument("--max-ambiguous", type=int, default=6,
                    help="max simultaneously-enumerated sites (2^N states)")
    ap.add_argument("--name-style",
                    choices=["short", "index", "charge", "ph"],
                    default="short",
                    help="'short' = fixed-width <stem><protomer><stereo> that "
                         "fits db2's 16-char field (ZINC -> C); others are "
                         "readable but may be truncated by db2")
    ap.add_argument("--max-stereo", type=int, default=DEFAULT_STEREO_CAP,
                    help="cap on stereoisomers per parent (--stereo enumerate)")
    ap.add_argument("--stereo-embed", action="store_true",
                    help="discard stereoisomers RDKit cannot embed in 3D "
                         "(slower, drops impossible ring configurations)")
    ap.add_argument("--stereo", choices=["carry", "enumerate", "keep", "strip"],
                    default="carry",
                    help="'carry' = keep the input's stereo and let "
                         "db2_converter enumerate anything undefined; "
                         "'enumerate' = expand undefined centres here, run "
                         "QupKake once per parent, propagate protomers across "
                         "the isomers; 'keep' = trust the SDF's 3D-derived "
                         "tags; 'strip' = flatten")
    ap.add_argument("--keep-amide-sites", action="store_true",
                    help="keep QupKake's amide N-H acid / amide N base sites "
                         "(dropped by default; see IMPLAUSIBLE_SITES)")
    ap.add_argument("--no-neutralize", action="store_true",
                    help="skip input neutralization (QupKake wants neutral input)")
    ap.add_argument("--keep-largest-fragment", action="store_true", default=True)
    ap.add_argument("-t", "--tautomerize", action="store_true",
                    help="pass -t to QupKake (xtb tautomer search; much slower)")
    ap.add_argument("--progress-interval", type=float, default=60,
                    help="seconds between QupKake stage/progress lines on stderr "
                         "(0 = off); a per-stage timing summary is always printed")
    ap.add_argument("-mp", "--nproc", type=int, default=1,
                    help="QupKake workers. Keep 1: QupKake also passes -P N to every "
                         "xtb, so N workers run N*N threads; parallelise with more "
                         "processes (array tasks) instead. Unset in QupKake itself, "
                         "this would be every core on the node.")
    ap.add_argument("--qupkake-exe", default="qupkake")
    ap.add_argument("--from-sdf", default=None,
                    help="reuse an existing QupKake output SDF; skip prediction")
    ap.add_argument("--max-conf", type=int, default=None,
                    help="emit a third column in protomers.smi with this "
                         "per-molecule conformer budget; db2_converter reads "
                         "'SMILES NAME CONF_NUM' and it overrides -n")
    ap.add_argument("--db2-cmd", default=None,
                    help="command template to run afterwards; {smi} is substituted, "
                         "e.g. 'build_ligand -i {smi} -m ccdc "
                         "--workingpath wp --outputpath out'")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    work = outdir / "qupkake_work"

    # ---- 1. read + standardize -------------------------------------------
    print(f"[1/5] reading {args.input}")
    raw = read_smi(Path(args.input))
    inputs: list[tuple[str, str]] = []
    input_mols: dict[str, Chem.Mol] = {}
    failed = []
    for smi, name in raw:
        mol = standardize(smi, neutralize=not args.no_neutralize,
                          largest_fragment=args.keep_largest_fragment)
        if mol is None:
            failed.append((smi, name))
            continue
        inputs.append((Chem.MolToSmiles(mol), name))
        input_mols[name] = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    print(f"      {len(inputs)} standardized, {len(failed)} rejected")

    # ---- 1b. stereoisomer expansion --------------------------------------
    stereo_sets: dict[str, list[Chem.Mol]] = {}
    stereo_overflow: list[str] = []
    for smi, name in inputs:
        if args.stereo == "enumerate":
            isomers, total = enumerate_stereoisomers(smi, args.max_stereo,
                                                     args.stereo_embed)
            if total > args.max_stereo:
                stereo_overflow.append(f"{name} ({total})")
            stereo_sets[name] = isomers or [Chem.MolFromSmiles(smi)]
        else:
            stereo_sets[name] = [Chem.MolFromSmiles(smi)]
    if args.stereo == "enumerate":
        n_iso = sum(len(v) for v in stereo_sets.values())
        multi = sum(1 for v in stereo_sets.values() if len(v) > 1)
        print(f"      {n_iso} stereoisomer(s) from {len(inputs)} parent(s); "
              f"{multi} parent(s) had undefined centres")
        print(f"      QupKake will still run {len(inputs)} job(s), "
              f"saving {n_iso - len(inputs)}")
        if stereo_overflow:
            print(f"  ! capped at {args.max_stereo} isomers: "
                  f"{', '.join(stereo_overflow[:5])}", file=sys.stderr)
        if args.tautomerize:
            print("  ! -t replaces the molecule with another tautomer, which "
                  "breaks the atom-order correspondence --stereo enumerate "
                  "relies on; protomers will fall back to the QupKake "
                  "skeleton.", file=sys.stderr)
    if failed:
        with open(outdir / "rejected_input.smi", "w") as fh:
            for smi, name in failed:
                fh.write(f"{smi} {name}\n")

    # ---- 2. QupKake -------------------------------------------------------
    if args.from_sdf:
        sdf_path = Path(args.from_sdf)
        print(f"[2/5] reusing {sdf_path}")
    else:
        print(f"[2/5] running QupKake on {len(inputs)} molecule(s)")
        if args.nproc and args.nproc > 1:
            print(f"  ! -mp {args.nproc}: QupKake gives each of its {args.nproc} workers a "
                  f"{args.nproc}-thread xtb, so {args.nproc * args.nproc} threads compete "
                  f"for {args.nproc} cores. -mp 1 and more processes is faster.",
                  file=sys.stderr)
        work.mkdir(parents=True, exist_ok=True)
        csv_path = work / "qupkake_input.csv"
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["smiles", "name"])
            w.writerows([[s, n] for s, n in inputs])
        sdf_path = run_qupkake(csv_path, work, "qupkake_output.sdf",
                               args.tautomerize, args.nproc, args.qupkake_exe,
                               n_mols=len(inputs),
                               progress_interval=args.progress_interval)

    grouped = parse_qupkake_sdf(sdf_path)
    n_records = sum(len(s) for _, s in grouped.values())
    print(f"      {n_records} site(s) across {len(grouped)} molecule(s)")
    dropped_sites = []
    if not args.keep_amide_sites:
        for name, (mol, sites) in grouped.items():
            bad = implausible_sites(mol)
            keep = [s for s in sites if (s.idx, s.kind) not in bad]
            dropped_sites += [(name, s.idx, s.kind, s.pka, bad[(s.idx, s.kind)])
                              for s in sites if (s.idx, s.kind) in bad]
            grouped[name] = (mol, keep)
        if dropped_sites:
            print(f"      dropped {len(dropped_sites)} implausible site(s): "
                  f"{dict(Counter(r for *_, r in dropped_sites))} (dropped_sites.tsv)")
    with open(outdir / "dropped_sites.tsv", "w") as fh:
        fh.write("name\tidx\tkind\tpka\treason\n")
        fh.writelines(f"{n}\t{i}\t{k}\t{p:.3f}\t{r}\n" for n, i, k, p, r in dropped_sites)

    qk_failed = audit_qupkake(sdf_path.parent.parent, sdf_path.name,
                              [n for _, n in inputs], grouped)
    if qk_failed is None:
        qk_failed = {}
        print("  ! no QupKake working tree next to the SDF: cannot tell a crashed "
              "molecule from one with no sites, so both pass through as "
              "no_sites_predicted", file=sys.stderr)
    elif qk_failed:
        print(f"  ! QupKake failed on {len(qk_failed)} molecule(s) "
              f"({Counter(qk_failed.values())}); excluded, see qupkake_failed.tsv. "
              f"Details in {sdf_path.parent.parent / 'logs'}", file=sys.stderr)
    with open(outdir / "qupkake_failed.tsv", "w") as fh:
        fh.write("name\treason\tsmiles\n")
        for smi, name in inputs:
            if name in qk_failed:
                fh.write(f"{name}\t{qk_failed[name]}\t{smi}\n")

    # ---- 3. enumerate -----------------------------------------------------
    print(f"[3/5] enumerating at pH {', '.join(str(p) for p in args.ph)}")
    by_parent: dict[str, list[Protomer]] = {}
    dropped_states = 0
    build_failures: list[str] = []
    fallback_parents: list[str] = []

    for smi, name in inputs:
        if name in qk_failed:
            continue
        isomers = stereo_sets[name]

        if name not in grouped:
            by_parent[name] = unchanged_protomers(name, isomers, args.ph)
            continue

        ref, sites = grouped[name]
        emitted, n_drop, n_fail, used_isomer = build_protomers(
            name, isomers, ref, sites, args.ph, args.margin, args.min_population,
            args.max_states, args.max_ambiguous, args.stereo)
        dropped_states += n_drop
        build_failures += [name] * n_fail
        if not used_isomer:
            fallback_parents.append(name)
        by_parent[name] = emitted

    # ---- 4. name + write --------------------------------------------------
    print("[4/5] labelling and writing")
    stems, effective_style = plan_names(list(by_parent), args.name_style)
    if effective_style != args.name_style:
        print(f"      names too long for {DB2_NAME_LIMIT} chars; "
              f"switched to compact serials (see name_map.tsv)")

    protomers: list[Protomer] = []
    for name, emitted in by_parent.items():
        n_iso = len(stereo_sets[name])
        for p in emitted:
            p.name = make_name(stems[name], p.protomer_index, p.stereo_index,
                               effective_style, p, n_iso)
            protomers.append(p)

    smi_out = outdir / "protomers.smi"
    with open(smi_out, "w") as fh:
        for p in protomers:
            if args.max_conf:
                fh.write(f"{p.smiles} {p.name} {args.max_conf}\n")
            else:
                fh.write(f"{p.smiles} {p.name}\n")

    with open(outdir / "protomers.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(PROTOMER_COLUMNS)
        w.writerows(protomer_row(p) for p in protomers)

    with open(outdir / "name_map.tsv", "w") as fh:
        fh.write("protomer_name\tparent_name\tnet_charge\tph_values"
                 "\tprotomer_index\tstereo_index\n")
        for p in protomers:
            fh.write(f"{p.name}\t{p.parent}\t{p.charge:+d}\t"
                     f"{';'.join(f'{v:g}' for v in sorted(set(p.ph_values)))}\t"
                     f"{p.protomer_index}\t{p.stereo_index}\n")

    # ---- 5. stats ---------------------------------------------------------
    stats = build_stats(protomers, grouped, inputs, args.ph, n_records)
    iso_counts = [len(v) for v in stereo_sets.values()] or [0]
    stats["stereo"] = {
        "mode": args.stereo,
        "parents": len(inputs),
        "total_stereoisomers": sum(iso_counts),
        "parents_with_undefined_centres": sum(1 for c in iso_counts if c > 1),
        "isomers_per_parent_max": max(iso_counts),
        "qupkake_jobs_run": len(inputs),
        "qupkake_jobs_saved": sum(iso_counts) - len(inputs),
        "capped_parents": len(stereo_overflow),
    }
    stats["input"]["molecules_rejected"] = len(failed)
    stats["input"]["qupkake_failed"] = dict(Counter(qk_failed.values()))
    stats["input"]["sites_dropped_as_implausible"] = dict(
        Counter(r for *_, r in dropped_sites))
    stats["input"]["states_failed_to_build"] = len(build_failures)
    stats["input"]["parents_on_qupkake_skeleton"] = len(set(fallback_parents))
    stats["protomers"]["states_dropped_by_filters"] = dropped_states
    stats["settings"] = {"ph": args.ph, "margin": args.margin,
                         "min_population": args.min_population,
                         "max_states": args.max_states,
                         "tautomerize": args.tautomerize,
                         "name_style": args.name_style,
                         "stereo": args.stereo,
                         "max_conf": args.max_conf}
    with open(outdir / "stats.json", "w") as fh:
        json.dump(stats, fh, indent=2)
    print_stats(stats)

    if build_failures:
        print(f"  ! {len(build_failures)} state(s) could not be built on "
              f"{len(set(build_failures))} molecule(s): "
              f"{', '.join(sorted(set(build_failures))[:5])}", file=sys.stderr)
    if fallback_parents:
        print(f"  ! {len(set(fallback_parents))} molecule(s) fell back to the "
              "QupKake skeleton (atom order differs)", file=sys.stderr)
    if stats["db2_naming"]["names_over_limit"]:
        print(f"  ! {stats['db2_naming']['names_over_limit']} name(s) exceed "
              f"{DB2_NAME_LIMIT} chars; db2 keeps only the LAST {DB2_NAME_LIMIT}.")
    if stats["db2_naming"]["truncation_collisions"]:
        print("  !! names collide after db2 truncation: "
              f"{stats['db2_naming']['truncation_collisions']}")

    print(f"  -> {smi_out}")

    # ---- optional handoff -------------------------------------------------
    if args.db2_cmd:
        cmd = args.db2_cmd.format(smi=str(smi_out))
        print(f"\n[db2] $ {cmd}")
        rc = subprocess.run(cmd, shell=True).returncode
        if rc != 0:
            print(f"  ! db2 command exited {rc}", file=sys.stderr)
            return rc
    return 0


if __name__ == "__main__":
    sys.exit(main())
