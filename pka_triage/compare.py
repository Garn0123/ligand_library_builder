#!/usr/bin/env python3
"""
compare.py -- MolGpKa (triage) against QupKake on the same parents.

    llb compare --parents parents --triage triage --qupkake protomers -o compare

--qupkake is an all-QupKake run over the SAME parents directory
(llb submit-qupkake parents protomers), so every parent has both predictions.
For a big set, run the triage with a huge --route-window first (everything
routed, so the QupKake run is over exactly the triage's parents) -- or just
point this at any QupKake run of the same parents.

What it answers:
  1. Do the two find the same ionizable sites, and how far apart are the pKas?
     sites.tsv + scatter.png, by site class (carboxylic acid, aliphatic amine,
     pyridine-like N, ...).
  2. Does the triage's decision hold? For each parent and target pH, the
     DOMINANT microstate (every site on whichever side of its pKa the pH puts
     it) from MolGpKa's pKas vs QupKake's. parents.tsv.
  3. What --route-window to use: for each window, the fraction of parents the
     triage would settle, and how many of THOSE settle on a different dominant
     state than QupKake gives (the misses routing exists to prevent).
     window.png + summary.json.

Both tools' site indices are mapped onto the parent as written in the shard
(substructure match), so nothing assumes the two kept the same atom order.
MolGpKa sites are taken as the triage used them (kept=1: amide-N bases
dropped). Agreement here is model against model -- neither is measured pKa.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "run_qupkake"))

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem.MolStandardize import rdMolStandardize  # noqa: E402

from qupkake_protomers import Site, parse_qupkake_sdf  # noqa: E402

RDLogger.DisableLog("rdApp.*")

WINDOWS = [0.25 * k for k in range(0, 17)]

# First match wins; the site atom is the pattern's first atom.
SITE_CLASSES = [
    ("carboxylic acid", "[OX2H1][CX3]=O"),
    ("phenol OH", "[OX2H1]c"),
    ("enol / oxime OH", "[OX2H1][$(C=*),$(N=*)]"),
    ("alcohol OH", "[OX2H1][CX4]"),
    ("sulfonamide NH", "[NX3;!H0]S(=O)=O"),
    ("amide / imide NH", "[NX3;!H0]C=[O,S]"),
    ("amidine / guanidine", "[NX2]=C[NX3]"),
    ("imidazole-like N", "[nX2;r5]"),
    ("azole NH", "[nX3;H1]"),
    ("N-substituted aromatic N", "[nX3;H0]"),
    ("pyridine-like N", "[nX2;r6]"),
    ("aniline N", "[NX3;!$(NC=[O,S,N]);!$(NS=O)]c"),
    ("aliphatic amine", "[NX3;!$(N[a]);!$(NC=[O,S,N]);!$(NS=O)]"),
    ("other", "[*]"),
]
_PATTERNS = [(name, Chem.MolFromSmarts(s)) for name, s in SITE_CLASSES]
_UNCHARGER = rdMolStandardize.Uncharger()


def site_class(mol: Chem.Mol, idx: int) -> str:
    for name, patt in _PATTERNS:
        if any(m[0] == idx for m in mol.GetSubstructMatches(patt)):
            return name
    return "other"


def context(mol: Chem.Mol, idx: int, radius: int = 2) -> str:
    """SMILES of the atoms within `radius` bonds of the site, site marked :1."""
    m = Chem.Mol(mol)
    m.GetAtomWithIdx(idx).SetAtomMapNum(1)
    for r in range(radius, 0, -1):
        bonds = list(Chem.FindAtomEnvironmentOfRadiusN(m, r, idx))
        if bonds:
            atoms = {idx} | {a for b in bonds for a in (m.GetBondWithIdx(b).GetBeginAtomIdx(),
                                                         m.GetBondWithIdx(b).GetEndAtomIdx())}
            return Chem.MolFragmentToSmiles(m, atomsToUse=sorted(atoms), bondsToUse=bonds,
                                            canonical=True)
    return Chem.MolToSmiles(m)


def mapping(src: Chem.Mol, parent: Chem.Mol) -> dict[int, int] | None:
    """src heavy-atom index -> parent atom index, by substructure match."""
    src = Chem.RemoveHs(src, sanitize=False)
    m = src.GetSubstructMatch(parent)          # m[parent_idx] = src_idx
    if len(m) != parent.GetNumAtoms():
        return None
    return {s: p for p, s in enumerate(m)}


def dominant(sites: list[Site], ph: float) -> frozenset:
    return frozenset((s.idx, s.kind) for s in sites if s.changed_fraction(ph) >= 0.5)


def charge(state: frozenset) -> int:
    return sum(1 if k == "basic" else -1 for _, k in state)


def read_shards(parents: Path) -> dict[str, str]:
    smiles = {}
    with open(parents / "shards.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for line in (parents / r["path"]).read_text().splitlines():
                if line.strip():
                    smi, pid = line.split()[:2]
                    smiles[pid] = smi
    return smiles


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parents", type=Path, required=True)
    ap.add_argument("--triage", type=Path, required=True)
    ap.add_argument("--qupkake", type=Path, required=True,
                    help="an all-QupKake run (llb submit-qupkake) over the same parents")
    ap.add_argument("-o", "--outdir", type=Path, default=Path("compare"))
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args(argv)

    settings = json.loads((args.triage / "settings.json").read_text())
    phs = settings["ph"]
    smiles = read_shards(args.parents)
    parent_mol = {p: Chem.MolFromSmiles(s) for p, s in smiles.items()}

    # ---- MolGpKa sites, mapped onto the parent ----------------------------------
    mg_sites: dict[str, list[Site]] = defaultdict(list)
    decisions = {}
    for d in sorted((args.triage / "shards").glob("shard_*")):
        with open(d / "decisions.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                decisions[r["parent_id"]] = r
        with open(d / "sites.tsv") as fh:
            raw = defaultdict(list)
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("kept", "1") == "1":
                    raw[r["parent_id"]].append(r)
        for pid, rows in raw.items():
            # the triage predicted on MolFromSmiles(MolToSmiles(uncharge(parent)))
            ref = Chem.MolFromSmiles(Chem.MolToSmiles(_UNCHARGER.uncharge(parent_mol[pid])))
            mp = mapping(ref, parent_mol[pid])
            if mp is None:
                print(f"  ! {pid}: MolGpKa atoms do not map onto the parent; skipped",
                      file=sys.stderr)
                continue
            mg_sites[pid] = [Site(mp[int(r["idx"])], r["kind"], float(r["pka"])) for r in rows]

    # ---- QupKake sites ---------------------------------------------------------------
    qk_sites: dict[str, list[Site]] = {}
    qk_failed: set[str] = set()
    qk_seen: set[str] = set()
    unmapped = 0
    for d in sorted(args.qupkake.glob("shard_*")):
        if not (d / "DONE").exists():
            continue
        with open(d / "protomers.csv") as fh:
            qk_seen |= {r["parent_name"] for r in csv.DictReader(fh)}
        f = d / "qupkake_failed.tsv"
        if f.exists():
            with open(f) as fh:
                qk_failed |= {r["name"] for r in csv.DictReader(fh, delimiter="\t")}
        sdf = d / "qupkake_work" / "output" / "qupkake_output.sdf"
        if not sdf.exists():
            continue
        for pid, (mol, sites) in parse_qupkake_sdf(sdf).items():
            if pid not in parent_mol:
                continue
            mp = mapping(mol, parent_mol[pid])
            if mp is None:
                unmapped += 1
                qk_failed.add(pid)
                continue
            qk_sites[pid] = [Site(mp[s.idx], s.kind, s.pka) for s in sites if s.idx in mp]
    # QupKake finished on a parent but listed no site: a genuine "no sites"
    for pid in qk_seen - set(qk_sites) - qk_failed:
        qk_sites[pid] = []

    both = sorted(p for p in parent_mol if p in qk_sites and p in decisions)
    args.outdir.mkdir(parents=True, exist_ok=True)

    # ---- 1. sites ----------------------------------------------------------------------
    site_rows, deltas = [], defaultdict(list)
    for pid in both:
        mol = parent_mol[pid]
        a = {(s.idx, s.kind): s.pka for s in mg_sites.get(pid, [])}
        b = {(s.idx, s.kind): s.pka for s in qk_sites[pid]}
        for key in sorted(set(a) | set(b)):
            cls = site_class(mol, key[0])
            pa, pb = a.get(key), b.get(key)
            site_rows.append({"parent_id": pid, "atom": key[0], "kind": key[1],
                              "class": cls,
                              "pka_molgpka": "" if pa is None else f"{pa:.2f}",
                              "pka_qupkake": "" if pb is None else f"{pb:.2f}",
                              "delta": "" if None in (pa, pb) else f"{pa - pb:+.2f}",
                              "found_by": "both" if None not in (pa, pb) else
                                          ("molgpka" if pb is None else "qupkake"),
                              "context": context(mol, key[0])})
            if None not in (pa, pb):
                deltas[cls].append(pa - pb)
    with open(args.outdir / "sites.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(site_rows[0]) if site_rows
                           else ["parent_id"])
        w.writeheader()
        w.writerows(site_rows)

    # ---- 2. dominant microstate per pH ----------------------------------------------
    parent_rows, gap = [], {}
    for pid in both:
        ms, qs = mg_sites.get(pid, []), qk_sites[pid]
        gap[pid] = min((abs(s.pka - ph) for s in ms for ph in phs), default=math.inf)
        row = {"parent_id": pid, "route": decisions[pid]["route"],
               "reason": decisions[pid]["reason"],
               "min_gap": "" if gap[pid] == math.inf else f"{gap[pid]:.2f}",
               "coupled": int(bool(decisions[pid]["coupled_pairs"]))}
        agree_all = True
        for ph in phs:
            dm, dq = dominant(ms, ph), dominant(qs, ph)
            same = dm == dq
            agree_all &= same
            row[f"agree_{ph:g}"] = int(same)
            row[f"charge_{ph:g}"] = f"{charge(dm):+d}/{charge(dq):+d}"
        row["agree_all"] = int(agree_all)
        mid = sorted(phs)[len(phs) // 2]
        row[f"qupkake_amide_anion_{mid:g}"] = int(any(
            k == "acidic" and site_class(parent_mol[pid], a) == "amide / imide NH"
            for a, k in dominant(qs, mid)))
        parent_rows.append(row)
    with open(args.outdir / "parents.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(parent_rows[0]) if parent_rows
                           else ["parent_id"])
        w.writeheader()
        w.writerows(parent_rows)

    # ---- 3. what each window would have done ---------------------------------------
    curve = []
    by_pid = {r["parent_id"]: r for r in parent_rows}
    n = len(parent_rows)
    for w_ in WINDOWS:
        settled = [p for p in both
                   if gap[p] > w_ and not (by_pid[p]["coupled"] and settings["coupled"] == "route")]
        missed = [p for p in settled if not by_pid[p]["agree_all"]]
        curve.append({"window": w_, "routed_fraction": round(1 - len(settled) / n, 3) if n else 0,
                      "settled": len(settled), "settled_disagree": len(missed),
                      "miss_rate_among_settled":
                          round(len(missed) / len(settled), 3) if settled else None})

    def stats(v):
        return {"n": len(v), "mae": round(sum(map(abs, v)) / len(v), 2),
                "mean_signed": round(sum(v) / len(v), 2)} if v else {"n": 0}
    found = Counter(r["found_by"] for r in site_rows)
    summary = {
        "parents_compared": n,
        "parents_qupkake_failed_or_unmapped": len(set(parent_mol) & qk_failed),
        "sites": dict(found),
        "pka_delta_molgpka_minus_qupkake": {"all": stats([d for v in deltas.values() for d in v]),
                                            **{c: stats(v) for c, v in sorted(deltas.items())}},
        "dominant_state_agreement": {f"{ph:g}": round(sum(r[f"agree_{ph:g}"] for r in parent_rows) / n, 3)
                                     for ph in phs} if n else {},
        "parents_where_qupkake_deprotonates_an_amide_nh": sum(
            v for r in parent_rows for k, v in r.items() if k.startswith("qupkake_amide_anion_")),
        "agree_at_every_ph": round(sum(r["agree_all"] for r in parent_rows) / n, 3) if n else None,
        "triage_as_run": {
            "route_window": settings["route_window"],
            "settled": sum(1 for r in parent_rows if r["route"] == "molgpka"),
            "settled_but_disagree": [r["parent_id"] for r in parent_rows
                                     if r["route"] == "molgpka" and not r["agree_all"]]},
        "by_window": curve,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))

    # ---- report ------------------------------------------------------------------------
    print(f"{n} parents compared ({summary['parents_qupkake_failed_or_unmapped']} QupKake "
          f"failures left out)")
    print(f"sites: {found.get('both', 0)} found by both, {found.get('molgpka', 0)} MolGpKa only, "
          f"{found.get('qupkake', 0)} QupKake only")
    print("pKa MolGpKa - QupKake, by site class:   n    MAE   mean")
    for c, s in summary["pka_delta_molgpka_minus_qupkake"].items():
        if s["n"]:
            print(f"  {c:<24} {s['n']:>5} {s['mae']:>6.2f} {s['mean_signed']:>+6.2f}")
    print("dominant microstate agrees with QupKake: " + ", ".join(
        f"pH {k} {v:.0%}" for k, v in summary["dominant_state_agreement"].items())
        + f"; all pH {summary['agree_at_every_ph']:.0%}")
    print(f"parents whose QupKake state at mid pH has an amide/imide NH deprotonated: "
          f"{summary['parents_where_qupkake_deprotonates_an_amide_nh']}  (see sites.tsv 'context')")
    t = summary["triage_as_run"]
    print(f"triage as run (window {t['route_window']}): settled {t['settled']}, of which "
          f"{len(t['settled_but_disagree'])} disagree with QupKake "
          f"{t['settled_but_disagree'][:5]}")
    print("window  routed  settled  settled-but-disagree")
    for c in curve:
        if c["window"] in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
            print(f"  {c['window']:<5} {c['routed_fraction']:>6.0%} {c['settled']:>8} "
                  f"{c['settled_disagree']:>8}")

    if not args.no_plots:
        plots(site_rows, curve, phs, settings["route_window"], args.outdir)
    print(f"-> {args.outdir}/summary.json, sites.tsv, parents.tsv"
          + ("" if args.no_plots else ", scatter.png, window.png"))
    return 0


def plots(site_rows, curve, phs, route_window, outdir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pairs = [r for r in site_rows if r["found_by"] == "both"]
    fig, ax = plt.subplots(figsize=(6.4, 6))
    lo, hi = -2, 16
    ax.axvspan(min(phs), max(phs), color="0.9", zorder=0)
    ax.axhspan(min(phs), max(phs), color="0.9", zorder=0)
    ax.plot([lo, hi], [lo, hi], color="0.3", lw=1)
    for off in (-1, 1):
        ax.plot([lo, hi], [lo + off, hi + off], color="0.6", lw=0.8, ls="--")
    classes = sorted({r["class"] for r in pairs})
    cmap = plt.get_cmap("tab10")
    for i, c in enumerate(classes):
        pts = [r for r in pairs if r["class"] == c]
        ax.scatter([float(r["pka_molgpka"]) for r in pts], [float(r["pka_qupkake"]) for r in pts],
                   s=22, color=cmap(i % 10), label=f"{c} ({len(pts)})", alpha=0.85)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("MolGpKa pKa")
    ax.set_ylabel("QupKake pKa")
    ax.set_title(f"Per-site pKa, {len(pairs)} sites found by both\n"
                 f"grey band = target pH {min(phs):g}-{max(phs):g}; dashed = +/-1")
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(outdir / "scatter.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4))
    ws = [c["window"] for c in curve]
    ax.plot(ws, [c["routed_fraction"] for c in curve], marker="o", label="fraction routed to QupKake")
    miss = [c["miss_rate_among_settled"] for c in curve]
    ax.plot([w for w, m in zip(ws, miss) if m is not None], [m for m in miss if m is not None],
            marker="s", label="settled but dominant state differs from QupKake")
    ax.axvline(route_window, color="0.5", ls=":", lw=1)
    ax.set_xlabel("--route-window (pKa units)")
    ax.set_ylabel("fraction")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8)
    ax.set_title("What each routing window would do")
    fig.tight_layout()
    fig.savefig(outdir / "window.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
