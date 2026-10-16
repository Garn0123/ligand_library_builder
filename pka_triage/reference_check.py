#!/usr/bin/env python3
"""
reference_check.py -- QupKake and MolGpKa on aromatic N bases with known pKas.

    llb qupkake pka_triage/reference/aromatic_n_bases.smi -o ref_qk --name-style index
    llb pka-reference ref_qk

On the drug-like 10k QupKake put pyrazine, oxazole, 1,2,3-triazole and
tetrazole ring N's at pKa 9-10 (measured pyrazine 0.6; the others are weaker
still), always in molecules that also carry a strong base, while MolGpKa gave
2.7-4.2. This asks whether QupKake is wrong on those rings in general, or only
next to a strong base:

  alone              the ring with a measured pKa (Williams compilation), or a
                     weak base with none (exp NA; expected < 4)
  tethered_guanidine the same ring two carbons from an N-ethylguanidine
  tethered_dimethylamine  ... from a dimethylamino group
  parent             the four 10k parents that showed it

For every molecule: the highest basic pKa each tool gives an aromatic N (the
ring N that would protonate first) and, when tethered, its strongest
non-aromatic base. QupKake sites are taken raw (the site filter has no rule for
aromatic N). Writes OUT/reference_check.tsv.

Reading it: QupKake ring-N values near the measured ones alone but at 8-10 when
tethered = an artefact of the strong base, and the fix is to prefer MolGpKa for
aromatic N in such molecules. Wrong alone too = QupKake should not decide
aromatic N sites at all.
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "run_qupkake"))

from rdkit import Chem, RDLogger  # noqa: E402

from molgpka import MolGpKa  # noqa: E402
from qupkake_protomers import parse_qupkake_sdf  # noqa: E402

RDLogger.DisableLog("rdApp.*")


def split(mol, sites):
    """-> (max basic pKa on an aromatic N, max basic pKa on any other N)."""
    ring, other = None, None
    for s in sites:
        if s.kind != "basic":
            continue
        a = mol.GetAtomWithIdx(s.idx)
        if a.GetAtomicNum() != 7:
            continue
        if a.GetIsAromatic():
            ring = s.pka if ring is None else max(ring, s.pka)
        else:
            other = s.pka if other is None else max(other, s.pka)
    return ring, other


def fmt(x):
    return "" if x is None else f"{x:.2f}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("qupkake_dir", type=Path, help="the -o of llb qupkake on the .smi")
    ap.add_argument("--set", type=Path, default=HERE / "reference" / "aromatic_n_bases.tsv")
    ap.add_argument("-o", "--out", type=Path, default=None,
                    help="default: QUPKAKE_DIR/reference_check.tsv")
    args = ap.parse_args(argv)

    sdf = args.qupkake_dir / "qupkake_work" / "output" / "qupkake_output.sdf"
    if not sdf.exists():
        raise SystemExit(f"{sdf} not found; run llb qupkake on {args.set.with_suffix('.smi')} first")
    qk = parse_qupkake_sdf(sdf)
    failed = set()
    f = args.qupkake_dir / "qupkake_failed.tsv"
    if f.exists():
        with open(f) as fh:
            failed = {r["name"] for r in csv.DictReader(fh, delimiter="\t")}
    model = MolGpKa()

    with open(args.set) as fh:
        ref = list(csv.DictReader(fh, delimiter="\t"))
    rows = []
    for r in ref:
        mol = Chem.MolFromSmiles(r["smiles"])
        mref, msites = model.predict(mol)
        mg_ring, mg_other = split(mref, msites)
        if r["name"] in qk:
            qmol, qsites = qk[r["name"]]
            qk_ring, qk_other = split(qmol, qsites)
            status = "ok"
        else:
            qk_ring = qk_other = None
            status = "qupkake_failed" if r["name"] in failed else "no_qupkake_site"
        rows.append({**r, "molgpka_ring_n": fmt(mg_ring), "qupkake_ring_n": fmt(qk_ring),
                     "molgpka_strong_base": fmt(mg_other), "qupkake_strong_base": fmt(qk_other),
                     "qupkake_status": status})

    out = args.out or args.qupkake_dir / "reference_check.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"{'name':<30} {'group':<23} {'exp':>5} {'MolGpKa':>8} {'QupKake':>8}   "
          f"strong base MG/QK")
    for r in rows:
        strong = (f"{r['molgpka_strong_base'] or '-'}/{r['qupkake_strong_base'] or '-'}"
                  if r["group"] != "alone" else "")
        note = "" if r["qupkake_status"] == "ok" else f"  ({r['qupkake_status']})"
        print(f"{r['name']:<30} {r['group']:<23} {r['exp_ring_n_pka']:>5} "
              f"{r['molgpka_ring_n'] or '-':>8} {r['qupkake_ring_n'] or '-':>8}   {strong}{note}")

    # measured rings: error of each tool
    print()
    for tool, col in (("MolGpKa", "molgpka_ring_n"), ("QupKake", "qupkake_ring_n")):
        d = [float(r[col]) - float(r["exp_ring_n_pka"]) for r in rows
             if r["exp_ring_n_pka"] != "NA" and r[col]]
        if d:
            print(f"{tool} vs measured ({len(d)} rings alone): MAE {statistics.mean(map(abs, d)):.2f}, "
                  f"mean {statistics.mean(d):+.2f}, worst {max(d, key=abs):+.2f}")
    # alone vs tethered, per ring, QupKake
    print("\nQupKake ring N, alone -> tethered (guanidine / dimethylamine):")
    by = {}
    for r in rows:
        if r["ring"] and r["qupkake_ring_n"]:
            by.setdefault(r["ring"], {}).setdefault(r["group"], []).append(float(r["qupkake_ring_n"]))
    for ring, g in sorted(by.items()):
        alone = g.get("alone")
        a = f"{max(alone):.1f}" if alone else "-"
        t = " / ".join(f"{max(g[k]):.1f}" if k in g else "-"
                       for k in ("tethered_guanidine", "tethered_dimethylamine"))
        p = f"   parent {max(g['parent']):.1f}" if "parent" in g else ""
        print(f"  {ring:<12} {a:>5} -> {t}{p}")
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
