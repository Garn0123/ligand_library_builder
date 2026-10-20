"""Test B: MolSanitizer's states and tautomers against a finished LLB library.

    python compare_library.py RUN/parents/parents.tsv RUN/library/library.tsv -o OUTDIR [-j 8]

Needs only an environment with MolSanitizer (pip install molsanitizer==0.7.4); LLB's
side is read from its files. For every parent LLB named:

  llb_7.4        net charges of LLB's names whose ph_values include 7.4 (invertomers and
                 stereoisomers collapse: this is a set)
  ms_7.4, ms_r1  MolSanitizer (neutralize, tautomerize, ionize) on the same neutral parent,
                 pH 7.4 range 0 and range 1 (LLB's 6.4-8.4)
  ms_tautomer    1 if MolSanitizer's tautomer layer changed the heavy-atom/H skeleton
  classes        which ionizable classes the parent carries (CLASSES below)

Writes OUTDIR/parents.tsv (one row per parent), OUTDIR/summary.json (agreement overall
and per class) and OUTDIR/review.tsv (a seeded random 100 of the disagreements, for a
person to judge against measured pKa where one exists).
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import time
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

CLASSES = {  # first-pass buckets for the summary; not exhaustive, overlaps allowed
    "aliphatic_amine": "[NX3;!$(N-[a,#6]=[O,S,N]);!$(N-S(=O)=O);!$(N-a);!$(N-[#7,#8])]",
    "piperazine": "N1CCNCC1",
    "amidine_guanidine": "[NX3;!$(N-C=O)]-[CX3](=[NX2;!$(N-C=O)])",
    "biguanide": "N-C(=N)-N-C(=N)-N",
    "aniline_like": "[NX3;H2,H1;!$(N-C=O)]-a",
    "pyridine_like": "[nX2;r6]",
    "imidazole_like": "[nX2;r5]:[c,n]:[nX3;r5]",
    "tetrazole": "c1nn[nH]n1",
    "carboxylic_acid": "C(=O)[OH]",
    "acyl_sulfonamide": "[#6]C(=O)[NH]S(=O)(=O)",
    "sulfonamide_NH": "[NX3;H1,H2;!$(N-C=O)]S(=O)(=O)",
    "phenol": "[OH]c",
    "imide": "C(=O)[NH]C(=O)",
}
_PATS = {k: Chem.MolFromSmarts(v) for k, v in CLASSES.items()}
_TAUT = _IONS = None


def _init():
    global _TAUT, _IONS
    from msani.moltransform.ionizer import Ionizer
    from msani.moltransform.tautomerizer import Tautomerizer
    _TAUT = Tautomerizer()
    _IONS = {"ms_7.4": Ionizer(pH=7.4, pH_range=0), "ms_r1": Ionizer(pH=7.4, pH_range=1)}


def skeleton(smi: str) -> str:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return ""
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)


def run_ms(item):
    pid, smi = item
    try:
        tauts = _TAUT.tautomerize(smiles=smi) or [smi]
        out = {"ms_tautomer": int({skeleton(t) for t in tauts} != {skeleton(smi)}),
               "ms_tautomers": " ".join(tauts)}
        for k, ion in _IONS.items():
            prods = sorted({s for t in tauts for s in (ion.ionize(smiles=t) or [t])})
            out[k] = sorted({Chem.GetFormalCharge(Chem.MolFromSmiles(s)) for s in prods})
            out[f"{k}_smiles"] = " ".join(prods)
        return pid, out
    except Exception as e:  # recorded, not fatal
        return pid, {"ms_error": repr(e)[:200]}


def fmt(cs):
    return ",".join(f"{c:+d}" if c else "0" for c in sorted(cs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("parents")
    ap.add_argument("library")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("-j", "--jobs", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20261020)
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    parents = {r["parent_id"]: r for r in csv.DictReader(open(a.parents), delimiter="\t")}
    llb, src = defaultdict(set), {}
    for r in csv.DictReader(open(a.library), delimiter="\t"):
        if "7.4" in r["ph_values"].replace(",", ";").split(";"):
            llb[r["parent_id"]].add(int(r["net_charge"]))
        src[r["parent_id"]] = r.get("site_source", "")
    todo = [(p, parents[p]["smiles"]) for p in src if p in parents]

    t0 = time.time()
    with Pool(a.jobs, initializer=_init) as pool:
        ms = dict(pool.imap_unordered(run_ms, todo, chunksize=20))
    seconds = time.time() - t0

    rows, by_class = [], defaultdict(Counter)
    for pid, smi in todo:
        m = ms[pid]
        mol = Chem.MolFromSmiles(smi)
        classes = [k for k, p in _PATS.items() if mol is not None and mol.HasSubstructMatch(p)]
        row = {"parent_id": pid, "smiles": smi, "site_source": src[pid], "classes": " ".join(classes),
               "llb_7.4": fmt(llb[pid])}
        if "ms_error" in m:
            row["ms_error"] = m["ms_error"]
        else:
            row.update({"ms_7.4": fmt(m["ms_7.4"]), "ms_r1": fmt(m["ms_r1"]),
                        "ms_tautomer": m["ms_tautomer"], "ms_tautomers": m["ms_tautomers"],
                        "ms_r1_smiles": m["ms_r1_smiles"]})
            row["agree_dominant"] = int(set(m["ms_7.4"]) <= llb[pid])   # MS's single state is among LLB's
            row["agree_set_r1"] = int(set(m["ms_r1"]) == llb[pid])
        rows.append(row)
        for c in classes or ["no_class"]:
            by_class[c]["n"] += 1
            by_class[c]["agree_dominant"] += row.get("agree_dominant", 0)
            by_class[c]["ms_tautomer"] += row.get("ms_tautomer", 0)

    cols = list(dict.fromkeys(k for r in rows for k in r))
    with open(out / "parents.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    dis = [r for r in rows if r.get("agree_dominant") == 0 or r.get("ms_tautomer") == 1]
    random.Random(a.seed).shuffle(dis)
    with open(out / "review.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols + ["verdict", "evidence"], delimiter="\t", restval="")
        w.writeheader()
        w.writerows(dis[:100])
    ok = [r for r in rows if "ms_error" not in r]
    summary = {
        "parents": len(rows), "ms_errors": len(rows) - len(ok), "ms_seconds": round(seconds, 1),
        "ms_ms_per_parent_wall": round(1000 * seconds / max(1, len(rows)), 1), "jobs": a.jobs,
        "agree_dominant": round(sum(r["agree_dominant"] for r in ok) / max(1, len(ok)), 4),
        "agree_set_r1": round(sum(r["agree_set_r1"] for r in ok) / max(1, len(ok)), 4),
        "ms_tautomer_changed": round(sum(r["ms_tautomer"] for r in ok) / max(1, len(ok)), 4),
        "llb_states_per_parent_7.4": round(sum(len(llb[r["parent_id"]]) for r in ok) / max(1, len(ok)), 3),
        "by_class": {c: {"n": v["n"], "agree_dominant": round(v["agree_dominant"] / v["n"], 3),
                         "ms_tautomer_changed": round(v["ms_tautomer"] / v["n"], 3)}
                     for c, v in sorted(by_class.items())},
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "by_class"}, indent=2))


if __name__ == "__main__":
    main()
