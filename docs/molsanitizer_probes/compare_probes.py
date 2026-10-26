"""Compare LLB's MolGpKa-path protomers with MolSanitizer on a small probe set.

Run in an environment with MolSanitizer (pip install molsanitizer==0.7.4):

    python compare_probes.py PROBES.tsv LLB_PARENTS_DIR LLB_TRIAGE_DIR -o OUT.tsv

PROBES.tsv: name, smiles, llb_rule, exp_charges_7.4 (comma-separated net charges
expected to be populated at pH 7.4), evidence.  LLB_PARENTS_DIR / LLB_TRIAGE_DIR are
`llb parents` and `llb triage` outputs for the same SMILES ("SMILES name" lines).

Scoring follows the MolSanitizer paper (section 3.3): over the set of net charges at
pH 7.4, correct = same set, incomplete = strict subset, exceed = strict superset,
incorrect = anything else. Charge sets, not SMILES, so that a different but
equivalent site choice (which N of an amidine) is not counted as a miss.

MolSanitizer runs as `msani --tautomers --protonation` does (neutralize, tautomerize,
ionize) through its Python API, in three settings:
  ms_7      pH 7,   range 0  (the CLI default; the CLI accepts integers only)
  ms_7.4    pH 7.4, range 0  (API only)
  ms_7.4r1  pH 7.4, range 1  (API only; LLB's 6.4-8.4 span)
LLB's charge set at 7.4 = net charges of the protomers whose ph_values include 7.4.
"""
from __future__ import annotations

import argparse
import csv
import glob
from collections import Counter

from rdkit import Chem, RDLogger

from msani.moltransform.ionizer import Ionizer
from msani.moltransform.tautomerizer import Tautomerizer

RDLogger.DisableLog("rdApp.*")


def charge(smi: str) -> int:
    return Chem.GetFormalCharge(Chem.MolFromSmiles(smi))


def score(pred: set, ref: set) -> str:
    if pred == ref:
        return "correct"
    if pred and pred < ref:
        return "incomplete"
    if pred > ref:
        return "exceed"
    return "incorrect"


def fmt(charges) -> str:
    return ",".join(f"{c:+d}" if c else "0" for c in sorted(charges))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("probes")
    ap.add_argument("parents")
    ap.add_argument("triage")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    probes = list(csv.DictReader(open(a.probes), delimiter="\t"))
    pid = {r["input_name"]: r for r in csv.DictReader(open(f"{a.parents}/parents.tsv"), delimiter="\t")}
    route = {r["parent_id"]: r for r in csv.DictReader(open(f"{a.triage}/decisions.tsv"), delimiter="\t")}
    llb = {}
    for f in glob.glob(f"{a.triage}/shards/*/protomers.csv"):
        for r in csv.DictReader(open(f)):
            llb.setdefault(r["parent_name"], []).append(r)

    taut = Tautomerizer()
    settings = {"ms_7": Ionizer(pH=7, pH_range=0), "ms_7.4": Ionizer(pH=7.4, pH_range=0),
                "ms_7.4r1": Ionizer(pH=7.4, pH_range=1)}

    rows, tally = [], {k: Counter() for k in ["llb", *settings]}
    for p in probes:
        ref = {int(x) for x in p["exp_charges_7.4"].split(",")}
        par = pid[p["name"]]
        prot = llb.get(par["parent_id"], [])
        llb_74 = {int(r["net_charge"]) for r in prot if "7.4" in r["ph_values"].split(";")}
        row = {"name": p["name"], "llb_rule": p["llb_rule"], "expected": fmt(ref),
               "llb_route": route[par["parent_id"]]["route"] + (f" ({route[par['parent_id']]['reason']})"
                                                                if route[par["parent_id"]]["reason"] else ""),
               "llb_tautomer_fix": par["tautomer_fix"], "llb_7.4": fmt(llb_74),
               "llb_score": score(llb_74, ref),
               "llb_smiles": " ".join(r["smiles"] for r in prot if "7.4" in r["ph_values"].split(";"))}
        tally["llb"][row["llb_score"]] += 1
        tauts = taut.tautomerize(smiles=p["smiles"]) or [p["smiles"]]
        row["ms_tautomers"] = " ".join(tauts)
        for key, ion in settings.items():
            out = sorted({s for t in tauts for s in (ion.ionize(smiles=t) or [t])})
            row[key] = fmt({charge(s) for s in out})
            row[f"{key}_score"] = score({charge(s) for s in out}, ref)
            tally[key][row[f"{key}_score"]] += 1
            if key == "ms_7.4r1":
                row["ms_7.4r1_smiles"] = " ".join(out)
        rows.append(row)

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    n = len(rows)
    for k, c in tally.items():
        print(f"{k:9s} " + "  ".join(f"{s} {c[s]}/{n}" for s in ["correct", "incomplete", "exceed", "incorrect"]))


if __name__ == "__main__":
    main()
