"""MolSanitizer vs LLB on pka_triage/reference/aromatic_n_bases.{smi,tsv}.

    python compare_aromatic_n.py REFERENCE_CHECK.tsv LLB_PARENTS_DIR LLB_TRIAGE_DIR -o OUT.tsv

REFERENCE_CHECK.tsv is `llb pka-reference` output (the drug-like 10k run's copy is
validation/2026-10_druglike10k/druglike_10k/ref_qk/reference_check.tsv on branch
data/validation-2026-10): measured ring-N pKa, MolGpKa's and QupKake's.

For the 'alone' rows with a measured pKa, the expected net-charge set at pH 7.4 holds
every charge with >= 10% population (one basic ring N). Columns:
  llb_molgpka   LLB's MolGpKa-path protomers at 7.4 (this run)
  llb_qupkake   what LLB's state model (margin 1.0) makes of QupKake's ring-N value
                in reference_check.tsv: +1 above 8.4, 0 below 6.4, both between
  llb_route     the triage route; for routed 'alone' rows production uses llb_qupkake
  ms_7.4 / ms_7.4r1   MolSanitizer (tautomerize + ionize) at pH 7.4, range 0 / 1
Tethered rows (ring + guanidine or dimethylamine) have no measurement; they show
whether each tool protonates the ring N beside a charged base.
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


def charge(smi):
    return Chem.GetFormalCharge(Chem.MolFromSmiles(smi))


def fmt(cs):
    return ",".join(f"{c:+d}" if c else "0" for c in sorted(cs))


def expected(pka, ph=7.4):
    f = 1 / (1 + 10 ** (ph - pka))
    return {c for c, keep in [(1, f >= 0.1), (0, f <= 0.9)] if keep}


def margin_set(pka, ph=7.4, margin=1.0):
    return {1} if pka > ph + margin else {0} if pka < ph - margin else {0, 1}


def score(pred, ref):
    if pred == ref:
        return "correct"
    if pred and pred < ref:
        return "incomplete"
    if pred > ref:
        return "exceed"
    return "incorrect"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference")
    ap.add_argument("parents")
    ap.add_argument("triage")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    ref = list(csv.DictReader(open(a.reference), delimiter="\t"))
    pid = {r["input_name"]: r["parent_id"] for r in csv.DictReader(open(f"{a.parents}/parents.tsv"), delimiter="\t")}
    route = {r["parent_id"]: r["route"] for r in csv.DictReader(open(f"{a.triage}/decisions.tsv"), delimiter="\t")}
    llb = {}
    for f in glob.glob(f"{a.triage}/shards/*/protomers.csv"):
        for r in csv.DictReader(open(f)):
            if "7.4" in r["ph_values"].split(";"):
                llb.setdefault(r["parent_name"], set()).add(int(r["net_charge"]))

    taut = Tautomerizer()
    ions = {"ms_7.4": Ionizer(pH=7.4, pH_range=0), "ms_7.4r1": Ionizer(pH=7.4, pH_range=1)}
    rows = []
    tally = {k: Counter() for k in ["llb_molgpka", "llb_production", *ions]}
    for r in ref:
        p = pid[r["name"]]
        row = {"name": r["name"], "group": r["group"], "exp_pka": r["exp_ring_n_pka"],
               "molgpka_ring_n": r["molgpka_ring_n"], "qupkake_ring_n": r["qupkake_ring_n"],
               "llb_route": route[p], "llb_molgpka": fmt(llb.get(p, set()))}
        qk = margin_set(float(r["qupkake_ring_n"])) if r["group"] == "alone" and r["qupkake_ring_n"] not in ("", "NA") else None
        row["llb_qupkake"] = fmt(qk) if qk is not None else ""
        tauts = taut.tautomerize(smiles=r["smiles"]) or [r["smiles"]]
        ms = {}
        for k, ion in ions.items():
            ms[k] = {charge(s) for t in tauts for s in (ion.ionize(smiles=t) or [t])}
            row[k] = fmt(ms[k])
        if r["group"] == "alone" and r["exp_ring_n_pka"] not in ("", "NA"):
            exp = expected(float(r["exp_ring_n_pka"]))
            prod = qk if route[p] == "qupkake" and qk is not None else llb.get(p, set())
            row["expected"] = fmt(exp)
            for k, pred in [("llb_molgpka", llb.get(p, set())), ("llb_production", prod), *ms.items()]:
                row[f"{k}_score"] = score(pred, exp)
                tally[k][row[f"{k}_score"]] += 1
        rows.append(row)

    cols = list(dict.fromkeys(k for row in rows for k in row))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    for k, c in tally.items():
        n = sum(c.values())
        print(f"{k:15s} " + "  ".join(f"{s} {c[s]}/{n}" for s in ["correct", "incomplete", "exceed", "incorrect"]))


if __name__ == "__main__":
    main()
