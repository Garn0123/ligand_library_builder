"""Tautomer probes: LLB's parents step against MolSanitizer's tautomerizer.

Two environments, so two passes, then a score:

    python compare_tautomers.py llb   tautomer_probes.tsv > llb.tsv     # llb_pka env, repo root on disk
    python compare_tautomers.py msani tautomer_probes.tsv > msani.tsv   # env with molsanitizer
    python compare_tautomers.py score tautomer_probes.tsv llb.tsv msani.tsv

llb   = qupkake_protomers.standardize + amino_tautomer, exactly what `llb parents` does.
msani = Tautomerizer() defaults (neutralize, RDKit layer + rule layer); it may return
        several tautomers. Scored as in the MolSanitizer paper: correct = only the
        expected form, exceed = expected form plus others, incorrect = expected absent.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def canon(smi: str) -> str:
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m is not None else f"UNPARSABLE:{smi}"


def main() -> None:
    mode, probes = sys.argv[1], list(csv.DictReader(open(sys.argv[2]), delimiter="\t"))
    if mode == "llb":
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "run_qupkake"))
        from qupkake_protomers import amino_tautomer, standardize
        for p in probes:
            m = standardize(p["input_smiles"])
            t = amino_tautomer(m) if m is not None else None
            print(p["name"], Chem.MolToSmiles(t if t is not None else m), sep="\t")
    elif mode == "msani":
        from msani.moltransform.tautomerizer import Tautomerizer
        taut = Tautomerizer()
        for p in probes:
            print(p["name"], " ".join(taut.tautomerize(smiles=p["input_smiles"]) or [p["input_smiles"]]), sep="\t")
    elif mode == "score":
        out = [dict(line.rstrip("\n").split("\t") for line in open(f)) for f in sys.argv[3:5]]
        tally = {"llb": {}, "msani": {}}
        w = csv.writer(sys.stdout, delimiter="\t")
        w.writerow(["name", "expected", "llb", "llb_score", "msani", "msani_score"])
        for p in probes:
            exp = canon(p["expected_smiles"])
            row = [p["name"], exp]
            for tool, got in zip(tally, out):
                forms = {canon(s) for s in got[p["name"]].split()}
                s = ("correct" if forms == {exp} else "exceed" if exp in forms else "incorrect")
                tally[tool][s] = tally[tool].get(s, 0) + 1
                row += [" ".join(sorted(forms)), s]
            w.writerow(row)
        for tool, c in tally.items():
            print(f"# {tool}: " + "  ".join(f"{k} {c.get(k, 0)}/{len(probes)}"
                                             for k in ["correct", "exceed", "incorrect"]), file=sys.stderr)


if __name__ == "__main__":
    main()
