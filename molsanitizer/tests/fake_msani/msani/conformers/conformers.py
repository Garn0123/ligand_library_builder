"""Fake gen_conf_chunk: db2 headers in the 0.7.4 format (db2/layout.serialize_db2).

Canned behaviour by SMILES: 'Br' -> an msani_error.err line; 'I' -> the process
dies (exit 139) on that molecule; 'Cl' -> wrong AMSOL charge; otherwise OK.
"""
import os
from rdkit import Chem


def gen_conf_chunk(df, args, input_file="0"):
    os.makedirs("db2", exist_ok=True)
    new = not os.path.exists("msani_timing.csv")
    with open("msani_timing.csv", "a") as timing:
        if new:
            timing.write("Name,Initial embedding,AMSOL,Torsional sampling,DB2 conversion,Total\n")
        for smi, name in zip(df["smiles"], df["ids"]):
            if "I" in smi:
                timing.flush()
                os._exit(139)
            if "Br" in smi:
                with open("msani_error.err", "a") as fh:
                    fh.write(f"{smi}\t{name}\tConformer generation: canned failure\n")
                continue
            q = Chem.GetFormalCharge(Chem.MolFromSmiles(smi)) + (1 if "Cl" in smi else 0)
            rec = ("M %16s %9s %3d %3d %6d %6d %6d %6d %6d %6d\n" % (name[-16:], "fake", 10, 10, 50, 12, 3, 6, 5, 0)
                   + "M %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n" % (q, -10.0, 1.0, -9.0, 300.0)
                   + "M %-76s\nM %-76s\nM %+10.4f\n" % (smi[-76:], "fake", 999.999)
                   + "A   1 C1    C.3    1  1   +0.0000     +0.000     +0.000     +0.000     0.000\nE\n")
            with open(f"db2/{name}.db2", "w") as fh:
                fh.write(rec * 2)                       # two rigid-fragment records
            timing.write(f"{name},0.1,0.2,0.3,0.05,0.65\n")
