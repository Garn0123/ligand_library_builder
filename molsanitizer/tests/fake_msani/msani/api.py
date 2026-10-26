from pandas import DataFrame

# input SMILES -> MolSanitizer-shaped output rows (suffix, smiles)
CANNED = {
    # a protonated tertiary amine whose N-H+ configurations come out as isomers
    "CCN(C)Cc1ccccc1": [(".1", "CC[N@H+](C)Cc1ccccc1"), (".2", "CC[N@@H+](C)Cc1ccccc1")],
    # two states, one of them with two stereoisomers
    "OC(=O)c1ccc(O)cc1": [("_1", "O=C([O-])c1ccc(O)cc1"), ("_2", "O=C(O)c1ccc(O)cc1")],
    "CC(N)C(=O)O": [("_1.1", "C[C@H]([NH3+])C(=O)[O-]"), ("_1.2", "C[C@@H]([NH3+])C(=O)[O-]"),
                    ("_2", "CC(N)C(=O)O")],
}


class Msani:
    def __init__(self, **kw):
        self.kw = kw

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def run(self, df, rejected_file=None, metal_error_file=None):
        rows = []
        for smi, pid in zip(df["smiles"], df["ids"]):
            if "CRASH2D" in smi or smi.startswith("[Xe]"):
                raise RuntimeError("canned 2D failure")
            if "[Na+]" in smi:
                with open(rejected_file, "a") as fh:
                    fh.write(f"{smi} {pid} organic salt removed whole\n")
                continue
            for suffix, out in CANNED.get(smi, [("", smi)]):
                rows.append({"smiles": out, "ids": pid + suffix})
        return DataFrame(rows, columns=["smiles", "ids"])
