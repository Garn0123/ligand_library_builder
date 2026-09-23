"""Build a fake QupKake output SDF exactly the way predict.py does, so the
pipeline can be tested without running xtb/torch."""
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools
import pandas as pd

# name -> (smiles, [(atom-match-smarts, pka_type, pka)])
CASES = {
    "benzoic_acid": ("OC(=O)c1ccccc1", [("[OX2H1][CX3]=O", "acidic", 4.20)]),
    "propranolol": ("CC(C)NCC(O)COc1cccc2ccccc12",
                    [("[NX3;H1]([CH](C)C)", "basic", 9.50)]),
    "histidine": ("N[C@@H](Cc1c[nH]cn1)C(=O)O",
                  [("[NX3;H2][CH]", "basic", 9.20),
                   ("[OX2H1][CX3]=O", "acidic", 2.20),
                   ("[nX2;H0]", "basic", 6.00)]),
    "borderline": ("c1ccc(cc1)C(=O)Nc1ccncc1",   # pyridine right at 7.4
                   [("n1ccccc1", "basic", 7.30)]),
}

rows = []
for name, (smi, sites) in CASES.items():
    # QupKake receives the CANONICAL standardized SMILES from the pipeline,
    # so site indices must be assigned against that, not the hand-written order
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(smi)))
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    AllChem.MMFFOptimizeMolecule(mol)
    for smarts, kind, pka in sites:
        patt = Chem.MolFromSmarts(smarts)
        match = mol.GetSubstructMatch(patt)
        assert match, f"{name}: {smarts} did not match"
        idx = match[0]
        rows.append({"ROMol": mol, "name": name, "idx": idx,
                     "pka_type": kind, "pka": pka})

df = pd.DataFrame(rows)
PandasTools.WriteSDF(df, "fake_qupkake_output.sdf", molColName="ROMol",
                     idName="name", properties=["idx", "pka_type", "pka"])
print(f"wrote {len(rows)} site records for {len(CASES)} molecules")
for name, (smi, _) in CASES.items():
    print(f"  {smi}  {name}")
