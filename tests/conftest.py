"""Shared fixtures for the library-build tests (run_qupkake/, pka_triage/).

    pytest                       # from the repo root; mol_compiler/ has its own suite

Most tests need only RDKit. Ones marked `molgpka` also need torch and a MolGpKa
checkout (MOLGPKA_DIR, as in config/hpc.env) and are skipped without it. The real
QupKake/xtb and db2_converter never run here: `fake_qupkake` writes a QupKake-shaped
SDF from the input CSV, with a basic site at a chosen pKa on every nitrogen, which
is the case that exposed the aromatic-N artefact.
"""
from __future__ import annotations

import os
import stat
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
for sub in ("run_qupkake", "pka_triage", "pipeline"):
    sys.path.insert(0, str(ROOT / sub))


def molgpka_available() -> bool:
    d = os.environ.get("MOLGPKA_DIR", "")
    if not d or not (Path(d) / "models" / "weight_acid.pth").exists():
        return False
    try:
        import torch  # noqa: F401
    except ImportError:
        return False
    return True


def pytest_collection_modifyitems(config, items):
    if molgpka_available():
        return
    skip = pytest.mark.skip(reason="needs torch and MOLGPKA_DIR (a MolGpKa checkout)")
    for item in items:
        if "molgpka" in item.keywords:
            item.add_marker(skip)


def pytest_configure(config):
    config.addinivalue_line("markers", "molgpka: needs torch and a MolGpKa checkout")


FAKE_QUPKAKE = r'''#!{python}
# QupKake stand-in: `qupkake file CSV -o OUT -r ROOT ...`. Embeds like QupKake's
# input step (dies on a molecule that does not embed), writes processed/NAME.pt
# and output/OUT with a basic site at pKa {pka} on every nitrogen.
import sys, csv, os
from rdkit import Chem
from rdkit.Chem import AllChem
a = sys.argv; root = a[a.index("-r") + 1]; out = a[a.index("-o") + 1]
rows = list(csv.DictReader(open(a[2])))
for d in ("processed", "raw", "output"):
    os.makedirs(f"{{root}}/{{d}}", exist_ok=True)
for r in rows:
    m = Chem.AddHs(Chem.MolFromSmiles(r["smiles"]))
    AllChem.EmbedMolecule(m); AllChem.MMFFOptimizeMolecule(m)
w = Chem.SDWriter(f"{{root}}/output/{{out}}")
for r in rows:
    open(f"{{root}}/processed/{{r['name']}}.pt", "w").close()
    m = Chem.AddHs(Chem.MolFromSmiles(r["smiles"])); AllChem.EmbedMolecule(m, randomSeed=1)
    for at in m.GetAtoms():
        if at.GetSymbol() == "N":
            mm = Chem.Mol(m); mm.SetProp("_Name", r["name"]); mm.SetProp("idx", str(at.GetIdx()))
            mm.SetProp("pka_type", "basic"); mm.SetProp("pka", "{pka}"); w.write(mm)
w.close()
'''


@pytest.fixture
def fake_qupkake(tmp_path):
    """Path to an executable QupKake stand-in (every N basic at pKa 9.0)."""
    exe = tmp_path / "fake_qupkake"
    exe.write_text(FAKE_QUPKAKE.format(python=sys.executable, pka="9.0"))
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    return exe
