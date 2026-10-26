"""naming/assign_names.py: MolSanitizer's NAME_state.stereo -> 16-character contract names.

The MolSanitizer outputs below are copied from real `msani` 0.7.4 runs
(--tautomers --protonation --pH 6 --pH_range 1 --stereoisomers), so the suffix
format is MolSanitizer's own, not a guess at it.
"""
import csv
import json

import pytest

pytest.importorskip("rdkit")
import assign_names                                              # noqa: E402
from assign_names import split_ms_name                           # noqa: E402

INPUTS = """\
CC(N)Cc1c[nH]cn1 ZINCaa000012wx34
CC(C)N1CCN(C)CC1C ZINCaa000013yz56
O=C1NC(=O)C(C)(c2ccccc2)N1 ZINCaa000011uv12
Cc1c[nH]cn1 ZINCaa000006kl12
CC(=O)CC(C)=O ZINCaa000009qr78
C1=CC=CC=1 weird_name_1.2
C.CC ZINCaa000016ef56
CCO ZINCaa000017gh78
CC(F)Cl ZINCaa000018ij90
CCC(C)CC ZINCaa000019kl12
OCC ZINCaa000020mn34
"""

MS_OUT = """\
C[C@H]([NH3+])Cc1c[nH]c[nH+]1 ZINCaa000012wx34_1.1
C[C@@H]([NH3+])Cc1c[nH]c[nH+]1 ZINCaa000012wx34_1.2
C[C@H]([NH3+])Cc1c[nH]cn1 ZINCaa000012wx34_2.1
C[C@@H]([NH3+])Cc1c[nH]cn1 ZINCaa000012wx34_2.2
CC(C)[N@@H+]1CCN(C)C[C@H]1C ZINCaa000013yz56_1.1
CC(C)[N@@H+]1CCN(C)C[C@@H]1C ZINCaa000013yz56_1.2
CC(C)[N@H+]1CCN(C)C[C@H]1C ZINCaa000013yz56_1.3
CC(C)[N@H+]1CCN(C)C[C@@H]1C ZINCaa000013yz56_1.4
C[C@]1(c2ccccc2)NC(=O)NC1=O ZINCaa000011uv12.1
C[C@@]1(c2ccccc2)NC(=O)NC1=O ZINCaa000011uv12.2
Cc1c[nH]c[nH+]1 ZINCaa000006kl12_1
Cc1c[nH]cn1 ZINCaa000006kl12_2
Cc1c[nH+]c[nH]1 ZINCaa000006kl12_3
Cc1cnc[nH]1 ZINCaa000006kl12_4
CC(=O)CC(C)=O ZINCaa000009qr78
C1=CC=CC=1 weird_name_1.2
C.CC ZINCaa000016ef56
CC(F)Cl ZINCaa000018ij90
CCO ZINCaa000020mn34
CCO ZINCaa000017gh78
CC ZINCzz999999zz99
"""


def test_split_ms_name_reads_state_and_stereo():
    known = {"ZINCaa000013yz56", "weird_name_1.2", "weird_name"}
    assert split_ms_name("ZINCaa000013yz56", known) == ("ZINCaa000013yz56", 0, 0)
    assert split_ms_name("ZINCaa000013yz56_2", known) == ("ZINCaa000013yz56", 2, 0)
    assert split_ms_name("ZINCaa000013yz56.3", known) == ("ZINCaa000013yz56", 0, 3)
    assert split_ms_name("ZINCaa000013yz56_2.3", known) == ("ZINCaa000013yz56", 2, 3)
    # an input name with '_' and '.' in it is matched whole, never cut
    assert split_ms_name("weird_name_1.2", known) == ("weird_name_1.2", 0, 0)
    assert split_ms_name("weird_name_1.2_3.1", known) == ("weird_name_1.2", 3, 1)
    assert split_ms_name("ZINCzz999999zz99", known) is None
    assert split_ms_name("ZINCaa000013yz56_x", known) is None


@pytest.fixture
def library(tmp_path):
    (tmp_path / "in.smi").write_text(INPUTS)
    (tmp_path / "in_clean.smi").write_text(MS_OUT)
    assign_names.main(["--inputs", str(tmp_path / "in.smi"), "--ms", str(tmp_path / "in_clean.smi"),
                       "-o", str(tmp_path / "lib"), "--ph", "6+-1"])
    lib = tmp_path / "lib"
    with open(lib / "library.tsv") as fh:
        rows = {r["ms_name"]: r for r in csv.DictReader(fh, delimiter="\t")}
    with open(lib / "dropped.tsv") as fh:
        dropped = {(r["input_name"] or r["ms_name"]): r["reason"]
                   for r in csv.DictReader(fh, delimiter="\t")}
    return rows, dropped, json.loads((lib / "stats.json").read_text()), lib


def test_state_and_stereo_become_the_last_two_characters(library):
    rows, *_ = library
    assert rows["ZINCaa000012wx34_1.1"]["name"] == "Caa000012wx3400C"
    assert rows["ZINCaa000012wx34_1.2"]["name"] == "Caa000012wx3401C"
    assert rows["ZINCaa000012wx34_2.1"]["name"] == "Caa000012wx3410C"
    assert rows["ZINCaa000012wx34_2.2"]["name"] == "Caa000012wx3411C"
    assert rows["ZINCaa000011uv12.2"]["name"] == "Caa000011uv1201C"     # stereo only
    assert rows["ZINCaa000006kl12_2"]["name"] == "Caa000006kl1210C"     # states only
    assert rows["ZINCaa000009qr78"]["name"] == "Caa000009qr7800C"       # neither
    r = rows["ZINCaa000012wx34_1.1"]
    assert (r["protomer_index"], r["stereo_index"], r["net_charge"]) == ("0", "0", "2")
    assert r["parent_id"] == r["zinc_id"] == "ZINCaa000012wx34"
    assert r["site_source"] == "molsanitizer" and r["ph_values"] == "6+-1"


def test_non_zinc_parent_gets_a_hashed_base_id(library):
    rows, *_ = library
    r = rows["weird_name_1.2"]
    assert r["parent_id"].startswith("LLBZ") and r["zinc_id"] == ""
    assert r["name"] == f"C{r['parent_id'][3:]}00C" and r["input_name"] == "weird_name_1.2"


def test_n_h_invertomers_point_at_their_partner(library):
    rows, *_ = library
    inv = {k: rows[k]["invertomer_of"] for k in rows if k.startswith("ZINCaa000013yz56")}
    # .1/.3 and .2/.4 differ only at the N-H+ centre; .1/.2 differ at a carbon
    assert inv == {"ZINCaa000013yz56_1.1": "", "ZINCaa000013yz56_1.2": "",
                   "ZINCaa000013yz56_1.3": "Caa000013yz5600C",
                   "ZINCaa000013yz56_1.4": "Caa000013yz5601C"}


def test_drops_are_recorded_with_a_reason(library):
    rows, dropped, stats, _ = library
    assert dropped["ZINCaa000016ef56"] == "multiple_fragments"
    assert dropped["ZINCaa000018ij90"] == "stereo_not_closed"       # MS left a centre open
    assert dropped["ZINCaa000020mn34"] == "duplicate_structure"     # same SMILES, later parent id
    assert rows["ZINCaa000017gh78"]["smiles"] == "CCO"
    assert dropped["ZINCzz999999zz99"] == "unknown_input_name"
    # _3 is _1's imidazolium with the + drawn on the other N: one species, one name
    assert "ZINCaa000006kl12_3" not in rows
    assert rows["ZINCaa000006kl12_4"]["name"] == "Caa000006kl1220C"   # a real tautomer
    assert dropped["ZINCaa000019kl12"] == "not_returned_by_molsanitizer"
    assert stats["inputs"] == 11 and stats["names"] == len(rows) == 16
    assert stats["n_h_invertomers"] == 2


def test_library_smi_is_what_pass_two_reads(library):
    rows, _, _, lib = library
    lines = (lib / "library.smi").read_text().splitlines()
    assert len(lines) == len(rows)
    for line in lines:
        smi, name = line.split()
        assert len(name) == 16 and name[0] == name[-1] == "C"


def test_refuses_to_overwrite(library):
    *_, lib = library
    with pytest.raises(SystemExit):
        assign_names.main(["--inputs", str(lib.parent / "in.smi"),
                           "--ms", str(lib.parent / "in_clean.smi"), "-o", str(lib)])
