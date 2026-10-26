"""Tests for the MolSanitizer build step (molsanitizer/), without MolSanitizer.

    pytest molsanitizer/tests                 # from the repo root
    python molsanitizer/tests/test_msani_build.py   # same tests, no pytest needed

naming/plan/collect need only the standard library. The end-to-end tests run
build_shard.py against tests/fake_msani (a canned stand-in for msani 0.7.4's
Python API and CLI) and need RDKit + pandas; they skip without them.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
MSB = HERE.parent
REPO = MSB.parent
FAKE = HERE / "fake_msani"
sys.path.insert(0, str(MSB))

import naming  # noqa: E402


def have_rdkit() -> bool:
    try:
        import pandas  # noqa: F401
        import rdkit  # noqa: F401
        return True
    except ImportError:
        return False


def skip(msg):
    try:
        import pytest
        pytest.skip(msg)
    except ImportError:
        raise _Skip(msg)


class _Skip(Exception):
    pass


# ---- naming ----------------------------------------------------------------------

def test_contract_names():
    assert naming.contract_name("h10000007sNS", 0, 0) == "Ch10000007sNS00C"
    assert naming.contract_name("h10000007sNS", 10, 61) == "Ch10000007sNSAzC"
    for bad in ((62, 0), (0, 62), (-1, 0)):
        try:
            naming.contract_name("h10000007sNS", *bad)
        except ValueError:
            continue
        raise AssertionError(f"{bad} accepted")


def test_parse_state_id():
    p = "ZINCh10000007sNS"
    assert naming.parse_state_id(p, p) == (0, 0)
    assert naming.parse_state_id(p, p + "_2") == (2, 0)
    assert naming.parse_state_id(p, p + ".3") == (0, 3)
    assert naming.parse_state_id(p, p + "_02.10") == (2, 10)        # zero-padded forms
    assert naming.parse_state_id(p, "ZINCh10000007sNT_1") is None    # another parent
    assert naming.parse_state_id(p, p + "_1_1") is None              # unknown format: stop


def test_assign_indices_dense_per_state():
    keys = [(1, 1), (1, 2), (2, 0), (3, 1), (3, 2), (3, 3)]
    assert naming.assign_indices(keys) == [(0, 0), (0, 1), (1, 0), (2, 0), (2, 1), (2, 2)]
    assert naming.assign_indices([(0, 0)]) == [(0, 0)]


def test_parent_ids_match_llb_ids():
    """naming.py copies run_qupkake/llb_ids.py; they must agree while both exist."""
    src = REPO / "run_qupkake" / "llb_ids.py"
    if not src.exists():
        skip("run_qupkake/llb_ids.py is gone; naming.py is now the one definition")
    import importlib.util
    spec = importlib.util.spec_from_file_location("llb_ids", src)
    llb_ids = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(llb_ids)
    names = ["ZINCh10000007sNS", "ZINC000012345678", "ZINC12345", "ZINC000012345678.0",
             "ZINC", "", "CHEMBL25", "my benzoic acid", "edge_piperazine", "ZINCabc"]
    for n in names:
        assert naming.parent_id(n) == llb_ids.parent_id(n), n
        pid = naming.parent_id(n)[0]
        if pid:
            assert naming.base_id_of(pid) == llb_ids.base_id_of(pid), n


# ---- plan + collect (stdlib) ---------------------------------------------------------

def plan(tmp: Path, lines: list[str], per_shard=2, extra=()) -> Path:
    inp = tmp / "H18.smi"
    inp.write_text("".join(l + "\n" for l in lines))
    run = tmp / "run"
    r = subprocess.run([sys.executable, str(MSB / "plan.py"), str(inp), "-o", str(run),
                        "--per-shard", str(per_shard), *extra], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return run


def test_plan_rejects_and_shards(tmp_path=None):
    tmp = Path(tmp_path or tempfile.mkdtemp())
    run = plan(tmp, ["CCO ZINCh10000007sNS", "CCN ZINC000012345678.0", "CCC ZINCh10000007sNS",
                     "justone", "c1ccccc1 benzene", "CO ZINC12345"])
    meta = json.loads((run / "run.json").read_text())
    assert meta["n_input_lines"] == 6 and meta["n_parents"] == 3 and meta["n_shards"] == 2
    rej = (run / "plan_rejected.tsv").read_text()
    assert "malformed_zinc_id" in rej and "duplicate_parent_id" in rej and "unparsable_line" in rej
    shard1 = (run / "shards" / "00001.tsv").read_text().splitlines()
    assert shard1[0].split("\t") == ["0", "CCO", "ZINCh10000007sNS", "ZINCh10000007sNS"]
    assert (run / "shards" / "00002.tsv").read_text().split("\t")[3].strip() == "ZINC000000012345"
    assert "PH: 7\n" in (run / "msani_configurations.yaml").read_text()
    # same plan again: no-op; different options: refused
    again = subprocess.run([sys.executable, str(MSB / "plan.py"), str(tmp / "H18.smi"), "-o", str(run),
                            "--per-shard", "2"], capture_output=True, text=True)
    assert again.returncode == 0 and "nothing to do" in again.stdout
    other = subprocess.run([sys.executable, str(MSB / "plan.py"), str(tmp / "H18.smi"), "-o", str(run),
                            "--per-shard", "2", "--pH", "7.4"], capture_output=True, text=True)
    assert other.returncode != 0


def test_ranges():
    from collect import ranges
    assert ranges([1, 2, 3, 7, 9, 10]) == "1-3,7,9-10"
    assert ranges([5]) == "5"


# ---- end to end against the fake MolSanitizer -------------------------------------------

def env_with_fake() -> dict:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(FAKE) + os.pathsep + env.get("PYTHONPATH", "")
    env["PATH"] = os.pathsep.join([str(FAKE / "bin"), str(Path(sys.executable).parent), env["PATH"]])
    return env


def build(run: Path, k: int, scratch: Path):
    return subprocess.run([sys.executable, str(MSB / "build_shard.py"), "--run", str(run), "--shard", str(k),
                           "--scratch", str(scratch)], capture_output=True, text=True, env=env_with_fake())


def collect(run: Path):
    return subprocess.run([sys.executable, str(MSB / "collect.py"), str(run)], capture_output=True, text=True)


def read(path: Path) -> list[dict]:
    import csv
    delim = "," if path.suffix == ".csv" else "\t"
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def test_end_to_end(tmp_path=None):
    if not have_rdkit():
        skip("needs RDKit + pandas")
    tmp = Path(tmp_path or tempfile.mkdtemp())
    lines = ["CCN(C)Cc1ccccc1 ZINCh10000007sNS",       # invertomer pair
             "OC(=O)c1ccc(O)cc1 acid_parent",          # two states, hashed name
             "CC(N)C(=O)O ZINCh10000008aaa",           # 3 names: state 0 has 2 isomers
             "BrCCO ZINCh10000009bbb",                 # 3D failure recorded by msani
             "CC(=O)[O-].[Na+] ZINCh1000000Asal",      # removed in 2D
             "C[Se]C ZINCh1000000Bsel",                # exotic for db2
             "ClCCO ZINCh1000000Ccll",                  # AMSOL charge mismatch
             "C(=O)(O)c1ccc(O)cc1 dup_parent",         # uncanned: one state, built
             "CCCI ZINCh1000000Diod",                  # msani process dies; retried alone
             "C1CC1C(=O)O(C ZINCh1000000Ebad",         # invalid SMILES
             "OC(=O)c1ccc(O)cc1 acid_again",           # same structures as acid_parent
             ]
    run = plan(tmp, lines, per_shard=4)
    scratch = tmp / "scratch"
    scratch.mkdir()
    for k in (1, 2):
        r = build(run, k, scratch)
        assert r.returncode == 0, r.stdout + r.stderr
    r = collect(run)
    assert r.returncode == 1 and "not done: 3" in r.stderr            # shard 3 missing
    r = build(run, 3, scratch)
    assert r.returncode == 0, r.stdout + r.stderr
    r = collect(run)
    assert r.returncode == 0, r.stdout + r.stderr

    ledger = {row["input_name"]: row for row in read(run / "parents_ledger.tsv")}
    assert len(ledger) == len(lines)
    assert ledger["ZINCh10000007sNS"]["status"] == "built"
    assert ledger["ZINCh10000009bbb"]["status"] == "failed" and "canned failure" in ledger["ZINCh10000009bbb"]["reason"]
    assert "removed_in_2d" in ledger["ZINCh1000000Asal"]["reason"]
    assert ledger["ZINCh1000000Ebad"]["reason"] == "invalid_smiles"
    assert "CHARGE_MISMATCH" in ledger["ZINCh1000000Ccll"]["reason"]
    assert "EXOTIC" in ledger["ZINCh1000000Bsel"]["reason"]
    assert "msani exit 139" in ledger["ZINCh1000000Diod"]["reason"]
    assert ledger["acid_again"]["status"] == "failed" and "DUPLICATE" in ledger["acid_again"]["reason"]

    lib = {row["name"]: row for row in read(run / "library" / "library.tsv")}
    assert all(naming.CONTRACT.match(n) for n in lib)
    amine = sorted(n for n, row in lib.items() if row["parent_id"] == "ZINCh10000007sNS")
    assert amine == ["Ch10000007sNS00C", "Ch10000007sNS01C"]
    assert lib["Ch10000007sNS01C"]["invertomer_of"] == "Ch10000007sNS00C"
    assert lib["Ch10000007sNS00C"]["invertomer_of"] == ""
    ala = sorted((row["protomer_index"], row["stereo_index"]) for row in lib.values()
                 if row["parent_id"] == "ZINCh10000008aaa")
    assert ala == [("0", "0"), ("0", "1"), ("1", "0")]
    assert lib["Ch10000008aaa00C"]["net_charge"] == "0" and lib["Ch10000008aaa10C"]["net_charge"] == "0"
    acid = [row for row in lib.values() if row["input_name"] == "acid_parent"]
    assert {row["net_charge"] for row in acid} == {"-1", "0"} and acid[0]["base_id"].startswith("Z")

    ver = {row["protomer_name"]: row for row in read(run / "db2_verification.csv")}
    assert set(ver) == set(lib)
    ok = [n for n, row in ver.items() if row["status"] == "OK"]
    for n in ok:
        assert (run / "db2_run" / "out" / ver[n]["db2_file"]).exists()
    assert len(ok) == 2 + 2 + 3 + 1    # amine pair, acid states, alanine names, dup_parent

    stats = json.loads((run / "collect.json").read_text())
    assert stats["input_lines"] == len(lines) and stats["n_h_invertomers"] == 1
    prov = json.loads((run / "provenance.json").read_text())
    assert prov["software"]["molsanitizer"] == "0.7.4-fake"


def test_resume_keeps_built_db2(tmp_path=None):
    if not have_rdkit():
        skip("needs RDKit + pandas")
    tmp = Path(tmp_path or tempfile.mkdtemp())
    run = plan(tmp, ["CC(N)C(=O)O ZINCh10000008aaa"], per_shard=5)
    scratch = tmp / "s"
    scratch.mkdir()
    assert build(run, 1, scratch).returncode == 0
    done = run / "tasks" / "00001" / "DONE.json"
    r = build(run, 1, scratch)
    assert "already done" in r.stdout
    # an interrupted shard: DONE gone, db2s present -> kept, not rebuilt
    done.unlink()
    gz = run / "db2_run" / "out" / "00001" / "Ch10000008aaa00C.db2.gz"
    before = gz.stat().st_mtime_ns
    r = build(run, 1, scratch)
    assert r.returncode == 0, r.stderr
    assert gz.stat().st_mtime_ns == before
    states = read(run / "tasks" / "00001" / "states.tsv")
    assert {s["detail"] for s in states} == {"kept from an earlier attempt"}


def test_check_db2_header_parsing(tmp_path=None):
    tmp = Path(tmp_path or tempfile.mkdtemp())
    from build_shard import check_db2
    rec = ("M %16s %9s %3d %3d %6d %6d %6d %6d %6d %6d\n" % ("Ch10000007sNS00C", "fake", 20, 21, 300, 40, 2, 7, 5, 0)
           + "M %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n" % (1.0002, -40.1, 2.0, -38.1, 310.0)
           + "M %-76s\nM %-76s\nM %+10.4f\n" % ("CC[NH+](C)Cc1ccccc1", "fake", 999.999))
    f = tmp / "x.db2"
    f.write_text(rec * 3)
    status, _, facts = check_db2(f, "Ch10000007sNS00C", 1, 0.05)
    assert status == "OK" and facts["db2_records"] == 3
    assert check_db2(f, "Ch10000007sNS00C", 0, 0.05)[0] == "CHARGE_MISMATCH"
    assert check_db2(f, "Ch10000007sNS01C", 1, 0.05)[0] == "RENAMED"
    f.write_text("")
    assert check_db2(f, "Ch10000007sNS00C", 1, 0.05)[0] == "EMPTY"


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"ok    {name}")
            except _Skip as e:
                print(f"skip  {name}: {e}")
            except Exception as e:  # noqa: BLE001
                failed += 1
                import traceback
                traceback.print_exc()
                print(f"FAIL  {name}: {e}")
    sys.exit(1 if failed else 0)
