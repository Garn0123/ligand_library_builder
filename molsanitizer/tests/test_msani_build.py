"""Tests for the MolSanitizer build step (molsanitizer/), without MolSanitizer.

    pytest molsanitizer/tests                 # from the repo root
    python molsanitizer/tests/test_msani_build.py   # same tests, no pytest needed

plan/collect need only the standard library. The end-to-end tests run pass 1,
the REAL naming/assign_names.py, and pass 2 against tests/fake_msani (a canned
stand-in for msani 0.7.4's Python API and CLI); they need RDKit + pandas and skip
without them. naming/ has its own tests (tests/test_names.py).
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

sys.path.insert(0, str(REPO / "naming"))
CONTRACT = __import__("re").compile(r"^C[0-9A-Za-z]{14}C$")


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
    shard1 = (run / "pass1" / "shards" / "00001.tsv").read_text().splitlines()
    assert shard1[0].split("\t") == ["0", "CCO", "ZINCh10000007sNS", "ZINCh10000007sNS"]
    assert (run / "pass1" / "shards" / "00002.tsv").read_text().split("\t")[3].strip() == "ZINC000000012345"
    assert "PH: 7\n" in (run / "msani_configurations.yaml").read_text()
    # same plan again: no-op; different options: refused
    again = subprocess.run([sys.executable, str(MSB / "plan.py"), str(tmp / "H18.smi"), "-o", str(run),
                            "--per-shard", "2"], capture_output=True, text=True)
    assert again.returncode == 0 and "nothing to do" in again.stdout
    other = subprocess.run([sys.executable, str(MSB / "plan.py"), str(tmp / "H18.smi"), "-o", str(run),
                            "--per-shard", "2", "--pH", "6"], capture_output=True, text=True)
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


def build(run: Path, npass: int, k: int, scratch: Path):
    return subprocess.run([sys.executable, str(MSB / "build_shard.py"), "--run", str(run), "--pass", str(npass),
                           "--shard", str(k), "--scratch", str(scratch)],
                          capture_output=True, text=True, env=env_with_fake())


def names_step(run: Path, per_shard: int):
    return subprocess.run([sys.executable, str(MSB / "names_step.py"), str(run), "--per-shard", str(per_shard)],
                          capture_output=True, text=True, env=env_with_fake())


def collect(run: Path):
    return subprocess.run([sys.executable, str(MSB / "collect.py"), str(run)], capture_output=True, text=True)


def read(path: Path) -> list[dict]:
    import csv
    delim = "," if path.suffix == ".csv" else "\t"
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def run_all(run: Path, scratch: Path, per_shard2: int = 3):
    n1 = json.loads((run / "run.json").read_text())["n_shards"]
    for k in range(1, n1 + 1):
        r = build(run, 1, k, scratch)
        assert r.returncode == 0, r.stdout + r.stderr
    r = names_step(run, per_shard2)
    assert r.returncode == 0, r.stdout + r.stderr
    n2 = json.loads((run / "pass2" / "plan.json").read_text())["n_shards"]
    for k in range(1, n2 + 1):
        r = build(run, 2, k, scratch)
        assert r.returncode == 0, r.stdout + r.stderr
    return n2


def test_end_to_end(tmp_path=None):
    if not have_rdkit():
        skip("needs RDKit + pandas")
    tmp = Path(tmp_path or tempfile.mkdtemp())
    lines = ["CCN(C)Cc1ccccc1 ZINCh10000007sNS",       # N-H+ invertomer pair
             "OC(=O)c1ccc(O)cc1 acid_parent",          # two states, hashed name
             "CC(N)C(=O)O ZINCh10000008aaa",           # 2 closed isomers + 1 state not stereo-closed
             "BrCCO ZINCh10000009bbb",                 # 3D failure recorded by msani
             "CC(=O)[O-].[Na+] ZINCh1000000Asal",      # removed in pass 1
             "C[Se]C ZINCh1000000Bsel",                # exotic for db2
             "ClCCO ZINCh1000000Ccll",                 # AMSOL charge mismatch
             "C(=O)(O)c1ccc(O)cc1 dup_parent",         # same structure as acid_parent's neutral state
             "CCCI ZINCh1000000Diod",                  # msani process dies; retried alone
             "C1CC1C(=O)O(C ZINCh1000000Ebad",         # invalid SMILES
             "OC(=O)c1ccc(O)cc1 acid_again",           # every state another parent's
             "CCO ZINC000012345678.0",                 # damaged ZINC id: plan-rejected
             ]
    run = plan(tmp, lines, per_shard=4)
    scratch = tmp / "scratch"
    scratch.mkdir()
    assert build(run, 1, 1, scratch).returncode == 0
    r = names_step(run, 3)
    assert r.returncode != 0 and "pass 1 not finished" in (r.stdout + r.stderr)
    assert collect(run).returncode == 1
    run_all(run, scratch)
    r = collect(run)
    assert r.returncode == 0, r.stdout + r.stderr

    ledger = {row["input_name"]: row for row in read(run / "parents_ledger.tsv")}
    assert len(ledger) == len(lines)
    assert ledger["ZINCh10000007sNS"]["status"] == "built"
    assert "canned failure" in ledger["ZINCh10000009bbb"]["reason"]
    assert ledger["ZINCh1000000Asal"]["reason"].startswith("pass1: removed_in_pass1")
    assert ledger["ZINCh1000000Ebad"]["reason"] == "pass1: invalid_smiles"
    assert "CHARGE_MISMATCH" in ledger["ZINCh1000000Ccll"]["reason"]
    assert "EXOTIC" in ledger["ZINCh1000000Bsel"]["reason"]
    assert "msani exit 139" in ledger["ZINCh1000000Diod"]["reason"]
    assert ledger["acid_again"]["status"] == "failed" and "duplicate_structure" in ledger["acid_again"]["reason"]
    assert ledger["ZINC000012345678.0"]["reason"] == "plan: malformed_zinc_id"
    for n in ("dup_parent", "acid_parent", "ZINCh10000008aaa"):
        assert ledger[n]["status"] == "built", ledger[n]

    lib = {row["name"]: row for row in read(run / "library" / "library.tsv")}
    assert all(CONTRACT.match(n) for n in lib)
    assert sorted(n for n, row in lib.items() if row["parent_id"] == "ZINCh10000007sNS") == \
        ["Ch10000007sNS00C", "Ch10000007sNS01C"]
    assert lib["Ch10000007sNS01C"]["invertomer_of"] == "Ch10000007sNS00C"
    assert lib["Ch10000007sNS00C"]["ms_name"] == "ZINCh10000007sNS.1"
    assert lib["Ch10000007sNS00C"]["ph_values"] == "7.4+-1"
    ala = [row for row in lib.values() if row["parent_id"] == "ZINCh10000008aaa"]
    assert sorted(r["stereo_index"] for r in ala) == ["0", "1"]          # the open state was dropped
    assert "stereo_not_closed" in (run / "library" / "dropped.tsv").read_text()

    ver = {row["protomer_name"]: row for row in read(run / "db2_verification.csv")}
    assert set(ver) == set(lib)
    ok = [n for n, row in ver.items() if row["status"] == "OK"]
    for n in ok:
        assert (run / "db2_run" / "out" / ver[n]["db2_file"]).exists()
    assert len(ok) == 2 + 1 + 2 + 1    # amine pair, acid anion, alanine isomers, dup_parent's neutral

    stats = json.loads((run / "collect.json").read_text())
    assert stats["input_lines"] == len(lines) and stats["n_h_invertomers"] == 1
    assert json.loads((run / "provenance.json").read_text())["software"]["molsanitizer"] == "0.7.4-fake"


def test_resume_keeps_built_db2(tmp_path=None):
    if not have_rdkit():
        skip("needs RDKit + pandas")
    tmp = Path(tmp_path or tempfile.mkdtemp())
    run = plan(tmp, ["CC(N)C(=O)O ZINCh10000008aaa"], per_shard=5)
    scratch = tmp / "s"
    scratch.mkdir()
    run_all(run, scratch, per_shard2=5)
    done = run / "pass2" / "tasks" / "00001" / "DONE.json"
    assert "already done" in build(run, 2, 1, scratch).stdout
    done.unlink()                       # an interrupted shard: db2s present -> kept, not rebuilt
    gz = run / "db2_run" / "out" / "00001" / "Ch10000008aaa00C.db2.gz"
    before = gz.stat().st_mtime_ns
    r = build(run, 2, 1, scratch)
    assert r.returncode == 0, r.stderr
    assert gz.stat().st_mtime_ns == before
    assert {s["detail"] for s in read(run / "pass2" / "tasks" / "00001" / "states.tsv")} == \
        {"kept from an earlier attempt"}
    assert collect(run).returncode == 0


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
