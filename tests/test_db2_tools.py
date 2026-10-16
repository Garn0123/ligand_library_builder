"""db2-side tooling: verify-db2 (shallow scan, several output paths, charge check),
llb timing's file discovery, and collapse-invertomers."""
import csv
import gzip

import pytest

import verify_db2_build
from timing_report import build_tables


def db2(path, name, charge):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write(f"M {name:<16} x\nM {charge} 0 0\nM CC\n")


@pytest.fixture
def run_tree(tmp_path):
    out = tmp_path / "out"
    db2(out / "0000" / "CAAAAAAAAAAAA00C.db2.gz", "CAAAAAAAAAAAA00C", 1.0)
    # a build's working copy synced mid-build: not a result
    db2(out / "0000" / "CAAAAAAAAAAAA00C" / "db2" / "output.0.db2.gz", "CAAAAAAAAAAAA00C", 1.0)
    db2(out / "0001" / "CAAAAAAAAAAAB00C.db2.gz", "CAAAAAAAAAAAB00C", 0.0)
    (out / "0001" / "s.0001.smi.conformator.faillist").write_text(
        "CC\tCAAAAAAAAAAAC00C\tconformator\t1generate\n")
    lib = tmp_path / "library.tsv"
    lib.write_text("name\tnet_charge\tinput_name\n"
                   "CAAAAAAAAAAAA00C\t1\tp1\nCAAAAAAAAAAAB00C\t0\tp2\n"
                   "CAAAAAAAAAAAC00C\t0\tp3\nCAAAAAAAAAAAD00C\t0\tp4\n")
    return tmp_path


def statuses(report):
    with open(report) as fh:
        return {r["protomer_name"] or r["db2_file"]: r["status"] for r in csv.DictReader(fh)}


def test_verify_ignores_working_copies_and_reads_faillists(run_tree):
    rep = run_tree / "rep.csv"
    rc = verify_db2_build.main([str(run_tree / "library.tsv"), "--outputpath",
                                str(run_tree / "out"), "-o", str(rep), "--threads", "2"])
    st = statuses(rep)
    assert st["CAAAAAAAAAAAA00C"] == "OK" and st["CAAAAAAAAAAAB00C"] == "OK"
    assert st["CAAAAAAAAAAAC00C"] == "BUILD_FAILED"
    assert st["CAAAAAAAAAAAD00C"] == "MISSING"
    assert "UNEXPECTED" not in st.values()
    assert rc == 1                                   # a silent loss fails the run


def test_verify_combines_output_paths(run_tree):
    db2(run_tree / "fallback" / "0000" / "CAAAAAAAAAAAC00C.db2.gz", "CAAAAAAAAAAAC00C", 0.0)
    db2(run_tree / "fallback" / "0000" / "CAAAAAAAAAAAD00C.db2.gz", "CAAAAAAAAAAAD00C", 0.0)
    rep = run_tree / "rep2.csv"
    rc = verify_db2_build.main([str(run_tree / "library.tsv"), "--outputpath",
                                str(run_tree / "out"), "--outputpath", str(run_tree / "fallback"),
                                "-o", str(rep)])
    assert set(statuses(rep).values()) == {"OK"} and rc == 0


def test_verify_catches_a_charge_mismatch(run_tree):
    db2(run_tree / "out" / "0001" / "CAAAAAAAAAAAB00C.db2.gz", "CAAAAAAAAAAAB00C", -1.0)
    rep = run_tree / "rep3.csv"
    verify_db2_build.main([str(run_tree / "library.tsv"), "--outputpath",
                           str(run_tree / "out"), "-o", str(rep)])
    assert statuses(rep)["CAAAAAAAAAAAB00C"] == "CHARGE_MISMATCH"


def test_timing_finds_build_tables_without_walking(tmp_path):
    (tmp_path / "out" / "0003").mkdir(parents=True)
    (tmp_path / "out" / "0003" / "llb_build.tsv").write_text("name\n")
    deep = tmp_path / "out" / "0003" / "CXXXXXXXXXXXX00C" / "llb_build.tsv"
    deep.parent.mkdir()
    deep.write_text("name\n")
    assert build_tables(tmp_path) == [tmp_path / "out" / "0003" / "llb_build.tsv"]


pd = pytest.importorskip("pandas")


def test_collapse_invertomers_keeps_best_score_per_group(tmp_path):
    import collapse_invertomers
    (tmp_path / "lib.tsv").write_text(
        "name\tinvertomer_of\nCA00C\t\nCA01C\tCA00C\nCB00C\t\nCC00C\t\nCC01C\tCC00C\n")
    (tmp_path / "s.tsv").write_text(
        "molecule_id\tt1\tt2\nCA00C\t-30\t-20\nCA01C\t-35\t-10\nCB00C\t-40\t\n"
        "CC00C\t\t-5\nCC01C\t-12\t-6\n")
    collapse_invertomers.main([str(tmp_path / "s.tsv"), "--library", str(tmp_path / "lib.tsv"),
                               "-o", str(tmp_path / "out.tsv")])
    out = pd.read_csv(tmp_path / "out.tsv", sep="\t").set_index("molecule_id")
    assert list(out.index) == ["CA00C", "CB00C", "CC00C"]
    assert out.loc["CA00C", "t1"] == -35 and out.loc["CA00C", "t2"] == -20
    assert out.loc["CC00C", "t1"] == -12 and out.loc["CC00C", "t2"] == -6
    prov = pd.read_csv(tmp_path / "out.invertomers.tsv", sep="\t")
    assert set(prov["species"]) == {"CA00C", "CC00C"}
