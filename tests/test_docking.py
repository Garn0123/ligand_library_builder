"""docking/ tools: collapse-invertomers."""
import pytest

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
