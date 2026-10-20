"""check_bcl.py finds BCL entries in db2_converter's config by name, and a
missing path among them fails the db2 preflight."""
import check_bcl


def test_finds_bcl_entries_by_level_name_or_value():
    config = {"all": {"AMSOLEXE": "/x/amsol", "BCLBASE": "/opt/bcl"},
              "bcl": {"EXE": "/opt/bcl/bin/bcl-apps-static.exe"},
              "conformator": {"CONF_EXE": "/home/u/bin/conformator"}}
    found = {(lvl, name) for lvl, name, _ in check_bcl.bcl_entries(config)}
    assert found == {("all", "BCLBASE"), ("bcl", "EXE")}


def test_skips_levels_that_are_not_mappings():
    assert check_bcl.bcl_entries({"version": "0.3", "bcl": {"X": "/b"}}) == [("bcl", "X", "/b")]
