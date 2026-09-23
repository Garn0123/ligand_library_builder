"""check_inputs.py -- the preflight integrity scanner.

Classifies every file in the tree and writes clean/dirty lists, so a damaged
library can be found and split in one pass instead of one crashed shard at a
time.
"""

import os

import db2gen


def test_scanner_classifies_and_splits_lists(tmp_path, run_script):
    src = tmp_path / "src"
    # One clean file, one corrupt (bad deflate body), one truncated (valid gzip
    # but no final E), and one plaintext-with-.gz-name (usable, just flagged).
    db2gen.write_gz(os.path.join(str(src), "ok.db2.gz"),
                    db2gen.records_text([("ZINC00000001", 2), ("ZINC00000002", 2)]))
    db2gen.write_corrupt_gz(
        os.path.join(str(src), "corrupt.db2.gz"),
        db2gen.records_text([("ZINC{:08d}".format(i), 2) for i in range(50)]))
    db2gen.write_gz(
        os.path.join(str(src), "trunc.db2.gz"),
        db2gen.make_record("ZINC00000003") + "M ZINC00000004 x\nA 0 C 0 0 0\n")
    db2gen.write_plain(os.path.join(str(src), "plain.db2.gz"),
                       db2gen.records_text([("ZINC00000005", 1)]))

    good = tmp_path / "good.txt"
    bad = tmp_path / "bad.txt"
    res = run_script("check_inputs.py", "-i", src,
                     "--good-list", good, "--bad-list", bad,
                     "--progress-interval", 0, check=False)

    assert res.returncode != 0                     # corrupt + truncated present
    assert "Traceback" not in res.stderr

    good_names = {os.path.basename(p) for p in good.read_text().split()}
    bad_names = {os.path.basename(p) for p in bad.read_text().split()}
    assert good_names == {"ok.db2.gz", "plain.db2.gz"}
    assert bad_names == {"corrupt.db2.gz", "trunc.db2.gz"}


def test_scanner_clean_tree_exits_zero(tmp_path, run_script):
    src = tmp_path / "src"
    db2gen.build_tree(str(src), {
        "a.db2.gz": [("ZINC00000001", 1)],
        "b.db2.gz": [("ZINC00000002", 2), ("ZINC00000003", 1)]})

    res = run_script("check_inputs.py", "-i", src, "--progress-interval", 0)
    assert res.returncode == 0
    assert "all files read cleanly" in res.stderr
