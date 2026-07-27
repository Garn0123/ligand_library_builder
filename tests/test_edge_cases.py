"""#6 — edge cases the pipeline promises to handle.

Truncated inputs must fail loudly (not silently drop a molecule), a file that
carries .gz but isn't gzipped must be read anyway and flagged, duplicate files
(hard/symlinks) must be deduped by inode, and an empty shard must still emit its
bookkeeping so phase 3 sees every shard.
"""

import os

import db2common as C  # noqa: F401  (kept for symmetry / future assertions)
import db2gen


def test_truncated_input_fails_loudly(tmp_path, run_script):
    src = tmp_path / "src"
    good = db2gen.make_record("ZINC00000001")
    partial = "M ZINC00000002 x\nA 0 C 0 0 0\n"     # no 'E' terminator
    db2gen.write_gz(os.path.join(str(src), "t.db2.gz"), good + partial)

    out = tmp_path / "out"
    res = run_script("01_chunk.py", "-i", src, "-o", out,
                     "--progress-interval", 0, check=False)

    assert res.returncode != 0
    assert "mid-record" in res.stderr or "truncat" in res.stderr.lower()


def test_plaintext_with_gz_name_is_read_and_flagged(tmp_path, run_script):
    src = tmp_path / "src"
    # Valid records, stored UNCOMPRESSED under a .db2.gz name.
    db2gen.write_plain(
        os.path.join(str(src), "p.db2.gz"),
        db2gen.records_text([("ZINC00000001", 1), ("ZINC00000002", 1)]))

    out = tmp_path / "out"
    res = run_script("01_chunk.py", "-i", src, "-o", out,
                     "--progress-interval", 0)

    assert res.returncode == 0
    assert "NOT gzipped" in res.stderr
    # The records were still read and chunked.
    ids = db2gen.chunk_ids(os.path.join(str(out), "chunk_00001.db2.gz"))
    assert ids == ["ZINC00000001", "ZINC00000002"]


def test_hardlink_duplicate_is_deduped(tmp_path, run_script):
    src = tmp_path / "src"
    a = os.path.join(str(src), "a.db2.gz")
    db2gen.write_gz(a, db2gen.records_text([("ZINC00000001", 1),
                                            ("ZINC00000002", 1)]))
    os.link(a, os.path.join(str(src), "b.db2.gz"))    # same inode

    out = tmp_path / "out"
    res = run_script("01_chunk.py", "-i", src, "-o", out,
                     "--progress-interval", 0)

    assert res.returncode == 0
    assert "duplicate" in res.stderr.lower()
    # Deduped to 2 molecules, not doubled to 4.
    all_ids = []
    for name in sorted(os.listdir(str(out))):
        if name.endswith(".db2.gz"):
            all_ids += db2gen.chunk_ids(os.path.join(str(out), name))
    assert sorted(all_ids) == ["ZINC00000001", "ZINC00000002"]


def test_empty_shard_still_emits_bookkeeping(tmp_path, run_script):
    src = tmp_path / "src"
    db2gen.build_tree(str(src), {"a.db2.gz": [("ZINC00000001", 1)]})  # one file

    work = tmp_path / "work"
    # More shards than files -> shards 1..3 are empty (the file goes to shard 0).
    run_script("make_shards.py", "-i", src, "-o", work, "-S", 4)
    run_script("p1_collect.py", "-w", work, "-s", 3, "-N", 2)

    counts = os.path.join(str(work), "counts", "shard_00003.tsv")
    man = os.path.join(str(work), "manifests", "shard_00003.tsv")
    assert os.path.exists(counts)
    assert os.path.exists(man)

    with open(counts) as fh:
        rows = [ln.strip() for ln in fh if ln.strip()]
    assert rows == ["0\t0\t0", "1\t0\t0"]     # one zero row per bin
    assert os.path.getsize(man) == 0          # empty manifest, but present
