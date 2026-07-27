"""#4 — labeling: duplicate suffixing and stale-label rejection.

02_label decides new ids for molecules whose catalog id legitimately repeats;
03_apply rewrites the chunks by POSITION and — crucially — verifies the id at
each position matches what labels.tsv expects before touching it. A mismatch
means the labels are stale (e.g. 01 was re-run with a different chunk size after
02), and the run must abort rather than silently corrupt data.
"""

import os
from collections import Counter

import db2common as C
import db2gen


def _make_src_with_dupes(src):
    """ZINC00000005 appears 3x, ZINC00000009 appears 2x, the rest are unique."""
    a = [("ZINC{:08d}".format(i), 1) for i in range(10)]   # includes ...05, ...09
    b = [("ZINC00000005", 1), ("ZINC00000009", 1), ("ZINC00000005", 1),
         ("ZINC00000015", 1), ("ZINC00000016", 1)]
    return db2gen.build_tree(str(src), {"a.db2.gz": a, "b.db2.gz": b})


def test_duplicate_labeling_happy_path(tmp_path, run_script):
    src = tmp_path / "src"
    _make_src_with_dupes(src)

    chunks = tmp_path / "chunks"
    labels = tmp_path / "labels.tsv"
    labelled = tmp_path / "labelled"

    run_script("01_chunk.py", "-i", src, "-o", chunks, "-n", 50,
               "--progress-interval", 0)
    run_script("02_label.py", "-m", os.path.join(str(chunks), "manifest.tsv"),
               "-o", labels, "--mode", "duplicates", "--progress-interval", 0)
    run_script("03_apply.py", "-c", chunks, "-L", labels, "-o", labelled,
               "--copy-unedited")

    ids = []
    for name in sorted(os.listdir(str(labelled))):
        if name.endswith(".db2.gz"):
            ids += db2gen.chunk_ids(os.path.join(str(labelled), name))
    counts = Counter(ids)

    # Duplicated ids are suffixed into distinct ids; the bare id is gone.
    assert "ZINC00000005" not in counts
    assert counts["ZINC00000005_1"] == 1
    assert counts["ZINC00000005_2"] == 1
    assert counts["ZINC00000005_3"] == 1
    assert counts["ZINC00000009_1"] == 1
    assert counts["ZINC00000009_2"] == 1
    # A never-duplicated id is left untouched.
    assert counts["ZINC00000000"] == 1
    # Net effect: every id is now unique, and the record count is unchanged.
    assert all(v == 1 for v in counts.values())
    assert len(ids) == 15


def test_stale_labels_are_rejected(tmp_path, run_script):
    src = tmp_path / "src"
    db2gen.build_tree(str(src), {
        "a.db2.gz": [("ZINC{:08d}".format(i), 1) for i in range(5)]})

    chunks = tmp_path / "chunks"
    run_script("01_chunk.py", "-i", src, "-o", chunks, "-n", 50,
               "--progress-interval", 0)

    # Hand-craft a labels.tsv that claims the WRONG original id at position 0.
    # (Deterministic: the real id there is ZINC00000000, never ZINC_WRONG.)
    labels = tmp_path / "labels.tsv"
    with open(str(labels), "w") as fh:
        fh.write("\t".join(C.LABELS_HEADER) + "\n")
        fh.write("chunk_00001.db2.gz\t0\tZINC_WRONG\tZINC_WRONG_1\n")

    out = tmp_path / "out"
    res = run_script("03_apply.py", "-c", chunks, "-L", labels, "-o", out,
                     check=False)

    assert res.returncode != 0                     # aborted, not silently applied
    # The mismatched chunk must not have been written (all-or-nothing).
    assert not os.path.exists(os.path.join(str(out), "chunk_00001.db2.gz"))
    assert not os.path.exists(os.path.join(str(out), "chunk_00001.db2.gz.tmp"))
    assert "MISMATCH" in res.stderr or "do not match" in res.stderr
