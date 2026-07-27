"""#3 — the serial and parallel pipelines produce the same library.

They assign molecules to chunks in different orders, so chunk *contents* differ,
but the guarantees that matter must hold for both: the id multiset is preserved
exactly (nothing lost or duplicated), and every chunk is a whole number of
records that starts with 'M' and ends with 'E'.
"""

import os
from collections import Counter

import db2common as C
import db2gen

FILES = {
    "s/a.db2.gz": [("ZINC{:08d}".format(i), (i % 3) + 1) for i in range(12)],
    "l/b.db2.gz": [("ZINC{:08d}".format(i), (i % 4) + 2) for i in range(12, 28)],
    "l/c.db2.gz": [("ZINC{:08d}".format(i), 1) for i in range(28, 40)],
}


def _multiset(chunk_dir):
    ids = Counter()
    for name in sorted(os.listdir(chunk_dir)):
        if name.endswith(".db2.gz"):
            for zid in db2gen.chunk_ids(os.path.join(chunk_dir, name)):
                ids[zid] += 1
    return ids


def _assert_whole_records(chunk_dir):
    found_any = False
    for name in sorted(os.listdir(chunk_dir)):
        if not name.endswith(".db2.gz"):
            continue
        found_any = True
        with C.open_gz_text(os.path.join(chunk_dir, name), "rt") as fh:
            recs = list(C.iter_records(fh))
        assert recs, "{} is empty".format(name)
        assert all(complete for _lines, complete in recs), name
        assert recs[0][0][0][:1] == "M", "{} does not start with M".format(name)
        assert recs[-1][0][-1].rstrip("\n") == "E", \
            "{} does not end with E".format(name)
    assert found_any, "no chunks written to {}".format(chunk_dir)


def test_serial_and_parallel_agree(tmp_path, run_script):
    src = tmp_path / "src"
    ids = db2gen.build_tree(str(src), FILES)

    serial_out = tmp_path / "serial"
    run_script("01_chunk.py", "-i", src, "-o", serial_out, "-n", 7,
               "--progress-interval", 0)

    work = tmp_path / "work"
    par_out = tmp_path / "par"
    shards, bins = 3, 6
    run_script("make_shards.py", "-i", src, "-o", work, "-S", shards)
    for s in range(shards):
        run_script("p1_collect.py", "-w", work, "-s", s, "-N", bins)
    for b in range(bins):
        run_script("p2_assemble.py", "-w", work, "-b", b, "-o", par_out)
    run_script("p3_finalize.py", "-w", work, "-o", par_out)

    _assert_whole_records(str(serial_out))
    _assert_whole_records(str(par_out))
    assert _multiset(str(serial_out)) == ids
    assert _multiset(str(par_out)) == ids
