"""Auto-sizing: make_shards --target-per-bin derives the chunk count.

The chunk count must be fixed before collect strides records, but the exact
molecule count isn't known until then. --target-per-bin estimates it from a
sample and writes work/bins.txt, which p1_collect reads when -N is omitted. On
a small fixture the sample covers every file, so the estimate (hence the bin
count) is exact and these assertions are deterministic.
"""

import os
from collections import Counter

import db2gen


def test_target_per_bin_writes_expected_bins(tmp_path, run_script):
    src = tmp_path / "src"
    ids = db2gen.build_tree(str(src), {
        "a.db2.gz": [("ZINC{:08d}".format(i), 2) for i in range(0, 10)],
        "b.db2.gz": [("ZINC{:08d}".format(i), 2) for i in range(10, 25)],
        "c.db2.gz": [("ZINC{:08d}".format(i), 2) for i in range(25, 40)],
    })
    total = sum(ids.values())          # 40

    work = tmp_path / "work"
    run_script("make_shards.py", "-i", src, "-o", work, "-S", 2,
               "--target-per-bin", 12)

    bins = int((work / "bins.txt").read_text().strip())
    assert bins == -(-total // 12)     # ceil(40/12) == 4


def test_autosized_pipeline_end_to_end(tmp_path, run_script):
    src = tmp_path / "src"
    ids = db2gen.build_tree(str(src), {
        "a.db2.gz": [("ZINC{:08d}".format(i), 2) for i in range(0, 20)],
        "b.db2.gz": [("ZINC{:08d}".format(i), 2) for i in range(20, 40)],
    })

    work = tmp_path / "work"
    out = tmp_path / "out"
    shards = 2
    run_script("make_shards.py", "-i", src, "-o", work, "-S", shards,
               "--target-per-bin", 10)          # 40 / 10 -> 4 bins
    bins = int((work / "bins.txt").read_text().strip())
    assert bins == 4

    for s in range(shards):
        run_script("p1_collect.py", "-w", work, "-s", s)   # NO -N: reads bins.txt
    for b in range(bins):
        run_script("p2_assemble.py", "-w", work, "-b", b, "-o", out)
    run_script("p3_finalize.py", "-w", work, "-o", out)

    chunks = [n for n in os.listdir(str(out)) if n.endswith(".db2.gz")]
    assert len(chunks) == bins

    got = Counter()
    for n in chunks:
        for zid in db2gen.chunk_ids(os.path.join(str(out), n)):
            got[zid] += 1
    assert got == ids                            # nothing lost or duplicated


def test_p1_errors_without_N_or_bins_txt(tmp_path, run_script):
    src = tmp_path / "src"
    db2gen.build_tree(str(src), {"a.db2.gz": [("ZINC00000001", 1)]})

    work = tmp_path / "work"
    run_script("make_shards.py", "-i", src, "-o", work, "-S", 1)   # no --target-per-bin

    res = run_script("p1_collect.py", "-w", work, "-s", 0, check=False)  # no -N
    assert res.returncode != 0
    assert "bins.txt" in res.stderr
