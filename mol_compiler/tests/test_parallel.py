"""#2 — the parallel global-index arithmetic (the crown-jewel test).

A record's position inside a finished chunk is NOT its position inside the shard
that produced it: chunk = part(shard 0) ++ part(shard 1) ++ ... So p3 computes a
global index as (sum of every earlier shard's contribution to that bin) + local
index. The README flags this arithmetic as "the thing most likely to break under
changes." This test runs the whole parallel pipeline and checks every single
manifest row against the id actually sitting at that byte position in the
assembled chunk, then confirms nothing was lost or duplicated versus the source.
"""

import os
from collections import Counter

import db2common as C
import db2gen

# Several files of different sizes across two "tranche" dirs, so make_shards has
# real packing to do and each bin ends up fed by multiple shards.
FILES = {
    "small/a.db2.gz": [("ZINC{:08d}".format(i), (i % 3) + 1) for i in range(0, 6)],
    "small/b.db2.gz": [("ZINC{:08d}".format(i), (i % 3) + 1) for i in range(6, 13)],
    "big/c.db2.gz": [("ZINC{:08d}".format(i), (i % 4) + 2) for i in range(13, 21)],
    "big/d.db2.gz": [("ZINC{:08d}".format(i), (i % 4) + 2) for i in range(21, 30)],
    "big/e.db2.gz": [("ZINC{:08d}".format(i), 1) for i in range(30, 40)],
}


def _run_parallel(src, work, out, run_script, shards, bins, mode="stride",
                  weight="count"):
    run_script("make_shards.py", "-i", src, "-o", work, "-S", shards)
    for s in range(shards):
        run_script("p1_collect.py", "-w", work, "-s", s, "-N", bins,
                   "--mode", mode, "--weight", weight)
    for b in range(bins):
        run_script("p2_assemble.py", "-w", work, "-b", b, "-o", out)
    run_script("p3_finalize.py", "-w", work, "-o", out)


def test_manifest_positions_match_assembled_bytes(tmp_path, run_script):
    src, work, out = tmp_path / "src", tmp_path / "work", tmp_path / "out"
    ids = db2gen.build_tree(str(src), FILES)

    _run_parallel(src, work, out, run_script, shards=3, bins=4)

    # Read every chunk's ids once, keyed by chunk filename.
    by_chunk = {}
    for name in os.listdir(str(out)):
        if name.endswith(".db2.gz"):
            by_chunk[name] = db2gen.chunk_ids(os.path.join(str(out), name))

    # Every manifest row must name the id that is byte-actually at that position.
    seen = Counter()
    rows = 0
    for chunk, idx, _src, _sidx, orig in C.read_tsv(
            os.path.join(str(out), "manifest.tsv"), C.MANIFEST_HEADER):
        assert by_chunk[chunk][int(idx)] == orig, \
            "manifest says {} at {}[{}] but bytes say {}".format(
                orig, chunk, idx, by_chunk[chunk][int(idx)])
        seen[orig] += 1
        rows += 1

    # Nothing lost, nothing duplicated, versus the source tree.
    assert seen == ids
    assert rows == sum(ids.values())


def test_greedy_mode_also_round_trips(tmp_path, run_script):
    # The greedy (cost-weighted) path uses the same index arithmetic; exercise
    # it too so a change there can't slip past.
    src, work, out = tmp_path / "src", tmp_path / "work", tmp_path / "out"
    ids = db2gen.build_tree(str(src), FILES)

    _run_parallel(src, work, out, run_script, shards=3, bins=4,
                  mode="greedy", weight="lines:C")

    by_chunk = {}
    for name in os.listdir(str(out)):
        if name.endswith(".db2.gz"):
            by_chunk[name] = db2gen.chunk_ids(os.path.join(str(out), name))

    seen = Counter()
    for chunk, idx, _src, _sidx, orig in C.read_tsv(
            os.path.join(str(out), "manifest.tsv"), C.MANIFEST_HEADER):
        assert by_chunk[chunk][int(idx)] == orig
        seen[orig] += 1
    assert seen == ids
