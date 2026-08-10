"""#5 — reindexing: keeping the manifest aligned after QC removes molecules.

QC tooling (db2tool subset / repair --strict) drops molecules from chunks, which
renumbers every molecule after the first deletion. 05_reindex applies the map
that tooling writes, and — like 03_apply — verifies the ids line up before
rewriting, because a stale map corrupts silently rather than erroring.
"""

import os

import db2common as C
import db2gen

MAP_HEADER = "chunk\torig_idx\tnew_idx\tname\tstatus\n"


def _chunked(tmp_path, run_script, n=6):
    """A one-chunk library plus its manifest."""
    src = tmp_path / "src"
    db2gen.build_tree(str(src), {
        "a.db2.gz": [("ZINC{:08d}".format(i), 1) for i in range(n)]})
    chunks = tmp_path / "chunks"
    run_script("01_chunk.py", "-i", src, "-o", chunks, "-n", 50,
               "--progress-interval", 0)
    return chunks, os.path.join(str(chunks), "manifest.tsv")


def _read_manifest(path):
    return list(C.read_tsv(path, C.MANIFEST_HEADER))


def _write_map(path, chunk, rows):
    """rows: (orig_idx, new_idx, name, status)"""
    with open(path, "w") as fh:
        fh.write(MAP_HEADER)
        for orig, new, name, status in rows:
            fh.write("{}\t{}\t{}\t{}\t{}\n".format(chunk, orig, new, name, status))


def test_drops_removed_rows_and_renumbers(tmp_path, run_script):
    chunks, manifest = _chunked(tmp_path, run_script)
    rows = _read_manifest(manifest)
    chunk = rows[0][0]

    # drop molecules 1 and 3; the rest renumber contiguously
    keep = [i for i in range(len(rows)) if i not in (1, 3)]
    mapping = [(i, -1, rows[i][4], "all_broken") if i in (1, 3)
               else (i, keep.index(i), rows[i][4], "kept")
               for i in range(len(rows))]
    mpath = tmp_path / "remap.tsv"
    _write_map(str(mpath), chunk, mapping)

    out = tmp_path / "out"
    run_script("05_reindex.py", "--map", mpath, "-m", manifest, "-o", out)

    new = _read_manifest(os.path.join(str(out), "manifest.tsv"))
    assert len(new) == len(rows) - 2
    # contiguous from 0, and each row still carries its own provenance
    assert [int(r[1]) for r in new] == list(range(len(new)))
    assert [r[4] for r in new] == [rows[i][4] for i in keep]


def test_rejects_a_stale_map(tmp_path, run_script):
    chunks, manifest = _chunked(tmp_path, run_script)
    rows = _read_manifest(manifest)
    chunk = rows[0][0]

    mapping = [(i, i, rows[i][4], "kept") for i in range(len(rows))]
    mapping[2] = (2, 2, "ZINC99999999", "kept")          # wrong id at position 2
    mpath = tmp_path / "stale.tsv"
    _write_map(str(mpath), chunk, mapping)

    out = tmp_path / "out"
    proc = run_script("05_reindex.py", "--map", mpath, "-m", manifest,
                      "-o", out, check=False)
    assert proc.returncode != 0
    assert "does not line up" in proc.stderr


def test_chunk_absent_from_map_passes_through(tmp_path, run_script):
    """Reindexing one chunk of a set must not disturb the others."""
    chunks, manifest = _chunked(tmp_path, run_script)
    rows = _read_manifest(manifest)

    mpath = tmp_path / "other.tsv"
    _write_map(str(mpath), "chunk_99999.db2.gz", [(0, -1, "ZINC00000000", "x")])

    out = tmp_path / "out"
    run_script("05_reindex.py", "--map", mpath, "-m", manifest, "-o", out)

    new = _read_manifest(os.path.join(str(out), "manifest.tsv"))
    assert new == rows          # untouched


def test_tolerates_truncated_ids_in_the_map(tmp_path, run_script):
    """db2tool reports the db2 header's 16-character id field.

    A ZINC22 base id fills that field exactly, so most tranches write it with
    the protomer suffix omitted, while zinc-22a writes the full id and loses two
    characters off the LEFT.  The map's name is therefore either the manifest id
    minus its suffix or its last 16 characters; both must be accepted, or the
    guard rejects most of a real library.
    """
    chunks, manifest = _chunked(tmp_path, run_script)
    rows = _read_manifest(manifest)
    chunk = rows[0][0]

    mapping = []
    for i, r in enumerate(rows):
        full = r[4]
        trunc = full[:-1] if i % 2 == 0 else full[1:]   # right-cut / left-cut
        mapping.append((i, i, trunc, "kept"))
    mpath = tmp_path / "trunc.tsv"
    _write_map(str(mpath), chunk, mapping)

    out = tmp_path / "out"
    run_script("05_reindex.py", "--map", mpath, "-m", manifest, "-o", out)
    assert _read_manifest(os.path.join(str(out), "manifest.tsv")) == rows


def test_id_match_accepts_both_truncation_directions():
    assert C.id_match("NC550000002ttm.0", "ZINC550000002ttm.0")   # left-cut
    assert C.id_match("ZINCa50000001eSu", "ZINCa50000001eSu.0")   # right-cut
    assert C.id_match("ZINC00000005", "ZINC00000005")
    assert not C.id_match("ZINC00000006", "ZINC00000005")
    assert not C.id_match("NO_ID", "ZINC00000005")
