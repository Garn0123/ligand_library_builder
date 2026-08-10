#!/usr/bin/env python3
"""
Stage 5 -- reindex a manifest after molecules were removed from the chunks.

QC tooling (db2tool subset, db2tool repair --strict) can drop molecules from a
chunk. That renumbers every molecule after the first deletion, so the manifest's
(chunk, chunk_idx) keys stop pointing at the right rows. Nothing errors when
that happens: the wrong provenance simply attaches to the wrong molecule.

This stage applies the map db2tool writes:

    chunk                 orig_idx  new_idx  name              status
    chunk_00001.db2.gz    0         0        ZINCa50000001eSk  kept
    chunk_00001.db2.gz    48        -1       ZINCa50000001eIH  all_broken

Rows whose new_idx is -1 were removed and are dropped from the manifest; the
rest get chunk_idx = new_idx. Provenance columns are carried through untouched.

    python3 05_reindex.py --map remap.tsv -m chunks/manifest.tsv -o clean

Maps from several chunks may be concatenated (keep one header) or passed as
repeated --map arguments. A chunk absent from the map is passed through
unchanged, so reindexing one chunk of a set is safe.

Like stage 3, this verifies before it rewrites: the map's id must match the
manifest's at the same (chunk, orig_idx), tolerating the db2 field truncation.
A mismatch means the map is stale -- generated against different chunks -- and
the run aborts rather than silently misattributing every downstream result.
"""

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root holds db2common.py
from db2common import read_tsv, id_match, MANIFEST_HEADER

MAP_HEADER = ["chunk", "orig_idx", "new_idx", "name", "status"]


def load_maps(paths):
    """map files -> {chunk: {orig_idx: (new_idx, name)}}"""
    by_chunk = defaultdict(dict)
    for path in paths:
        for chunk, orig, new, name, _status in read_tsv(path, MAP_HEADER):
            by_chunk[chunk][int(orig)] = (int(new), name)
    return by_chunk


def reindex(manifest, maps, out_fh, strict_ids=True):
    """Stream the manifest, applying the maps. Returns (kept, dropped, passed)."""
    kept = dropped = passed = 0
    mismatches = []

    out_fh.write("\t".join(MANIFEST_HEADER) + "\n")

    for row in read_tsv(manifest, MANIFEST_HEADER):
        chunk, chunk_idx = row[0], int(row[1])
        edits = maps.get(chunk)

        if edits is None:
            # chunk was not part of this QC pass; carry it through untouched
            passed += 1
            out_fh.write("\t".join(row) + "\n")
            continue

        entry = edits.get(chunk_idx)
        if entry is None:
            # The map covers this chunk but not this row -- it describes a
            # different set of molecules than the manifest does.
            mismatches.append((chunk, chunk_idx, row[4], "(absent from map)"))
            continue

        new_idx, name = entry
        # Stage 3's guard, for the same reason: ids are position-addressed, so a
        # stale map corrupts silently.  id_match tolerates the truncated form
        # the db2 header carries.
        if strict_ids and not id_match(name, row[4]):
            mismatches.append((chunk, chunk_idx, row[4], name))
            continue

        if new_idx < 0:
            dropped += 1
            continue

        row[1] = str(new_idx)
        kept += 1
        out_fh.write("\t".join(row) + "\n")

    if mismatches:
        sys.stderr.write(
            "\nERROR: the map does not line up with the manifest.\n"
            "       This means they were generated from different chunks; applying\n"
            "       it would attach the wrong provenance to every later molecule.\n")
        for chunk, idx, want, got in mismatches[:10]:
            sys.stderr.write("  {} idx {}: manifest id {!r}, map says {!r}\n"
                             .format(chunk, idx, want, got))
        if len(mismatches) > 10:
            sys.stderr.write("  ... and {} more\n".format(len(mismatches) - 10))
        raise SystemExit(1)

    return kept, dropped, passed


def main():
    ap = argparse.ArgumentParser(
        description="Reindex a manifest after QC removed molecules from chunks.")
    ap.add_argument("--map", action="append", required=True, metavar="FILE",
                    help="db2tool map tsv; repeat for several chunks")
    ap.add_argument("-m", "--manifest", required=True,
                    help="manifest.tsv to reindex")
    ap.add_argument("-o", "--output-dir", required=True,
                    help="directory to write the new manifest.tsv into")
    ap.add_argument("--no-verify-ids", action="store_true",
                    help="skip the id cross-check (not recommended; the check "
                         "is what catches a stale map)")
    args = ap.parse_args()

    maps = load_maps(args.map)
    if not maps:
        raise SystemExit("no rows in the map file(s)")

    os.makedirs(args.output_dir, exist_ok=True)
    dst = os.path.join(args.output_dir, "manifest.tsv")
    if os.path.abspath(dst) == os.path.abspath(args.manifest):
        raise SystemExit("refusing to overwrite the input manifest; "
                         "choose a different -o")

    with open(dst, "w") as fh:
        kept, dropped, passed = reindex(args.manifest, maps, fh,
                                        strict_ids=not args.no_verify_ids)

    sys.stderr.write(
        "reindexed {} row(s), dropped {}, passed through {} (chunks not in the map)\n"
        "-> {}\n".format(kept, dropped, passed, dst))


if __name__ == "__main__":
    main()
