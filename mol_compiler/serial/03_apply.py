#!/usr/bin/env python3
"""
Stage 3 -- apply labels.tsv to the chunk files.

Records are located by their position within a chunk, not by id, so
duplicate ids are never ambiguous. Before rewriting a record this stage
verifies the id at that position matches what labels.tsv expects; a mismatch
means the labels are stale (stage 1 re-run with a different chunk size, say)
and the run aborts rather than corrupting data.

    python3 03_apply.py -c chunks -L labels.tsv -o chunks_labelled

One chunk at a time, for a job array:

    python3 03_apply.py -c chunks -L labels.tsv -o chunks_labelled \\
        --only chunk_00007.db2.gz
"""

import argparse
import os
import sys
import zlib
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root holds db2common.py
from db2common import (iter_records, extract_id, relabel, open_gz_text,
                       read_tsv, id_match, LABELS_HEADER)


def load_labels(path, only=None):
    """labels.tsv -> {chunk: {chunk_idx: (orig, new)}}"""
    by_chunk = defaultdict(dict)
    for chunk, idx, orig, new in read_tsv(path, LABELS_HEADER):
        if only and chunk != only:
            continue
        by_chunk[chunk][int(idx)] = (orig, new)
    return by_chunk


def process_chunk(src, dst, edits, compresslevel):
    """Rewrite one chunk applying edits. Returns (records, applied, mismatches)."""
    applied = 0
    records = 0
    mismatches = []
    with open_gz_text(src, "rt") as fin, open_gz_text(dst, "wt", compresslevel) as fout:
        for lines, complete in iter_records(fin):
            if not complete:
                mismatches.append((records, "TRUNCATED", "record has no E terminator"))
                fout.writelines(lines)
                break
            edit = edits.get(records)
            if edit is not None:
                orig, new = edit
                found = extract_id(lines)
                if not id_match(found, orig):
                    mismatches.append((records, orig, found))
                else:
                    # Relabel the id that is actually in the record. For ZINC22
                    # the header carries a truncated form (a suffix of `orig`),
                    # so rewrite `found` -> found + the same suffix rather than
                    # `orig` -> `new` (which wouldn't be found in the record).
                    suffix = new[len(orig):] if new.startswith(orig) else None
                    if suffix is not None and found != orig:
                        lines = relabel(lines, found, found + suffix)
                    else:
                        lines = relabel(lines, orig, new)
                    applied += 1
            fout.writelines(lines)
            records += 1
    return records, applied, mismatches


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--chunk-dir", required=True)
    ap.add_argument("-L", "--labels", required=True)
    ap.add_argument("-o", "--output-dir", required=True)
    ap.add_argument("--only", help="process just this chunk filename")
    ap.add_argument("-l", "--compresslevel", type=int, default=6)
    ap.add_argument("--copy-unedited", action="store_true",
                    help="also copy chunks that have no edits (default: skip "
                         "them, leaving originals in place)")
    args = ap.parse_args()

    by_chunk = load_labels(args.labels, args.only)
    if args.only and not by_chunk:
        sys.stderr.write("no edits for {} -- nothing to do\n".format(args.only))
        return

    chunks = sorted(f for f in os.listdir(args.chunk_dir) if f.endswith(".db2.gz"))
    if args.only:
        chunks = [c for c in chunks if c == args.only]
        if not chunks:
            sys.exit("chunk not found in {}: {}".format(args.chunk_dir, args.only))

    targets = [c for c in chunks if c in by_chunk or args.copy_unedited]
    if not targets:
        sys.stderr.write("no chunks need edits\n")
        return

    os.makedirs(args.output_dir, exist_ok=True)
    sys.stderr.write("{} chunk(s) to process\n".format(len(targets)))

    total_records = total_applied = 0
    all_mismatches = []

    for name in targets:
        src = os.path.join(args.chunk_dir, name)
        dst = os.path.join(args.output_dir, name)
        tmp = dst + ".tmp"
        edits = by_chunk.get(name, {})
        try:
            recs, applied, mism = process_chunk(src, tmp, edits, args.compresslevel)
        except (OSError, EOFError, zlib.error) as exc:
            if os.path.exists(tmp):
                os.unlink(tmp)
            sys.exit("failed reading {}: {}".format(src, exc))

        missing = set(edits) - set(range(recs))
        if missing:
            mism.extend((i, edits[i][0], "PAST END OF CHUNK") for i in sorted(missing))

        if mism:
            os.unlink(tmp)
            all_mismatches.append((name, mism))
            sys.stderr.write("  {}: {} MISMATCH(es), not written\n".format(
                name, len(mism)))
            continue

        os.replace(tmp, dst)
        total_records += recs
        total_applied += applied
        sys.stderr.write("  {}: {} records, {} relabelled\n".format(
            name, recs, applied))

    expected = sum(len(v) for v in by_chunk.values())
    sys.stderr.write("\nrecords read:  {}\n".format(total_records))
    sys.stderr.write("relabelled:    {} of {} expected\n".format(
        total_applied, expected))
    sys.stderr.write("output: {}\n".format(args.output_dir))

    if all_mismatches:
        sys.stderr.write("\nERROR: labels do not match the chunk contents.\n"
                         "Stage 1 was probably re-run with different settings "
                         "after stage 2.\n")
        for name, mism in all_mismatches[:5]:
            sys.stderr.write("  {}:\n".format(name))
            for idx, want, got in mism[:5]:
                sys.stderr.write("    record {}: expected {}, found {}\n".format(
                    idx, want, got))
        sys.exit(1)

    if total_applied != expected:
        sys.exit("ERROR: applied {} edits but labels.tsv had {}".format(
            total_applied, expected))


if __name__ == "__main__":
    main()
