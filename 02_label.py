#!/usr/bin/env python3
"""
Stage 2 -- decide new ids from the manifest. Touches no db2 data.

Reads manifest.tsv from stage 1 and writes labels.tsv, containing only the
records whose id needs to change. Because this stage sees the whole id
population as plain text, it can label just the genuine duplicates and leave
unique molecules with their catalog id untouched.

    python3 02_label.py -m chunks/manifest.tsv -o labels.tsv --mode duplicates

Modes:
    duplicates  only ids occurring 2+ times are suffixed (default)
    occurrence  every id gets the nth-seen suffix
    serial      every id gets a global running index

Output: chunk, chunk_idx, original_id, new_id
"""

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from db2common import read_tsv, Progress, MANIFEST_HEADER, LABELS_HEADER


def load_dup_ids(path):
    """Read a precomputed duplicate-id list (one id per line)."""
    with open(path) as fh:
        return set(ln.strip() for ln in fh if ln.strip())


def count_dup_ids(manifest):
    """First pass: which ids occur more than once.

    Holds one entry per unique id. For very large trees pass --dup-ids
    instead, computed with sort/uniq on disk.
    """
    seen = set()
    dups = set()
    for chunk, idx, src, sidx, mol_id in read_tsv(manifest, MANIFEST_HEADER):
        if mol_id in seen:
            dups.add(mol_id)
        else:
            seen.add(mol_id)
    return dups


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-m", "--manifest", required=True)
    ap.add_argument("-o", "--output", required=True, help="labels.tsv to write")
    ap.add_argument("--mode", choices=["duplicates", "occurrence", "serial"],
                    default="duplicates")
    ap.add_argument("--dup-ids", metavar="FILE",
                    help="precomputed duplicate ids (one per line); avoids "
                         "holding every unique id in memory")
    ap.add_argument("--keep-first", action="store_true",
                    help="in duplicates mode, leave the first occurrence "
                         "unsuffixed and number only the later ones")
    ap.add_argument("--sep", default="_", help="suffix separator (default: _)")
    ap.add_argument("--width", type=int, default=0,
                    help="zero-pad suffix to this width (0 = no padding)")
    ap.add_argument("--progress-interval", type=float, default=30.0, metavar="SEC",
                    help="seconds between progress lines; 0 disables (default: 30)")
    ap.add_argument("--max-id-len", type=int, default=0,
                    help="warn if a new id exceeds this length (0 = no check)")
    args = ap.parse_args()

    dups = None
    if args.mode == "duplicates":
        if args.dup_ids:
            dups = load_dup_ids(args.dup_ids)
            sys.stderr.write("loaded {} duplicate id(s)\n".format(len(dups)))
        else:
            sys.stderr.write("pass 1: finding duplicate ids...\n")
            dups = count_dup_ids(args.manifest)
            sys.stderr.write("  {} id(s) occur more than once\n".format(len(dups)))

    occ = defaultdict(int)
    prog = Progress(args.progress_interval, unit="records")
    written = total = no_id = 0
    overlong = []

    def suffix(n):
        return "{:0{w}d}".format(n, w=args.width) if args.width else str(n)

    with open(args.output, "w") as out:
        out.write("\t".join(LABELS_HEADER) + "\n")
        for chunk, idx, src, sidx, mol_id in read_tsv(args.manifest, MANIFEST_HEADER):
            total += 1
            if total % 100000 == 0:
                prog.tick(total)
            if mol_id == "NO_ID":
                no_id += 1
                continue

            if args.mode == "serial":
                new_id = mol_id + args.sep + suffix(total)
            elif args.mode == "occurrence":
                occ[mol_id] += 1
                new_id = mol_id + args.sep + suffix(occ[mol_id])
            else:  # duplicates
                if mol_id not in dups:
                    continue
                occ[mol_id] += 1
                n = occ[mol_id]
                if args.keep_first and n == 1:
                    continue
                new_id = mol_id + args.sep + suffix(n)

            if args.max_id_len and len(new_id) > args.max_id_len and len(overlong) < 20:
                overlong.append(new_id)
            out.write("{}\t{}\t{}\t{}\n".format(chunk, idx, mol_id, new_id))
            written += 1

    sys.stderr.write("\nmanifest records: {}\n".format(total))
    sys.stderr.write("records to relabel: {} ({:.1f}%)\n".format(
        written, 100.0 * written / total if total else 0.0))
    sys.stderr.write("labels written to {}\n".format(args.output))
    if no_id:
        sys.stderr.write("skipped {} record(s) with no ZINC id\n".format(no_id))
    if overlong:
        sys.stderr.write("\nWARNING: new ids exceed --max-id-len, e.g.:\n")
        for i in overlong[:5]:
            sys.stderr.write("  {} ({} chars)\n".format(i, len(i)))


if __name__ == "__main__":
    main()
