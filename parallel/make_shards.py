#!/usr/bin/env python3
"""
Plan A -- split the input tree into size-balanced shards for a SLURM array.

Sharding by file COUNT would leave tasks wildly uneven, because tranche files
differ hugely in size. This packs shards by cumulative bytes instead, so array
tasks finish at roughly the same time.

    python3 make_shards.py -i /path/to/3D -o work --shards 200

With --target-per-bin it also estimates the molecule count (from a small sample)
and writes the chunk count to work/bins.txt, so p1_collect / submit.slurm can
size the run without you guessing:

    python3 make_shards.py -i /path/to/3D -o work --shards 200 --target-per-bin 50000

Writes work/shards/shard_00000.txt ..., work/plan.tsv, and (with
--target-per-bin) work/bins.txt
"""

import argparse
import os
import sys
import zlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root holds db2common.py
from db2common import (find_inputs, dedupe_by_inode, iter_records_bytes,
                       iter_sources, INPUT_SUFFIXES)


def count_mols(path):
    """Count E-terminated records in one source (archive or bare file).

    None if the whole source is unreadable (e.g. a corrupt archive); a single
    bad member is skipped so the rest of the archive still counts.
    """
    errs = []
    n = 0
    for _label, fh in iter_sources([path], lambda p, e: errs.append(e)):
        try:
            for _rec, complete in iter_records_bytes(fh):
                if complete:
                    n += 1
        except (OSError, EOFError, zlib.error):
            continue
    if errs and n == 0:
        return None
    return n


def estimate_bins(sized, total_bytes, target_per_bin, sample, work_dir):
    """Estimate the molecule count from a sample and derive the bin count.

    Samples files evenly across the size-sorted list (so the calibration spans
    small and large molecules), counts their records, and scales
    molecules-per-compressed-byte up to the whole library. Writes work/bins.txt
    and returns the bin count.
    """
    k = min(sample, len(sized))
    step = len(sized) / float(k)
    idx = sorted(set(int(i * step) for i in range(k)))
    smol = sbytes = used = 0
    for j in idx:
        nbytes, path = sized[j]
        c = count_mols(path)
        if c is None:
            continue
        smol += c
        sbytes += nbytes
        used += 1
    if used == 0 or sbytes == 0:
        sys.exit("could not read any sample files to estimate the count; "
                 "pass an explicit bin count to p1_collect with -N")

    mols_per_byte = smol / float(sbytes)
    est_total = int(round(mols_per_byte * total_bytes))
    bins = max(1, -(-est_total // target_per_bin))   # ceil(est / target)

    with open(os.path.join(work_dir, "bins.txt"), "w") as fh:
        fh.write("{}\n".format(bins))

    sys.stderr.write("\nmolecule-count estimate ({} sampled file(s)):\n".format(used))
    sys.stderr.write("  {:,} mol across {:.2f} GB sampled\n".format(
        smol, sbytes / 1e9))
    sys.stderr.write("  estimated total: {:,} molecules\n".format(est_total))
    sys.stderr.write("  target {:,}/chunk -> BINS = {}  (written to {}/bins.txt)\n".format(
        target_per_bin, bins, work_dir))
    sys.stderr.write("  NOTE: an estimate (~+/-10-20%). Override with p1_collect "
                     "-N, or use an exact count pass if you need a hard ceiling.\n")
    return bins


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input-dir", default=".")
    ap.add_argument("-o", "--work-dir", required=True)
    ap.add_argument("-S", "--shards", type=int, required=True)
    ap.add_argument("--suffix", action="append", metavar="EXT",
                    help="filename suffix to collect; repeatable (default .db2.gz)")
    ap.add_argument("--file-list", help="use this list instead of walking the tree")
    ap.add_argument("--keep-duplicate-paths", action="store_true")
    ap.add_argument("--target-per-bin", type=int, default=0, metavar="N",
                    help="target molecules per finished chunk. Estimates the "
                         "total count from a sample and writes the derived bin "
                         "count to work/bins.txt (0 = off, the default)")
    ap.add_argument("--sample", type=int, default=300, metavar="K",
                    help="files to parse when estimating the count for "
                         "--target-per-bin (default 300)")
    args = ap.parse_args()

    if args.target_per_bin < 0:
        sys.exit("--target-per-bin must be >= 0")

    if args.file_list:
        with open(args.file_list) as fh:
            inputs = [ln.strip() for ln in fh if ln.strip()]
    else:
        sys.stderr.write("walking {}...\n".format(args.input_dir))
        inputs = find_inputs(args.input_dir, tuple(args.suffix or INPUT_SUFFIXES))
    if not inputs:
        sys.exit("no input files found")
    sys.stderr.write("found {} files\n".format(len(inputs)))

    if not args.keep_duplicate_paths:
        inputs, skipped = dedupe_by_inode(inputs)
        if skipped:
            sys.stderr.write("skipped {} symlink/hardlink duplicate(s)\n".format(
                len(skipped)))

    sized = []
    total = 0
    for p in inputs:
        try:
            n = os.path.getsize(p)
        except OSError:
            n = 0
        sized.append((n, p))
        total += n
    sys.stderr.write("total input: {:.1f} GB\n".format(total / 1e9))

    # Longest-processing-time first: biggest files into the lightest shard.
    sized.sort(reverse=True)
    loads = [0] * args.shards
    members = [[] for _ in range(args.shards)]
    import heapq
    heap = [(0, i) for i in range(args.shards)]
    heapq.heapify(heap)
    for n, p in sized:
        load, i = heapq.heappop(heap)
        members[i].append(p)
        loads[i] = load + n
        heapq.heappush(heap, (loads[i], i))

    sdir = os.path.join(args.work_dir, "shards")
    os.makedirs(sdir, exist_ok=True)
    with open(os.path.join(args.work_dir, "plan.tsv"), "w") as plan:
        plan.write("shard\tfiles\tbytes\n")
        for i in range(args.shards):
            # restore path order within a shard for reproducibility
            members[i].sort()
            with open(os.path.join(sdir, "shard_{:05d}.txt".format(i)), "w") as fh:
                fh.write("\n".join(members[i]) + ("\n" if members[i] else ""))
            plan.write("{}\t{}\t{}\n".format(i, len(members[i]), loads[i]))

    live = [l for l in loads if l]
    sys.stderr.write("wrote {} shard lists to {}\n".format(args.shards, sdir))
    if live:
        sys.stderr.write("  shard bytes: min {:.2f} GB  max {:.2f} GB  "
                         "spread {:.2f}x\n".format(
                             min(live) / 1e9, max(live) / 1e9,
                             max(live) / min(live)))
    if len(live) < args.shards:
        sys.stderr.write("  NOTE: {} shard(s) are empty -- you have fewer files "
                         "than shards\n".format(args.shards - len(live)))

    if args.target_per_bin:
        estimate_bins(sized, total, args.target_per_bin, args.sample,
                      args.work_dir)

    sys.stderr.write("\nSLURM array range: 0-{}\n".format(args.shards - 1))


if __name__ == "__main__":
    main()
