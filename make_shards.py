#!/usr/bin/env python3
"""
Plan A -- split the input tree into size-balanced shards for a SLURM array.

Sharding by file COUNT would leave tasks wildly uneven, because tranche files
differ hugely in size. This packs shards by cumulative bytes instead, so array
tasks finish at roughly the same time.

    python3 make_shards.py -i /path/to/3D -o work --shards 200

Writes work/shards/shard_00000.txt ... and work/plan.tsv
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from db2common import find_inputs, dedupe_by_inode


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
    args = ap.parse_args()

    if args.file_list:
        with open(args.file_list) as fh:
            inputs = [ln.strip() for ln in fh if ln.strip()]
    else:
        sys.stderr.write("walking {}...\n".format(args.input_dir))
        inputs = find_inputs(args.input_dir, tuple(args.suffix or [".db2.gz"]))
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
    sys.stderr.write("\nSLURM array range: 0-{}\n".format(args.shards - 1))


if __name__ == "__main__":
    main()
