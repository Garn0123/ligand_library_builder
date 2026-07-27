#!/usr/bin/env python3
"""
Phase 1 -- one SLURM array task per shard. Embarrassingly parallel.

Each task reads its own files and writes a PART of every global bin. Because
every shard contributes an even slice of its own records to every bin, the
bins come out globally balanced without any cross-task communication -- even
though each shard may hold only one region of the size distribution.

    python3 p1_collect.py -w work -s $SLURM_ARRAY_TASK_ID -N 170

Writes:
    work/parts/bin_00001/part_00007.db2.gz
    work/manifests/shard_00007.tsv    bin, local_idx, source_file, source_idx, id
    work/counts/shard_00007.tsv       bin, records, weight
"""

import argparse
import gzip
import heapq
import os
import resource
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root holds db2common.py
from db2common import (iter_records_bytes, extract_id_bytes, make_weight_fn_bytes,
                       is_gzip)


def check_fd_limit(n):
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    need = n + 32
    if soft >= need:
        return
    if hard >= need:
        resource.setrlimit(resource.RLIMIT_NOFILE, (need, hard))
        sys.stderr.write("raised open-file limit to {}\n".format(need))
        return
    sys.exit("need {} open files for {} bins, hard limit is {}. Lower -N or "
             "raise ulimit -n.".format(need, n, hard))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-w", "--work-dir", required=True)
    ap.add_argument("-s", "--shard", type=int, required=True)
    ap.add_argument("-N", "--bins", type=int, required=True,
                    help="number of global output chunks")
    ap.add_argument("--mode", choices=["stride", "greedy"], default="stride")
    ap.add_argument("--weight", default="count", metavar="SPEC",
                    help="cost metric for greedy: count, bytes, lines:X")
    ap.add_argument("-l", "--compresslevel", type=int, default=1,
                    help="gzip level for parts. These are intermediates that "
                         "get concatenated, so 1 is usually right (default: 1)")
    ap.add_argument("--progress-interval", type=float, default=60.0)
    args = ap.parse_args()

    check_fd_limit(args.bins)

    shard_list = os.path.join(args.work_dir, "shards",
                              "shard_{:05d}.txt".format(args.shard))
    if not os.path.exists(shard_list):
        sys.exit("no shard list: {}".format(shard_list))
    with open(shard_list) as fh:
        files = [ln.strip() for ln in fh if ln.strip()]

    tag = "{:05d}".format(args.shard)
    sys.stderr.write("shard {} : {} files, {} bins\n".format(
        tag, len(files), args.bins))
    if not files:
        # still emit empty bookkeeping so phase 3 sees every shard
        for sub in ("manifests", "counts"):
            os.makedirs(os.path.join(args.work_dir, sub), exist_ok=True)
        open(os.path.join(args.work_dir, "manifests",
                          "shard_{}.tsv".format(tag)), "w").close()
        with open(os.path.join(args.work_dir, "counts",
                               "shard_{}.tsv".format(tag)), "w") as fh:
            for b in range(args.bins):
                fh.write("{}\t0\t0\n".format(b))
        return

    weight_fn = make_weight_fn_bytes(args.weight)
    for sub in ("parts", "manifests", "counts"):
        os.makedirs(os.path.join(args.work_dir, sub), exist_ok=True)

    writers = []
    for b in range(args.bins):
        d = os.path.join(args.work_dir, "parts", "bin_{:05d}".format(b + 1))
        os.makedirs(d, exist_ok=True)
        writers.append(gzip.open(os.path.join(d, "part_{}.db2.gz".format(tag)),
                                 "wb", compresslevel=args.compresslevel))

    counts = [0] * args.bins
    weights = [0.0] * args.bins
    heap = [(0.0, i) for i in range(args.bins)]
    heapq.heapify(heap)

    man_path = os.path.join(args.work_dir, "manifests", "shard_{}.tsv".format(tag))
    man = open(man_path, "w")

    seq = 0
    truncated = []
    no_id = 0
    start = time.time()
    last = start

    for fi, path in enumerate(files):
        opener = gzip.open if is_gzip(path) else open
        try:
            with opener(path, "rb") as fh:
                for src_idx, (rec, complete) in enumerate(iter_records_bytes(fh)):
                    if not complete:
                        truncated.append(path)
                        break
                    mol_id = extract_id_bytes(rec)
                    if mol_id == "NO_ID":
                        no_id += 1
                    if args.mode == "stride":
                        # offset by shard so bin 0 is not always fed the
                        # first (smallest) record of every shard
                        b = (seq + args.shard) % args.bins
                    else:
                        _w, b = heapq.heappop(heap)
                    w = weight_fn(rec)
                    writers[b].write(rec)
                    man.write("{}\t{}\t{}\t{}\t{}\n".format(
                        b, counts[b], path, src_idx, mol_id))
                    counts[b] += 1
                    weights[b] += w
                    if args.mode == "greedy":
                        heapq.heappush(heap, (weights[b], b))
                    seq += 1
        except (OSError, EOFError) as exc:
            sys.stderr.write("  UNREADABLE {}: {}\n".format(path, exc))
            continue

        now = time.time()
        if args.progress_interval > 0 and now - last >= args.progress_interval:
            last = now
            el = now - start
            sys.stderr.write("  {}/{} files  {} molecules  {:,.0f} mol/s\n".format(
                fi + 1, len(files), seq, seq / el if el else 0))
            sys.stderr.flush()

    for w in writers:
        w.close()
    man.close()

    with open(os.path.join(args.work_dir, "counts",
                           "shard_{}.tsv".format(tag)), "w") as fh:
        for b in range(args.bins):
            fh.write("{}\t{}\t{:.0f}\n".format(b, counts[b], weights[b]))

    el = time.time() - start
    sys.stderr.write("shard {} done: {} molecules in {:.0f}s ({:,.0f} mol/s)\n".format(
        tag, seq, el, seq / el if el else 0))
    if no_id:
        sys.stderr.write("  {} record(s) had no ZINC id\n".format(no_id))
    if truncated:
        sys.stderr.write("  WARNING: {} truncated file(s):\n".format(len(truncated)))
        for p in truncated[:10]:
            sys.stderr.write("    {}\n".format(p))
        sys.exit(1)


if __name__ == "__main__":
    main()
