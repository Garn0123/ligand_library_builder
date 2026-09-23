#!/usr/bin/env python3
"""
Phase 1 -- one SLURM array task per shard. Embarrassingly parallel.

Each task reads its own sources and writes a PART of every global bin. Because
every shard contributes an even slice of its own records to every bin, the
bins come out globally balanced without any cross-task communication -- even
though each shard may hold only one region of the size distribution.

Sources may be old bare `.db2.gz` files or new ZINC22 `.db2.tgz` archives (read
member-by-member in place via iter_sources, never extracted). For a tar member
the manifest source_file is `<archive>::<member>` and the id comes from the
member filename (the full `ZINC...`); for a bare file it comes from the header.

    python3 p1_collect.py -w work -s $SLURM_ARRAY_TASK_ID -N 170
    python3 p1_collect.py -w work -s $SLURM_ARRAY_TASK_ID       # -N from work/bins.txt

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
import zlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root holds db2common.py
from db2common import (iter_records_bytes, extract_id_bytes, make_weight_fn_bytes,
                       iter_sources, id_from_name)


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
    ap.add_argument("-N", "--bins", type=int, default=None,
                    help="number of global output chunks. If omitted, read from "
                         "<work-dir>/bins.txt (written by make_shards "
                         "--target-per-bin)")
    ap.add_argument("--mode", choices=["stride", "greedy"], default="stride")
    ap.add_argument("--weight", default="count", metavar="SPEC",
                    help="cost metric for greedy: count, bytes, lines:X")
    ap.add_argument("-l", "--compresslevel", type=int, default=1,
                    help="gzip level for parts. These are intermediates that "
                         "get concatenated, so 1 is usually right (default: 1)")
    ap.add_argument("--skip-corrupt", action="store_true",
                    help="continue past corrupt/truncated input files instead of "
                         "failing the task: keep whatever reads cleanly, drop the "
                         "rest, and exit 0 so assembly can proceed. The dropped "
                         "files are still listed. Use to finish a run over a "
                         "library with a few known-bad files.")
    ap.add_argument("--progress-interval", type=float, default=60.0)
    args = ap.parse_args()

    if args.bins is None:
        bpath = os.path.join(args.work_dir, "bins.txt")
        if not os.path.exists(bpath):
            sys.exit("no -N/--bins given and no {}; run make_shards with "
                     "--target-per-bin, or pass -N.".format(bpath))
        with open(bpath) as fh:
            args.bins = int(fh.read().strip())
        sys.stderr.write("bins: {} (from {})\n".format(args.bins, bpath))
    if args.bins < 1:
        sys.exit("--bins must be >= 1")

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
    members = 0
    truncated = []
    corrupt = []
    no_id = 0
    start = time.time()
    last = start

    def on_error(p, e):
        # Container-level failure (corrupt archive, unopenable file). zlib.error
        # is NOT an OSError: a valid gzip header over a damaged deflate body
        # raises it deep in gzip.read and would otherwise crash the whole task.
        corrupt.append((p, str(e)))
        sys.stderr.write("  UNREADABLE {}: {}\n".format(p, e))

    for label, fh in iter_sources(files, on_error):
        members += 1
        try:
            for src_idx, (rec, complete) in enumerate(iter_records_bytes(fh)):
                if not complete:
                    truncated.append(label)
                    break
                # New ZINC22 tar members carry the full id in the member name;
                # the record header only has a truncated form. Old bare files
                # carry it in the header.
                if "::" in label:
                    mol_id = id_from_name(label)
                    if mol_id == "NO_ID":
                        mol_id = extract_id_bytes(rec)
                else:
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
                    b, counts[b], label, src_idx, mol_id))
                counts[b] += 1
                weights[b] += w
                if args.mode == "greedy":
                    heapq.heappush(heap, (weights[b], b))
                seq += 1
        except (OSError, EOFError, zlib.error) as exc:
            # Malformed record inside one member: skip the member, keep the run.
            corrupt.append((label, str(exc)))
            sys.stderr.write("  UNREADABLE {}: {}\n".format(label, exc))
            continue

        now = time.time()
        if args.progress_interval > 0 and now - last >= args.progress_interval:
            last = now
            el = now - start
            sys.stderr.write("  {} members  {} molecules  {:,.0f} mol/s\n".format(
                members, seq, seq / el if el else 0))
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
    if corrupt:
        sys.stderr.write("  {}: {} corrupt/unreadable file(s) (bad gzip data):\n".format(
            "WARNING" if args.skip_corrupt else "ERROR", len(corrupt)))
        for p, e in corrupt[:10]:
            sys.stderr.write("    {}: {}\n".format(p, e))
    if truncated or corrupt:
        n = len(truncated) + len(corrupt)
        if args.skip_corrupt:
            sys.stderr.write("  --skip-corrupt: kept what read cleanly and dropped "
                             "the {} damaged file(s) above; exiting 0.\n".format(n))
        else:
            sys.stderr.write("  Repair or exclude these, then resubmit this shard "
                             "(or rerun with --skip-corrupt to drop them).\n")
            sys.exit(1)


if __name__ == "__main__":
    main()
