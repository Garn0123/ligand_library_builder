#!/usr/bin/env python3
"""
Stage 4 -- rebalance existing chunks without re-reading the source tree.

Redistributes molecules across a new set of chunks at the RECORD level.
Shuffling whole chunk files accomplishes nothing when input and output chunk
sizes match -- each output would just be one input -- so this streams every
record and reassigns it individually.

    python3 04_rebalance.py -c chunks -m chunks/manifest.tsv \\
        -o rebalanced --chunks 10 --mode greedy --weight lines:C

Modes:
    stride  record i goes to bin i % N. Deterministic systematic sampling
            across the input order; needs no cost model.
    greedy  each record goes to whichever bin currently has the least
            accumulated weight. Balances an explicit cost metric.

Writes a fresh manifest.tsv with the same schema as stage 1, carrying the
original provenance through, so stages 2 and 3 work on the output unchanged.

Cost: one decompress+recompress pass over the chunks. No source tree walk.
"""

import argparse
import heapq
import os
import resource
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root holds db2common.py
from db2common import (iter_records, extract_id, open_gz_text, read_tsv,
                       make_weight_fn, Progress, MANIFEST_HEADER)


def check_fd_limit(n):
    """N output chunks means N simultaneously open gzip writers."""
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    need = n + 16  # headroom for stdio, inputs, manifest
    if soft >= need:
        return
    if hard >= need:
        resource.setrlimit(resource.RLIMIT_NOFILE, (need, hard))
        sys.stderr.write("raised open-file limit to {}\n".format(need))
        return
    sys.exit(
        "need {} open files for --chunks {} but the hard limit is {}.\n"
        "Either lower --chunks, or raise the limit (ulimit -n).".format(
            need, n, hard))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--chunk-dir", required=True)
    ap.add_argument("-o", "--output-dir", required=True)
    ap.add_argument("-m", "--manifest",
                    help="stage 1 manifest.tsv, to carry original provenance "
                         "through. Without it, provenance is recorded as the "
                         "old chunk name and position.")
    ap.add_argument("-N", "--chunks", type=int, required=True,
                    help="number of output chunks")
    ap.add_argument("--mode", choices=["stride", "greedy"], default="stride")
    ap.add_argument("--weight", default="count", metavar="SPEC",
                    help="cost metric for greedy: count, bytes, or lines:X "
                         "(default: count)")
    ap.add_argument("-p", "--prefix", default="chunk")
    ap.add_argument("-l", "--compresslevel", type=int, default=6)
    ap.add_argument("--progress-interval", type=float, default=30.0, metavar="SEC")
    args = ap.parse_args()

    if args.chunks < 1:
        sys.exit("--chunks must be >= 1")
    check_fd_limit(args.chunks)

    chunks = sorted(f for f in os.listdir(args.chunk_dir)
                    if f.endswith(".db2.gz"))
    if not chunks:
        sys.exit("no .db2.gz chunks found in {}".format(args.chunk_dir))
    sys.stderr.write("{} input chunk(s), redistributing into {}\n".format(
        len(chunks), args.chunks))

    weight_fn = make_weight_fn(args.weight)
    os.makedirs(args.output_dir, exist_ok=True)

    digits = max(5, len(str(args.chunks)))
    names = ["{}_{:0{d}d}.db2.gz".format(args.prefix, i + 1, d=digits)
             for i in range(args.chunks)]
    writers = [open_gz_text(os.path.join(args.output_dir, n), "wt",
                            args.compresslevel) for n in names]
    counts = [0] * args.chunks
    weights = [0.0] * args.chunks
    heap = [(0.0, i) for i in range(args.chunks)]
    heapq.heapify(heap)

    man_rows = read_tsv(args.manifest, MANIFEST_HEADER) if args.manifest else None
    man_fh = open(os.path.join(args.output_dir, "manifest.tsv"), "w")
    man_fh.write("\t".join(MANIFEST_HEADER) + "\n")

    try:
        total_bytes = sum(os.path.getsize(os.path.join(args.chunk_dir, c))
                          for c in chunks)
    except OSError:
        total_bytes = 0
    prog = Progress(args.progress_interval, total=len(chunks), unit="chunks",
                    total_bytes=total_bytes)

    seq = 0
    done_bytes = 0
    misaligned = []

    for ci, cname in enumerate(chunks):
        path = os.path.join(args.chunk_dir, cname)
        with open_gz_text(path, "rt") as fh:
            for idx, (lines, complete) in enumerate(iter_records(fh)):
                if not complete:
                    sys.exit("{}: record {} has no E terminator; "
                             "refusing to rebalance a damaged chunk".format(
                                 cname, idx))

                src_file, src_idx, mol_id = cname, str(idx), extract_id(lines)
                if man_rows is not None:
                    try:
                        m_chunk, m_idx, m_src, m_sidx, m_id = next(man_rows)
                    except StopIteration:
                        sys.exit("manifest ran out at record {} -- it does not "
                                 "match this chunk directory".format(seq))
                    if m_chunk != cname or int(m_idx) != idx:
                        misaligned.append((seq, cname, idx, m_chunk, m_idx))
                        if len(misaligned) > 5:
                            break
                    else:
                        src_file, src_idx, mol_id = m_src, m_sidx, m_id

                if args.mode == "stride":
                    b = seq % args.chunks
                else:
                    w, b = heapq.heappop(heap)

                rec_w = weight_fn(lines)
                writers[b].writelines(lines)
                man_fh.write("{}\t{}\t{}\t{}\t{}\n".format(
                    names[b], counts[b], src_file, src_idx, mol_id))
                counts[b] += 1
                weights[b] += rec_w
                if args.mode == "greedy":
                    heapq.heappush(heap, (weights[b], b))

                seq += 1
                if seq % 10000 == 0:
                    prog.tick(ci, seq, done_bytes=done_bytes)

        done_bytes += os.path.getsize(path)
        prog.tick(ci + 1, seq, done_bytes=done_bytes)
        if misaligned and len(misaligned) > 5:
            break

    for w in writers:
        w.close()
    man_fh.close()

    if misaligned:
        sys.stderr.write("\nERROR: manifest does not line up with the chunks.\n"
                         "Was it produced by a different stage 1 run?\n")
        for seqn, cn, i, mc, mi in misaligned[:5]:
            sys.stderr.write("  record {}: chunk {} idx {} vs manifest "
                             "{} idx {}\n".format(seqn, cn, i, mc, mi))
        sys.exit(1)

    if man_rows is not None:
        leftover = sum(1 for _ in man_rows)
        if leftover:
            sys.exit("\nERROR: manifest has {} more rows than the chunks "
                     "contain".format(leftover))

    prog.tick(len(chunks), seq, force=True, done_bytes=done_bytes)
    sys.stderr.write("\nmolecules redistributed: {}\n".format(seq))
    sys.stderr.write("output: {} ({} chunks)\n".format(args.output_dir,
                                                       args.chunks))
    lo, hi = min(counts), max(counts)
    sys.stderr.write("  molecules  min {}  max {}\n".format(lo, hi))
    if args.weight != "count":
        wlo, whi = min(weights), max(weights)
        sys.stderr.write("  weight ({})  min {:,.0f}  max {:,.0f}  "
                         "spread {:.2f}x\n".format(
                             args.weight, wlo, whi,
                             whi / wlo if wlo else float("inf")))
    sys.stderr.write("manifest: {}/manifest.tsv\n".format(args.output_dir))
    sys.stderr.write("\nNOTE: positions have changed. Any labels.tsv from a "
                     "previous run is now stale -- re-run 02_label.py against "
                     "the new manifest.\n")


if __name__ == "__main__":
    main()
