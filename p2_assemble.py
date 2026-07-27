#!/usr/bin/env python3
"""
Phase 2 -- one SLURM array task per bin. Embarrassingly parallel.

Concatenates that bin's parts into a finished chunk. Gzip members concatenate
into a valid stream, so this is a pure byte copy: no decompression, no
recompression, no parsing. It runs at disk speed.

    python3 p2_assemble.py -w work -b $SLURM_ARRAY_TASK_ID -o chunks

Parts are joined in shard order, which is the order phase 3 assumes when it
computes global record indices.
"""

import argparse
import glob
import os
import shutil
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-w", "--work-dir", required=True)
    ap.add_argument("-b", "--bin", type=int, required=True,
                    help="0-based bin index (SLURM_ARRAY_TASK_ID)")
    ap.add_argument("-o", "--output-dir", required=True)
    ap.add_argument("-p", "--prefix", default="chunk")
    ap.add_argument("--keep-parts", action="store_true",
                    help="do not delete parts after a successful join")
    ap.add_argument("--bufsize", type=int, default=16 << 20)
    args = ap.parse_args()

    bindir = os.path.join(args.work_dir, "parts",
                          "bin_{:05d}".format(args.bin + 1))
    if not os.path.isdir(bindir):
        sys.exit("no such bin directory: {}".format(bindir))

    parts = sorted(glob.glob(os.path.join(bindir, "part_*.db2.gz")))
    if not parts:
        sys.exit("no parts in {}".format(bindir))

    os.makedirs(args.output_dir, exist_ok=True)
    out = os.path.join(args.output_dir,
                       "{}_{:05d}.db2.gz".format(args.prefix, args.bin + 1))
    tmp = out + ".tmp"

    written = 0
    with open(tmp, "wb") as fout:
        for p in parts:
            n = os.path.getsize(p)
            if n == 0:
                continue          # a shard that contributed nothing to this bin
            with open(p, "rb") as fin:
                shutil.copyfileobj(fin, fout, args.bufsize)
            written += n
    os.replace(tmp, out)

    sys.stderr.write("bin {}: joined {} part(s), {:.1f} MB -> {}\n".format(
        args.bin + 1, len(parts), written / 1e6, out))

    if not args.keep_parts:
        for p in parts:
            os.unlink(p)
        try:
            os.rmdir(bindir)
        except OSError:
            pass
        sys.stderr.write("  parts removed\n")


if __name__ == "__main__":
    main()
