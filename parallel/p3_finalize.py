#!/usr/bin/env python3
"""
Phase 3 -- stitch the per-shard manifests into one global manifest.

A record's position inside a finished chunk is not its position inside the
shard that produced it: chunk = part(shard 0) + part(shard 1) + ... So the
global index is the sum of every earlier shard's contribution to that bin,
plus the record's local index. The counts files give those offsets.

    python3 p3_finalize.py -w work -o chunks

Writes chunks/manifest.tsv with the same schema as the serial pipeline, so
02_label.py and 03_apply.py work on the result unchanged.
"""

import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root holds db2common.py
from db2common import MANIFEST_HEADER


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-w", "--work-dir", required=True)
    ap.add_argument("-o", "--output-dir", required=True)
    ap.add_argument("-p", "--prefix", default="chunk")
    ap.add_argument("--keep-intermediates", action="store_true")
    args = ap.parse_args()

    cfiles = sorted(glob.glob(os.path.join(args.work_dir, "counts",
                                           "shard_*.tsv")))
    if not cfiles:
        sys.exit("no counts files in {}/counts".format(args.work_dir))

    shards = []
    counts = {}          # shard -> {bin: n}
    weights = {}
    nbins = 0
    for cf in cfiles:
        s = os.path.basename(cf)[len("shard_"):-len(".tsv")]
        shards.append(s)
        counts[s] = {}
        weights[s] = {}
        with open(cf) as fh:
            for line in fh:
                b, n, w = line.split("\t")
                counts[s][int(b)] = int(n)
                weights[s][int(b)] = float(w)
                nbins = max(nbins, int(b) + 1)
    sys.stderr.write("{} shards, {} bins\n".format(len(shards), nbins))

    # offset[bin][shard] = records written to that bin by all earlier shards
    offset = {}
    bin_totals = [0] * nbins
    bin_weights = [0.0] * nbins
    for b in range(nbins):
        run = 0
        for s in shards:
            offset[(b, s)] = run
            run += counts[s].get(b, 0)
            bin_weights[b] += weights[s].get(b, 0.0)
        bin_totals[b] = run

    total = sum(bin_totals)
    sys.stderr.write("total molecules: {}\n".format(total))

    os.makedirs(args.output_dir, exist_ok=True)
    tmpdir = os.path.join(args.work_dir, "manifest_parts")
    os.makedirs(tmpdir, exist_ok=True)

    # One pass over shard manifests, routing rows to per-bin files. Shards are
    # visited in order, so each bin file comes out already sorted by index.
    outs = [open(os.path.join(tmpdir, "bin_{:05d}.tsv".format(b + 1)), "w")
            for b in range(nbins)]
    seen = 0
    try:
        for s in shards:
            mpath = os.path.join(args.work_dir, "manifests",
                                 "shard_{}.tsv".format(s))
            if not os.path.exists(mpath):
                sys.exit("missing manifest for shard {}".format(s))
            with open(mpath) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    b, local, src, sidx, mol_id = line.split("\t")
                    b = int(b)
                    gidx = offset[(b, s)] + int(local)
                    outs[b].write("{}_{:05d}.db2.gz\t{}\t{}\t{}\t{}\n".format(
                        args.prefix, b + 1, gidx, src, sidx, mol_id))
                    seen += 1
    finally:
        for o in outs:
            o.close()

    if seen != total:
        sys.exit("ERROR: manifests hold {} rows but counts say {}".format(
            seen, total))

    final = os.path.join(args.output_dir, "manifest.tsv")
    with open(final, "w") as out:
        out.write("\t".join(MANIFEST_HEADER) + "\n")
        for b in range(nbins):
            p = os.path.join(tmpdir, "bin_{:05d}.tsv".format(b + 1))
            with open(p) as fh:
                for line in fh:
                    out.write(line)
    sys.stderr.write("manifest: {}\n".format(final))

    lo, hi = min(bin_totals), max(bin_totals)
    sys.stderr.write("\nchunk balance:\n")
    sys.stderr.write("  molecules  min {}  max {}  spread {:.3f}x\n".format(
        lo, hi, float(hi) / lo if lo else float("inf")))
    wlo, whi = min(bin_weights), max(bin_weights)
    if whi > 0 and wlo > 0 and whi != float(hi):
        sys.stderr.write("  weight     min {:,.0f}  max {:,.0f}  "
                         "spread {:.3f}x\n".format(wlo, whi, whi / wlo))

    # cross-check against what phase 2 actually produced
    missing = []
    for b in range(nbins):
        p = os.path.join(args.output_dir,
                         "{}_{:05d}.db2.gz".format(args.prefix, b + 1))
        if not os.path.exists(p):
            missing.append(os.path.basename(p))
    if missing:
        sys.stderr.write("\nWARNING: {} chunk file(s) missing -- did every "
                         "phase 2 task succeed?\n".format(len(missing)))
        for m in missing[:10]:
            sys.stderr.write("  {}\n".format(m))
        sys.exit(1)

    if not args.keep_intermediates:
        for p in glob.glob(os.path.join(tmpdir, "bin_*.tsv")):
            os.unlink(p)
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass


if __name__ == "__main__":
    main()
