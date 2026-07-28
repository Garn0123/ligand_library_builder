#!/usr/bin/env python3
"""
Stage 1 -- chunk a tree of .db2.gz files into fixed-size pieces.

Splits only on record boundaries, never mid-molecule. Writes a record-level
manifest that stages 2 and 3 use to address individual molecules by their
position within a chunk (position, not id -- ids may repeat).

    python3 01_chunk.py -i . -o chunks -n 50000

Outputs:
    chunks/chunk_00001.db2.gz, ...
    chunks/manifest.tsv        chunk, chunk_idx, source_file, source_idx, original_id
    chunks/chunk_sources.tsv   chunk, source_file, molecules   (human summary)
"""

import argparse
import os
import random
import sys
import zlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root holds db2common.py
from db2common import (find_inputs, dedupe_by_inode, iter_records, extract_id,
                       open_gz_text, is_gzip, Progress, make_weight_fn,
                       MANIFEST_HEADER)


class ChunkWriter:
    """Rotates gzip output every `chunk_size` records."""

    def __init__(self, outdir, prefix, chunk_size, compresslevel=6, digits=5,
                 target_weight=0.0):
        self.outdir = outdir
        self.prefix = prefix
        self.chunk_size = chunk_size
        self.compresslevel = compresslevel
        self.digits = digits
        self.target_weight = target_weight
        self.index = 0
        self.fh = None
        self.count_in_chunk = 0
        self.weight_in_chunk = 0.0
        self.total = 0
        self.stats = []   # (chunk_name, molecules, weight) per finished chunk

    @property
    def current_name(self):
        return "{}_{:0{d}d}.db2.gz".format(self.prefix, self.index, d=self.digits)

    def write_record(self, lines, weight=1.0):
        """Write one record. Returns (chunk_name, index_within_chunk).

        Rotates when EITHER the molecule cap or the weight target is hit,
        so a weight-balanced run still respects --chunk-size as a ceiling.
        """
        if self.fh is None:
            self.index += 1
            self.fh = open_gz_text(os.path.join(self.outdir, self.current_name),
                                   "wt", self.compresslevel)
            self.count_in_chunk = 0
            self.weight_in_chunk = 0.0
        self.fh.writelines(lines)
        idx = self.count_in_chunk
        self.count_in_chunk += 1
        self.weight_in_chunk += weight
        self.total += 1
        if (self.count_in_chunk >= self.chunk_size
                or (self.target_weight
                    and self.weight_in_chunk >= self.target_weight)):
            self._finish()
        return self.current_name, idx

    def _finish(self):
        self.stats.append((self.current_name, self.count_in_chunk,
                           self.weight_in_chunk))
        self.fh.close()
        self.fh = None

    def close(self):
        if self.fh is not None:
            self._finish()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input-dir", default=".")
    ap.add_argument("-o", "--output-dir", required=True)
    ap.add_argument("-n", "--chunk-size", type=int, default=50000)
    ap.add_argument("-p", "--prefix", default="chunk")
    ap.add_argument("-l", "--compresslevel", type=int, default=6,
                    help="gzip level 1-9; 1 is much faster (default: 6)")
    ap.add_argument("--file-list", help="read input paths from this file")
    ap.add_argument("--suffix", action="append", metavar="EXT",
                    help="filename suffix to collect; repeatable. "
                         "default: .db2.gz")
    ap.add_argument("--keep-duplicate-paths", action="store_true",
                    help="do not skip symlinks/hardlinks to files already in the set")
    ap.add_argument("--order", choices=["path", "shuffle"], default="path",
                    help="input file order. 'path' is sorted, which in a tranche "
                         "tree means sorted by molecular size. 'shuffle' mixes "
                         "tranches so every chunk is a representative sample "
                         "(default: path)")
    ap.add_argument("--seed", type=int, default=0,
                    help="seed for --order shuffle (default: 0)")
    ap.add_argument("--weight", default="count", metavar="SPEC",
                    help="per-molecule cost: count, bytes, or lines:X "
                         "(e.g. lines:C to weight by conformers). default: count")
    ap.add_argument("--target-weight", type=float, default=0.0,
                    help="close a chunk once its summed weight reaches this. "
                         "--chunk-size still applies as a hard ceiling. "
                         "Use --dry-run to see total weight first.")
    ap.add_argument("--progress-interval", type=float, default=30.0,
                    metavar="SEC",
                    help="seconds between progress lines; 0 disables "
                         "(default: 30)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if args.chunk_size < 1:
        sys.exit("chunk-size must be >= 1")

    if args.file_list:
        with open(args.file_list) as fh:
            inputs = [ln.strip() for ln in fh if ln.strip()]
    else:
        inputs = find_inputs(args.input_dir, tuple(args.suffix or [".db2.gz"]))
    if not inputs:
        sys.exit("no .db2.gz files found")

    sys.stderr.write("found {} input files\n".format(len(inputs)))
    weight_fn = make_weight_fn(args.weight)

    if not args.keep_duplicate_paths:
        inputs, skipped = dedupe_by_inode(inputs)
        if skipped:
            sys.stderr.write("skipping {} symlink/hardlink duplicate(s):\n".format(
                len(skipped)))
            for dup, orig in skipped[:20]:
                sys.stderr.write("  {} -> same file as {}\n".format(dup, orig))
            if len(skipped) > 20:
                sys.stderr.write("  ...and {} more\n".format(len(skipped) - 20))
            sys.stderr.write("processing {} unique files\n".format(len(inputs)))

    if args.order == "shuffle":
        random.Random(args.seed).shuffle(inputs)
        sys.stderr.write("input order shuffled (seed {})\n".format(args.seed))

    writer = None
    man_fh = src_fh = None
    if not args.dry_run:
        os.makedirs(args.output_dir, exist_ok=True)
        writer = ChunkWriter(args.output_dir, args.prefix, args.chunk_size,
                             args.compresslevel,
                             target_weight=args.target_weight)
        man_fh = open(os.path.join(args.output_dir, "manifest.tsv"), "w")
        man_fh.write("\t".join(MANIFEST_HEADER) + "\n")
        src_fh = open(os.path.join(args.output_dir, "chunk_sources.tsv"), "w")
        src_fh.write("chunk\tsource_file\tmolecules\n")

    total = zinc_lines = no_id = 0
    total_weight = 0.0
    bad_files, truncated, plaintext = [], [], []
    try:
        total_bytes = sum(os.path.getsize(p) for p in inputs)
    except OSError:
        total_bytes = 0
    done_bytes = 0
    prog = Progress(args.progress_interval, total=len(inputs), unit="files",
                    total_bytes=total_bytes)

    for n, path in enumerate(inputs, 1):
        run_chunk, run_count = None, 0
        if not is_gzip(path):
            plaintext.append(path)
        try:
            with open_gz_text(path, "rt") as fh:
                for src_idx, (lines, complete) in enumerate(iter_records(fh)):
                    if not complete:
                        truncated.append(path)
                        break
                    zinc_lines += sum(1 for ln in lines if "ZINC" in ln)
                    w = weight_fn(lines)
                    total_weight += w
                    mol_id = extract_id(lines)
                    if mol_id == "NO_ID":
                        no_id += 1
                    if writer is not None:
                        chunk, idx = writer.write_record(lines, w)
                        man_fh.write("{}\t{}\t{}\t{}\t{}\n".format(
                            chunk, idx, path, src_idx, mol_id))
                        if chunk != run_chunk:
                            if run_chunk is not None:
                                src_fh.write("{}\t{}\t{}\n".format(
                                    run_chunk, path, run_count))
                            run_chunk, run_count = chunk, 0
                        run_count += 1
                    total += 1
                    if total % 10000 == 0:
                        prog.tick(n - 1, total, done_bytes=done_bytes)
        except (OSError, EOFError, zlib.error) as exc:
            # zlib.error (bad deflate body) is not an OSError; catch it too so a
            # corrupt .db2.gz is reported, not an uncaught traceback.
            bad_files.append((path, str(exc)))
            continue

        if src_fh is not None and run_chunk is not None:
            src_fh.write("{}\t{}\t{}\n".format(run_chunk, path, run_count))

        try:
            done_bytes += os.path.getsize(path)
        except OSError:
            pass
        prog.tick(n, total, done_bytes=done_bytes)

    if writer is not None:
        writer.close()
    for fh in (man_fh, src_fh):
        if fh is not None:
            fh.close()

    prog.tick(len(inputs), total, force=True, done_bytes=done_bytes)
    sys.stderr.write("\nmolecules (E-terminated records): {}\n".format(total))
    sys.stderr.write("lines containing 'ZINC':          {}\n".format(zinc_lines))
    if total != zinc_lines:
        sys.stderr.write("NOTE: counts differ -- 'grep -c ZINC' is not a reliable\n"
                         "      molecule count here; the record count is.\n")
    if writer is not None:
        sys.stderr.write("wrote {} chunks to {}\n".format(writer.index,
                                                          args.output_dir))
        if not args.target_weight:
            # only meaningful for pure count-based chunking
            expected = -(-total // args.chunk_size)
            assert writer.index == expected or total == 0, "chunk count mismatch"
        sys.stderr.write("manifest: {}/manifest.tsv\n".format(args.output_dir))
        if writer.stats:
            mols = [m for _, m, _ in writer.stats]
            wts = [w for _, _, w in writer.stats]
            sys.stderr.write("\nchunk balance ({}):\n".format(args.weight))
            sys.stderr.write("  molecules  min {}  max {}\n".format(
                min(mols), max(mols)))
            if args.weight != "count":
                lo, hi = min(wts), max(wts)
                sys.stderr.write("  weight     min {:,.0f}  max {:,.0f}"
                                 "  spread {:.2f}x\n".format(
                                     lo, hi, hi / lo if lo else float("inf")))
    else:
        sys.stderr.write("dry run: would write {} chunks of {} molecules\n".format(
            -(-total // args.chunk_size), args.chunk_size))
        if args.weight != "count":
            sys.stderr.write("total weight ({}): {:,.0f}\n".format(
                args.weight, total_weight))
            for k in (10, 50, 100, 500):
                sys.stderr.write("  --target-weight {:,.0f}  -> ~{} chunks\n"
                                 .format(total_weight / k, k).replace(",", ""))

    if plaintext:
        sys.stderr.write("\nNOTE: {} input(s) are NOT gzipped despite their name; "
                         "read as plain text:\n".format(len(plaintext)))
        for p in plaintext[:10]:
            sys.stderr.write("  {}\n".format(p))
        if len(plaintext) > 10:
            sys.stderr.write("  ...and {} more\n".format(len(plaintext) - 10))
    if no_id:
        sys.stderr.write("\nWARNING: {} record(s) have no ZINC token on their "
                         "M lines\n".format(no_id))
    if bad_files:
        sys.stderr.write("\nWARNING: {} unreadable file(s):\n".format(len(bad_files)))
        for p, e in bad_files[:20]:
            sys.stderr.write("  {}: {}\n".format(p, e))
    if truncated:
        sys.stderr.write("\nWARNING: {} file(s) ended mid-record:\n".format(
            len(truncated)))
        for p in truncated[:20]:
            sys.stderr.write("  {}\n".format(p))

    if bad_files or truncated:
        sys.exit(1)


if __name__ == "__main__":
    main()
