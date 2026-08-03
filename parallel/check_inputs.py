#!/usr/bin/env python3
"""
Preflight -- scan a db2 tree for unreadable, corrupt, or truncated files.

A single damaged .db2.gz (a partial download, a bad rsync) makes a collect task
raise deep inside gzip/zlib. Because that error surfaces mid-stream, the planner
(make_shards.py) can't see it -- it only stats file sizes. This tool actually
decompresses every file at the record level and classifies it:

    OK         reads cleanly to a record boundary
    CORRUPT    gzip/zlib error or otherwise unreadable (bad compressed data)
    TRUNCATED  ends mid-record (a valid gzip, but no final E terminator)
    NOT_GZIP   readable, but not gzip despite the .gz name (read as plain text)

    python3 check_inputs.py -i /path/to/3D --good-list good.txt --bad-list bad.txt

Exit status is nonzero if any file is CORRUPT or TRUNCATED, so it drops into a
shell `&&` chain. Feed --good-list into make_shards.py / 01_chunk.py with
--file-list to run the pipeline over just the clean files while you repair or
re-fetch the bad ones. A full scan costs roughly one decompression pass over the
library; use --jobs to spread it across cores.
"""

import argparse
import gzip
import os
import sys
import time
import zlib
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root holds db2common.py
from db2common import (find_inputs, dedupe_by_inode, is_gzip, iter_records_bytes,
                       iter_sources, TAR_EXTS, INPUT_SUFFIXES)


def check_archive(path):
    """Fully stream every member of a .db2.tgz. Returns (path, status, detail)."""
    container = []
    members = 0
    bad = []
    for label, fh in iter_sources([path], lambda p, e: container.append(e)):
        members += 1
        try:
            last_complete = True
            for _rec, complete in iter_records_bytes(fh):
                last_complete = complete
            if not last_complete:
                bad.append(label.split("::")[-1])
        except (OSError, EOFError, zlib.error) as exc:
            bad.append("{} ({})".format(label.split("::")[-1], exc))
    if container:
        return (path, "CORRUPT", "archive: {}".format(container[0]))
    if bad:
        return (path, "TRUNCATED", "{} bad member(s), e.g. {}".format(
            len(bad), bad[0]))
    return (path, "OK", "{} members".format(members))


def check_file(path):
    """Fully stream one source. Returns (path, status, detail)."""
    if path.endswith(TAR_EXTS):
        return check_archive(path)
    gz = is_gzip(path)
    opener = gzip.open if gz else open
    try:
        n = 0
        last_complete = True
        with opener(path, "rb") as fh:
            for _rec, complete in iter_records_bytes(fh):
                last_complete = complete
                n += 1
        if not last_complete:
            return (path, "TRUNCATED", "{} records, last has no E".format(n))
        if not gz:
            return (path, "NOT_GZIP", "{} records (plain text)".format(n))
        return (path, "OK", str(n))
    except (OSError, EOFError, zlib.error) as exc:
        # zlib.error (bad deflate body) is not an OSError -- catch it explicitly.
        return (path, "CORRUPT", str(exc))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input-dir", default=".")
    ap.add_argument("--file-list", help="check these paths instead of walking the tree")
    ap.add_argument("--suffix", action="append", metavar="EXT",
                    help="filename suffix to collect; repeatable (default .db2.gz)")
    ap.add_argument("--good-list", metavar="FILE",
                    help="write OK + NOT_GZIP paths here (feed to --file-list)")
    ap.add_argument("--bad-list", metavar="FILE",
                    help="write CORRUPT + TRUNCATED paths here")
    ap.add_argument("-j", "--jobs", type=int, default=1,
                    help="parallel worker processes (default 1)")
    ap.add_argument("--keep-duplicate-paths", action="store_true",
                    help="do not dedupe symlink/hardlink duplicates by inode")
    ap.add_argument("--progress-interval", type=float, default=30.0, metavar="SEC",
                    help="seconds between progress lines; 0 disables (default 30)")
    args = ap.parse_args()

    if args.file_list:
        with open(args.file_list) as fh:
            inputs = [ln.strip() for ln in fh if ln.strip()]
    else:
        sys.stderr.write("walking {}...\n".format(args.input_dir))
        inputs = find_inputs(args.input_dir, tuple(args.suffix or INPUT_SUFFIXES))
    if not inputs:
        sys.exit("no input files found")
    if not args.keep_duplicate_paths:
        inputs, _skipped = dedupe_by_inode(inputs)
    sys.stderr.write("checking {} file(s) with {} worker(s)...\n".format(
        len(inputs), args.jobs))

    results = []
    start = last = time.time()

    def note(done):
        nonlocal last
        now = time.time()
        if args.progress_interval > 0 and now - last >= args.progress_interval:
            last = now
            el = now - start
            sys.stderr.write("  {}/{} checked  ({:,.0f}/s)\n".format(
                done, len(inputs), done / el if el else 0))
            sys.stderr.flush()

    if args.jobs > 1:
        with Pool(args.jobs) as pool:
            for i, r in enumerate(pool.imap_unordered(check_file, inputs), 1):
                results.append(r)
                note(i)
    else:
        for i, path in enumerate(inputs, 1):
            results.append(check_file(path))
            note(i)

    by_status = {}
    for path, status, detail in results:
        by_status.setdefault(status, []).append((path, detail))
    ok = sorted(by_status.get("OK", []))
    not_gzip = sorted(by_status.get("NOT_GZIP", []))
    corrupt = sorted(by_status.get("CORRUPT", []))
    truncated = sorted(by_status.get("TRUNCATED", []))

    sys.stderr.write("\nchecked {} file(s): {} OK, {} not-gzip, {} CORRUPT, "
                     "{} TRUNCATED\n".format(len(inputs), len(ok), len(not_gzip),
                                             len(corrupt), len(truncated)))

    if args.good_list:
        with open(args.good_list, "w") as fh:
            for path, _ in ok + not_gzip:
                fh.write(path + "\n")
        sys.stderr.write("good-list ({} files): {}\n".format(
            len(ok) + len(not_gzip), args.good_list))
    if args.bad_list:
        with open(args.bad_list, "w") as fh:
            for path, _ in corrupt + truncated:
                fh.write(path + "\n")
        sys.stderr.write("bad-list ({} files): {}\n".format(
            len(corrupt) + len(truncated), args.bad_list))

    for label, items in (("CORRUPT", corrupt), ("TRUNCATED", truncated)):
        if items:
            sys.stderr.write("\n{} {} file(s):\n".format(len(items), label))
            for path, detail in items[:50]:
                sys.stderr.write("  {}: {}\n".format(path, detail))
            if len(items) > 50:
                sys.stderr.write("  ...and {} more\n".format(len(items) - 50))
    if not_gzip:
        sys.stderr.write("\nNOTE: {} file(s) are not gzip despite their .gz name; "
                         "the pipeline reads them as plain text.\n".format(
                             len(not_gzip)))

    if corrupt or truncated:
        sys.exit(1)
    sys.stderr.write("\nall files read cleanly.\n")


if __name__ == "__main__":
    main()
