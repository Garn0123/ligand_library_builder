#!/usr/bin/env python3
"""timing_summary.py -- median and sum of MolSanitizer's per-step seconds (msani_timing.csv or timing.tsv).

    python3 timing_summary.py msani_timing.csv
    python3 timing_summary.py RUN/pass2/tasks/*/timing.tsv

"Initial embedding" includes starting msani's per-molecule embedding process (and,
without the forkserver preload, its import of msani.conformers); compare it with
"Torsional sampling" to see where the time goes. stdlib only.
"""
import csv
import statistics
import sys

COLS = {"Initial embedding": "embedding", "AMSOL": "amsol", "Torsional sampling": "sampling",
        "DB2 conversion": "db2_conversion", "Total": "msani_total"}

rows = []
for path in sys.argv[1:]:
    with open(path, newline="") as fh:
        rows += list(csv.DictReader(fh, delimiter="\t" if path.endswith(".tsv") else ","))
print(f"   {len(rows)} names from {len(sys.argv) - 1} file(s)")
for csv_col, tsv_col in COLS.items():
    xs = []
    for r in rows:
        v = r.get(csv_col, r.get(tsv_col, ""))
        try:
            xs.append(float(v))
        except (TypeError, ValueError):
            pass
    if xs:
        print(f"   {csv_col:<20} median {statistics.median(xs):7.2f} s   sum {sum(xs):9.1f} s")
