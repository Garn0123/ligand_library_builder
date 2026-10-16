#!/usr/bin/env python3
"""
timing_report.py -- what did a run cost, and what came out of it?

    llb timing triage protomers_qk db2_run
    llb timing db2_run --library library/library.tsv      # db2 time by heavy atoms

Reads whatever each directory holds and reports it:

  triage dir      shards/*/DONE (parents, elapsed_s)       ms per parent
  QupKake dir     shard_*/DONE  (parents, elapsed_s)       s per parent, the slowest
                                                           shard (for --time), CPU-h
  db2 run dir     out/*/llb_build.tsv, one row per         s per molecule (median,
                  molecule, written by the build launcher  p90, max), CPU-h, failures
                  as each finishes                         by code, conformers
                                                           produced, db2 size

and, for a db2 run, projects CPU-hours and disk per 100k molecules from the
medians and means. Shards written before timing was recorded (no elapsed_s, no
llb_build.tsv) are counted as "untimed", never guessed.

stdlib only.
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from collections import Counter
from pathlib import Path

N_OUT_BINS = [(0, 0, "0"), (1, 24, "1-24"), (25, 99, "25-99"), (100, 299, "100-299"),
              (300, 10 ** 9, "300+")]


def read_kv(p: Path) -> dict:
    return dict(l.split("=", 1) for l in p.read_text().splitlines() if "=" in l)


def pct(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, int(q * (len(v) - 1) + 0.5))] if v else float("nan")


def shard_report(label: str, dones: list[Path], unit_ms: bool) -> None:
    timed, untimed = [], 0
    for d in dones:
        kv = read_kv(d)
        if "elapsed_s" in kv and "parents" in kv:
            timed.append((float(kv["elapsed_s"]), int(kv["parents"]), d.parent.name))
        else:
            untimed += 1
    print(f"\n== {label}: {len(dones)} finished shard(s)"
          + (f", {untimed} untimed (finished before timing was recorded)" if untimed else ""))
    if not timed:
        return
    sec = sum(t for t, _, _ in timed)
    n = sum(p for _, p, _ in timed)
    per = [t / p for t, p, _ in timed if p]
    slow = max(timed)
    if unit_ms:
        print(f"  {n:,} parents in {sec:,.1f} s: {1000 * sec / n:.1f} ms per parent")
    else:
        print(f"  {n:,} parents, {sec / 3600:,.2f} CPU-h: {sec / n:,.1f} s per parent overall; "
              f"per shard median {statistics.median(per):,.1f} s, p90 {pct(per, 0.9):,.1f} s")
        print(f"  slowest shard {slow[2]}: {slow[0] / 60:,.1f} min for {slow[1]} parents "
              f"(size --time from this, with headroom)")


def build_tables(run: Path) -> list[Path]:
    """llb_build.tsv files: RUN/out/NNNN/ (submit_db2.sh), or an --outputpath
    itself / its subdirectories. Listed, not walked: a recursive walk crawls
    every molecule's working directory and on NFS takes hours."""
    return sorted({*run.glob("llb_build.tsv"), *run.glob("*/llb_build.tsv"),
                   *run.glob("out/*/llb_build.tsv")})


def db2_report(run: Path, library: Path | None) -> None:
    rows = {}
    files = build_tables(run)
    for f in files:
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                rows[r["name"]] = r                  # a rerun's later row wins
    print(f"\n== db2 build {run}: {len(rows):,} molecule(s) timed, from {len(files)} shard file(s)")
    if not rows:
        print("  no llb_build.tsv: built before per-molecule timing, or by plain build_ligand. "
              "Fallback: db2_converter's ALL SUMMARY table in out/NNNN/*.log")
        return
    el = [float(r["elapsed_s"]) for r in rows.values()]
    built = [r for r in rows.values() if int(r["db2_bytes"] or 0) > 0]
    fails = Counter(code for r in rows.values() for code in r["failed"].split("_")
                    if code and code != "-")
    n_out = [int(r["n_out"]) for r in rows.values() if r["n_out"] not in ("", "None", "NA")]
    sizes = [int(r["db2_bytes"]) for r in built]
    sec = sum(el)
    print(f"  time per molecule: median {statistics.median(el):,.1f} s, p90 {pct(el, 0.9):,.1f} s, "
          f"max {max(el):,.1f} s; total {sec / 3600:,.2f} CPU-h")
    print(f"  built {len(built):,} of {len(rows):,} ({100 * len(built) / len(rows):.1f}%)"
          + ("; failures: " + ", ".join(f"{k} {v}" for k, v in fails.most_common()) if fails else ""))
    if n_out:
        bins = Counter()
        for v in n_out:
            for lo, hi, lab in N_OUT_BINS:
                if lo <= v <= hi:
                    bins[lab] += 1
        print("  conformers produced (n_out): " +
              ", ".join(f"{lab} {bins[lab]}" for _, _, lab in N_OUT_BINS) +
              f"  -> {100 * (bins['0'] + bins['1-24']) / len(n_out):.1f}% under 25")
    if sizes:
        print(f"  db2.gz size: median {statistics.median(sizes) / 1024:,.0f} KB, "
              f"mean {statistics.mean(sizes) / 1024:,.0f} KB, total "
              f"{sum(sizes) / 1024 ** 3:,.2f} GB")
        print(f"  per 100k molecules at these rates: {statistics.mean(el) * 1e5 / 3600:,.0f} CPU-h, "
              f"{statistics.mean(sizes) * 1e5 / 1024 ** 3:,.0f} GB of db2")
    if library is not None:
        with open(library) as fh:
            ha = {r["name"]: int(r["heavy_atoms"]) for r in csv.DictReader(fh, delimiter="\t")}
        by = {}
        for name, r in rows.items():
            if name in ha:
                by.setdefault(ha[name], []).append(float(r["elapsed_s"]))
        print("  by heavy atoms: " + "  ".join(
            f"H{h}: {statistics.median(v):,.0f}s (n={len(v)})" for h, v in sorted(by.items())))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", type=Path,
                    help="triage -o, a QupKake OUT dir, and/or a db2 RUNDIR, in any order")
    ap.add_argument("--library", type=Path, default=None,
                    help="library.tsv, to break db2 time down by heavy atoms")
    args = ap.parse_args(argv)
    for d in args.dirs:
        if not d.is_dir():
            print(f"\n== {d}: not a directory", file=sys.stderr)
            continue
        if (d / "settings.json").exists() and (d / "shards").is_dir():
            shard_report(f"triage {d}", sorted((d / "shards").glob("shard_*/DONE")), True)
        elif any(d.glob("shard_*/DONE")):
            shard_report(f"QupKake {d}", sorted(d.glob("shard_*/DONE")), False)
        elif (d / "out").is_dir() or build_tables(d):
            db2_report(d, args.library)
        else:
            print(f"\n== {d}: no triage, QupKake or db2 run found here", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
