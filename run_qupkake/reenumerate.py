#!/usr/bin/env python3
"""
reenumerate.py -- rebuild a finished QupKake array's protomers from its SDFs.

    llb reenumerate PARENTS QK_IN QK_OUT [qupkake_protomers.py options]

QupKake's expensive part (xtb + the GNN) ends in
QK_IN/shard_NNNNN/qupkake_work/output/qupkake_output.sdf. Everything after that
-- site filtering, Henderson-Hasselbalch states, stereo -- is cheap and lives in
qupkake_protomers.py. When that code changes (e.g. the amide-site filter,
2026-10-08), this re-runs it for every DONE shard with --from-sdf, so no xtb
runs again. Seconds per shard.
QK_OUT/shard_NNNNN/qupkake_work is a symlink to QK_IN's, so keep QK_IN.

PARENTS is the parents directory QK_IN was run on (for a triage run:
triage/routed). Each shard is rebuilt with the pH / margin / min population /
max states recorded in its DONE, and QK_OUT/shard_NNNNN/DONE carries the old
DONE plus what was re-run, so llb merge and llb names accept QK_OUT exactly as
they accepted QK_IN. QK_IN is never written to.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import io
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import qupkake_protomers  # noqa: E402


def read_done(p: Path) -> dict:
    return dict(l.split("=", 1) for l in p.read_text().splitlines() if "=" in l)


def commit() -> str:
    try:
        return subprocess.run(["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True).stdout.strip() or "unknown"
    except OSError:
        return "unknown"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("parents", type=Path, help="parents directory QK_IN was run on")
    ap.add_argument("qk_in", type=Path)
    ap.add_argument("qk_out", type=Path)
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="show qupkake_protomers.py's own output per shard")
    args, extra = ap.parse_known_args(argv)
    if args.qk_out.resolve() == args.qk_in.resolve():
        raise SystemExit("QK_OUT must be a new directory; QK_IN is kept as it was")

    with open(args.parents / "shards.tsv") as fh:
        manifest = {int(r["index"]): r for r in csv.DictReader(fh, delimiter="\t")}
    rev = commit()
    done_n = skipped = 0
    drops = 0
    for d in sorted(args.qk_in.glob("shard_*")):
        if not (d / "DONE").exists():
            skipped += 1
            continue
        old = read_done(d / "DONE")
        idx = int(old.get("shard_index", d.name.split("_")[1]))
        r = manifest.get(idx)
        if r is None or old.get("shard_sha256") != r["sha256"]:
            raise SystemExit(f"{d}: its DONE does not match shard {idx} of {args.parents}; "
                             f"is PARENTS the directory this array was run on?")
        sdf = d / "qupkake_work" / "output" / "qupkake_output.sdf"
        if not sdf.exists():
            raise SystemExit(f"{sdf} missing; this shard cannot be rebuilt without QupKake")
        out = args.qk_out / d.name
        if (out / "DONE").exists():
            continue
        cli = [str(args.parents / r["path"]), "-o", str(out), "--from-sdf", str(sdf),
               "--stereo", "carry", "--name-style", "index",
               "--ph", *old.get("ph", "6.4 7.4 8.4").split(),
               "--margin", old.get("margin", "1.0"),
               "--min-population", old.get("min_population", "0.01"),
               "--max-states", old.get("max_states", "8"), *extra]
        sink = contextlib.nullcontext() if args.verbose else contextlib.redirect_stdout(io.StringIO())
        with sink:
            rc = qupkake_protomers.main(cli)
        if rc:
            raise SystemExit(f"{d.name}: qupkake_protomers.py exited {rc}")
        # point at the QupKake working tree it was built from, so llb compare
        # and anything else reading qupkake_work/ sees the same SDF
        link = out / "qupkake_work"
        if not link.exists():
            link.symlink_to((d / "qupkake_work").resolve(), target_is_directory=True)
        with open(out / "dropped_sites.tsv") as fh:
            drops += sum(1 for _ in fh) - 1
        (out / "DONE").write_text(
            (d / "DONE").read_text().rstrip("\n") + "\n"
            + f"reenumerated_from={d.resolve()}\nreenumerated_commit={rev}\n"
            + f"reenumerated_args={' '.join(extra)}\n"
            + f"reenumerated={time.strftime('%Y-%m-%dT%H:%M:%S')}\n")
        done_n += 1
    print(f"{done_n} shard(s) rebuilt from their QupKake SDFs -> {args.qk_out}  "
          f"({drops} implausible site(s) dropped; {skipped} shard(s) without DONE skipped)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
