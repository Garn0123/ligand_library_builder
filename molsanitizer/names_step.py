#!/usr/bin/env python3
"""
names_step.py -- between the passes: `llb names` over every pass-1 file, then
pass 2's shards.

    python names_step.py RUN [--per-shard 40]        # in an env with RDKit (the msani env)

1. Refuses unless every pass-1 shard is done (a missing shard would silently
   leave its parents out of the library).
2. Runs naming/assign_names.py, the one naming step (UPSTREAM.md section 2):
       --inputs <the run's input files>  --ms RUN/pass1/out/*_clean.smi
       -o RUN/library  --ph <the run's pH>
   It writes library/library.smi, library.tsv, dropped.tsv, stats.json and
   checks stereo closure, one name per structure, invariants 1-3.
3. Cuts library.smi into RUN/pass2/shards/NNNNN.smi of --per-shard names, with
   pass2/shards.tsv and pass2/plan.json (the library.smi sha256 they came from).

Re-running with the same library is a no-op. A different library under existing
pass-2 shards is refused: start a new run directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ASSIGN = HERE.parent / "naming" / "assign_names.py"
sys.path.insert(0, str(HERE))
from collect import ranges, shard_status  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ph_label(o: dict) -> str:
    return f"{o['pH']:g}+-{o['pH_range']:g}" if o["pH_range"] else f"{o['pH']:g}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--per-shard", type=int, default=40, help="names per pass-2 array task")
    args = ap.parse_args(argv)
    run = args.run.resolve()
    meta = json.loads((run / "run.json").read_text())

    shards, done, missing = shard_status(run, 1)
    if missing:
        raise SystemExit(f"pass 1 not finished: shards {ranges(missing)} (msani.sh submit RUN 1 --tasks missing)")

    lib = run / "library"
    p2 = run / "pass2"
    if not (lib / "library.smi").exists():
        ms = sorted((run / "pass1" / "out").glob("*_clean.smi"))
        cmd = [sys.executable, str(ASSIGN), "--inputs", *[i["path"] for i in meta["inputs"]],
               "--ms", *map(str, ms), "-o", str(lib), "--ph", ph_label(meta["options"])]
        print(f"$ {ASSIGN.name} --inputs <{len(meta['inputs'])} files> --ms <{len(ms)} files> "
              f"-o {lib} --ph {ph_label(meta['options'])}")
        r = subprocess.run(cmd)
        if r.returncode:
            return r.returncode
    lib_sha = sha256(lib / "library.smi")

    if (p2 / "plan.json").exists():
        old = json.loads((p2 / "plan.json").read_text())
        if old["library_smi_sha256"] != lib_sha or old["per_shard"] != args.per_shard:
            raise SystemExit(f"{p2} was planned from a different library or shard size; use a new run")
        print(f"{p2}: same plan already written ({old['n_shards']} shards)")
        return 0

    names = (lib / "library.smi").read_text().splitlines()
    (p2 / "shards").mkdir(parents=True, exist_ok=True)
    rows = []
    for k in range(0, len(names), args.per_shard):
        idx = len(rows) + 1
        path = p2 / "shards" / f"{idx:05d}.smi"
        chunk = names[k:k + args.per_shard]
        path.write_text("".join(l + "\n" for l in chunk))
        rows.append((idx, len(chunk), sha256(path)))
    with open(p2 / "shards.tsv", "w") as fh:
        fh.write("shard\tn\tsha256\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    plan = {"library_smi_sha256": lib_sha, "per_shard": args.per_shard, "n_shards": len(rows),
            "n_names": len(names), "assign_names": str(ASSIGN)}
    (p2 / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    print(f"{len(names):,} names -> {len(rows):,} pass-2 shards of <= {args.per_shard}")
    print(f"next: msani.sh submit {run} 2")
    return 0


if __name__ == "__main__":
    sys.exit(main())
