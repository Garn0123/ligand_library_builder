#!/usr/bin/env python3
"""
plan.py -- `SMILES NAME` files -> a MolSanitizer build run directory.

Writes, into RUN (refuses a RUN that already holds a different plan):

  run.json                 options (the chemistry), inputs + sha256, shard count,
                           ligand_library_builder commit; every task reads it
  msani_configurations.yaml  MolSanitizer's defaults file for this run
                           (tasks export MSANI_CONFIG to it, so a personal
                           ~/.config/msani file cannot change a build)
  shards/NNNNN.tsv         gidx, smiles, input_name, parent_id   (no header)
  shards.tsv               shard, n, first_gidx, last_gidx, sha256
  plan_rejected.tsv        input lines that never reach a shard, with the reason
                           (unparsable line, damaged ZINC id, repeated parent id)

gidx numbers every input line 0.. in the order given, across all files, so
collect.py can show that each one ends up built or failed-with-reason.

stdlib only (runs on a login node with any python3 >= 3.8).

    python3 plan.py samples/H18.smi samples/H19.smi ... -o RUN --per-shard 40
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from naming import parent_id  # noqa: E402

# The chemistry. Defaults are MolSanitizer 0.7.4's own (the settings its paper
# benchmarks), except where noted. Every value lands in run.json.
DEFAULTS = {
    "removesalts": True,        # paper's recommended prep; ZINC22 2D is one fragment anyway
    "neutralize": True,         # MolSanitizer turns it on with tautomers/protonation
    "tautomers": True,
    "extended_tautomers": False,
    "protonation": True,
    "pH": 7.0,
    "pH_range": 0.0,
    "metal": "strict",
    "max_stereoisomers": 8,
    "stereo_timeout": 60,       # s per molecule
    "numconfs": 2000,
    "energywindow": 25.0,
    "rmsd": 0.5,
    "nringconfs": 1,
    "embed_timeout": 2,         # min, before the Open Babel fallback
    "random_seed": 42,
    "mol_timeout": 1800,        # s, our cap on one molecule's 3D build when retried alone
    "charge_tol": 0.05,         # |AMSOL total charge - formal charge|, as llb verify-db2
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def open_text(path: Path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path)


def git_state(repo: Path) -> dict:
    def git(*a):
        r = subprocess.run(["git", "-C", str(repo), *a], capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else ""
    return {"commit": git("rev-parse", "HEAD"), "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
            "dirty": bool(git("status", "--porcelain", "--untracked-files=no"))}


def msani_yaml(o: dict) -> str:
    # Keys of msani/msani_configurations.yaml (0.7.4). Batch-mode keys are unused:
    # submission is ours, not msani_batch's.
    return (
        "# Written by ligand_library_builder molsanitizer/plan.py for this run.\n"
        "# Tasks export MSANI_CONFIG=<this file>; every value is also passed on the\n"
        "# command line, so this only shuts out a personal defaults file.\n"
        f"EMBED_METHOD: 'rdkit'\nCORINA: 'corina'\nENERGY_WINDOW: {o['energywindow']}\n"
        f"NUMCONFS: {o['numconfs']}\nMAX_STEREOISOMERS: {o['max_stereoisomers']}\n"
        # PH stays at MolSanitizer's 7/0 whatever the run's pH: the file is only
        # read by the 3D call, and a non-default pH there switches --protonation
        # on (parsers.Sanitycheck), re-protonating the named states. The run's pH
        # goes to the 2D step through the Python API (build_shard.py).
        f"TIMEOUT: {o['embed_timeout']}\nPH: 7\nPH_RANGE: 0\n"
        "SLURM_ACCOUNT: null\nSLURM_PARTITION: null\nLINES_PER_JOB: 200\nTIME_LIMIT: 72\n"
        "MAX_ARRAY_SIZE: 1000\nMAX_JOBS: 2000\nMAX_LIMIT_PROJECT: 5000\n"
        "WHOLE_NODE: false\nWHOLE_NODE_CORES: 32\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", type=Path, help="`SMILES NAME` files (.smi or .smi.gz)")
    ap.add_argument("-o", "--run", type=Path, required=True)
    ap.add_argument("--per-shard", type=int, default=40, help="parents per array task")
    ap.add_argument("--label", default="", help="free text kept in run.json")
    for k, v in DEFAULTS.items():
        flag = "--" + k.replace("_", "-")
        if isinstance(v, bool):
            ap.add_argument(flag, dest=k, action=argparse.BooleanOptionalAction, default=v)
        else:
            ap.add_argument(flag, dest=k, type=type(v), default=v)
    args = ap.parse_args(argv)
    opts = {k: getattr(args, k) for k in DEFAULTS}

    run = args.run.resolve()
    inputs = []
    for p in args.inputs:
        if not p.is_file():
            raise SystemExit(f"no such input: {p}")
        inputs.append({"path": str(p.resolve()), "sha256": sha256(p)})

    if (run / "run.json").exists():
        old = json.loads((run / "run.json").read_text())
        same_inputs = [(i["path"], i["sha256"]) for i in old["inputs"]] == [(i["path"], i["sha256"]) for i in inputs]
        if not same_inputs or old["options"] != opts or old["per_shard"] != args.per_shard:
            raise SystemExit(f"{run}/run.json holds a different plan; use a new run directory")
        print(f"{run}: same plan already written; nothing to do")
        return 0

    rows, rejected, seen = [], [], {}
    gidx = 0
    for inp in inputs:
        inp["first_gidx"] = gidx
        with open_text(Path(inp["path"])) as fh:
            for lineno, line in enumerate(fh, 1):
                if not line.strip():
                    continue
                parts = line.split()
                where = f"{Path(inp['path']).name}:{lineno}"
                if len(parts) < 2:
                    rejected.append((gidx, where, parts[0] if parts else "", "", "", "unparsable_line"))
                else:
                    smi, name = parts[0], parts[1]
                    pid, kind = parent_id(name)
                    if pid is None:
                        rejected.append((gidx, where, smi, name, "", kind))
                    elif pid in seen:
                        rejected.append((gidx, where, smi, name, pid, f"duplicate_parent_id (first: gidx {seen[pid]})"))
                    else:
                        seen[pid] = gidx
                        rows.append((gidx, smi, name, pid))
                gidx += 1
        inp["n_lines"] = gidx - inp["first_gidx"]

    (run / "shards").mkdir(parents=True, exist_ok=True)
    shard_rows = []
    for k in range(0, len(rows), args.per_shard):
        chunk = rows[k:k + args.per_shard]
        idx = len(shard_rows) + 1                       # SLURM array index, 1-based
        path = run / "shards" / f"{idx:05d}.tsv"
        path.write_text("".join("\t".join(map(str, r)) + "\n" for r in chunk))
        shard_rows.append((idx, len(chunk), chunk[0][0], chunk[-1][0], sha256(path)))
    with open(run / "shards.tsv", "w") as fh:
        fh.write("shard\tn\tfirst_gidx\tlast_gidx\tsha256\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in shard_rows)
    with open(run / "plan_rejected.tsv", "w") as fh:
        fh.write("gidx\twhere\tsmiles\tinput_name\tparent_id\treason\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rejected)
    (run / "msani_configurations.yaml").write_text(msani_yaml(opts))

    meta = {
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "host": os.uname().nodename,
        "label": args.label,
        "inputs": inputs,
        "n_input_lines": gidx,
        "n_parents": len(rows),
        "n_plan_rejected": len(rejected),
        "per_shard": args.per_shard,
        "n_shards": len(shard_rows),
        "options": opts,
        "llb": git_state(HERE.parent),
    }
    tmp = run / "run.json.tmp"
    tmp.write_text(json.dumps(meta, indent=2) + "\n")
    tmp.replace(run / "run.json")
    print(f"{gidx:,} input lines -> {len(rows):,} parents in {len(shard_rows):,} shards "
          f"of <= {args.per_shard}; {len(rejected):,} rejected at plan (plan_rejected.tsv)")
    if meta["llb"]["dirty"]:
        print("NOTE: the ligand_library_builder checkout has uncommitted changes (run.json llb.dirty)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
