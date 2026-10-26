#!/usr/bin/env python3
"""
collect.py -- every finished shard of a MolSanitizer build -> one verified library.

    python3 collect.py RUN                    # the whole run; refuses if a shard is missing
    python3 collect.py RUN --status           # which shards are done / missing (array ranges)
    python3 collect.py RUN --drap-dir DRAP    # also run DRAP's stage4/validate_names.py

Checks, all over the WHOLE library (no shard can see the others):
  * every input line (gidx 0..N-1) is accounted for exactly once: built, or
    failed with a reason (plan_rejected.tsv + each shard's parents.tsv);
  * names are 16 characters, C...C, unique, and so are their db2 tails name[-16:];
  * a base id belongs to one parent;
  * one structure, one name: a SMILES built under two parents keeps the first
    (lowest gidx); the later name is marked DUPLICATE and is not docked;
  * every shard ran the same MolSanitizer / RDKit / options.

Writes into RUN (the layout rna_pipeline pipeline/library_docking/prepare.py reads):
  db2_verification.csv   protomer_name, parent, expected_charge, status, detail, ...
                         status OK = docked; everything else is not
  db2_run/out/NNNNN/     (written by the tasks) <name>.db2.gz
  library/library.tsv    LLB's library.tsv columns (DRAP collapse_invertomers, join)
  library/library.smi    `SMILES NAME`, for validate_names.py
  library/dropped.tsv    named states that are not in the library, and why
  parents_ledger.tsv     one row per input line: status, reason, counts, seconds
  timing.tsv             one row per name: MolSanitizer's step times, ours
  collect.json           counts, reasons, attrition per input file, timing summary
  provenance.json        versions, options, llb commit, hosts, SLURM jobs, time span

stdlib only.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from build_shard import LIBRARY_COLS, PARENT_COLS  # noqa: E402
from naming import CONTRACT  # noqa: E402

VERIFY_COLS = ["protomer_name", "parent", "expected_charge", "status", "detail",
               "db2_file", "db2_name", "amsol_charge"]


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_rows(path: Path, cols: list[str], rows, delimiter="\t") -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter=delimiter, lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def ranges(nums: list[int]) -> str:
    """[1,2,3,7,9,10] -> '1-3,7,9-10' (sbatch --array syntax)."""
    out, nums = [], sorted(nums)
    i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(str(nums[i]) if i == j else f"{nums[i]}-{nums[j]}")
        i = j + 1
    return ",".join(out)


def shard_status(run: Path):
    shards = read_tsv(run / "shards.tsv")
    done, missing = [], []
    for r in shards:
        k = int(r["shard"])
        d = run / "tasks" / f"{k:05d}" / "DONE.json"
        (done if d.exists() and json.loads(d.read_text()).get("input_sha256") == r["sha256"]
         else missing).append(k)
    return shards, done, missing


def pct(xs: list[float], q: float) -> float:
    xs = sorted(xs)
    return round(xs[min(len(xs) - 1, int(q * len(xs)))], 2) if xs else 0.0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--status", action="store_true", help="report done/missing shards and stop")
    ap.add_argument("--drap-dir", type=Path, default=None,
                    help="DRAP checkout; runs stage4/validate_names.py on library.smi")
    args = ap.parse_args(argv)
    run = args.run.resolve()
    meta = json.loads((run / "run.json").read_text())

    shards, done, missing = shard_status(run)
    if args.status:
        print(f"{len(done)}/{len(shards)} shards done")
        if missing:
            print(f"missing: {ranges(missing)}")
        return 0
    if missing:
        print(f"{len(missing)} shard(s) not done: {ranges(missing)}", file=sys.stderr)
        print("resubmit them (msani.sh submit RUN --tasks missing), then collect again", file=sys.stderr)
        return 1

    # ---- read every shard --------------------------------------------------------
    states, parents, timing, dones = [], [], [], []
    for k in done:
        td = run / "tasks" / f"{k:05d}"
        dones.append(json.loads((td / "DONE.json").read_text()))
        for s in read_tsv(td / "states.tsv"):
            s["shard"] = k
            states.append(s)
        parents += read_tsv(td / "parents.tsv")
        timing += [dict(t, shard=k) for t in read_tsv(td / "timing.tsv")]

    problems = []

    # ---- one software stack -----------------------------------------------------
    stacks = Counter(json.dumps({k: v for k, v in d["versions"].items() if k != "msani_path"},
                                sort_keys=True) for d in dones)
    if len(stacks) > 1:
        problems.append(f"shards ran {len(stacks)} different software stacks: {list(stacks)[:2]}")
    if len({d["options_hash"] for d in dones}) > 1:
        problems.append("shards ran with different options")

    # ---- every input line once --------------------------------------------------
    rejected = read_tsv(run / "plan_rejected.tsv")
    ledger = {}
    for r in rejected:
        ledger[int(r["gidx"])] = {"gidx": r["gidx"], "input_name": r["input_name"], "parent_id": r["parent_id"],
                                  "input_smiles": r["smiles"], "status": "failed",
                                  "reason": f"plan: {r['reason']}", "n_states": 0, "n_names": 0, "n_ok": 0}
    dup_gidx = []
    for p in parents:
        g = int(p["gidx"])
        if g in ledger:
            dup_gidx.append(g)
        ledger[g] = p
    n_in = meta["n_input_lines"]
    absent = [g for g in range(n_in) if g not in ledger]
    if dup_gidx or absent or len(ledger) != n_in:
        problems.append(f"ledger: {len(absent)} input line(s) unaccounted (e.g. {absent[:5]}), "
                        f"{len(dup_gidx)} accounted twice")

    # ---- one structure, one name (across shards) ----------------------------------
    owner = {}
    for s in sorted(states, key=lambda s: (int(s["gidx"]), s["name"])):
        if s["status"] in ("DUPLICATE", "EXOTIC"):
            continue
        if s["smiles"] in owner and owner[s["smiles"]] != s["name"]:
            s["status"], s["detail"] = "DUPLICATE", f"same structure as {owner[s['smiles']]}"
        else:
            owner[s["smiles"]] = s["name"]

    library = [s for s in states if s["status"] not in ("DUPLICATE", "EXOTIC")]
    dropped = [s for s in states if s["status"] in ("DUPLICATE", "EXOTIC")]

    # a parent whose every built name turned out to be another parent's structure
    # is not built after all
    ok_by_gidx = Counter(int(s["gidx"]) for s in library if s["status"] == "OK")
    for p in parents:
        g = int(p["gidx"])
        p["n_ok"] = ok_by_gidx.get(g, 0)
        if p["status"] == "built" and not p["n_ok"]:
            why = sorted({f"{s['status']}: {s['detail']}" for s in states if int(s["gidx"]) == g})
            p["status"], p["reason"] = "failed", ("no state built -- " + " | ".join(why))[:500]

    # ---- the naming contract, whole library ----------------------------------------
    names = [s["name"] for s in library]
    bad = [n for n in names if not CONTRACT.match(n)]
    if bad:
        problems.append(f"{len(bad)} name(s) are not C + 14 base62 + C: {bad[:3]}")
    for label, keys in (("names", names), ("db2 tails name[-16:]", [n[-16:] for n in names])):
        dup = [k for k, c in Counter(keys).items() if c > 1]
        if dup:
            problems.append(f"{len(dup)} duplicate {label}, e.g. {dup[:3]}")
    owners = defaultdict(set)
    for s in library:
        owners[s["base_id"]].add(s["parent_id"])
    shared = [b for b, v in owners.items() if len(v) > 1]
    if shared:
        problems.append(f"{len(shared)} base id(s) map to more than one parent: {shared[:3]}")

    # invertomer links must point inside the library: re-point a group whose first
    # member was dropped to its first remaining member
    in_lib = {s["name"] for s in library}
    groups = defaultdict(list)
    for s in library:
        groups[s["invertomer_of"] or s["name"]].append(s)
    for head, members in groups.items():
        if head not in in_lib:
            new = min(m["name"] for m in members)
            for m in members:
                m["invertomer_of"] = "" if m["name"] == new else new

    # ---- write -----------------------------------------------------------------------
    if problems:
        print("collect REFUSED -- nothing written:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 1

    libdir = run / "library"
    libdir.mkdir(exist_ok=True)
    library.sort(key=lambda s: s["name"])
    write_rows(libdir / "library.tsv", LIBRARY_COLS, library)
    (libdir / "library.smi").write_text("".join(f"{s['smiles']} {s['name']}\n" for s in library))
    write_rows(libdir / "dropped.tsv", ["name", "parent_id", "status", "detail", "smiles"], dropped)
    verify = [{"protomer_name": s["name"], "parent": s["parent_id"], "expected_charge": s["net_charge"],
               "status": s["status"], "detail": s["detail"],
               "db2_file": f"{int(s['shard']):05d}/{s['name']}.db2.gz" if s["status"] == "OK" else "",
               "db2_name": s["name"] if s["status"] == "OK" else "", "amsol_charge": s["amsol_charge"]}
              for s in library]
    write_rows(run / "db2_verification.csv", VERIFY_COLS, verify, delimiter=",")
    ledger_rows = [ledger[g] for g in sorted(ledger)]
    write_rows(run / "parents_ledger.tsv", PARENT_COLS, ledger_rows)
    tcols = ["name", "shard", "how", "seconds_wall", "embedding", "amsol", "sampling",
             "db2_conversion", "msani_total"]
    write_rows(run / "timing.tsv", tcols, sorted(timing, key=lambda t: t["name"]))

    # missing db2 for an OK name would make prepare.py stop; check now
    lost = [v["protomer_name"] for v in verify if v["status"] == "OK"
            and not (run / "db2_run" / "out" / v["db2_file"]).exists()]
    if lost:
        print(f"WARNING: {len(lost)} OK name(s) have no db2.gz on disk, e.g. {lost[:3]}", file=sys.stderr)

    # ---- summary ----------------------------------------------------------------------
    def secs(t):
        try:
            return float(t.get("seconds_wall") or t.get("msani_total") or "nan")
        except ValueError:
            return float("nan")
    per_name = [x for x in (secs(t) for t in timing if t.get("how") in ("batch", "alone")) if x == x]
    by_file = {}
    for inp in meta["inputs"]:
        lo, n = inp["first_gidx"], inp["n_lines"]
        rows = [ledger[g] for g in range(lo, lo + n) if g in ledger]
        by_file[Path(inp["path"]).name] = {
            "lines": n, "built": sum(r["status"] == "built" for r in rows),
            "names_ok": sum(int(r.get("n_ok") or 0) for r in rows)}
    task_seconds = [d["seconds"] for d in dones]
    stats = {
        "input_lines": n_in,
        "parents_built": sum(r["status"] == "built" for r in ledger_rows),
        "parents_failed": sum(r["status"] == "failed" for r in ledger_rows),
        "parent_fail_reasons": dict(Counter(r["reason"].split(":")[0].split(" --")[0]
                                            for r in ledger_rows if r["status"] == "failed").most_common()),
        "names": len(library),
        "names_by_status": dict(Counter(s["status"] for s in library).most_common()),
        "dropped_states": dict(Counter(s["status"] for s in dropped)),
        "names_per_built_parent": round(
            sum(s["status"] == "OK" for s in library) / max(1, sum(r["status"] == "built" for r in ledger_rows)), 3),
        "n_h_invertomers": sum(1 for s in library if s["invertomer_of"]),
        "net_charge": dict(sorted(Counter(int(s["net_charge"]) for s in library).items())),
        "by_input_file": by_file,
        "db2_gz_bytes_total": sum(int(s["db2_bytes"] or 0) for s in library),
        "db2_gz_bytes_median": (statistics.median([int(s["db2_bytes"]) for s in library if s["db2_bytes"]])
                                if any(s["db2_bytes"] for s in library) else None),
        "db2_records_median": (statistics.median([int(s["db2_records"]) for s in library if s["db2_records"]])
                               if any(s["db2_records"] for s in library) else None),
        "timing": {
            "names_timed": len(per_name),
            "seconds_per_name_median": round(statistics.median(per_name), 2) if per_name else None,
            "seconds_per_name_p90": pct(per_name, 0.9),
            "seconds_per_name_max": round(max(per_name), 1) if per_name else None,
            "task_seconds_median": round(statistics.median(task_seconds), 1) if task_seconds else None,
            "task_seconds_max": round(max(task_seconds), 1) if task_seconds else None,
            "cpu_hours_total": round(sum(task_seconds) / 3600, 2),
            "cpu_hours_per_1000_parents": round(sum(task_seconds) / 3600 / max(1, len(parents)) * 1000, 2),
        },
    }
    (run / "collect.json").write_text(json.dumps(stats, indent=2) + "\n")
    prov = {
        "run": {k: meta[k] for k in ("created", "label", "inputs", "per_shard", "n_shards", "options", "llb")},
        "software": json.loads(next(iter(stacks))) if stacks else {},
        "msani_path": sorted({d["versions"].get("msani_path", "") for d in dones}),
        "hosts": len({d["host"] for d in dones}),
        "cpus": dict(Counter(d["cpu"] for d in dones)),
        "slurm_jobs": sorted({d["slurm_job"] for d in dones if d["slurm_job"]}),
        "first_start": min(d["started"] for d in dones) if dones else "",
        "last_finish": max(d["finished"] for d in dones) if dones else "",
    }
    (run / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n")

    print(f"{n_in:,} input lines: {stats['parents_built']:,} parents built, "
          f"{stats['parents_failed']:,} failed (reasons in parents_ledger.tsv)")
    print(f"{len(library):,} names: {stats['names_by_status']}; "
          f"{stats['names_per_built_parent']} OK names per built parent")
    print(f"timing: median {stats['timing']['seconds_per_name_median']} s/name, "
          f"{stats['timing']['cpu_hours_total']} CPU-h total, "
          f"{stats['timing']['cpu_hours_per_1000_parents']} CPU-h per 1,000 parents")
    print(f"wrote {run}/db2_verification.csv, library/, parents_ledger.tsv, collect.json, provenance.json")

    if args.drap_dir:
        v = args.drap_dir / "stage4" / "validate_names.py"
        print(f"\n$ python3 {v} --smi {libdir / 'library.smi'}")
        return subprocess.run([sys.executable, str(v), "--smi", str(libdir / "library.smi")]).returncode
    return 0


if __name__ == "__main__":
    sys.exit(main())
