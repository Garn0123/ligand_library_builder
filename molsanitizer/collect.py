#!/usr/bin/env python3
"""
collect.py -- pass 1 + `llb names` + pass 2 -> one verified, accounted-for library.

    python3 collect.py RUN                    # the whole run; refuses if anything is missing
    python3 collect.py RUN --status           # pass 1 / names / pass 2: done and missing shards
    python3 collect.py RUN --drap-dir DRAP    # also run DRAP's stage4/validate_names.py

Checks, over the WHOLE run:
  * every input line (gidx 0..N-1) is accounted for exactly once: built (>= 1
    name with a verified db2), or failed with the reason from whichever step
    lost it -- plan (plan_rejected.tsv), pass 1 (parents.tsv), naming
    (library/dropped.tsv: stereo not closed, duplicate structure, ...) or
    pass 2 (states.tsv: BUILD_FAILED, CHARGE_MISMATCH, EXOTIC, ...);
  * every library name was built in pass 2 exactly once, and nothing else was;
  * every shard of a pass ran the same MolSanitizer / RDKit / options, and both
    passes the same software.

Writes into RUN (the layout rna_pipeline pipeline/library_docking/prepare.py reads):
  db2_verification.csv   protomer_name, parent, expected_charge, status, detail, ...
                         status OK = docked; anything else is not
  db2_run/out/NNNNN/     (written by pass 2) <name>.db2.gz
  library/               (written by `llb names`) library.tsv / .smi / dropped.tsv / stats.json
  parents_ledger.tsv     one row per input line: status, reason, names, OK names, seconds
  timing.tsv             one row per name: MolSanitizer's step times, ours
  collect.json           counts, reasons, attrition per input file, timing per pass
  provenance.json        versions, options, LLB commit, hosts, SLURM jobs, time span

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

VERIFY_COLS = ["protomer_name", "parent", "expected_charge", "status", "detail",
               "db2_file", "db2_name", "amsol_charge"]
LEDGER_COLS = ["gidx", "input_name", "parent_id", "input_smiles", "status", "reason",
               "n_pass1_rows", "n_names", "n_ok", "seconds_pass1", "seconds_pass2"]


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


def shard_status(run: Path, npass: int):
    """(shard rows, done, missing) for a pass; pass 2 has no rows before naming."""
    f = run / f"pass{npass}" / "shards.tsv"
    if not f.exists():
        return [], [], []
    shards = read_tsv(f)
    done, missing = [], []
    for r in shards:
        k = int(r["shard"])
        d = run / f"pass{npass}" / "tasks" / f"{k:05d}" / "DONE.json"
        (done if d.exists() and json.loads(d.read_text()).get("input_sha256") == r["sha256"]
         else missing).append(k)
    return shards, done, missing


def pct(xs: list[float], q: float) -> float:
    xs = sorted(xs)
    return round(xs[min(len(xs) - 1, int(q * len(xs)))], 2) if xs else 0.0


def reason_key(reason: str) -> str:
    """'pass1: removed_in_pass1: organic salt' -> 'pass1: removed_in_pass1' (first cause only)."""
    first = reason.split(" | ")[0]
    parts = [x.strip() for x in first.split(":")]
    return ": ".join(parts[:2]).split(" x")[0] if len(parts) > 1 else first


def stack(d: dict) -> str:
    return json.dumps({k: v for k, v in d["versions"].items() if k != "msani_path"}, sort_keys=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--status", action="store_true", help="report done/missing shards and stop")
    ap.add_argument("--drap-dir", type=Path, default=None,
                    help="DRAP checkout; runs stage4/validate_names.py on library/library.smi")
    args = ap.parse_args(argv)
    run = args.run.resolve()
    meta = json.loads((run / "run.json").read_text())
    s1, d1, m1 = shard_status(run, 1)
    s2, d2, m2 = shard_status(run, 2)
    named = (run / "pass2" / "plan.json").exists()

    if args.status:
        print(f"pass 1: {len(d1)}/{len(s1)} shards done" + (f"; missing {ranges(m1)}" if m1 else ""))
        print("names:  " + ("done" if named else "not yet (msani.sh names RUN, once pass 1 is done)"))
        if named:
            print(f"pass 2: {len(d2)}/{len(s2)} shards done" + (f"; missing {ranges(m2)}" if m2 else ""))
        return 0
    if m1 or not named or m2:
        print("not finished: " + "; ".join(x for x in (
            f"pass 1 shards {ranges(m1)}" if m1 else "", "" if named else "names not run",
            f"pass 2 shards {ranges(m2)}" if m2 else "") if x), file=sys.stderr)
        return 1

    problems = []
    dones = {1: [], 2: []}
    p1 = {}
    for k in d1:
        td = run / "pass1" / "tasks" / f"{k:05d}"
        dones[1].append(json.loads((td / "DONE.json").read_text()))
        for p in read_tsv(td / "parents.tsv"):
            p1[int(p["gidx"])] = p
    states, timing = {}, []
    for k in d2:
        td = run / "pass2" / "tasks" / f"{k:05d}"
        dones[2].append(json.loads((td / "DONE.json").read_text()))
        for s in read_tsv(td / "states.tsv"):
            if s["name"] in states:
                problems.append(f"{s['name']} built in two pass-2 shards")
            s["shard"] = k
            states[s["name"]] = s
        timing += [dict(t, shard=k) for t in read_tsv(td / "timing.tsv")]

    stacks = Counter(stack(d) for n in (1, 2) for d in dones[n])
    if len(stacks) > 1:
        problems.append(f"shards ran {len(stacks)} different software stacks: {list(stacks)[:2]}")
    if len({d["options_hash"] for n in (1, 2) for d in dones[n]}) > 1:
        problems.append("shards ran with different options")

    library = read_tsv(run / "library" / "library.tsv")
    dropped = read_tsv(run / "library" / "dropped.tsv")
    lib_names = {r["name"] for r in library}
    if lib_names != set(states):
        problems.append(f"library vs pass 2: {len(lib_names - set(states))} name(s) never built, "
                        f"{len(set(states) - lib_names)} built but not in the library")

    # ---- one ledger row per input line --------------------------------------------
    by_input = defaultdict(list)
    for r in library:
        by_input[r["input_name"]].append(r["name"])
    why_dropped = defaultdict(Counter)
    for d in dropped:
        if d["input_name"]:
            why_dropped[d["input_name"]][d["reason"]] += 1
    p2_seconds = defaultdict(float)
    for t in timing:
        try:
            p2_seconds[t["name"]] += float(t.get("seconds_wall") or t.get("msani_total") or 0)
        except ValueError:
            pass

    ledger = {}
    for r in read_tsv(run / "plan_rejected.tsv"):
        ledger[int(r["gidx"])] = {"gidx": r["gidx"], "input_name": r["input_name"], "parent_id": r["parent_id"],
                                  "input_smiles": r["smiles"], "status": "failed", "reason": f"plan: {r['reason']}"}
    for g, p in p1.items():
        if g in ledger:
            problems.append(f"gidx {g} both plan-rejected and in pass 1")
        row = {"gidx": g, "input_name": p["input_name"], "parent_id": p["parent_id"],
               "input_smiles": p["input_smiles"], "n_pass1_rows": p["n_rows"], "seconds_pass1": p["seconds_2d"]}
        names = by_input.get(p["input_name"], [])
        ok = [n for n in names if states.get(n, {}).get("status") == "OK"]
        row.update(n_names=len(names), n_ok=len(ok),
                   seconds_pass2=round(sum(p2_seconds[n] for n in names), 3))
        if p["status"] != "ok":
            row.update(status="failed", reason=f"pass1: {p['reason']}")
        elif ok:
            row.update(status="built", reason="")
        else:
            why = [f"names: {k} x{v}" for k, v in sorted(why_dropped[p["input_name"]].items())]
            why += sorted({f"pass2: {states[n]['status']}: {states[n]['detail']}" for n in names if n in states})
            row.update(status="failed", reason=(" | ".join(why) or "no reason recorded")[:500])
        ledger[g] = row
    n_in = meta["n_input_lines"]
    absent = [g for g in range(n_in) if g not in ledger]
    if absent or len(ledger) != n_in:
        problems.append(f"ledger: {len(absent)} input line(s) unaccounted, e.g. gidx {absent[:5]}")

    if problems:
        print("collect REFUSED -- nothing written:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        return 1

    # ---- write -----------------------------------------------------------------------
    verify = []
    for r in sorted(library, key=lambda r: r["name"]):
        s = states[r["name"]]
        ok = s["status"] == "OK"
        verify.append({"protomer_name": r["name"], "parent": r["parent_id"], "expected_charge": r["net_charge"],
                       "status": s["status"], "detail": s["detail"],
                       "db2_file": f"{int(s['shard']):05d}/{r['name']}.db2.gz" if ok else "",
                       "db2_name": r["name"] if ok else "", "amsol_charge": s["amsol_charge"]})
        if str(s["net_charge"]) != str(r["net_charge"]):
            print(f"WARNING: {r['name']} charge {s['net_charge']} in pass 2, {r['net_charge']} in library.tsv",
                  file=sys.stderr)
    write_rows(run / "db2_verification.csv", VERIFY_COLS, verify, delimiter=",")
    ledger_rows = [ledger[g] for g in sorted(ledger)]
    write_rows(run / "parents_ledger.tsv", LEDGER_COLS, ledger_rows)
    tcols = ["name", "shard", "how", "seconds_wall", "embedding", "amsol", "sampling",
             "db2_conversion", "msani_total"]
    write_rows(run / "timing.tsv", tcols, sorted(timing, key=lambda t: t["name"]))
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
    per_parent_p1 = [float(p["seconds_2d"]) for p in p1.values()]
    by_file = {}
    for inp in meta["inputs"]:
        lo, n = inp["first_gidx"], inp["n_lines"]
        rows = [ledger[g] for g in range(lo, lo + n)]
        by_file[Path(inp["path"]).name] = {"lines": n, "built": sum(r["status"] == "built" for r in rows),
                                           "names_ok": sum(int(r.get("n_ok") or 0) for r in rows)}
    ok_states = [s for s in states.values() if s["status"] == "OK"]
    n_built = sum(r["status"] == "built" for r in ledger_rows)
    cpu_h = {n: round(sum(d["seconds"] for d in dones[n]) / 3600, 2) for n in (1, 2)}
    task_s = {n: [d["seconds"] for d in dones[n]] for n in (1, 2)}
    names_stats = json.loads((run / "library" / "stats.json").read_text())
    stats = {
        "input_lines": n_in,
        "parents_built": n_built,
        "parents_failed": n_in - n_built,
        "parent_fail_reasons": dict(Counter(reason_key(r["reason"])
                                            for r in ledger_rows if r["status"] == "failed").most_common()),
        "names": len(library),
        "names_by_status": dict(Counter(s["status"] for s in states.values()).most_common()),
        "naming_dropped": names_stats.get("dropped", {}),
        "ok_names_per_built_parent": round(len(ok_states) / max(1, n_built), 3),
        "n_h_invertomers": sum(1 for r in library if r["invertomer_of"]),
        "net_charge": dict(sorted(Counter(int(r["net_charge"]) for r in library).items())),
        "by_input_file": by_file,
        "db2_gz_bytes_total": sum(int(s["db2_bytes"] or 0) for s in ok_states),
        "db2_gz_bytes_median": statistics.median([int(s["db2_bytes"] or 0) for s in ok_states]) if ok_states else None,
        "db2_records_median": statistics.median([int(s["db2_records"] or 0) for s in ok_states]) if ok_states else None,
        "timing": {
            "pass1_seconds_per_parent_median": round(statistics.median(per_parent_p1), 3) if per_parent_p1 else None,
            "pass1_seconds_per_parent_p90": pct(per_parent_p1, 0.9),
            "pass2_names_timed": len(per_name),
            "pass2_seconds_per_name_median": round(statistics.median(per_name), 2) if per_name else None,
            "pass2_seconds_per_name_p90": pct(per_name, 0.9),
            "pass2_seconds_per_name_max": round(max(per_name), 1) if per_name else None,
            "task_seconds_median": {n: round(statistics.median(task_s[n]), 1) if task_s[n] else None for n in (1, 2)},
            "task_seconds_max": {n: round(max(task_s[n]), 1) if task_s[n] else None for n in (1, 2)},
            "cpu_hours": cpu_h,
            "cpu_hours_per_1000_parents": round((cpu_h[1] + cpu_h[2]) / max(1, len(p1)) * 1000, 2),
        },
    }
    (run / "collect.json").write_text(json.dumps(stats, indent=2) + "\n")
    alld = dones[1] + dones[2]
    prov = {
        "run": {k: meta[k] for k in ("created", "label", "inputs", "per_shard", "n_shards", "options", "llb")},
        "pass2_plan": json.loads((run / "pass2" / "plan.json").read_text()),
        "software": json.loads(next(iter(stacks))) if stacks else {},
        "msani_path": sorted({d["versions"].get("msani_path", "") for d in alld}),
        "hosts": len({d["host"] for d in alld}),
        "cpus": dict(Counter(d["cpu"] for d in alld)),
        "slurm_jobs": {n: sorted({d["slurm_job"] for d in dones[n] if d["slurm_job"]}) for n in (1, 2)},
        "first_start": min(d["started"] for d in alld) if alld else "",
        "last_finish": max(d["finished"] for d in alld) if alld else "",
    }
    (run / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n")

    t = stats["timing"]
    print(f"{n_in:,} input lines: {n_built:,} parents built, {n_in - n_built:,} failed "
          f"(reasons in parents_ledger.tsv)")
    print(f"{len(library):,} names: {stats['names_by_status']}; "
          f"{stats['ok_names_per_built_parent']} OK names per built parent")
    print(f"timing: pass 1 median {t['pass1_seconds_per_parent_median']} s/parent, "
          f"pass 2 median {t['pass2_seconds_per_name_median']} s/name; CPU-h {cpu_h}, "
          f"{t['cpu_hours_per_1000_parents']} CPU-h per 1,000 parents")
    print(f"wrote {run}/db2_verification.csv, parents_ledger.tsv, timing.tsv, collect.json, provenance.json")

    if args.drap_dir:
        v = args.drap_dir / "stage4" / "validate_names.py"
        smi = run / "library" / "library.smi"
        print(f"\n$ python3 {v} --smi {smi}")
        return subprocess.run([sys.executable, str(v), "--smi", str(smi)]).returncode
    return 0


if __name__ == "__main__":
    sys.exit(main())
