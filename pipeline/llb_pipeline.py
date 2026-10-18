#!/usr/bin/env python3
"""
llb_pipeline.py -- SMILES to verified db2s, stage by stage, with every loss recorded.

    llb pipeline start SMILES... -o RUN [--validate] [--set key=value ...] [--here]
    llb pipeline status RUN
    llb pipeline advance RUN          (what each control job runs; safe by hand)
    llb pipeline resume RUN [--set key=value ...]
    llb pipeline cancel RUN

Stages, in order (docs/PIPELINE_DESIGN.md):
    parents  triage  qupkake  merge  names  checks  db2  verify
Inline stages run inside a control job; qupkake, names and db2 are SLURM jobs, each
followed by a control job (`--dependency=afterany`) that checks the outputs, applies
the stage's stop threshold and submits the next stage. A task that failed or timed
out is resubmitted once at twice the time limit; a threshold that is exceeded stops
the run (the control job exits 3, so SLURM mails LLB_MAIL). `resume` re-checks a
stopped stage after you have looked, optionally with a changed threshold.

Records in RUN: pipeline.json (inputs, settings, versions; written once),
pipeline/state.json, pipeline/logs/, ledger.tsv, attrition.tsv, REPORT.md.

Stdlib only; each step runs through `llb`, in its own environment.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import datetime as dt
import fcntl
import hashlib
import json
import math
import os
import re
import shlex
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

LLB_ROOT = Path(__file__).resolve().parent.parent
LLB = str(LLB_ROOT / "bin" / "llb")

STAGES = ["parents", "triage", "qupkake", "merge", "names", "checks", "db2", "verify"]
JOB_STAGES = {"qupkake", "names", "db2"}

DEFAULTS = {
    "validate": False,
    "shard_size": 10, "preflight": 10,
    "qupkake_time": "00:30:00",
    "names_cpus": 16, "names_mem": "32G", "names_time": "03:00:00",
    "db2_ncpus": 8, "db2_names_per_task": 25, "db2_time": "01:00:00", "db2_mmff": "restrained",
    "control_time": "04:00:00", "control_mem": "8G",
    "max_array_tasks": 4500,
    # stop thresholds: the loss fraction above which a stage halts the run instead
    # of recording and continuing (about twice what the drug-like 10k measured)
    "stop": {"parents": 0.01, "triage": 0.01, "qupkake": 0.10, "names": 0.02,
             "checks": 0.005, "db2": 0.05, "end_to_end": 0.05},
}
# parents rejected for these reasons are losses; duplicates map onto a kept parent
PARENT_LOSSES = {"no_name", "malformed_zinc_id", "unparsable", "id_conflict", "hash_collision"}
CONTRACT = re.compile(r"^C[0-9A-Za-z]{14}C$")


class Stop(Exception):
    """A stage cannot continue: (phase, message). phase: run | collect | gate."""

    def __init__(self, phase: str, message: str):
        super().__init__(message)
        self.phase, self.message = phase, message


# ---------------------------------------------------------------------------
# the outside world: llb steps, sbatch, squeue, sacct (replaced in tests)
# ---------------------------------------------------------------------------
class Runner:
    def __init__(self, cwd: Path, log: Path):
        self.cwd, self.log = cwd, log

    def llb(self, *args: str) -> int:
        """Run an llb step in the run directory, output appended to the stage log."""
        self.log.parent.mkdir(parents=True, exist_ok=True)
        with open(self.log, "a") as fh:
            fh.write(f"\n$ llb {shlex.join(args)}\n")
            fh.flush()
            return subprocess.run([LLB, *args], cwd=self.cwd, stdout=fh,
                                  stderr=subprocess.STDOUT).returncode

    def llb_capture(self, *args: str) -> tuple[int, str]:
        res = subprocess.run([LLB, *args], cwd=self.cwd, capture_output=True, text=True)
        self.log.parent.mkdir(parents=True, exist_ok=True)
        with open(self.log, "a") as fh:
            fh.write(f"\n$ llb {shlex.join(args)}\n{res.stdout}{res.stderr}")
        return res.returncode, res.stdout + res.stderr

    def sbatch(self, args: list[str]) -> str:
        res = subprocess.run(["sbatch", "--parsable", *args], cwd=self.cwd,
                             capture_output=True, text=True)
        if res.returncode != 0:
            raise Stop("run", f"sbatch failed: {res.stderr.strip()}")
        return res.stdout.strip().split(";")[0]

    def active(self, ids: list[str]) -> bool:
        if not ids:
            return False
        res = subprocess.run(["squeue", "-h", "-j", ",".join(ids), "-o", "%i"],
                             capture_output=True, text=True)
        return bool(res.stdout.strip())

    def failed_tasks(self, ids: list[str]) -> int:
        """Array tasks / jobs that ended in anything but COMPLETED."""
        if not ids:
            return 0
        res = subprocess.run(["sacct", "-n", "-X", "-P", "-j", ",".join(ids),
                              "--format=JobID,State"], capture_output=True, text=True)
        return sum(1 for line in res.stdout.splitlines()
                   if line.strip() and not line.split("|")[1].startswith("COMPLETED"))

    def failed_task_indices(self, ids: list[str]) -> list[int]:
        """Array task indices (JOBID_N) that ended in anything but COMPLETED."""
        if not ids:
            return []
        res = subprocess.run(["sacct", "-n", "-X", "-P", "-j", ",".join(ids),
                              "--format=JobID,State"], capture_output=True, text=True)
        out = []
        for line in res.stdout.splitlines():
            jid, _, state = line.partition("|")
            if "_" in jid and not state.startswith("COMPLETED"):
                idx = jid.split("_", 1)[1]
                if idx.isdigit():
                    out.append(int(idx))
        return sorted(set(out))

    def cancel(self, ids: list[str]) -> None:
        if ids:
            subprocess.run(["scancel", *ids])

    def config(self) -> dict:
        keys = "LLB_ACCOUNT LLB_PARTITION LLB_MAIL LLB_ARRAY_PARALLEL DB2C_SRC LLB_CONFIG"
        script = (f'source "{LLB_ROOT}/common/env.sh"; llb_load_config >/dev/null 2>&1; '
                  + "; ".join(f'printf "%s\\t%s\\n" {k} "${{{k}:-}}"' for k in keys.split()))
        res = subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                             env={**os.environ, "LLB_ROOT": str(LLB_ROOT)})
        return dict(l.split("\t", 1) for l in res.stdout.splitlines() if "\t" in l)


SUBMITTED = re.compile(r"Submitted batch job (\d+)")


# ---------------------------------------------------------------------------
# the run directory
# ---------------------------------------------------------------------------
class Run:
    def __init__(self, path: Path, runner: Runner | None = None):
        self.path = path.resolve()
        self.pdir = self.path / "pipeline"
        self.spec = json.loads((self.path / "pipeline.json").read_text())
        self.settings = self.spec["settings"]
        self.state = json.loads((self.pdir / "state.json").read_text())
        self.runner = runner or Runner(self.path, self.pdir / "logs" / "steps.log")
        # thresholds and times changed by `resume --set` apply from then on
        for o in self.state.get("overrides", []):
            set_setting(self.settings, o["key"], o["value"])

    @property
    def name(self) -> str:
        return self.path.name

    def save(self) -> None:
        tmp = self.pdir / "state.json.tmp"
        tmp.write_text(json.dumps(self.state, indent=2))
        tmp.replace(self.pdir / "state.json")

    @contextlib.contextmanager
    def lock(self):
        with open(self.pdir / "lock", "w") as fh:
            try:
                fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise SystemExit("another `llb pipeline advance` is running on this run")
            yield

    def stage(self, name: str) -> dict:
        return self.state["stages"][name]

    def f(self, *parts) -> Path:
        return self.path.joinpath(*parts)


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def now() -> str:
    return dt.datetime.now().isoformat(timespec="seconds")


def code_version() -> str:
    """The checkout's commit now (code can be pulled mid-run), '+dirty' if modified."""
    head = subprocess.run(["git", "-C", str(LLB_ROOT), "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True).stdout.strip() or "unknown"
    dirty = subprocess.run(["git", "-C", str(LLB_ROOT), "status", "--porcelain", "--untracked-files=no"],
                           capture_output=True, text=True).stdout.strip()
    return head + ("+dirty" if dirty else "")


def read_tsv(p: Path) -> list[dict]:
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_csv(p: Path) -> list[dict]:
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def double_time(t: str) -> str:
    """'HH:MM:SS' or 'D-HH:MM:SS' -> twice as long, same format."""
    days, _, rest = t.rpartition("-")
    h, m, s = (int(x) for x in rest.split(":"))
    total = 2 * ((int(days or 0) * 24 + h) * 3600 + m * 60 + s)
    d, rem = divmod(total, 86400)
    hh, rem = divmod(rem, 3600)
    mm, ss = divmod(rem, 60)
    return (f"{d}-" if d else "") + f"{hh:02d}:{mm:02d}:{ss:02d}"


def set_setting(settings: dict, key: str, value) -> None:
    """key 'stop.db2' or 'db2_time'; value cast to the default's type."""
    if key.startswith("stop."):
        k = key.split(".", 1)[1]
        if k not in DEFAULTS["stop"]:
            raise SystemExit(f"no stop threshold '{k}' (have {sorted(DEFAULTS['stop'])})")
        settings["stop"][k] = float(value)
        return
    if key not in DEFAULTS or key == "stop":
        raise SystemExit(f"unknown setting '{key}' (have {sorted(k for k in DEFAULTS if k != 'stop')})")
    kind = type(DEFAULTS[key])
    settings[key] = (str(value).lower() in ("1", "true", "yes")) if kind is bool else kind(value)


def frac(lost: int, total: int) -> float:
    return lost / total if total else 0.0


def gate(run: Run, key: str, lost: int, total: int, what: str) -> dict:
    """Record a stage's loss; raise Stop when it is above the threshold."""
    f, limit = frac(lost, total), run.settings["stop"][key]
    result = {"lost": lost, "total": total, "fraction": round(f, 5), "threshold": limit, "what": what}
    if f > limit:
        raise Stop("gate", f"{what}: {lost:,} of {total:,} ({f:.2%}) is above the "
                           f"{limit:.2%} stop threshold for {key}. Look at it, then "
                           f"`llb pipeline resume {run.path}` (add --set stop.{key}=X to accept more)")
    return result


def ranges(idx: list[int]) -> str:
    """[3, 6, 7, 8, 17] -> '3,6-8,17' (a SLURM --array list)."""
    out, idx = [], sorted(set(idx))
    for i in idx:
        if out and i == out[-1][1] + 1:
            out[-1][1] = i
        else:
            out.append([i, i])
    return ",".join(str(a) if a == b else f"{a}-{b}" for a, b in out)


def qupkake_input(run: Run) -> str | None:
    """Parents directory QupKake runs on: all parents (--validate) or the routed set."""
    if run.settings["validate"]:
        return "parents"
    return "triage/routed" if run.f("triage", "routed", "shards.tsv").exists() else None


def shard_ids(run: Run, parents_dir: str) -> list[str]:
    return [f"shard_{int(r['index']):05d}" for r in read_tsv(run.f(parents_dir, "shards.tsv"))]


def parents_in(run: Run, parents_dir: str) -> int:
    return sum(int(r["n_parents"]) for r in read_tsv(run.f(parents_dir, "shards.tsv")))


# ---------------------------------------------------------------------------
# stages: inline ones run here; job ones have submit() and collect()
# ---------------------------------------------------------------------------
def do_parents(run: Run) -> None:
    s = run.settings
    if run.f("parents", "shards.tsv.meta").exists():       # finished before an interruption
        return
    if run.f("parents", "shards.tsv").exists():
        raise Stop("run", "parents/ is half-written (a control job ended during llb parents); "
                          "remove parents/, then `llb pipeline resume`")
    if run.runner.llb("parents", *run.spec["inputs_abs"], "-o", "parents",
                      "--shard-size", str(s["shard_size"]), "--preflight", str(s["preflight"])):
        raise Stop("run", "llb parents failed (pipeline/logs/steps.log)")


def check_parents(run: Run) -> dict:
    kept = read_tsv(run.f("parents", "parents.tsv"))
    rejected = read_tsv(run.f("parents", "rejected.tsv"))
    lost = sum(1 for r in rejected if r["reason"] in PARENT_LOSSES)
    return gate(run, "parents", lost, len(kept) + len(rejected), "inputs rejected")


def do_triage(run: Run) -> None:
    if run.f("triage", "stats.json").exists():              # finished before an interruption
        return
    if run.f("triage", "routed", "shards.tsv").exists():
        raise Stop("run", "triage/routed/ was written but stats.json was not (a control job "
                          "ended during the finish step); remove triage/routed/, then resume")
    if run.runner.llb("triage", "parents", "-o", "triage"):
        raise Stop("run", "llb triage failed (pipeline/logs/steps.log)")


def check_triage(run: Run) -> dict:
    dec = read_tsv(run.f("triage", "decisions.tsv"))
    lost = sum(1 for d in dec if d["reason"] in ("molgpka_failed", "unparsable"))
    return gate(run, "triage", lost, len(dec), "parents MolGpKa failed on")


def submit_qupkake(run: Run, retry: bool) -> list[str]:
    src = qupkake_input(run)
    if src is None:
        return []
    n = len(shard_ids(run, src))
    if n > run.settings["max_array_tasks"]:
        raise Stop("run", f"{n:,} QupKake shards is more than max_array_tasks "
                          f"({run.settings['max_array_tasks']:,}); waves are not built yet -- "
                          f"split the input into smaller libraries")
    t = run.settings["qupkake_time"]
    # a retry submits only the shards with no DONE, not the whole array again
    extra = ["--only-missing"] if retry else []
    rc, out = run.runner.llb_capture("submit-qupkake", src, "protomers_qk", *extra,
                                     f"--time={double_time(t) if retry else t}")
    ids = SUBMITTED.findall(out)
    if retry and not rc and "nothing to submit" in out:
        return []
    if rc or not ids:
        raise Stop("run", "llb submit-qupkake failed (pipeline/logs/steps.log)")
    return ids


def collect_qupkake(run: Run) -> list[str]:
    src = qupkake_input(run)
    if src is None:
        return []
    return [sid for sid in shard_ids(run, src) if not run.f("protomers_qk", sid, "DONE").exists()]


def gate_qupkake(run: Run) -> dict:
    src = qupkake_input(run)
    if src is None:
        return {"note": "nothing routed to QupKake"}
    failed = sum(len(read_tsv(run.f("protomers_qk", sid, "qupkake_failed.tsv")))
                 for sid in shard_ids(run, src))
    return gate(run, "qupkake", failed, parents_in(run, src),
                "routed parents QupKake failed on (MolGpKa's states used instead)")


def do_merge(run: Run) -> None:
    src = qupkake_input(run)
    r = run.runner
    if src is None:
        if r.llb("merge", "--parents", "parents", "--triage", "triage", "-o", "merged"):
            raise Stop("run", "llb merge failed (pipeline/logs/steps.log)")
        return
    if r.llb("reenumerate", src, "protomers_qk", "protomers_qk_h", "--molgpka-triage", "triage"):
        raise Stop("run", "llb reenumerate failed (pipeline/logs/steps.log)")
    extra = ["--qupkake-parents", "parents"] if run.settings["validate"] else []
    if r.llb("merge", "--parents", "parents", "--triage", "triage",
             "--qupkake", "protomers_qk_h", *extra, "-o", "merged"):
        raise Stop("run", "llb merge refused or failed (pipeline/logs/steps.log)")
    if run.settings["validate"]:
        r.llb("compare", "--parents", "parents", "--triage", "triage",
              "--qupkake", "protomers_qk_h", "-o", "compare")


def check_merge(run: Run) -> dict:
    if qupkake_input(run) is None:
        return {"note": "nothing routed; MolGpKa states only"}
    swaps = sum(len({row["name"] for row in read_tsv(p)})
                for p in run.f("protomers_qk_h").glob("shard_*/hybrid_sites.tsv"))
    return {"aromatic_n_swapped_parents": swaps}


def submit_names(run: Run, retry: bool) -> list[str]:
    s, c = run.settings, run.state["cluster"]
    t = double_time(s["names_time"]) if retry else s["names_time"]
    return [run.runner.sbatch([
        "-A", c["LLB_ACCOUNT"], "-p", c["LLB_PARTITION"], "-c", str(s["names_cpus"]),
        "--mem", s["names_mem"], "-t", t, "-J", f"llb-{run.name}-names",
        "-o", str(run.pdir / "logs" / "names.%j.log"),
        "--wrap", f"{LLB} names --shards-dir merged --parents parents -o library"])]


def collect_names(run: Run) -> list[str]:
    st = run.stage("names")
    # written by this attempt (5 min slack: compute-node clocks and NFS times drift)
    done = run.f("library", "stats.json").exists() and \
        run.f("library", "stats.json").stat().st_mtime >= st.get("submitted_epoch", 0) - 300
    return [] if done and not run.runner.failed_tasks(st.get("attempt", [])) else ["names"]


def gate_names(run: Run) -> dict:
    kept = {r["parent_id"] for r in read_tsv(run.f("parents", "parents.tsv"))}
    named = {r["parent_id"] for r in read_tsv(run.f("library", "library.tsv"))}
    return gate(run, "names", len(kept - named), len(kept), "parents with no name")


def do_checks(run: Run) -> None:
    r = run.runner
    if r.llb("validate-names", "--smi", "library/library.smi"):
        raise Stop("run", "llb validate-names failed: the name contract is broken")
    r.llb("check-stereo", "library/library.smi")      # exit 1 when names would split


def check_checks(run: Run) -> dict:
    split = {x["name"] for x in read_tsv(run.f("library", "library.stereo_split.tsv"))}
    crashed = {x["name"] for x in read_tsv(run.f("library", "library.stereo_crashed.tsv"))}
    if not run.f("library", "library.stereo_crashed.tsv").exists():
        raise Stop("run", "llb check-stereo wrote no results (pipeline/logs/steps.log)")
    excluded = split | crashed
    lines = run.f("library", "library.smi").read_text().splitlines()
    with open(run.f("library", "build.smi"), "w") as fh:
        fh.writelines(l + "\n" for l in lines if l.split()[1] not in excluded)
    with open(run.f("library", "excluded.tsv"), "w") as fh:
        fh.write("name\treason\n")
        fh.writelines(f"{n}\tcheckstereo_would_split\n" for n in sorted(split))
        fh.writelines(f"{n}\tstereo_enumeration_crash\n" for n in sorted(crashed - split))
    res = gate(run, "checks", len(excluded), len(lines), "names excluded by check-stereo")
    res.update(split=len(split), crashed=len(crashed))
    return res


def submit_db2(run: Run, retry: bool) -> list[str]:
    s = run.settings
    n = len(run.f("library", "build.smi").read_text().splitlines())
    nshards = max(1, math.ceil(n / s["db2_names_per_task"]))
    tasks = math.ceil(nshards / s["db2_ncpus"])
    if tasks > s["max_array_tasks"]:
        raise Stop("run", f"{tasks:,} db2 tasks is more than max_array_tasks; split the library")
    t = double_time(s["db2_time"]) if retry else s["db2_time"]
    extra = []
    if retry:   # only the tasks that failed or timed out, not the whole array again
        failed = run.runner.failed_task_indices(run.stage("db2").get("attempt", []))
        if failed:
            extra = ["--tasks", ranges(failed)]
    rc, out = run.runner.llb_capture("submit-db2", "db2_run", "library/build.smi",
                                     "--ncpus", str(s["db2_ncpus"]), "--nshards", str(nshards),
                                     "--mmff", s["db2_mmff"], "--time", t, *extra)
    ids = SUBMITTED.findall(out)
    if rc or not ids:
        raise Stop("run", "llb submit-db2 failed (pipeline/logs/steps.log)")
    return ids


def collect_db2(run: Run) -> list[str]:
    n = run.runner.failed_tasks(run.stage("db2").get("attempt", []))
    return [f"{n} task(s)"] if n else []


def do_verify(run: Run) -> None:
    r = run.runner
    r.llb("verify-db2", "library/library.tsv", "--outputpath", "db2_run/out",
          "-o", "db2_verification.csv")                 # exit 1 for excluded names too
    rc, timing = r.llb_capture("timing", "triage", "protomers_qk", "db2_run",
                               "--library", "library/library.tsv")
    run.f("timing.txt").write_text(timing)


def check_verify(run: Run) -> dict:
    rows = read_csv(run.f("db2_verification.csv"))
    if not rows:
        raise Stop("run", "llb verify-db2 wrote no report (pipeline/logs/steps.log)")
    excluded = {x["name"] for x in read_tsv(run.f("library", "excluded.tsv"))}
    bad = Counter(x["status"] for x in rows
                  if x["status"] in ("CHARGE_MISMATCH", "NAME_MISMATCH", "RENAMED", "UNREADABLE"))
    silent = [x["protomer_name"] for x in rows
              if x["status"] == "MISSING" and x["protomer_name"] not in excluded]
    if bad or silent:
        raise Stop("gate", f"verify-db2 found {dict(bad)} and {len(silent)} name(s) missing with "
                           f"no reason (e.g. {silent[:3]}): see db2_verification.csv")
    ledger, attrition = build_ledger(run)
    named = {x["parent_id"] for x in read_tsv(run.f("library", "library.tsv"))}
    built = {x["parent_id"] for x in ledger if x["fate"] == "built"}
    res = {"db2": gate(run, "db2", len(named - built), len(named),
                       "named parents with nothing built")}
    kept = [x for x in ledger if x["parent_status"] == "kept"]
    res["end_to_end"] = gate(run, "end_to_end", sum(1 for x in kept if x["fate"] != "built"),
                             len(kept), "input parents with no db2")
    return res


# ---------------------------------------------------------------------------
# ledger, attrition, report
# ---------------------------------------------------------------------------
LEDGER_COLUMNS = ["input_name", "input_smiles", "source", "heavy_atoms", "parent_id",
                  "parent_status", "tautomer_fix", "route", "route_reason", "site_source",
                  "aromatic_n_swapped", "n_protomers", "n_names", "n_stereo_dropped",
                  "n_excluded_checks", "n_built", "n_build_failed", "build_fail_reasons",
                  "fate", "loss_reason"]


def build_ledger(run: Run) -> tuple[list[dict], list[dict]]:
    """One row per input molecule (kept or rejected parent), written to ledger.tsv;
    stage-by-stage counts by heavy-atom bin to attrition.tsv."""
    dec = {d["parent_id"]: d for d in read_tsv(run.f("triage", "decisions.tsv"))}
    lib = read_tsv(run.f("library", "library.tsv"))
    by_parent = defaultdict(list)
    for x in lib:
        by_parent[x["parent_id"]].append(x)
    dropped = Counter(d["parent_id"] for d in read_tsv(run.f("library", "dropped.tsv"))
                      if d.get("reason") == "stereo_cannot_embed")
    excluded = {x["name"] for x in read_tsv(run.f("library", "excluded.tsv"))}
    ver = {x["protomer_name"]: x for x in read_csv(run.f("db2_verification.csv"))}
    swapped = {row["name"] for p in run.f("protomers_qk_h").glob("shard_*/hybrid_sites.tsv")
               for row in read_tsv(p)}
    ledger = []
    for p in read_tsv(run.f("parents", "parents.tsv")):
        pid, names = p["parent_id"], by_parent.get(p["parent_id"], [])
        status = [ver.get(n["name"], {}).get("status", "") for n in names]
        fails = Counter(ver[n["name"]]["detail"] for n in names
                        if ver.get(n["name"], {}).get("status") == "BUILD_FAILED")
        n_built = status.count("OK")
        if n_built:
            fate, why = "built", ""
        elif not names:
            fate, why = "lost_at_names", "no buildable stereoisomer" if dropped[pid] else "no name"
        elif all(n["name"] in excluded for n in names):
            fate, why = "lost_at_checks", "every name excluded by check-stereo"
        else:
            fate, why = "lost_at_db2", ";".join(sorted(fails)) or "not built"
        d = dec.get(pid, {})
        ledger.append({
            "input_name": p["input_name"], "input_smiles": p["input_smiles"],
            "source": p["source"], "heavy_atoms": p["heavy_atoms"], "parent_id": pid,
            "parent_status": "kept", "tautomer_fix": p.get("tautomer_fix", ""),
            "route": d.get("route", ""), "route_reason": d.get("reason", ""),
            "site_source": names[0].get("site_source", "") if names else "",
            "aromatic_n_swapped": int(pid in swapped),
            "n_protomers": len({n["protomer_index"] for n in names}), "n_names": len(names),
            "n_stereo_dropped": dropped[pid],
            "n_excluded_checks": sum(n["name"] in excluded for n in names),
            "n_built": n_built, "n_build_failed": sum(fails.values()),
            "build_fail_reasons": ";".join(f"{k}:{v}" for k, v in sorted(fails.items())),
            "fate": fate, "loss_reason": why})
    for rj in read_tsv(run.f("parents", "rejected.tsv")):
        ledger.append({**{c: "" for c in LEDGER_COLUMNS}, "input_name": rj["input_name"],
                       "parent_status": rj["reason"], "fate":
                       "duplicate" if rj["reason"].startswith("duplicate") else "lost_at_parents",
                       "loss_reason": rj["detail"][:200]})
    with open(run.f("ledger.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LEDGER_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(ledger)

    order = ["inputs", "parents", "named", "passed_checks", "built"]
    bins = defaultdict(Counter)
    for x in ledger:
        b = x["heavy_atoms"] or "NA"
        if x["fate"] == "duplicate":
            continue
        bins[b]["inputs"] += 1
        if x["parent_status"] == "kept":
            bins[b]["parents"] += 1
            bins[b]["named"] += int(x["n_names"] > 0)
            bins[b]["passed_checks"] += int(x["n_names"] > x["n_excluded_checks"])
            bins[b]["built"] += int(x["fate"] == "built")
    attrition = [{"heavy_atoms": b, **{k: bins[b][k] for k in order}}
                 for b in sorted(bins, key=lambda v: (v == "NA", int(v) if v.isdigit() else 0))]
    attrition.append({"heavy_atoms": "all", **{k: sum(r[k] for r in attrition) for k in order}})
    with open(run.f("attrition.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["heavy_atoms", *order], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(attrition)
    return ledger, attrition


def mol2_commands() -> str:
    return ("llb mol2 db2_run/out --library library/library.tsv -o library/library.mol2\n"
            "llb mol2-library library/library.mol2 --library library/library.tsv "
            "-o mol2_lib --shard-size 1000")


def write_report(run: Run) -> None:
    att = read_tsv(run.f("attrition.tsv"))
    lines = [f"# Library run {run.name}", "",
             f"Started {run.state['started']}, finished {run.state.get('finished', '')}. "
             f"LLB at start: {run.spec['versions'].get('llb_commit', '?')[:12]} "
             f"(per stage below).", "",
             "## Attrition (molecules, by heavy-atom bin)", "",
             "| " + " | ".join(att[0].keys()) + " |" if att else "",
             "| " + " | ".join("---" for _ in att[0]) + " |" if att else ""]
    lines += ["| " + " | ".join(r.values()) + " |" for r in att]
    lines += ["", "## Stages", "",
              "Code: the LLB commit each stage ran with (code pulled mid-run shows up here).", "",
              "| stage | status | code | jobs | retries | result |",
              "| --- | --- | --- | --- | --- | --- |"]
    for st in STAGES:
        s = run.stage(st)
        lines.append(f"| {st} | {s['status']} | {' '.join(s.get('commits', []))} | "
                     f"{' '.join(s.get('jobs', []))} | {s.get('retries', 0)} | "
                     f"{json.dumps(s.get('result', {}))} |")
    if run.state.get("overrides"):
        lines += ["", "## Overrides after a stop", ""]
        lines += [f"- {o['at']}: {o['key']} = {o['value']}" for o in run.state["overrides"]]
    lines += ["", "## Settings", "", "```json", json.dumps(run.spec["settings"], indent=2), "```",
              "", "## Versions", "", "```json", json.dumps(run.spec["versions"], indent=2), "```",
              "", "## Cost", "", "```", run.f("timing.txt").read_text() if run.f("timing.txt").exists()
              else "(no timing)", "```",
              "", "## Next: the mol2 library", "", "```bash", mol2_commands(), "```", ""]
    run.f("REPORT.md").write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# the state machine
# ---------------------------------------------------------------------------
INLINE = {"parents": (do_parents, check_parents), "triage": (do_triage, check_triage),
          "merge": (do_merge, check_merge), "checks": (do_checks, check_checks),
          "verify": (do_verify, check_verify)}
SUBMIT = {"qupkake": submit_qupkake, "names": submit_names, "db2": submit_db2}
COLLECT = {"qupkake": collect_qupkake, "names": collect_names, "db2": collect_db2}
GATE = {"qupkake": gate_qupkake, "names": gate_names, "db2": lambda run: {}}


def submit_control(run: Run, after: list[str], last: bool) -> str:
    c, s = run.state["cluster"], run.settings
    args = ["-A", c["LLB_ACCOUNT"], "-p", c["LLB_PARTITION"], "-c", "1",
            "--mem", s["control_mem"], "-t", s["control_time"],
            "-J", f"llb-{run.name}-control", "-o", str(run.pdir / "logs" / "control.%j.log")]
    if after:
        args.append(f"--dependency=afterany:{':'.join(after)}")
    if c.get("LLB_MAIL"):
        # TIME_LIMIT too: a control job that runs out of time ends the chain
        args += ["--mail-user", c["LLB_MAIL"],
                 f"--mail-type={'END,FAIL,TIME_LIMIT' if last else 'FAIL,TIME_LIMIT'}"]
    args += ["--wrap", f"{LLB} pipeline advance {shlex.quote(str(run.path))}"]
    return run.runner.sbatch(args)


def stop_stage(run: Run, stage: str, err: Stop) -> int:
    s = run.stage(stage)
    s.update(status="stopped", phase=err.phase, message=err.message, ended=now())
    run.save()
    print(f"STOPPED at {stage}: {err.message}", file=sys.stderr)
    return 3


def advance(run: Run) -> int:
    """Move the run forward as far as it can go now; idempotent."""
    commit = code_version()
    for stage in STAGES:
        s = run.stage(stage)
        if s["status"] == "done":
            continue
        if s["status"] == "stopped":
            print(f"stopped at {stage}: {s['message']}", file=sys.stderr)
            return 3
        try:
            if s["status"] == "check":                  # resumed after a threshold stop
                if commit not in s.setdefault("commits", []):
                    s["commits"].append(commit)
                check = INLINE[stage][1] if stage in INLINE else GATE[stage]
                s.update(result=check(run), status="done", ended=now())
                run.save()
                continue
            if s["status"] == "running":
                if run.runner.active(s["jobs"]):
                    print(f"{stage}: jobs {' '.join(s['jobs'])} still queued or running")
                    return 0
                missing = COLLECT[stage](run)
                if missing:
                    if s.get("retries", 0) >= 1:
                        raise Stop("collect", f"{len(missing)} unfinished after a retry at twice "
                                              f"the time, e.g. {missing[:5]}")
                    epoch = dt.datetime.now().timestamp()
                    ids = SUBMIT[stage](run, retry=True)
                    if not ids:
                        raise Stop("collect", f"{len(missing)} unfinished but nothing could be "
                                              f"resubmitted, e.g. {missing[:5]}")
                    s.update(jobs=s["jobs"] + ids, attempt=ids, retries=1, submitted_epoch=epoch)
                    s.setdefault("commits", []).append(commit)
                    s.setdefault("control", []).append(submit_control(run, ids, stage == "db2"))
                    run.save()
                    print(f"{stage}: {len(missing)} unfinished; resubmitted once at twice the time")
                    return 0
                s.update(result=GATE[stage](run), status="done", ended=now())
                if commit not in s.setdefault("commits", []):
                    s["commits"].append(commit)
                run.save()
                continue
            # pending, or in_progress: a control job ended during this inline stage
            # (time limit, node failure); each do_* step skips or refuses safely
            s["started"] = s.get("started") or now()
            if stage in INLINE:
                do, check = INLINE[stage]
                if s["status"] == "in_progress":
                    print(f"{stage}: the last control job ended during it; picking it up")
                s["status"] = "in_progress"
                if commit not in s.setdefault("commits", []):
                    s["commits"].append(commit)
                run.save()
                do(run)
                s.update(result=check(run), status="done", ended=now())
                run.save()
                continue
            epoch = dt.datetime.now().timestamp()
            ids = SUBMIT[stage](run, retry=False)
            if not ids:          # e.g. nothing routed to QupKake
                s.update(status="done", result=GATE[stage](run), ended=now())
                run.save()
                continue
            s.update(status="running", jobs=ids, attempt=ids, retries=0, submitted_epoch=epoch)
            if commit not in s.setdefault("commits", []):
                s["commits"].append(commit)
            s["control"] = [submit_control(run, ids, stage == "db2")]
            run.save()
            print(f"{stage}: submitted {' '.join(ids)}; a control job continues when it ends")
            return 0
        except Stop as err:
            return stop_stage(run, stage, err)
    run.state["finished"] = run.state.get("finished") or now()
    run.save()
    write_report(run)
    print(f"done: {run.f('REPORT.md')}\nmol2 library next:\n{mol2_commands()}")
    return 0


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------
def versions(runner: Runner, cluster: dict) -> dict:
    out = {}
    git = subprocess.run(["git", "-C", str(LLB_ROOT), "rev-parse", "HEAD"],
                         capture_output=True, text=True)
    out["llb_commit"] = git.stdout.strip()
    if cluster.get("DB2C_SRC"):
        g = subprocess.run(["git", "-C", cluster["DB2C_SRC"], "rev-parse", "HEAD"],
                           capture_output=True, text=True)
        out["db2_converter_commit"] = g.stdout.strip()
    for env in ("PREP", "PKA", "QUPKAKE", "DB2C"):
        rc, text = runner.llb_capture("env", env, "python", str(LLB_ROOT / "common" / "versions.py"))
        # the JSON line, wherever it is: `module load` banners can print after it
        found = None
        for line in reversed(text.splitlines()):
            line = line.strip()
            if line.startswith("{"):
                try:
                    found = json.loads(line)
                    break
                except ValueError:
                    continue
        out[env.lower()] = found if found is not None else {"error": text.strip()[-300:]}
    return out


def create_run(run_dir: Path, smiles: list[str], settings: dict, cluster: dict,
               vers: dict, runner: Runner | None = None) -> Run:
    """Write pipeline.json and an all-pending state; return the Run."""
    inputs = [Path(p).resolve() for p in smiles]
    for p in inputs:
        if not p.exists():
            raise SystemExit(f"no such input: {p}")
    (run_dir / "pipeline" / "logs").mkdir(parents=True, exist_ok=True)
    spec = {"created": now(), "inputs": [str(p) for p in smiles],
            "inputs_abs": [str(p) for p in inputs],
            "input_sha256": {str(p): sha256(p) for p in inputs},
            "settings": settings, "versions": vers,
            "cluster": {k: v for k, v in cluster.items() if k != "LLB_MAIL"}}
    (run_dir / "pipeline.json").write_text(json.dumps(spec, indent=2))
    state = {"started": now(), "cluster": cluster,
             "stages": {st: {"status": "pending"} for st in STAGES}}
    (run_dir / "pipeline" / "state.json").write_text(json.dumps(state, indent=2))
    return Run(run_dir, runner)


def cmd_start(a) -> int:
    run_dir = Path(a.outdir).resolve()
    if (run_dir / "pipeline.json").exists():
        raise SystemExit(f"{run_dir} already holds a run; `llb pipeline status {run_dir}` "
                         f"or choose a new -o")
    settings = json.loads(json.dumps(DEFAULTS))
    settings["validate"] = a.validate
    for kv in a.set or []:
        k, _, v = kv.partition("=")
        set_setting(settings, k, v)
    (run_dir / "pipeline" / "logs").mkdir(parents=True, exist_ok=True)
    runner = Runner(run_dir, run_dir / "pipeline" / "logs" / "steps.log")
    cluster = runner.config()
    if not cluster.get("LLB_ACCOUNT"):
        raise SystemExit("LLB_ACCOUNT is not set in your hpc.env")
    run = create_run(run_dir, a.smiles, settings, cluster, versions(runner, cluster))
    if a.here:
        with run.lock():
            return advance(run)
    jid = submit_control(run, [], last=False)
    run.state["first_control"] = jid
    run.save()
    print(f"run {run_dir}: control job {jid} runs parents and triage, then submits the rest.\n"
          f"follow it with: llb pipeline status {run_dir}")
    return 0


def cmd_status(a) -> int:
    run = Run(Path(a.run))
    print(f"run {run.path}  (started {run.state['started']})")
    print(f"{'stage':<9} {'status':<9} {'jobs':<24} result")
    for st in STAGES:
        s = run.stage(st)
        res = s.get("message") if s["status"] == "stopped" else json.dumps(s.get("result", ""))
        print(f"{st:<9} {s['status']:<9} {' '.join(s.get('jobs', []))[:24]:<24} {res}")
    if run.state.get("finished"):
        print(f"\nfinished {run.state['finished']}: {run.f('REPORT.md')}\n"
              f"mol2 library next:\n{mol2_commands()}")
    return 0


def cmd_advance(a) -> int:
    run = Run(Path(a.run))
    with run.lock():
        return advance(run)


def cmd_resume(a) -> int:
    run = Run(Path(a.run))
    with run.lock():
        for kv in a.set or []:
            k, _, v = kv.partition("=")
            set_setting(run.settings, k, v)          # validates the key
            run.state.setdefault("overrides", []).append({"at": now(), "key": k, "value": v})
        stopped = [st for st in STAGES if run.stage(st)["status"] == "stopped"]
        if not stopped:
            print("nothing is stopped")
            return 0
        st = stopped[0]
        s = run.stage(st)
        phase = s.get("phase")
        if phase == "gate":            # outputs are there; re-check with today's thresholds
            s["status"] = "check"
        elif phase == "collect":       # resubmit the unfinished work once more
            s.update(status="running", retries=0, attempt=[])
        else:                          # the step itself failed: rerun it
            s["status"] = "pending"
        s.pop("message", None)
        run.save()
        print(f"resuming {st} (stopped at {phase})")
    if a.here:
        return cmd_advance(a)
    jid = submit_control(run, [], last=False)
    print(f"control job {jid} picks the run up")
    return 0


def cmd_cancel(a) -> int:
    run = Run(Path(a.run))
    ids = [j for st in STAGES for j in run.stage(st).get("jobs", []) + run.stage(st).get("control", [])]
    run.runner.cancel(ids)
    print(f"cancelled {len(ids)} job(s)")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("start", help="create a run and submit its first control job")
    p.add_argument("smiles", nargs="+", help="'SMILES NAME' files (plain or .gz)")
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--validate", action="store_true",
                   help="QupKake on every parent too, plus llb compare (new chemical space)")
    p.add_argument("--set", action="append", metavar="KEY=VALUE",
                   help="e.g. qupkake_time=00:45:00, stop.db2=0.08")
    p.add_argument("--here", action="store_true",
                   help="run the first stages in this shell instead of a control job")
    p.set_defaults(fn=cmd_start)
    for name, fn, hlp in (("status", cmd_status, "where the run is"),
                          ("advance", cmd_advance, "one step forward (control jobs run this)"),
                          ("cancel", cmd_cancel, "scancel every job of the run")):
        q = sub.add_parser(name, help=hlp)
        q.add_argument("run")
        q.set_defaults(fn=fn)
    r = sub.add_parser("resume", help="continue after a stop")
    r.add_argument("run")
    r.add_argument("--set", action="append", metavar="KEY=VALUE")
    r.add_argument("--here", action="store_true")
    r.set_defaults(fn=cmd_resume)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
