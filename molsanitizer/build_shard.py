#!/usr/bin/env python3
"""
build_shard.py -- one array task of a MolSanitizer build, pass 1 or pass 2.

Runs inside the MolSanitizer env (see environment/). The build (UPSTREAM.md
section 2) is two MolSanitizer passes with LLB's naming step between them,
because a db2 keeps only the last 16 characters of a name and a ZINC22 id is
already 16:

  pass 1  chemistry only, SMILES out. MolSanitizer's Python API on one parent at
          a time, named by its INPUT name: salts, neutralize, tautomers,
          protonation at the run's pH (the API takes 7.4; the CLI only
          integers), stereoisomers. Writes MolSanitizer's own `INPUT_state.stereo`
          rows to pass1/out/NNNNN_clean.smi, plus a reason for every parent that
          returns nothing.
  names   `llb names` over every pass-1 file at once (msani.sh names), which
          writes library/library.smi with contract names.
  pass 2  3D + db2 only, on contract names: the msani CLI with --no-neutralize
          --no-stereoisomers --metal off -3d -f db2, so it neither changes a
          state nor renames it. One call per shard; molecules that call leaves
          without a db2 or an error (it died or timed out) are retried one at a
          time under mol_timeout, so one crash costs one molecule.

Every db2 is checked before it is kept (the guarantee `llb verify-db2` gave):
each record's M-line name is the requested name, its AMSOL total charge is the
formal charge within charge_tol, and it has conformations.

Outputs (K = shard):
  pass 1: RUN/pass1/out/K_clean.smi           MolSanitizer rows, `SMILES MS_NAME`
          RUN/pass1/tasks/K/parents.tsv       every shard line: ok / failed + reason
  pass 2: RUN/db2_run/out/K/<name>.db2.gz     one per built name (prepare.py's layout)
          RUN/pass2/tasks/K/states.tsv        every name: status, detail, db2 facts
          RUN/pass2/tasks/K/timing.tsv        MolSanitizer's step times + ours
  both:   RUN/passN/tasks/K/DONE.json         last; versions, host, seconds, input sha256

Resumable: a shard whose DONE.json matches its input and options is skipped. An
interrupted pass-2 shard keeps every db2.gz already written (each is written to a
temporary file and renamed, so a present one is whole).

    python build_shard.py --run RUN --pass 1 --shard 17
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

PARENT_COLS = ["gidx", "input_name", "parent_id", "input_smiles", "status", "reason",
               "n_rows", "seconds_2d"]
STATE_COLS = ["name", "net_charge", "status", "detail", "db2_records", "amsol_charge",
              "db2_bytes", "smiles"]
BATCH = int(os.environ.get("LLB_MSANI_BATCH", "20"))   # names per msani call in pass 2
TIMING_COLS = ["name", "how", "seconds_wall", "embedding", "amsol", "sampling",
               "db2_conversion", "msani_total"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def options_hash(opts: dict) -> str:
    return hashlib.sha256(json.dumps(opts, sort_keys=True).encode()).hexdigest()[:16]


def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def versions() -> dict:
    import numpy
    import pandas
    import rdkit
    import msani
    v = {"molsanitizer": msani.__version__, "rdkit": rdkit.__version__,
         "python": platform.python_version(), "numpy": numpy.__version__,
         "pandas": pandas.__version__, "msani_path": str(Path(msani.__file__).parent)}
    obabel = shutil.which("obabel")
    if obabel:
        r = subprocess.run([obabel, "-V"], capture_output=True, text=True)
        v["openbabel"] = r.stdout.strip()
    return v


def cpu_model() -> str:
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


# ---- pass 1: chemistry ----------------------------------------------------------

def pass1(rows, opts, work: Path, out_smi: Path) -> list[dict]:
    """Every parent through MolSanitizer's 2D steps; its rows -> out_smi."""
    from pandas import DataFrame
    from rdkit import Chem, RDLogger
    from msani.api import Msani
    RDLogger.DisableLog("rdApp.*")

    rejected = work / "rejected_2d.txt"
    metal_err = work / "metal_error_2d.smi"
    proc = Msani(removesalts=opts["removesalts"], neutralize=opts["neutralize"],
                 tautomers=opts["tautomers"], extended_tautomers=opts["extended_tautomers"],
                 protonation=opts["protonation"], pH=opts["pH"], pH_range=opts["pH_range"],
                 stereoisomers=True, max_stereoisomers=opts["max_stereoisomers"],
                 stereo_timeout=opts["stereo_timeout"], randomSeed=opts["random_seed"],
                 metal=opts["metal"], numcores=1)
    parents, lines = [], []
    with proc:
        for gidx, smi, input_name, pid in rows:
            t0 = time.perf_counter()
            p = {"gidx": gidx, "input_name": input_name, "parent_id": pid,
                 "input_smiles": smi, "status": "ok", "reason": "", "n_rows": 0}
            if Chem.MolFromSmiles(smi) is None:
                p.update(status="failed", reason="invalid_smiles")
            else:
                try:
                    # the INPUT name, not the parent id: `llb names` matches
                    # MolSanitizer's suffixes against the names in the sample files
                    df = proc.run(DataFrame({"smiles": [smi], "ids": [input_name]}),
                                  rejected_file=str(rejected), metal_error_file=str(metal_err))
                except Exception as e:  # one parent's crash is that parent's reason
                    p.update(status="failed", reason=f"msani_2d_error: {type(e).__name__}: {e}"[:300])
                    df = None
                if df is not None and df.empty:
                    p.update(status="failed",
                             reason="removed_in_pass1: " + grep_reason(input_name, rejected, metal_err))
                elif df is not None:
                    got = [(str(s), str(i)) for s, i in zip(df["smiles"], df["ids"])]
                    lines += [f"{s} {i}\n" for s, i in got]
                    p["n_rows"] = len(got)
            p["seconds_2d"] = round(time.perf_counter() - t0, 3)
            parents.append(p)
    tmp = out_smi.with_name(out_smi.name + ".tmp")
    tmp.write_text("".join(lines))
    tmp.replace(out_smi)
    return parents


def grep_reason(name: str, *files: Path) -> str:
    for f in files:
        if f.exists():
            for line in open(f):
                if name in line.split():
                    return line.strip().split(name, 1)[-1].strip()[:200] or f.name
    return "no reason recorded by MolSanitizer"


# ---- pass 2: 3D + db2 -------------------------------------------------------------

def msani_flags(opts: dict) -> list[str]:
    """The 3D-only msani flags: no 2D step may change or rename a named state."""
    return ["--prefix", "todo", "--no-neutralize", "--no-stereoisomers", "--metal", "off",
            "-3d", "-f", "db2", "--numconfs", str(opts["numconfs"]),
            "--energywindow", str(opts["energywindow"]), "--rmsd", str(opts["rmsd"]),
            "--nringconfs", str(opts["nringconfs"]), "--timeout", str(opts["embed_timeout"]),
            "--randomSeed", str(opts["random_seed"]), "-j", "1", "--timing"]


def run_msani3d(states, opts, work: Path, on_chunk, chunk: int) -> tuple[int | str, float]:
    """msani3d.py on `states` in a fresh `work`; on_chunk(a, b) after each finished chunk.

    Each chunk has a deadline of mol_timeout x its size; past it the process is killed
    and "timeout" returned. Returns (exit code or "timeout", wall seconds).
    """
    import queue
    import threading
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    todo = work / "todo.smi"
    todo.write_text("".join(f"{s['smiles']} {s['name']}\n" for s in states))
    cmd = [sys.executable, str(Path(__file__).resolve().parent / "msani3d.py"),
           "--todo", str(todo), "--chunk", str(chunk), "--", *msani_flags(opts)]
    t0 = time.perf_counter()
    log = open(work / "msani.stdout", "w")
    proc = subprocess.Popen(cmd, cwd=work, stdout=subprocess.PIPE, stderr=log, text=True)
    lines: "queue.Queue[str | None]" = queue.Queue()
    threading.Thread(target=lambda: ([lines.put(l) for l in proc.stdout], lines.put(None)),
                     daemon=True).start()
    start, rc = 0, None
    while start < len(states):
        n = min(chunk, len(states) - start)
        deadline = time.monotonic() + opts["mol_timeout"] * n
        line = ""
        while not line.startswith("CHUNK "):
            try:
                line = lines.get(timeout=max(0.0, deadline - time.monotonic()))
            except queue.Empty:
                proc.kill()
                proc.wait()
                rc = "timeout"
                break
            if line is None:
                break
            log.write(line)
        if rc == "timeout" or line is None:
            break
        _, a, b, _ = line.split()
        on_chunk(int(a), int(b))
        start = int(b)
    if rc is None:
        rc = proc.wait()
    log.close()
    return rc, time.perf_counter() - t0


def read_errors(work: Path) -> dict[str, str]:
    errs = {}
    f = work / "msani_error.err"
    if f.exists():
        for line in open(f):
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                errs[parts[1]] = (parts[2] if len(parts) > 2 else "failed, no reason given")[:300]
    return errs


def read_timing(work: Path) -> dict[str, dict]:
    f = work / "msani_timing.csv"
    out = {}
    if f.exists():
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh):
                out[r["Name"]] = {"embedding": r.get("Initial embedding", ""), "amsol": r.get("AMSOL", ""),
                                  "sampling": r.get("Torsional sampling", ""),
                                  "db2_conversion": r.get("DB2 conversion", ""), "msani_total": r.get("Total", "")}
    return out


def check_db2(path: Path, name: str, charge: int, tol: float) -> tuple[str, str, dict]:
    """(status, detail, header facts) for one db2 written by MolSanitizer."""
    heads, charges, nconfs = [], [], []
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            lines = fh.read().splitlines()
    else:
        lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        t = line.split()
        if len(t) >= 11 and t[0] == "M" and not t[1].lstrip("+-").replace(".", "").isdigit():
            heads.append(t[1])
            nconfs.append(int(t[6]))
            try:
                charges.append(float(lines[i + 1].split()[1]))
            except (IndexError, ValueError):
                charges.append(None)
    facts = {"db2_records": len(heads), "amsol_charge": charges[0] if charges else ""}
    if not heads:
        return "EMPTY", "no db2 record header", facts
    if any(h != name for h in heads):
        return "RENAMED", f"M line holds {sorted(set(heads) - {name})[:2]}", facts
    if any(q is None for q in charges):
        return "UNREADABLE", "no parsable charge line", facts
    if any(abs(q - charge) > tol for q in charges):
        return "CHARGE_MISMATCH", f"AMSOL {charges} vs formal {charge}", facts
    if any(n <= 0 for n in nconfs):
        return "EMPTY", "a record has no conformations", facts
    return "OK", "", facts


def keep_db2(src: Path, dest: Path) -> int:
    tmp = dest.with_name(dest.name + ".tmp")
    with open(src, "rb") as fi, gzip.open(tmp, "wb", compresslevel=6) as fo:
        shutil.copyfileobj(fi, fo)
    tmp.replace(dest)
    return dest.stat().st_size


def keep_log(work: Path, taskdir: Path, tag: str) -> None:
    for f in ("todo.log", "msani.stdout", "msani_error.err"):
        if (work / f).exists():
            shutil.copy2(work / f, taskdir / f"{tag}.{f}")


def pass2(states, opts, scratch: Path, outdir: Path, taskdir: Path) -> dict:
    """Build every state's db2; set each state's status. Returns per-name timing."""
    from pandas import DataFrame
    from rdkit import Chem, RDLogger
    from msani.filtering.filters import Filters
    RDLogger.DisableLog("rdApp.*")

    timing: dict[str, dict] = {}
    for s in states:
        s["net_charge"] = Chem.GetFormalCharge(Chem.MolFromSmiles(s["smiles"]))
    # what pass 2's gen3d would drop without a word (msani filters.remove_exotic_chem_to_db2)
    df = DataFrame({"smiles": [s["smiles"] for s in states], "ids": [s["name"] for s in states]})
    df["mol"] = df["smiles"].apply(Chem.MolFromSmiles)
    kept = set(Filters.remove_exotic_chem_to_db2(df)["ids"]) if len(df) else set()
    for s in states:
        if s["name"] not in kept:
            s.update(status="EXOTIC", detail="not db2-compatible (MolSanitizer remove_exotic_chem_to_db2)")
        elif (outdir / f"{s['name']}.db2.gz").exists():           # resumed shard: already kept
            gz = outdir / f"{s['name']}.db2.gz"
            status, detail, facts = check_db2(gz, s["name"], int(s["net_charge"]), opts["charge_tol"])
            if status == "OK":
                s.update(status="OK", detail="kept from an earlier attempt", db2_bytes=gz.stat().st_size, **facts)
                timing[s["name"]] = {"how": "earlier_attempt"}
            else:
                gz.unlink()                                         # rebuild it
    todo = [s for s in states if s["status"] == "pending"]

    def alone(s):
        rc, wall = run_msani3d([s], opts, scratch / "alone", lambda a, b: None, 1)
        harvest(scratch / "alone", [s], "alone", rc, wall)
        if rc == "timeout":
            s.update(status="BUILD_FAILED", detail=f"timeout after {opts['mol_timeout']} s alone")

    def harvest(work: Path, batch, how: str, rc, wall: float):
        errs, times = read_errors(work), read_timing(work)
        for s in batch:
            f = work / "db2" / f"{s['name']}.db2"
            t = {"how": how, **times.get(s["name"], {})}
            if how == "alone":
                t["seconds_wall"] = round(wall, 3)
            timing[s["name"]] = t
            if f.exists() and f.stat().st_size:
                status, detail, facts = check_db2(f, s["name"], int(s["net_charge"]), opts["charge_tol"])
                s.update(status=status, detail=detail, **facts)
                if status == "OK":
                    s["db2_bytes"] = keep_db2(f, outdir / f"{s['name']}.db2.gz")
            elif s["name"] in errs:
                s.update(status="BUILD_FAILED", detail=errs[s["name"]])
            elif how == "alone":
                s.update(status="BUILD_FAILED",
                         detail=f"msani exit {rc}, no db2 and no error recorded" if rc != 0
                         else "no db2 and no error recorded (see msani.stdout)")

    # One msani3d.py process for the shard (msani's conformer module preloaded into its
    # forkserver once, not imported per molecule); each chunk of BATCH names is checked
    # and kept as soon as it is reported, so a task killed at its time limit keeps all but
    # the chunk in flight. If the process dies or a chunk overruns, that chunk's unbuilt
    # names are retried alone and a new process takes the rest.
    t0 = time.perf_counter()
    rest = list(todo)
    while rest:
        def on_chunk(a, b, rest=rest):
            harvest(scratch / "batch", rest[a:b], "batch", 0, 0.0)
            ok = sum(s["status"] == "OK" for s in todo)
            print(f"  {len(todo) - len(rest) + b}/{len(todo)} names done, {ok} OK, "
                  f"{time.perf_counter() - t0:.0f} s", flush=True)
        rc, _ = run_msani3d(rest, opts, scratch / "batch", on_chunk, BATCH)
        keep_log(scratch / "batch", taskdir, "batch")
        if rc == 0:
            break
        tried = [s for s in rest if s["name"] in timing]
        in_flight = [s for s in rest if s["name"] not in timing][:BATCH]
        print(f"  msani3d stopped ({rc}); retrying {len(in_flight)} name(s) alone", flush=True)
        for s in in_flight:
            alone(s)
        rest = [s for s in rest if s not in tried and s not in in_flight]
    for s in [s for s in todo if s["status"] == "pending"]:
        alone(s)
    return timing


# ---- main ----------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--pass", dest="npass", type=int, choices=(1, 2), required=True)
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--scratch", type=Path, default=None,
                    help="node-local work dir (default: $SLURM_TMPDIR, $TMPDIR, /tmp)")
    args = ap.parse_args(argv)

    run = args.run.resolve()
    opts = json.loads((run / "run.json").read_text())["options"]
    key = f"{args.shard:05d}"
    pdir = run / f"pass{args.npass}"
    shard = pdir / "shards" / (f"{key}.tsv" if args.npass == 1 else f"{key}.smi")
    taskdir = pdir / "tasks" / key
    taskdir.mkdir(parents=True, exist_ok=True)

    in_sha, o_hash = sha256(shard), options_hash(opts)
    done = taskdir / "DONE.json"
    if done.exists():
        d = json.loads(done.read_text())
        if d.get("input_sha256") == in_sha and d.get("options_hash") == o_hash:
            print(f"pass {args.npass} shard {key}: already done ({d['counts']})")
            return 0
        done.unlink()
    if args.npass == 2 and (run / "pack.json").exists() and \
            json.loads((run / "pack.json").read_text()).get("loose_deleted"):
        raise SystemExit(f"{run} is packed and its loose db2s deleted (pack.py); pass 2 would "
                         "rebuild every name. Build into a new run instead.")

    os.environ["MSANI_CONFIG"] = str(run / "msani_configurations.yaml")
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = "1"
    started, t_all = time.strftime("%Y-%m-%dT%H:%M:%S%z"), time.perf_counter()
    base = args.scratch or Path(os.environ.get("SLURM_TMPDIR") or os.environ.get("TMPDIR") or "/tmp")
    scratch = Path(tempfile.mkdtemp(prefix=f"msani_p{args.npass}_{key}_", dir=base))
    try:
        if args.npass == 1:
            rows = []
            for line in open(shard):
                g, smi, name, pid = line.rstrip("\n").split("\t")
                rows.append((int(g), smi, name, pid))
            (pdir / "out").mkdir(exist_ok=True)
            parents = pass1(rows, opts, scratch, pdir / "out" / f"{key}_clean.smi")
            write_tsv(taskdir / "parents.tsv", PARENT_COLS, parents)
            counts = {"parents": len(parents), "ok": sum(p["status"] == "ok" for p in parents),
                      "rows": sum(p["n_rows"] for p in parents)}
            extra = {"seconds_2d": round(sum(p["seconds_2d"] for p in parents), 1)}
        else:
            states = []
            for line in open(shard):
                smi, name = line.split()[:2]
                states.append({"name": name, "smiles": smi, "status": "pending", "detail": ""})
            outdir = run / "db2_run" / "out" / key
            outdir.mkdir(parents=True, exist_ok=True)
            wanted = {f"{s['name']}.db2.gz" for s in states}
            for f in outdir.iterdir():              # strays from a different earlier plan
                if f.name not in wanted:
                    f.unlink()
            timing = pass2(states, opts, scratch, outdir, taskdir)
            write_tsv(taskdir / "states.tsv", STATE_COLS, states)
            write_tsv(taskdir / "timing.tsv", TIMING_COLS,
                      [{"name": n, **t} for n, t in sorted(timing.items())])
            counts = {"names": len(states), "ok": sum(s["status"] == "OK" for s in states)}
            extra = {}
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    d = {"pass": args.npass, "shard": args.shard, "input_sha256": in_sha, "options_hash": o_hash,
         "counts": counts, "versions": versions(), "host": os.uname().nodename, "cpu": cpu_model(),
         "slurm_job": os.environ.get("SLURM_ARRAY_JOB_ID") or os.environ.get("SLURM_JOB_ID", ""),
         "slurm_task": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
         "started": started, "finished": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
         "seconds": round(time.perf_counter() - t_all, 1), **extra}
    tmp = done.with_suffix(".tmp")
    tmp.write_text(json.dumps(d, indent=2) + "\n")
    tmp.replace(done)
    print(f"pass {args.npass} shard {key}: {counts}, {d['seconds']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
