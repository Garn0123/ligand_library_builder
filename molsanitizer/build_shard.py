#!/usr/bin/env python3
"""
build_shard.py -- one array task: a shard of parents -> named states -> db2.gz.

Runs inside the MolSanitizer env (see environment/). Three steps per shard:

  2D   MolSanitizer's Python API on one parent at a time: salts, neutralize,
       tautomers, protonation at the run's pH, stereoisomers. One parent per
       call so each failure has its own reason and each parent its own seconds.
       (The API takes a non-integer pH; the msani CLI does not.)
  name each state -> a 16-character contract name (naming.py).
  3D   the msani CLI on the named states: --no-neutralize --no-stereoisomers
       -3d -f db2, so it neither changes a state nor renames it. One call for the
       shard; if that call dies or leaves molecules unaccounted, each remaining
       molecule is retried alone under --mol-timeout, so one crash costs one
       molecule, not the shard.

Every db2 is checked before it is kept (the guarantee `llb verify-db2` gave):
each record's M-line name is the requested name, and its AMSOL total charge is
the state's formal charge within charge_tol.

Outputs (RUN from run.json, K = shard number):
  RUN/db2_run/out/K/<name>.db2.gz    one per built name (prepare.py's layout)
  RUN/tasks/K/states.tsv             one row per named state: library.tsv columns
                                     + status/detail/db2 header fields
  RUN/tasks/K/parents.tsv            one row per shard line: built or failed + reason
  RUN/tasks/K/timing.tsv             per name: MolSanitizer's own step times + ours
  RUN/tasks/K/DONE.json              last; versions, host, seconds, input sha256

Resumable: a shard whose DONE.json matches its input and options is skipped; an
interrupted one is redone, keeping every db2.gz already written for an unchanged
name (each is written to a temporary file and renamed, so a present one is whole).

    python build_shard.py --run RUN --shard 17
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from naming import (assign_indices, base_id, base_id_of, contract_name,  # noqa: E402
                    invertomer_key, is_zinc, parse_state_id)

LIBRARY_COLS = ["name", "base_id", "parent_id", "zinc_id", "input_name", "protomer_index",
                "stereo_index", "net_charge", "heavy_atoms", "ph_values", "population_estimate",
                "invertomer_of", "note", "smiles", "site_source"]
STATE_COLS = LIBRARY_COLS + ["gidx", "msani_id", "status", "detail", "db2_records",
                             "amsol_charge", "db2_bytes"]
PARENT_COLS = ["gidx", "input_name", "parent_id", "input_smiles", "status", "reason",
               "n_states", "n_names", "n_ok", "seconds_2d", "seconds_3d"]
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


# ---- 2D ------------------------------------------------------------------------

def prep_2d(rows, opts, work: Path):
    """Per parent: (parent row, [state dicts]). States carry MolSanitizer's id."""
    from pandas import DataFrame
    from rdkit import Chem, RDLogger
    from msani.api import Msani
    from msani.filtering.filters import Filters
    RDLogger.DisableLog("rdApp.*")

    rejected = work / "rejected_2d.txt"
    metal_err = work / "metal_error_2d.smi"
    proc = Msani(removesalts=opts["removesalts"], neutralize=opts["neutralize"],
                 tautomers=opts["tautomers"], extended_tautomers=opts["extended_tautomers"],
                 protonation=opts["protonation"], pH=opts["pH"], pH_range=opts["pH_range"],
                 stereoisomers=True, max_stereoisomers=opts["max_stereoisomers"],
                 stereo_timeout=opts["stereo_timeout"], randomSeed=opts["random_seed"],
                 metal=opts["metal"], numcores=1)
    out = []
    with proc:
        for gidx, smi, input_name, pid in rows:
            t0 = time.perf_counter()
            prow = {"gidx": gidx, "input_name": input_name, "parent_id": pid,
                    "input_smiles": smi, "status": "", "reason": ""}
            states = []
            if Chem.MolFromSmiles(smi) is None:
                prow.update(status="failed", reason="invalid_smiles")
            else:
                try:
                    df = proc.run(DataFrame({"smiles": [smi], "ids": [pid]}),
                                  rejected_file=str(rejected), metal_error_file=str(metal_err))
                except Exception as e:  # one parent's crash is that parent's reason
                    prow.update(status="failed", reason=f"msani_2d_error: {type(e).__name__}: {e}"[:300])
                    df = None
                if df is not None and df.empty:
                    prow.update(status="failed", reason="removed_in_2d: " + grep_reason(pid, rejected, metal_err))
                elif df is not None:
                    kept = Filters.remove_exotic_chem_to_db2(df.copy())
                    kept_ids = set(kept["ids"])
                    for _, r in df.iterrows():
                        states.append({"msani_id": str(r["ids"]), "smiles": str(r["smiles"]),
                                       "exotic": r["ids"] not in kept_ids})
            prow["seconds_2d"] = round(time.perf_counter() - t0, 3)
            out.append((prow, states))
    return out


def grep_reason(pid: str, *files: Path) -> str:
    for f in files:
        if f.exists():
            for line in open(f):
                if pid in line.split():
                    return line.strip().split(pid, 1)[-1].strip()[:200] or f.name
    return "no reason recorded by MolSanitizer"


# ---- names ---------------------------------------------------------------------

def name_states(prepped, opts):
    """Contract names, invertomer links, in-shard duplicate structures."""
    from rdkit import Chem
    owner: dict[str, str] = {}
    ph = (f"{opts['pH'] - opts['pH_range']:g}-{opts['pH'] + opts['pH_range']:g}"
          if opts["pH_range"] else f"{opts['pH']:g}")
    for prow, states in prepped:
        if prow["status"]:
            continue
        pid = prow["parent_id"]
        keys = [parse_state_id(pid, s["msani_id"]) for s in states]
        if any(k is None for k in keys):
            bad = [s["msani_id"] for s, k in zip(states, keys) if k is None]
            prow.update(status="failed", reason=f"unparsed_msani_id: {bad[:3]}")
            states.clear()
            continue
        idx = assign_indices(keys)
        if max(i for i, _ in idx) > 61 or max(j for _, j in idx) > 61:
            prow.update(status="failed", reason="more than 62 states or isomers")
            states.clear()
            continue
        base = base_id_of(pid)
        first_of: dict[tuple[int, str], str] = {}
        for s, (i, j) in zip(states, idx):
            name = contract_name(base, i, j)
            mol = Chem.MolFromSmiles(s["smiles"])
            unassigned = sum(1 for _, t in Chem.FindMolChiralCenters(
                mol, includeUnassigned=True, useLegacyImplementation=False) if t == "?")
            inv = first_of.setdefault((i, invertomer_key(s["smiles"])), name)
            s.update({
                "name": name, "base_id": base_id(name), "parent_id": pid,
                "zinc_id": pid if is_zinc(pid) else "", "input_name": prow["input_name"],
                "protomer_index": i, "stereo_index": j,
                "net_charge": Chem.GetFormalCharge(mol), "heavy_atoms": mol.GetNumHeavyAtoms(),
                "ph_values": ph, "population_estimate": "",
                "invertomer_of": inv if inv != name else "",
                "note": f"unassigned_stereo={unassigned}" if unassigned else "",
                "site_source": "molsanitizer", "gidx": prow["gidx"],
                "status": "pending", "detail": ""})
            canon = Chem.MolToSmiles(mol)
            if s.pop("exotic"):
                s.update(status="EXOTIC", detail="not db2-compatible (MolSanitizer remove_exotic_chem_to_db2)")
            elif canon in owner:
                s.update(status="DUPLICATE", detail=f"same structure as {owner[canon]}")
            else:
                owner[canon] = name
    return prepped


# ---- 3D ------------------------------------------------------------------------

def msani_cmd(todo: Path, opts: dict) -> list[str]:
    exe = shutil.which("msani") or str(Path(sys.executable).parent / "msani")
    return [exe, "-i", str(todo), "--prefix", "todo",
            "--no-neutralize", "--no-stereoisomers", "--metal", "off",
            "-3d", "-f", "db2", "--numconfs", str(opts["numconfs"]),
            "--energywindow", str(opts["energywindow"]), "--rmsd", str(opts["rmsd"]),
            "--nringconfs", str(opts["nringconfs"]), "--timeout", str(opts["embed_timeout"]),
            "--randomSeed", str(opts["random_seed"]), "-j", "1", "--timing"]


def run_msani(states, opts, work: Path, timeout: float) -> tuple[int | str, float]:
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    todo = work / "todo.smi"
    todo.write_text("".join(f"{s['smiles']} {s['name']}\n" for s in states))
    t0 = time.perf_counter()
    with open(work / "msani.stdout", "w") as log:
        try:
            r = subprocess.run(msani_cmd(todo, opts), cwd=work, stdout=log, stderr=subprocess.STDOUT,
                               timeout=timeout)
            rc: int | str = r.returncode
        except subprocess.TimeoutExpired:
            rc = "timeout"
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


def build_3d(states, opts, scratch: Path, outdir: Path, taskdir: Path, timing: dict):
    """Run msani on every pending state; set each state's status."""
    todo = [s for s in states if s["status"] == "pending"]
    for s in todo:                                  # resumed shard: already kept
        if (outdir / f"{s['name']}.db2.gz").exists():
            s.update(status="OK", detail="kept from an earlier attempt")
            timing[s["name"]] = {"how": "earlier_attempt"}
    todo = [s for s in todo if s["status"] == "pending"]

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

    if todo:
        # one call for the shard; generous ceiling, the per-molecule retry is the real guard
        rc, wall = run_msani(todo, opts, scratch / "batch", timeout=opts["mol_timeout"] * len(todo))
        harvest(scratch / "batch", todo, "batch", rc, wall)
        keep_log(scratch / "batch", taskdir, "batch")
    for s in [s for s in todo if s["status"] == "pending"]:
        rc, wall = run_msani([s], opts, scratch / "alone", timeout=opts["mol_timeout"])
        harvest(scratch / "alone", [s], "alone", rc, wall)
        if rc == "timeout":
            s.update(status="BUILD_FAILED", detail=f"timeout after {opts['mol_timeout']} s alone")


def keep_log(work: Path, taskdir: Path, tag: str) -> None:
    for f in ("todo.log", "msani.stdout", "msani_error.err"):
        if (work / f).exists():
            shutil.copy2(work / f, taskdir / f"{tag}.{f}")


# ---- main ----------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--scratch", type=Path, default=None,
                    help="node-local work dir (default: $SLURM_TMPDIR, $TMPDIR, /tmp)")
    args = ap.parse_args(argv)

    run = args.run.resolve()
    meta = json.loads((run / "run.json").read_text())
    opts = meta["options"]
    key = f"{args.shard:05d}"
    shard = run / "shards" / f"{key}.tsv"
    taskdir = run / "tasks" / key
    outdir = run / "db2_run" / "out" / key
    taskdir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    in_sha, o_hash = sha256(shard), options_hash(opts)
    done = taskdir / "DONE.json"
    if done.exists():
        d = json.loads(done.read_text())
        if d.get("input_sha256") == in_sha and d.get("options_hash") == o_hash:
            print(f"shard {key}: already done ({d['counts']})")
            return 0
        done.unlink()

    os.environ["MSANI_CONFIG"] = str(run / "msani_configurations.yaml")
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = "1"
    started, t_all = time.strftime("%Y-%m-%dT%H:%M:%S%z"), time.perf_counter()
    rows = []
    for line in open(shard):
        g, smi, name, pid = line.rstrip("\n").split("\t")
        rows.append((int(g), smi, name, pid))

    base = args.scratch or Path(os.environ.get("SLURM_TMPDIR") or os.environ.get("TMPDIR") or "/tmp")
    scratch = Path(tempfile.mkdtemp(prefix=f"msani_{key}_", dir=base))
    try:
        prepped = name_states(prep_2d(rows, opts, scratch), opts)
        states = [s for _, ss in prepped for s in ss]

        # A resumed shard keeps old db2s only for names whose structure is unchanged.
        prev = taskdir / "states.partial.tsv"
        if prev.exists():
            with open(prev, newline="") as fh:
                old = {r["name"]: r["smiles"] for r in csv.DictReader(fh, delimiter="\t")}
            for s in states:
                if s["name"] in old and old[s["name"]] != s["smiles"]:
                    (outdir / f"{s['name']}.db2.gz").unlink(missing_ok=True)
        write_tsv(prev, STATE_COLS, states)
        wanted = {f"{s['name']}.db2.gz" for s in states}
        for f in outdir.iterdir():                  # strays from a different earlier plan
            if f.name not in wanted:
                f.unlink()

        timing: dict[str, dict] = {}
        build_3d(states, opts, scratch, outdir, taskdir, timing)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    parents = []
    for prow, ss in prepped:
        ok = sum(1 for s in ss if s["status"] == "OK")
        secs = 0.0
        for s in ss:
            t = timing.get(s.get("name", ""), {})
            try:
                secs += float(t.get("seconds_wall") or t.get("msani_total") or 0)
            except ValueError:
                pass
        if not prow["status"]:
            if ok:
                prow["status"] = "built"
            else:
                prow["status"] = "failed"
                why = sorted({f"{s['status']}: {s['detail']}" for s in ss})
                prow["reason"] = ("no state built -- " + " | ".join(why))[:500] if why else "no states"
        prow.update(n_states=len({s.get("protomer_index") for s in ss}), n_names=len(ss),
                    n_ok=ok, seconds_3d=round(secs, 3))
        parents.append(prow)

    write_tsv(taskdir / "states.tsv", STATE_COLS, states)
    write_tsv(taskdir / "parents.tsv", PARENT_COLS, parents)
    write_tsv(taskdir / "timing.tsv", TIMING_COLS,
              [{"name": n, **t} for n, t in sorted(timing.items())])
    prev.unlink(missing_ok=True)
    counts = {"parents": len(parents), "built": sum(p["status"] == "built" for p in parents),
              "names": len(states), "ok": sum(s["status"] == "OK" for s in states)}
    d = {"shard": args.shard, "input_sha256": in_sha, "options_hash": o_hash, "counts": counts,
         "versions": versions(), "host": os.uname().nodename, "cpu": cpu_model(),
         "slurm_job": os.environ.get("SLURM_ARRAY_JOB_ID") or os.environ.get("SLURM_JOB_ID", ""),
         "slurm_task": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
         "started": started, "finished": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
         "seconds": round(time.perf_counter() - t_all, 1),
         "seconds_2d": round(sum(p["seconds_2d"] for p in parents), 1)}
    tmp = done.with_suffix(".tmp")
    tmp.write_text(json.dumps(d, indent=2) + "\n")
    tmp.replace(done)
    print(f"shard {key}: {counts['built']}/{counts['parents']} parents built, "
          f"{counts['ok']}/{counts['names']} names OK, {d['seconds']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
