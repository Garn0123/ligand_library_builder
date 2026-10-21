#!/usr/bin/env python3
"""llb doctor -- is this machine set up to run LLB? One line per check, and on
failure the fix.

    llb doctor [--stage LIST] [--node] [--smoke] [--method M] [--pbfilter on|off|auto]
               [--json FILE]

Host checks run here (stdlib only, login node OK). Each env stage is ONE
subprocess, `with_env.sh TOOL python common/doctor_env.py STAGE`, so a setup
hook or a conda env that does not activate is a FAIL line, not a crash.
Exit 0 = no FAIL, 1 = a FAIL, 64 = usage error.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STAGES = ("host", "prep", "pka", "qupkake", "db2", "dock")
TOOL = {"prep": "PREP", "pka": "PKA", "qupkake": "QUPKAKE", "db2": "DB2C", "dock": "DOCK"}
STAGE_TIMEOUT = int(os.environ.get("LLB_DOCTOR_TIMEOUT", "180"))
SMOKE_TIMEOUT = 1800
# two small molecules, one rotatable bond and an amide between them
SMOKE_SMI = "Oc1ccccc1 llb_smoke_phenol\nCC(=O)Nc1ccc(O)cc1 llb_smoke_paracetamol\n"


class Doctor:
    def __init__(self):
        self.checks = []

    def add(self, stage, check, status, detail="", fix=""):
        self.checks.append({"stage": stage, "check": check, "status": status,
                            "detail": " ".join(str(detail).split()),
                            "fix": fix if status == "FAIL" else ""})
        print(f"{status:<5} {stage:<8} {check:<18} {self.checks[-1]['detail']}", flush=True)
        if status == "FAIL" and fix:
            print(f"      fix: {fix}", flush=True)

    def count(self, status):
        return sum(c["status"] == status for c in self.checks)


def resolve_pbfilter(method, pbfilter):
    """auto = on when the method uses bcl, rdkit or confgenx (submit_db2.sh's rule)."""
    if pbfilter != "auto":
        return pbfilter
    return "on" if re.search(r"bcl|rdkit|confgenx", method) else "off"


def default_stages(cfg):
    out = ["host", "prep", "pka", "qupkake", "db2"]
    if any(cfg.get(k) for k in ("DOCK_SETUP", "DOCK_ENV", "DB2TOOL_EXE")):
        out.append("dock")
    return out


def config_path():
    """Where llb_load_config looks first (common/env.sh); for messages only."""
    if os.environ.get("LLB_CONFIG"):
        return os.environ["LLB_CONFIG"]
    site = ROOT.parent / "hpc.env"
    return str(site if site.is_file() else ROOT / "config" / "hpc.env")


def load_config():
    """(env with hpc.env applied, None) or (None, error text), via common/env.sh."""
    script = ('set -a; source "$LLB_ROOT/common/env.sh" && llb_load_config || exit 2; '
              'exec "$0" -c "import json, os; print(json.dumps(dict(os.environ)))"')
    env = dict(os.environ, LLB_ROOT=str(ROOT))
    try:
        p = subprocess.run(["bash", "-c", script, sys.executable], env=env, timeout=60,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                           stdin=subprocess.DEVNULL)
    except subprocess.TimeoutExpired:
        return None, "sourcing the config took over 60 s"
    if p.returncode != 0:
        fatal = [x for x in p.stderr.splitlines() if x.startswith("FATAL: ")]
        return None, (fatal[0][7:] if fatal else tail(p.stderr)) or f"llb_load_config exited {p.returncode}"
    try:
        return json.loads(p.stdout.strip().splitlines()[-1]), None
    except (ValueError, IndexError):
        return None, f"cannot read the loaded config: {tail(p.stdout + p.stderr)}"


def tail(text, n=3):
    lines = [x.strip() for x in (text or "").strip().splitlines() if x.strip()]
    return " | ".join(lines[-n:])


def run(cmd, timeout=30, env=None):
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                           timeout=timeout, env=env, stdin=subprocess.DEVNULL)
    except (OSError, subprocess.TimeoutExpired) as e:
        return None, str(e)
    return p.returncode, p.stdout


def zbh_license_problem(path):
    """None if PATH is a usable ZBH license directory (what ~/ZBH is: one INI-style
    file per tool, `[license]` then `unicon=<key>` / `conformator=<key>`), else
    why not. The keys themselves are never returned."""
    if not os.path.exists(path):
        return "does not exist"
    if not os.path.isdir(path):
        return "not a directory (ZBH_LICENSE is the directory that plays ~/ZBH's role)"
    if not os.access(path, os.R_OK | os.X_OK):
        return "directory not readable"
    found = set()
    for name in sorted(os.listdir(path)):
        f = os.path.join(path, name)
        if not os.path.isfile(f):
            continue
        try:
            text = open(f, encoding="utf-8", errors="replace").read()
        except OSError:
            continue
        section = None
        for line in text.splitlines():
            line = line.strip()
            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1].strip().lower()
            elif section == "license" and "=" in line:
                key, value = (x.strip() for x in line.split("=", 1))
                if key.lower() in ("unicon", "conformator") and value:
                    found.add(key.lower())
    missing = [t for t in ("unicon", "conformator") if t not in found]
    if missing:
        return f"no [license] {' or '.join(t + '=' for t in missing)} key in any file there"
    return None


# ---- host ---------------------------------------------------------------------

def host_checks(doc, cfg, err):
    cpath = (cfg or {}).get("LLB_CONFIG") or config_path()
    if cfg is None:
        doc.add("host", "config", "FAIL", err,
                f"cp {ROOT}/config/hpc.env.example {cpath} and edit it "
                f"({ROOT.parent}/hpc.env for the site layout), or set LLB_CONFIG (UPSTREAM.md §0)")
    else:
        user_cfg = cfg.get("LLB_USER_CONFIG", "")
        doc.add("host", "config", "PASS", cpath + (f" + {user_cfg}" if user_cfg else ""))
        env = cfg
        acct, part = cfg.get("LLB_ACCOUNT", ""), cfg.get("LLB_PARTITION", "")
        if not acct:
            doc.add("host", "account", "FAIL", "LLB_ACCOUNT is empty", f"set LLB_ACCOUNT in {cpath}")
        elif not shutil.which("sacctmgr", path=env.get("PATH")):
            doc.add("host", "account", "WARN", f"{acct}: sacctmgr not on PATH, not verified")
        else:
            user = env.get("USER") or os.environ.get("USER", "")
            rc, out = run(["sacctmgr", "-nP", "show", "assoc", f"user={user}", "format=account"],
                          env=env)
            have = sorted({x.strip() for x in out.splitlines() if x.strip()}) if rc == 0 else []
            if acct in have:
                doc.add("host", "account", "PASS", f"{acct} (user {user})")
            else:
                doc.add("host", "account", "FAIL",
                        f"{acct} not among {user}'s accounts: {', '.join(have) or tail(out) or 'none'}",
                        f"set LLB_ACCOUNT in {cpath} to one of your SLURM accounts")
        if not part:
            doc.add("host", "partition", "FAIL", "LLB_PARTITION is empty",
                    f"set LLB_PARTITION in {cpath}")
        elif not shutil.which("sinfo", path=env.get("PATH")):
            doc.add("host", "partition", "WARN", f"{part}: sinfo not on PATH, not verified")
        else:
            rc, out = run(["sinfo", "-h", "-p", part], env=env)
            if rc == 0 and out.strip():
                doc.add("host", "partition", "PASS", part)
            else:
                doc.add("host", "partition", "FAIL", f"sinfo does not know partition {part}"
                        + (f": {tail(out)}" if out.strip() else ""),
                        f"set LLB_PARTITION in {cpath} (sinfo -s lists them)")
        envs = [k for k in ("QUPKAKE_ENV", "DB2C_ENV", "PREP_ENV", "PKA_ENV", "DOCK_ENV")
                if cfg.get(k)]
        base = cfg.get("CONDA_BASE", "")
        if not envs:
            doc.add("host", "conda", "PASS", "no *_ENV set, conda not needed")
        elif base and os.path.isfile(os.path.join(base, "etc/profile.d/conda.sh")):
            doc.add("host", "conda", "PASS", base)
        else:
            doc.add("host", "conda", "FAIL",
                    f"no etc/profile.d/conda.sh under CONDA_BASE='{base}' ({', '.join(envs)} set)",
                    f"set CONDA_BASE in {cpath}")
        drap = cfg.get("DRAP_DIR", "")
        if not drap:
            v = ROOT / "run_qupkake" / "validate_names.py"
            if v.is_file():
                doc.add("host", "drap", "PASS", f"{v} (bundled copy)")
            else:
                doc.add("host", "drap", "FAIL", f"{v} missing",
                        "restore it: git checkout -- run_qupkake/validate_names.py")
        else:
            v = os.path.join(drap, "stage4", "validate_names.py")
            if os.path.isfile(v):
                doc.add("host", "drap", "PASS", v)
            else:
                doc.add("host", "drap", "FAIL", f"no stage4/validate_names.py under DRAP_DIR='{drap}'",
                        f"set DRAP_DIR in {cpath} to a DRAP checkout, or empty it to use "
                        "the bundled run_qupkake/validate_names.py")
        lic = cfg.get("ZBH_LICENSE", "")
        if lic:
            site = cfg.get("LLB_SITE", "")
            user = cfg.get("USER") or os.environ.get("USER", "")
            problem = zbh_license_problem(lic)
            if not problem:
                doc.add("host", "zbh_license", "PASS", f"{lic}: unicon and conformator keys")
            elif (site and lic == f"{site}/licenses/{user}/ZBH" and not os.path.exists(lic)
                  and not cfg.get("LLB_IMAGE") and not cfg.get("ZBH_DIR")):
                # only the site default, and no shim reads it (ZBH still uses ~/ZBH)
                doc.add("host", "zbh_license", "WARN",
                        f"{lic} (the site default) does not exist; copy your ~/ZBH directory there")
            else:
                doc.add("host", "zbh_license", "FAIL", f"ZBH_LICENSE={lic}: {problem}",
                        f"copy your ~/ZBH directory to {lic}, or set ZBH_LICENSE in "
                        f"{cfg.get('LLB_USER_CONFIG') or cpath} (per-user: users/$USER.env beside hpc.env)")

    try:
        with tempfile.NamedTemporaryFile(dir=os.getcwd(), prefix=".llb_doctor."):
            pass
        doc.add("host", "rundir", "PASS", os.getcwd())
    except OSError as e:
        doc.add("host", "rundir", "FAIL", f"{os.getcwd()}: {e}",
                "run llb from a writable run directory (bin/llb header)")
    path = (cfg or {}).get("PATH") or os.environ.get("PATH")
    for name, exe, why in (("git", "git", "provenance"),
                           ("wget", (cfg or {}).get("WGET") or "wget", "llb fetch"),
                           ("parallel", (cfg or {}).get("PARALLEL") or "parallel", "llb fetch")):
        p = shutil.which(exe, path=path)
        if p:
            doc.add("host", name, "PASS", p)
        else:
            doc.add("host", name, "WARN", f"{exe} not on PATH (needed for {why})")


# ---- the toolchain image -----------------------------------------------------------

def image_tools(cfg):
    """Stages whose tool runs in LLB_IMAGE (its <TOOL>_ENV empty; never dock)."""
    if not cfg.get("LLB_IMAGE"):
        return set()
    return {s for s, t in TOOL.items() if s != "dock" and not cfg.get(f"{t}_ENV")}


def image_checks(doc, cfg):
    """Runtime, image file, toolchain version. False = env stages cannot run here."""
    img, cpath = cfg["LLB_IMAGE"], cfg.get("LLB_CONFIG", config_path())
    p = subprocess.run(["bash", "-c", 'source "$LLB_ROOT/common/env.sh" && llb_container'],
                       env=dict(cfg, LLB_ROOT=str(ROOT)), stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, text=True, stdin=subprocess.DEVNULL)
    rt = p.stdout.strip()
    if p.returncode or not rt:
        doc.add("host", "container", "WARN",
                "no apptainer or singularity here (a login node?); image stages not checked: "
                "run `llb doctor --node`")
        return False
    rc, out = run([rt, "--version"])
    doc.add("host", "container", "PASS", f"{rt}: {tail(out, 1)}")
    if not os.path.isfile(img):
        doc.add("host", "image", "FAIL", f"LLB_IMAGE={img} does not exist",
                f"pull it (container/README.md) or fix LLB_IMAGE in {cpath}")
        return False
    want = (ROOT / "container" / "TOOLCHAIN_VERSION").read_text().strip()
    rc, out = run([rt, "inspect", "--labels", img], timeout=60)
    got = re.search(r"org\.llb\.toolchain\.version\W+(\S+)", out or "")
    if not got:
        doc.add("host", "image", "FAIL", f"{img}: no toolchain version label ({tail(out, 1)})",
                "rebuild or re-pull the image (container/README.md)")
        return False
    if got.group(1) != want:
        doc.add("host", "image", "FAIL", f"{img} is toolchain {got.group(1)}; this checkout wants {want}",
                f"pull llb-toolchain:{want} and point LLB_IMAGE at it in {cpath}")
        return True
    doc.add("host", "image", "PASS", f"{img} (toolchain {want})")
    return True


# ---- env stages -------------------------------------------------------------------

def env_stage(doc, stage, cfg, method, pbfilter):
    tool = TOOL[stage]
    py = "python3" if stage == "dock" else "python"
    cmd = [str(ROOT / "common" / "with_env.sh"), tool, py, str(ROOT / "common" / "doctor_env.py"),
           stage]
    if stage == "db2":
        cmd += ["--method", method, "--pbfilter", pbfilter]
    env = dict(cfg, LLB_ROOT=str(ROOT))
    cpath = cfg.get("LLB_CONFIG", config_path())
    fix = f"check {tool}_SETUP and {tool}_ENV (and CONDA_BASE) in {cpath}"
    if stage == "prep":
        fix += "; empty PREP_ENV means the QupKake env"
    try:
        p = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             stdin=subprocess.DEVNULL, text=True, start_new_session=True)
    except OSError as e:
        doc.add(stage, "activate", "FAIL", str(e), fix)
        return
    try:
        out, err = p.communicate(timeout=STAGE_TIMEOUT)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(p.pid, signal.SIGKILL)
        except OSError:
            pass
        out, err = p.communicate()
        doc.add(stage, "activate", "FAIL",
                f"no answer in {STAGE_TIMEOUT} s (a setup hook waiting on input?): {tail(err)}", fix)
        return
    rows, banner = [], ""
    for line in out.splitlines():
        if line.startswith("{"):
            try:
                rows.append(json.loads(line))
                continue
            except ValueError:
                pass
        if line.startswith(f"[{tool}]") and not banner:
            banner = line[len(tool) + 2:].strip()
    if not rows:
        doc.add(stage, "activate", "FAIL", tail(err) or f"exit {p.returncode}, no output", fix)
        return
    doc.add(stage, "activate", "PASS", banner or tool)
    for r in rows:
        doc.add(stage, r.get("check", "?"), r.get("status", "FAIL"), r.get("detail", ""),
                r.get("fix", ""))
    if p.returncode != 0:
        doc.add(stage, "doctor_env", "FAIL", f"exit {p.returncode}: {tail(err)}",
                "a crash in common/doctor_env.py; report it with this line")


# ---- --node, --smoke ---------------------------------------------------------------

def run_on_node(argv, cfg):
    rest = [a for a in argv if a != "--node"]
    acct, part = cfg.get("LLB_ACCOUNT", ""), cfg.get("LLB_PARTITION", "")
    if not (acct and part):
        print(f"doctor --node: set LLB_ACCOUNT and LLB_PARTITION in {cfg.get('LLB_CONFIG')}",
              file=sys.stderr)
        return 1
    cmd = ["srun", "-A", acct, "-p", part, "-c", "1", "--mem=4G", "-t", "00:15:00",
           str(ROOT / "bin" / "llb"), "doctor"] + rest
    print("doctor: " + " ".join(cmd), flush=True)
    try:
        return subprocess.call(cmd, env=dict(cfg, LLB_ROOT=str(ROOT)))
    except OSError as e:
        print(f"doctor --node: cannot run srun: {e}", file=sys.stderr)
        return 1


def smoke(doc, cfg, method, pbfilter):
    wd = Path(tempfile.mkdtemp(prefix="llb_doctor_smoke.", dir=os.getcwd()))
    (wd / "in.smi").write_text(SMOKE_SMI)
    names = [x.split()[1] for x in SMOKE_SMI.splitlines()]
    # db2_array.sbatch's flags (restrained MMFF), with fewer conformers
    cmd = [str(ROOT / "bin" / "llb"), "build-ligand", "-i", str(wd / "in.smi"), "-n", "50",
           "--keep_max_conf", "-m", method, "--checkstereo", "--reseth", "--rotateh", "--dock38",
           "--MMFFopt"] + (["--PBfilter"] if pbfilter == "on" else []) \
        + ["--workingpath", str(wd), "--outputpath", str(wd)]
    log = wd / "build.log"
    try:
        with open(log, "w") as fh:
            rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT, timeout=SMOKE_TIMEOUT,
                                 env=dict(cfg, LLB_ROOT=str(ROOT)), stdin=subprocess.DEVNULL)
    except (OSError, subprocess.TimeoutExpired) as e:
        rc = str(e)
    missing = [n for n in names if not list(wd.rglob(f"{n}.db2.gz"))]
    if not missing:
        doc.add("smoke", "build_ligand", "PASS",
                f"{len(names)} db2.gz built ({method}, pbfilter {pbfilter})")
        shutil.rmtree(wd, ignore_errors=True)
    else:
        doc.add("smoke", "build_ligand", "FAIL",
                f"no db2.gz for {', '.join(missing)} (exit {rc}); log {log}",
                f"read {log}; then DB2C_SETUP / DB2C_ENV in {cfg.get('LLB_CONFIG')} (UPSTREAM.md §5)")


# ---- main ------------------------------------------------------------------------

class Parser(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        print(f"llb doctor: {message}", file=sys.stderr)
        sys.exit(64)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    ap = Parser(prog="llb doctor", description=__doc__.split("\n")[0])
    ap.add_argument("--stage", help=f"comma list of {','.join(STAGES)}")
    ap.add_argument("--node", action="store_true", help="run on a compute node via srun")
    ap.add_argument("--smoke", action="store_true", help="then build 2 molecules")
    ap.add_argument("--method", default="conformator_bcl")
    ap.add_argument("--pbfilter", choices=("on", "off", "auto"), default="auto")
    ap.add_argument("--json", metavar="FILE")
    a = ap.parse_args(argv)
    stages = None
    if a.stage:
        stages = [s.strip() for s in a.stage.split(",") if s.strip()]
        bad = [s for s in stages if s not in STAGES]
        if bad or not stages:
            ap.error(f"unknown stage {','.join(bad) or a.stage!r}; choose from {','.join(STAGES)}")
    pbfilter = resolve_pbfilter(a.method, a.pbfilter)

    cfg, err = load_config()
    doc = Doctor()
    if a.node:
        if cfg is None:
            host_checks(doc, None, err)
            return finish(doc, a, None)
        return run_on_node(argv, cfg)
    if stages is None:
        stages = default_stages(cfg or {})
    if "host" in stages or cfg is None:
        host_checks(doc, cfg, err)
    skipped = [s for s in stages if s != "host"] if cfg is None else []
    if cfg is not None:
        in_image = image_tools(cfg) if not os.environ.get("LLB_IN_IMAGE") else set()
        if in_image and not image_checks(doc, cfg):
            stages = [s for s in stages if s not in in_image]
        for s in stages:
            if s != "host":
                env_stage(doc, s, cfg, a.method, pbfilter)
    if skipped:
        print(f"doctor: not run without a config: {','.join(skipped)}")
    if a.smoke:
        if doc.count("FAIL"):
            print("doctor: --smoke skipped, fix the FAILs first")
        else:
            smoke(doc, cfg, a.method, pbfilter)
    return finish(doc, a, cfg)


def finish(doc, a, cfg):
    cpath = (cfg or {}).get("LLB_CONFIG") or config_path()
    if a.json:
        Path(a.json).write_text(json.dumps({
            "llb_root": str(ROOT), "config": cpath, "host": socket.gethostname(),
            "when": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
            "checks": doc.checks}, indent=1) + "\n")
    print(f"doctor: {doc.count('PASS')} passed, {doc.count('WARN')} warnings, "
          f"{doc.count('FAIL')} failed (config: {cpath})")
    return 1 if doc.count("FAIL") else 0


if __name__ == "__main__":
    sys.exit(main())
