#!/usr/bin/env python3
"""doctor_env.py STAGE [--method M --pbfilter on|off] -- `llb doctor`'s checks
inside one tool's environment (run by common/doctor.py through with_env.sh).

Prints one JSON object per check, {"stage","check","status","detail","fix"}.
Only the stdlib is imported at top level, and each check catches everything,
so a missing package is a FAIL line instead of a traceback that hides the rest.
Anything a package or tool prints goes to stderr: stdout carries only the JSON.
"""
from __future__ import annotations

import argparse
import datetime as dt
import importlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys

ROOT = os.environ.get("LLB_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG = os.environ.get("LLB_CONFIG") or "hpc.env"
RDKIT_PIN = "2024.03.5"
XTB_PIN = "6.4.1"
ZBH_URL = "https://software.zbh.uni-hamburg.de"
LICENSE_RE = re.compile(r"Your license is valid until:\s*([A-Za-z]+, \d{1,2} [A-Za-z]+ \d{4})")
LICENSE_WARN_DAYS = 30

_out = sys.stdout
STAGE = ""


def emit(check, status, detail="", fix=""):
    _out.write(json.dumps({"stage": STAGE, "check": check, "status": status,
                           "detail": str(detail), "fix": fix if status == "FAIL" else ""}) + "\n")
    _out.flush()


def run(cmd, timeout=60):
    """(returncode, stdout+stderr); returncode None when it could not run or timed out."""
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           stdin=subprocess.DEVNULL, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, f"timed out after {timeout} s"
    except OSError as e:
        return None, str(e)
    return p.returncode, p.stdout.decode("utf-8", "replace")


def tail(text, n=3):
    lines = [x.strip() for x in text.strip().splitlines() if x.strip()]
    return " | ".join(lines[-n:])


def resolve_exe(var, default):
    """Like llb_require_exe in common/env.sh: VAR's path, else DEFAULT on PATH."""
    want = os.environ.get(var, "")
    exe = want or default
    if "/" in exe:
        return exe if os.path.isfile(exe) and os.access(exe, os.X_OK) else None
    return shutil.which(exe)


def check_exe(name, var, default, fix_hint=""):
    exe = resolve_exe(var, default)
    if exe:
        emit(name, "PASS", exe)
    elif os.environ.get(var):
        emit(name, "FAIL", f"{var}={os.environ[var]}: not executable",
             f"set {var} in {CFG} to the {default} executable{fix_hint}")
    else:
        emit(name, "FAIL", f"{default} not on PATH in this env",
             f"set {var} in {CFG}{fix_hint}")
    return exe


def check_import(module, fix, name=None):
    try:
        m = importlib.import_module(module)
    except BaseException as e:            # noqa: BLE001 -- a SystemExit or segfault-adjacent error is a FAIL too
        emit(name or module, "FAIL", f"import {module}: {type(e).__name__}: {e}", fix)
        return None
    emit(name or module, "PASS", f"{getattr(m, '__version__', '')} {getattr(m, '__file__', '')}".strip())
    return m


def check_rdkit(fix):
    try:
        import rdkit
    except BaseException as e:            # noqa: BLE001
        emit("rdkit", "FAIL", f"import rdkit: {type(e).__name__}: {e}", fix)
        return
    v = rdkit.__version__
    if v == RDKIT_PIN:
        emit("rdkit", "PASS", v)
    else:
        emit("rdkit", "WARN", f"{v}, not the pinned {RDKIT_PIN}")


def check_on_path(name, exe, status_missing, fix):
    p = shutil.which(exe)
    if p:
        emit(name, "PASS", p)
    else:
        emit(name, status_missing, f"{exe} not on PATH", fix)


def license_status(text, today):
    """(status, detail) from unicon's banner; past or missing is FAIL, < 30 days WARN."""
    m = LICENSE_RE.search(text)
    if not m:
        return "FAIL", f"no 'Your license is valid until:' line in unicon's output ({tail(text) or 'empty'})"
    try:
        until = dt.datetime.strptime(m.group(1), "%A, %d %B %Y").date()
    except ValueError:
        return "FAIL", f"cannot parse license date '{m.group(1)}'"
    days = (until - today).days
    if days < 0:
        return "FAIL", f"expired on {m.group(1)} ({-days} days ago)"
    if days < LICENSE_WARN_DAYS:
        return "WARN", f"valid until {m.group(1)} ({days} days)"
    return "PASS", f"valid until {m.group(1)} ({days} days)"


def xtb_version(text):
    m = re.search(r"xtb version (\d+\.\d+\.\d+)", text)
    return m.group(1) if m else None


# ---- stages -------------------------------------------------------------------

def stage_prep(a):
    fix = f"PREP_ENV in {CFG} must name an env with RDKit (empty = the QupKake env)"
    check_rdkit(fix)
    sys.path.insert(0, os.path.join(ROOT, "run_qupkake"))
    for mod in ("prepare_parents", "assign_names"):
        check_import(mod, fix)


def stage_pka(a):
    fix = f"PKA_ENV in {CFG}: create it from environments/llb_pka.yml"
    check_rdkit(fix)
    check_import("torch", fix)
    d = os.environ.get("MOLGPKA_DIR", "")
    wfix = (f"set MOLGPKA_DIR in {CFG} to a clone of github.com/Xundrug/MolGpKa "
            "at 4dc8352 (config/hpc.env.example, pKa triage)")
    if not d:
        emit("molgpka_weights", "FAIL", "MOLGPKA_DIR is not set", wfix)
        return
    sys.path.insert(0, os.path.join(ROOT, "pka_triage"))
    try:
        import molgpka
    except BaseException as e:            # noqa: BLE001
        emit("molgpka_weights", "FAIL", f"cannot import pka_triage/molgpka.py: {e}", fix)
        return
    bad = []
    for fname, want in molgpka.KNOWN_WEIGHTS.items():
        p = os.path.join(d, "models", fname)
        if not os.path.exists(p):
            bad.append(f"{p} missing")
        elif molgpka._sha256(molgpka.Path(p)) != want:
            bad.append(f"{p} sha256 differs from 4dc8352")
    if bad:
        emit("molgpka_weights", "FAIL", "; ".join(bad), wfix)
    else:
        emit("molgpka_weights", "PASS", f"{d}/models: {', '.join(molgpka.KNOWN_WEIGHTS)} sha256 ok")


def stage_qupkake(a):
    fix = f"QUPKAKE_ENV in {CFG} (run_qupkake/README.md §1)"
    check_exe("qupkake_cli", "QUPKAKE_EXE", "qupkake")
    qk = check_import("qupkake", fix)
    check_rdkit(fix)
    xfix = f"build xtb {XTB_PIN} (run_qupkake/README.md §1) and set XTBPATH in {CFG}"
    xtb = os.environ.get("XTBPATH", "")
    label = f"$XTBPATH={xtb}"
    if not xtb:
        xtb = getattr(qk, "XTB_LOCATION", "") if qk else ""
        if not xtb:
            spec = importlib.util.find_spec("qupkake")
            if spec and spec.origin:
                xtb = os.path.join(os.path.dirname(spec.origin), "xtb-641", "bin", "xtb")
        label = f"QupKake's bundled xtb {xtb}" if xtb else ""
    if not xtb:
        emit("xtb", "FAIL", "XTBPATH unset and QupKake not found, so no xtb", xfix)
        return
    if not (os.path.isfile(xtb) and os.access(xtb, os.X_OK)):
        emit("xtb", "FAIL", f"{label}: not executable", xfix)
        return
    rc, out = run([xtb, "--version"])
    v = xtb_version(out)
    if v == XTB_PIN:
        emit("xtb", "PASS", f"{label}: xtb {v}")
    else:
        emit("xtb", "FAIL", f"{label}: `xtb --version` says {v or tail(out) or 'nothing'}, "
                            f"want {XTB_PIN}", xfix)


def stage_db2(a):
    fix = f"DB2C_ENV / DB2C_SETUP in {CFG} (UPSTREAM.md §5)"
    check_exe("build_ligand", "BUILD_LIGAND_EXE", "build_ligand")
    check_rdkit(fix)
    d2c = check_import("db2_converter", fix + "; or set DB2C_SRC")
    sys.path.insert(0, os.path.join(ROOT, "run_qupkake"))
    check_import("build_ligand_relaxed", fix, name="build_launcher")
    if a.pbfilter == "on":
        check_import("posebusters", f"install posebusters into DB2C_ENV ({CFG}), "
                                    "or run with --pbfilter off")
    if "bcl" in a.method:
        rc, out = run([sys.executable, os.path.join(ROOT, "run_qupkake", "check_bcl.py")])
        if rc == 0:
            emit("check_bcl", "PASS", tail(out, 2).replace("  bcl  ", ""))
        else:
            emit("check_bcl", "FAIL", f"method {a.method} needs BCL: {tail(out)}",
                 f"DB2C_SETUP in {CFG} must make BCL findable (BCL_ROOT, PATH; UPSTREAM.md §5)")
    afix = (f"DB2C_SETUP in {CFG} must load what AMSOL links against "
            "(e.g. the module giving libgfortran.so.5)")
    amsol = ""
    if d2c is not None:
        try:
            from db2_converter.config import config
            amsol = str(config["all"]["AMSOLEXE"])
        except BaseException as e:        # noqa: BLE001
            emit("amsol", "FAIL", f"no AMSOLEXE in db2_converter.config: {e}", fix)
    else:
        emit("amsol", "FAIL", "not checked: db2_converter not importable", fix)
    if amsol:
        if not os.path.exists(amsol):
            emit("amsol", "FAIL", f"AMSOLEXE={amsol} does not exist", fix)
        elif not shutil.which("ldd"):
            emit("amsol", "WARN", f"{amsol}: no ldd here, libraries not checked")
        else:
            rc, out = run(["ldd", amsol])
            missing = [x.strip() for x in out.splitlines() if "not found" in x]
            if missing:
                emit("amsol", "FAIL", f"{amsol}: {'; '.join(missing)}", afix)
            else:
                emit("amsol", "PASS", amsol)
    zbh = os.environ.get("ZBH_DIR", "")
    if zbh:
        # ZBH_DIR set: the binaries must be there, whatever PATH says
        zfix = f"ZBH_DIR={zbh} in {CFG} must hold the unicon and conformator binaries ({ZBH_URL})"
        found = {}
        for name in ("conformator", "unicon"):
            p = os.path.join(zbh, name)
            if os.path.isfile(p) and os.access(p, os.X_OK):
                emit(name, "PASS", p)
                found[name] = p
            else:
                emit(name, "FAIL", f"{p}: {'not executable' if os.path.exists(p) else 'missing'}",
                     zfix)
        unicon = found.get("unicon")
    else:
        if "conformator" in a.method:
            check_on_path("conformator", "conformator", "FAIL",
                          f"put conformator on PATH via DB2C_SETUP in {CFG}, or set ZBH_DIR "
                          f"({ZBH_URL})")
        unicon = shutil.which("unicon")
        if unicon:
            emit("unicon", "PASS", unicon)
        else:
            emit("unicon", "FAIL", "unicon not on PATH",
                 f"put unicon on PATH via DB2C_SETUP in {CFG}, or set ZBH_DIR ({ZBH_URL})")
    lfix = f"renew or install the UNICON license ({ZBH_URL}); see ZBH_LICENSE / DB2C_SETUP in {CFG}"
    if unicon:
        rc, out = run([unicon])
        status, detail = license_status(out, dt.date.today())
        emit("unicon_license", status, detail, lfix)
    else:
        emit("unicon_license", "FAIL", "not checked: no unicon", lfix)
    check_on_path("singularity", "singularity", "WARN", "")


def stage_dock(a):
    check_exe("db2tool", "DB2TOOL_EXE", "db2tool", " (or load it in DOCK_SETUP)")


STAGES = {"prep": stage_prep, "pka": stage_pka, "qupkake": stage_qupkake,
          "db2": stage_db2, "dock": stage_dock}


def main(argv=None):
    global _out, STAGE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=sorted(STAGES))
    ap.add_argument("--method", default="conformator_bcl")
    ap.add_argument("--pbfilter", choices=("on", "off"), default="on")
    a = ap.parse_args(argv)
    STAGE = a.stage
    # JSON on a private copy of stdout; fd 1 (imports, child tools) goes to stderr
    sys.stdout.flush()
    _out = os.fdopen(os.dup(1), "w")
    os.dup2(2, 1)
    try:
        STAGES[a.stage](a)
    except BaseException as e:            # noqa: BLE001 -- keep the lines already printed
        emit("doctor_env", "FAIL", f"{type(e).__name__}: {e}",
             "a bug in common/doctor_env.py; report it with this line")
    return 0


if __name__ == "__main__":
    sys.exit(main())
