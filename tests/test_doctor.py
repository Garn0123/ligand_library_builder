"""llb doctor against a fake machine: a temp hpc.env, a no-op conda, and fake
sacctmgr/sinfo/unicon/xtb on PATH. No SLURM and no chemistry tools run."""
import datetime as dt
import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "common"))

import doctor  # noqa: E402
import doctor_env  # noqa: E402

DOCTOR = ROOT / "common" / "doctor.py"
DOCTOR_ENV = ROOT / "common" / "doctor_env.py"


def exe(path, body):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/bash\n" + body + "\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def banner(days):
    d = dt.date.today() + dt.timedelta(days=days)
    return f"UNICON 1.4.2\nYour license is valid until: {d.strftime('%A, %d %B %Y')}\n"


@pytest.fixture
def machine(tmp_path):
    """A fake cluster that passes every host check; returns a dict to tweak."""
    bin_ = tmp_path / "bin"
    exe(bin_ / "sacctmgr", "printf 'other\\nmaom99\\n'")
    exe(bin_ / "sinfo", '[[ "$3" == standard ]] && echo "standard up 14-00:00:00 10 idle"; true')
    for name in ("git", "wget", "parallel"):
        exe(bin_ / name, "true")
    conda = tmp_path / "conda"
    (conda / "etc/profile.d").mkdir(parents=True)
    (conda / "etc/profile.d/conda.sh").write_text("conda() { :; }\n")
    drap = tmp_path / "drap"
    (drap / "stage4").mkdir(parents=True)
    (drap / "stage4/validate_names.py").write_text("")
    run = tmp_path / "run"
    run.mkdir()
    keys = {"CONDA_BASE": conda, "LLB_ACCOUNT": "maom99", "LLB_PARTITION": "standard",
            "QUPKAKE_ENV": "qupkake", "DB2C_ENV": "db2_converter", "PKA_ENV": "llb_pka",
            "PREP_ENV": "", "DRAP_DIR": drap, "DOCK_SETUP": "", "DOCK_ENV": "",
            "DB2TOOL_EXE": ""}
    return {"tmp": tmp_path, "bin": bin_, "run": run, "keys": keys,
            "path": f"{bin_}:{Path(sys.executable).parent}:/usr/bin:/bin"}


def llb_doctor(m, *args):
    cfg = m["tmp"] / "hpc.env"
    cfg.write_text("".join(f'{k}="{v}"\n' for k, v in m["keys"].items()))
    env = {"PATH": m["path"], "HOME": str(m["tmp"]), "USER": "tester", "LLB_CONFIG": str(cfg),
           "LLB_DOCTOR_TIMEOUT": "60"}
    p = subprocess.run([sys.executable, str(DOCTOR), *args], cwd=m["run"], env=env,
                       capture_output=True, text=True, timeout=300)
    return p.returncode, p.stdout, p.stderr


def lines(out):
    """{(stage, check): status} from doctor's check lines."""
    got = {}
    for line in out.splitlines():
        parts = line.split(None, 3)
        if len(parts) >= 3 and parts[0] in ("PASS", "WARN", "FAIL"):
            got[(parts[1], parts[2])] = parts[0]
    return got


def run_env(m, stage, *args, extra=None):
    env = {"PATH": m["path"], "HOME": str(m["tmp"]), "LLB_ROOT": str(ROOT),
           "LLB_CONFIG": "hpc.env", **(extra or {})}
    p = subprocess.run([sys.executable, str(DOCTOR_ENV), stage, *args], env=env,
                       capture_output=True, text=True, timeout=120)
    assert p.returncode == 0, p.stderr
    rows = [json.loads(x) for x in p.stdout.splitlines()]   # stdout is JSON only
    return {r["check"]: r for r in rows}


# ---- pure pieces ------------------------------------------------------------------

def test_license_date_past_near_far():
    today = dt.date(2026, 10, 20)
    assert doctor_env.license_status(
        "Your license is valid until: Friday, 13 November 2026", today) == \
        ("WARN", "valid until Friday, 13 November 2026 (24 days)")
    assert doctor_env.license_status(
        "Your license is valid until: Monday, 19 October 2026\n", today)[0] == "FAIL"
    assert doctor_env.license_status(
        "x\nYour license is valid until: Wednesday, 20 October 2027\n", today)[0] == "PASS"
    status, detail = doctor_env.license_status("License expired.", today)
    assert status == "FAIL" and "valid until" in detail


@pytest.mark.parametrize("method, given, want", [
    ("conformator_bcl", "auto", "on"), ("conformator", "auto", "off"),
    ("rdkit", "auto", "on"), ("confgenx", "auto", "on"), ("conformator_bcl", "off", "off"),
    ("conformator", "on", "on")])
def test_pbfilter_auto_rule(method, given, want):
    assert doctor.resolve_pbfilter(method, given) == want


def test_dock_only_when_configured():
    assert "dock" not in doctor.default_stages({"DOCK_SETUP": "", "DB2TOOL_EXE": ""})
    for k in ("DOCK_SETUP", "DOCK_ENV", "DB2TOOL_EXE"):
        assert doctor.default_stages({k: "x"})[-1] == "dock"


def test_xtb_version_parse():
    assert doctor_env.xtb_version("      * xtb version 6.4.1 (060f3b5) compiled by") == "6.4.1"
    assert doctor_env.xtb_version("normal termination") is None


# ---- doctor.py end to end ---------------------------------------------------------------

def test_host_all_pass_format_and_exit_0(machine):
    rc, out, err = llb_doctor(machine, "--stage", "host")
    assert rc == 0, out + err
    assert set(lines(out).values()) == {"PASS"}
    assert "PASS  host     account            maom99 (user tester)" in out
    assert out.strip().splitlines()[-1] == \
        f"doctor: 9 passed, 0 warnings, 0 failed (config: {machine['tmp'] / 'hpc.env'})"


def test_wrong_account_and_partition_fail_with_fix(machine):
    machine["keys"].update(LLB_ACCOUNT="nobody", LLB_PARTITION="gpu")
    rc, out, _ = llb_doctor(machine, "--stage", "host")
    assert rc == 1
    got = lines(out)
    assert got[("host", "account")] == got[("host", "partition")] == "FAIL"
    assert "      fix: set LLB_ACCOUNT in " in out
    assert "      fix: set LLB_PARTITION in " in out


def test_missing_slurm_tools_and_wget_only_warn(machine):
    for name in ("sacctmgr", "sinfo", "wget", "parallel"):
        (machine["bin"] / name).unlink()
    machine["path"] = f"{machine['bin']}:/usr/bin:/bin"
    machine["keys"].update(WGET="llb_no_wget", PARALLEL="llb_no_parallel")  # CI has wget
    rc, out, _ = llb_doctor(machine, "--stage", "host")
    got = lines(out)
    assert got[("host", "account")] == got[("host", "partition")] == "WARN"
    assert got[("host", "wget")] == "WARN"
    assert rc == 0


def test_no_config_fails_cleanly(machine, tmp_path):
    env = {"PATH": machine["path"], "HOME": str(tmp_path), "LLB_CONFIG": str(tmp_path / "nope.env")}
    p = subprocess.run([sys.executable, str(DOCTOR)], cwd=machine["run"], env=env,
                       capture_output=True, text=True, timeout=60)
    assert p.returncode == 1
    assert lines(p.stdout)[("host", "config")] == "FAIL"
    assert "fix: cp " in p.stdout and "not run without a config: prep,pka,qupkake,db2" in p.stdout
    assert "Traceback" not in p.stderr


def test_usage_errors_exit_64(machine):
    assert llb_doctor(machine, "--stage", "host,bogus")[0] == 64
    assert llb_doctor(machine, "--pbfilter", "maybe")[0] == 64
    assert llb_doctor(machine, "--frobnicate")[0] == 64


def test_activation_failure_is_fail_activate(machine):
    machine["keys"]["PKA_ENV"] = ""
    rc, out, _ = llb_doctor(machine, "--stage", "pka")
    assert rc == 1
    assert lines(out) == {("pka", "activate"): "FAIL"}
    assert "PKA_ENV is empty" in out and "fix: check PKA_SETUP and PKA_ENV" in out


def test_failing_setup_hook_is_fail_activate(machine):
    machine["keys"]["DB2C_SETUP"] = "echo module purge failed >&2; false"
    rc, out, _ = llb_doctor(machine, "--stage", "db2")
    assert lines(out) == {("db2", "activate"): "FAIL"} and rc == 1


def test_default_stages_add_dock_only_when_configured(machine):
    for k in ("QUPKAKE_ENV", "DB2C_ENV", "PKA_ENV"):
        machine["keys"][k] = ""           # every env stage fails to activate, fast
    machine["keys"]["CONDA_BASE"] = machine["tmp"] / "noconda"
    _, out, _ = llb_doctor(machine)
    stages = {s for s, _ in lines(out)}
    assert stages == {"host", "prep", "pka", "qupkake", "db2"}
    machine["keys"]["DB2TOOL_EXE"] = "/no/db2tool"
    _, out, _ = llb_doctor(machine)
    assert lines(out)[("dock", "db2tool")] == "FAIL"


def test_dock_stage_runs_and_checks_db2tool(machine):
    machine["keys"]["DB2TOOL_EXE"] = str(exe(machine["tmp"] / "dock/db2tool", "true"))
    rc, out, _ = llb_doctor(machine, "--stage", "dock")
    assert lines(out) == {("dock", "activate"): "PASS", ("dock", "db2tool"): "PASS"}
    assert rc == 0


def test_json_shape(machine):
    rc, out, _ = llb_doctor(machine, "--stage", "host", "--json", "d.json")
    d = json.loads((machine["run"] / "d.json").read_text())
    assert set(d) == {"llb_root", "config", "host", "when", "checks"}
    assert d["llb_root"] == str(ROOT) and d["config"] == str(machine["tmp"] / "hpc.env")
    assert all(set(c) == {"stage", "check", "status", "detail", "fix"} for c in d["checks"])
    assert len(d["checks"]) == len(lines(out))


def test_node_reruns_under_srun_and_keeps_exit_code(machine):
    log = machine["tmp"] / "srun.args"
    exe(machine["bin"] / "srun", f'printf "%s\\n" "$@" > {log}; exit 3')
    rc, _, _ = llb_doctor(machine, "--node", "--stage", "host")
    assert rc == 3
    assert log.read_text().split() == [
        "-A", "maom99", "-p", "standard", "-c", "1", "--mem=4G", "-t", "00:15:00",
        str(ROOT / "bin/llb"), "doctor", "--stage", "host"]


FAKE_PYTHON = '''#!/bin/bash
# build_ligand stand-in: one NAME.db2.gz per input line; anything else -> real python
if [[ "$1" == *build_ligand_relaxed.py ]]; then
    shift; while [[ $# -gt 0 ]]; do
        case "$1" in -i) in="$2"; shift ;; --outputpath) out="$2"; shift ;; esac; shift
    done
    while read -r smi name; do echo db2 > "$out/$name.db2.gz"; done < "$in"
    exit 0
fi
exec {python} "$@"
'''


def test_smoke_builds_two_molecules(machine):
    fake = machine["tmp"] / "fakepy"
    exe(fake / "python", FAKE_PYTHON.format(python=sys.executable)[len("#!/bin/bash\n"):])
    machine["path"] = f"{fake}:{machine['path']}"
    rc, out, _ = llb_doctor(machine, "--stage", "host", "--smoke", "--method", "conformator")
    assert lines(out)[("smoke", "build_ligand")] == "PASS", out
    assert "pbfilter off" in out and rc == 0
    assert not list(machine["run"].glob("llb_doctor_smoke.*"))


def test_smoke_skipped_after_a_fail(machine):
    machine["keys"]["LLB_ACCOUNT"] = ""
    rc, out, _ = llb_doctor(machine, "--stage", "host", "--smoke")
    assert rc == 1 and "--smoke skipped" in out and ("smoke", "build_ligand") not in lines(out)


# ---- doctor_env.py in a (fake) env ---------------------------------------------------------

@pytest.mark.parametrize("days, want", [(-3, "FAIL"), (10, "WARN"), (300, "PASS")])
def test_db2_unicon_license(machine, days, want):
    exe(machine["bin"] / "unicon", f"cat <<'EOF'\n{banner(days)}EOF")
    got = run_env(machine, "db2", "--method", "conformator", "--pbfilter", "off")
    assert got["unicon"]["status"] == "PASS"
    assert got["unicon_license"]["status"] == want
    if want == "FAIL":
        assert "software.zbh.uni-hamburg.de" in got["unicon_license"]["fix"]
    assert "check_bcl" not in got and "posebusters" not in got   # not asked for


def test_db2_missing_tools_are_fail_lines_not_a_crash(machine):
    got = run_env(machine, "db2", "--method", "conformator_bcl", "--pbfilter", "on",
                  extra={"BUILD_LIGAND_EXE": "/no/build_ligand"})
    assert got["build_ligand"]["status"] == "FAIL"
    assert "BUILD_LIGAND_EXE" in got["build_ligand"]["fix"]
    assert got["unicon_license"]["status"] == "FAIL"
    assert got["singularity"]["status"] == "WARN"
    assert {"check_bcl", "posebusters", "conformator", "amsol"} <= set(got)
    assert all(r["fix"] for r in got.values() if r["status"] == "FAIL")


@pytest.mark.parametrize("version, want", [("6.4.1", "PASS"), ("6.6.1", "FAIL")])
def test_qupkake_xtb_version(machine, version, want):
    xtb = exe(machine["tmp"] / "xtb/bin/xtb", f'echo "      * xtb version {version} (abc) compiled"')
    got = run_env(machine, "qupkake", extra={"XTBPATH": str(xtb)})
    assert got["xtb"]["status"] == want
    assert got["qupkake_cli"]["status"] == "FAIL"
    assert "QUPKAKE_EXE" in got["qupkake_cli"]["fix"]


def test_qupkake_xtb_not_executable(machine):
    xtb = machine["tmp"] / "xtb"
    xtb.write_text("")
    got = run_env(machine, "qupkake", extra={"XTBPATH": str(xtb)})
    assert got["xtb"]["detail"] == f"$XTBPATH={xtb}: not executable"
    assert "set XTBPATH in hpc.env" in got["xtb"]["fix"]


def test_pka_without_molgpka_dir(machine):
    got = run_env(machine, "pka", extra={"MOLGPKA_DIR": ""})
    assert got["molgpka_weights"]["status"] == "FAIL"
    assert "MOLGPKA_DIR" in got["molgpka_weights"]["fix"]
