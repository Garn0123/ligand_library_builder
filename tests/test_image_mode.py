"""Image mode (with_env.sh -> apptainer exec ... micromamba run -n ENV) and the ZBH
shims. A fake `apptainer` records its arguments; a fake `unicon` shows what it
sees under $HOME/ZBH. No container runtime, no ZBH binary."""
import os
import stat
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SHIM_DIR = ROOT / "container" / "shims"


def exe(path, body):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/bash\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


@pytest.fixture
def fake(tmp_path):
    bin_ = tmp_path / "bin"
    # records argv and the env it was given, one line each
    exe(bin_ / "apptainer", 'printf "%s\\n" "$@" > "$FAKE_LOG"; '
        'echo "IN_IMAGE=${LLB_IN_IMAGE:-}" >> "$FAKE_LOG"; '
        'echo "PREPEND=${APPTAINERENV_PREPEND_PATH:-}" >> "$FAKE_LOG"\n')
    image = tmp_path / "llb.sif"
    image.write_text("not really an image")
    return tmp_path, bin_, image


def with_env(tmp_path, bin_, cfg_text, *cmd):
    cfg = tmp_path / "hpc.env"
    cfg.write_text(cfg_text)
    env = {"PATH": f"{bin_}:/usr/bin:/bin", "HOME": str(tmp_path), "USER": "tester",
           "LLB_CONFIG": str(cfg), "LLB_ROOT": str(ROOT), "CONDA_BASE": "/no/conda",
           "FAKE_LOG": str(tmp_path / "argv.txt"), "TMPDIR": str(tmp_path)}
    p = subprocess.run([str(ROOT / "common" / "with_env.sh"), *cmd], env=env, cwd=tmp_path,
                       capture_output=True, text=True, timeout=30)
    log = (tmp_path / "argv.txt").read_text().splitlines() if (tmp_path / "argv.txt").exists() else []
    return p, log


def test_empty_env_runs_the_tool_in_the_image(fake):
    tmp_path, bin_, image = fake
    p, log = with_env(tmp_path, bin_, f"LLB_IMAGE={image}\nZINC_DATA={tmp_path}\n",
                      "QUPKAKE", "python", "-V")
    assert p.returncode == 0, p.stderr
    assert "[QUPKAKE] env=qupkake image=" in p.stdout
    assert log[0] == "exec"
    i = log.index(str(image))
    assert log[i + 1:i + 6] == ["micromamba", "run", "-n", "qupkake", "python"]
    binds = [log[k + 1] for k, a in enumerate(log[:i]) if a == "--bind"]
    assert str(ROOT) in binds and str(tmp_path) in binds
    assert "IN_IMAGE=1" in log
    assert f"PREPEND={SHIM_DIR}" in log


@pytest.mark.parametrize("tool,env", [("DB2C", "db2_converter"), ("PKA", "llb_pka"),
                                      ("PREP", "llb_pka")])
def test_image_env_per_tool(fake, tool, env):
    tmp_path, bin_, image = fake
    p, log = with_env(tmp_path, bin_, f"LLB_IMAGE={image}\n", tool, "true")
    assert p.returncode == 0, p.stderr
    i = log.index(str(image))
    assert log[i + 1:i + 4] == ["micromamba", "run", "-n"] and log[i + 4] == env


def test_a_set_env_keeps_that_tool_native(fake):
    tmp_path, bin_, image = fake
    p, log = with_env(tmp_path, bin_, f"LLB_IMAGE={image}\nQUPKAKE_ENV=myqk\n", "QUPKAKE", "true")
    assert log == []                                   # never reached the container
    assert "activating QUPKAKE failed" in p.stderr or "conda not found" in p.stderr


def test_no_runtime_says_what_to_do(fake):
    tmp_path, bin_, image = fake
    (bin_ / "apptainer").unlink()
    p, _ = with_env(tmp_path, bin_, f"LLB_IMAGE={image}\nLLB_CONTAINER=apptainer\n", "PKA", "true")
    assert p.returncode == 2
    assert "no apptainer or singularity here" in p.stderr and "llb doctor --node" in p.stderr


def test_missing_image_fails(fake):
    tmp_path, bin_, _ = fake
    p, log = with_env(tmp_path, bin_, f"LLB_IMAGE={tmp_path}/nope.sif\n", "PKA", "true")
    assert p.returncode == 2 and "does not exist" in p.stderr and log == []


# ---- the ZBH shims ------------------------------------------------------------------

def run_shim(tmp_path, tool, **env):
    base = {"PATH": "/usr/bin:/bin", "HOME": str(tmp_path / "realhome"), "TMPDIR": str(tmp_path)}
    return subprocess.run([str(SHIM_DIR / tool), "--version"], env={**base, **env},
                          capture_output=True, text=True, timeout=30)


@pytest.fixture
def zbh(tmp_path):
    zdir = tmp_path / "zbh"
    for tool in ("unicon", "conformator"):
        exe(zdir / tool, f'echo "{tool} args=$*"; echo "home=$HOME"; cat "$HOME/ZBH/{tool}.lic"\n')
    lic = tmp_path / "licenses" / "ZBH"
    lic.mkdir(parents=True)
    (lic / "unicon.lic").write_text("[license]\nunicon=KEY1\n")
    (lic / "conformator.lic").write_text("[license]\nconformator=KEY2\n")
    return zdir, lic


@pytest.mark.parametrize("tool,key", [("unicon", "KEY1"), ("conformator", "KEY2")])
def test_shim_runs_zbh_dir_binary_with_the_license(tmp_path, zbh, tool, key):
    zdir, lic = zbh
    p = run_shim(tmp_path, tool, ZBH_DIR=str(zdir), ZBH_LICENSE=str(lic))
    assert p.returncode == 0, p.stderr
    assert f"{tool} args=--version" in p.stdout and f"{tool}={key}" in p.stdout
    home = next(x for x in p.stdout.splitlines() if x.startswith("home="))[5:]
    assert home.startswith(str(tmp_path / "llb-zbh."))
    assert not os.path.exists(home)                    # cleaned up


def test_shim_without_license_keeps_home(tmp_path, zbh):
    zdir, _ = zbh
    p = run_shim(tmp_path, "unicon", ZBH_DIR=str(zdir))
    assert f"home={tmp_path / 'realhome'}" in p.stdout


def test_shim_finds_the_next_one_on_path_and_passes_the_exit_code(tmp_path):
    other = tmp_path / "other"
    exe(other / "unicon", "exit 7\n")
    p = run_shim(tmp_path, "unicon", PATH=f"{SHIM_DIR}:{other}:/usr/bin:/bin")
    assert p.returncode == 7


def test_shim_with_no_binary_says_so(tmp_path):
    p = run_shim(tmp_path, "unicon", PATH=f"{SHIM_DIR}:/usr/bin:/bin")
    assert p.returncode == 127 and "set ZBH_DIR" in p.stderr


# ---- llb doctor's image checks ----------------------------------------------------

def doctor(tmp_path, bin_, cfg_text, label):
    exe(bin_ / "apptainer",
        'case "$1" in\n'
        '  --version) echo "apptainer version 1.3.6" ;;\n'
        f'  inspect) echo "org.llb.toolchain.version: {label}" ;;\n'
        '  *) exit 0 ;;\n'
        'esac\n')
    cfg = tmp_path / "hpc.env"
    cfg.write_text("LLB_ACCOUNT=a\nLLB_PARTITION=p\n" + cfg_text)
    env = {"PATH": f"{bin_}:/usr/bin:/bin", "HOME": str(tmp_path), "USER": "tester",
           "LLB_CONFIG": str(cfg), "CONDA_BASE": "/no/conda", "TMPDIR": str(tmp_path)}
    return subprocess.run(["python3", str(ROOT / "common" / "doctor.py"), "--stage", "pka"],
                          env=env, cwd=tmp_path, capture_output=True, text=True, timeout=60)


def test_doctor_image_version_matches(fake):
    tmp_path, bin_, image = fake
    want = (ROOT / "container" / "TOOLCHAIN_VERSION").read_text().strip()
    p = doctor(tmp_path, bin_, f"LLB_IMAGE={image}\n", want)
    assert f"PASS  host     image" in p.stdout and f"toolchain {want}" in p.stdout


def test_doctor_image_version_mismatch_fails(fake):
    tmp_path, bin_, image = fake
    p = doctor(tmp_path, bin_, f"LLB_IMAGE={image}\n", "1999.01.1")
    assert p.returncode == 1
    assert "FAIL  host     image" in p.stdout and "1999.01.1" in p.stdout


def test_doctor_without_runtime_skips_image_stages(fake):
    tmp_path, bin_, image = fake
    p = doctor(tmp_path, bin_, f"LLB_IMAGE={image}\nLLB_CONTAINER=nosuchruntime\n", "x")
    assert "WARN  host     container" in p.stdout and "llb doctor --node" in p.stdout
    assert " pka " not in p.stdout                     # not attempted here
