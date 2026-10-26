"""common/env.sh llb_load_config: where hpc.env is found, the users/$USER.env
overlay, LLB_SITE, and the FATAL message. A fake site directory with a copy of
the checkout's env.sh; nothing is activated."""
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SHOW = ("LLB_CONFIG", "LLB_SITE", "LLB_USER_CONFIG", "LLB_ACCOUNT", "LLB_MAIL",
        "ZBH_DIR", "ZBH_LICENSE")


@pytest.fixture
def site(tmp_path):
    """$LLB_SITE/ligand_library_builder: just enough of a checkout for llb_find_root."""
    site = tmp_path / "site"
    llb = site / "ligand_library_builder"
    (llb / "common").mkdir(parents=True)
    (llb / "config").mkdir()
    shutil.copy(ROOT / "common" / "env.sh", llb / "common" / "env.sh")
    shutil.copy(ROOT / "config" / "hpc.env.example", llb / "config" / "hpc.env.example")
    return site


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def load(site, **env):
    """Run llb_load_config under `set -euo pipefail`; (rc, {var: value}, stderr)."""
    script = ('set -euo pipefail; source "$LLB_ROOT/common/env.sh"; llb_load_config; '
              + "; ".join(f'echo "{k}=${{{k}:-}}"' for k in SHOW))
    base = {"PATH": "/usr/bin:/bin", "HOME": str(site), "USER": "tester",
            "LLB_ROOT": str(site / "ligand_library_builder"), "CONDA_BASE": "/no/conda"}
    p = subprocess.run(["bash", "-c", script], env={**base, **env}, cwd=site,
                       capture_output=True, text=True, timeout=30)
    got = dict(x.split("=", 1) for x in p.stdout.splitlines() if "=" in x)
    return p.returncode, got, p.stderr


def test_site_hpc_env_beside_the_checkout_wins(site):
    write(site / "hpc.env", "LLB_ACCOUNT=site_acct\nZBH_DIR=$LLB_SITE/zbh\n")
    write(site / "ligand_library_builder/config/hpc.env", "LLB_ACCOUNT=old_acct\n")
    rc, got, err = load(site)
    assert rc == 0, err
    assert got["LLB_ACCOUNT"] == "site_acct"
    assert Path(got["LLB_CONFIG"]).resolve() == (site / "hpc.env").resolve()
    assert got["LLB_SITE"] == str(site)
    assert got["ZBH_DIR"] == f"{site}/zbh"                    # hpc.env can use $LLB_SITE
    assert got["ZBH_LICENSE"] == f"{site}/licenses/tester/ZBH"   # the site default
    assert got["LLB_USER_CONFIG"] == ""


def test_old_config_hpc_env_still_works(site):
    cfg = write(site / "ligand_library_builder/config/hpc.env", "LLB_ACCOUNT=old_acct\n")
    rc, got, err = load(site)
    assert rc == 0, err
    assert got["LLB_ACCOUNT"] == "old_acct" and got["LLB_CONFIG"] == str(cfg)
    assert got["LLB_SITE"] == "" and got["ZBH_LICENSE"] == ""


def test_llb_config_overrides_both(site, tmp_path):
    write(site / "hpc.env", "LLB_ACCOUNT=site_acct\n")
    own = write(tmp_path / "elsewhere/hpc.env", "LLB_ACCOUNT=own_acct\n")
    rc, got, _ = load(site, LLB_CONFIG=str(own))
    assert got["LLB_ACCOUNT"] == "own_acct" and got["LLB_CONFIG"] == str(own)
    assert got["LLB_SITE"] == "" and rc == 0


def test_inherited_llb_config_keeps_llb_site(site):
    """A SLURM job inherits LLB_CONFIG=$LLB_SITE/hpc.env and must see the same site."""
    write(site / "hpc.env", "LLB_ACCOUNT=site_acct\n")
    _, first, _ = load(site)
    rc, again, _ = load(site, LLB_CONFIG=first["LLB_CONFIG"])
    assert rc == 0 and again["LLB_SITE"] == first["LLB_SITE"] == str(site)


def test_user_file_overrides_the_site(site):
    write(site / "hpc.env", "LLB_ACCOUNT=site_acct\nLLB_MAIL=\nZBH_LICENSE=\n")
    user = write(site / "users/tester.env",
                 "LLB_ACCOUNT=my_acct\nLLB_MAIL=me@example.org\nZBH_LICENSE=/home/tester/ZBH\n")
    write(site / "users/someone_else.env", "LLB_ACCOUNT=not_mine\n")
    rc, got, err = load(site)
    assert rc == 0, err
    assert got["LLB_ACCOUNT"] == "my_acct" and got["LLB_MAIL"] == "me@example.org"
    assert got["ZBH_LICENSE"] == "/home/tester/ZBH"                # not the site default
    assert Path(got["LLB_USER_CONFIG"]).resolve() == user.resolve()


def test_user_file_beside_an_llb_config_file(site, tmp_path):
    own = write(tmp_path / "cfg/hpc.env", "LLB_ACCOUNT=own_acct\n")
    write(tmp_path / "cfg/users/tester.env", "LLB_ACCOUNT=my_acct\n")
    rc, got, _ = load(site, LLB_CONFIG=str(own))
    assert rc == 0 and got["LLB_ACCOUNT"] == "my_acct"
    assert got["LLB_USER_CONFIG"] == str(tmp_path / "cfg/users/tester.env")


def test_no_config_fatal_lists_every_place_looked(site):
    rc, _, err = load(site)
    assert rc == 2
    first = err.splitlines()[0]
    assert first.startswith("FATAL: no hpc.env; looked for: ")
    assert str(site / "hpc.env") in first
    assert str(site / "ligand_library_builder/config/hpc.env") in first
    assert "cp " in err and "hpc.env.example" in err


def test_missing_llb_config_is_not_silently_skipped(site, tmp_path):
    write(site / "hpc.env", "LLB_ACCOUNT=site_acct\n")
    rc, _, err = load(site, LLB_CONFIG=str(tmp_path / "nope.env"))
    assert rc == 2
    assert err.splitlines()[0] == f"FATAL: no hpc.env; looked for: LLB_CONFIG={tmp_path / 'nope.env'}"
