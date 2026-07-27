"""Shared test fixtures and path setup.

Makes the pipeline package (``db2common`` and the numbered scripts) and the test
helpers (``db2gen``) importable, and provides ``run_script`` for driving the CLI
stages as subprocesses — the honest way to test their real argparse/exit
behavior without entangling sys.argv/sys.exit into the test process.
"""

import subprocess
import sys
from pathlib import Path

import pytest

PIPE_DIR = Path(__file__).resolve().parent.parent
TESTS_DIR = Path(__file__).resolve().parent

for _p in (str(PIPE_DIR), str(TESTS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# The stage scripts live under serial/ and parallel/; db2common.py stays at the
# repo root. Callers pass a bare basename and we locate it.
_SCRIPT_DIRS = ("serial", "parallel")


def _resolve_script(name):
    if "/" in name:                       # explicit path given
        return PIPE_DIR / name
    for sub in _SCRIPT_DIRS:
        cand = PIPE_DIR / sub / name
        if cand.exists():
            return cand
    return PIPE_DIR / name                 # fallback: repo root


def _run_script(name, *args, check=True, cwd=None):
    """Run a pipeline script with the current interpreter.

    ``-u`` keeps stderr unbuffered (matching how submit.slurm runs them). With
    ``check=True`` (default) a nonzero exit raises AssertionError with the
    captured output; pass ``check=False`` to assert on a failure yourself.
    """
    cmd = [sys.executable, "-u", str(_resolve_script(name))] + [str(a) for a in args]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    if check and res.returncode != 0:
        raise AssertionError(
            "{} exited {}\n--- stdout ---\n{}\n--- stderr ---\n{}".format(
                name, res.returncode, res.stdout, res.stderr))
    return res


@pytest.fixture
def run_script():
    return _run_script
