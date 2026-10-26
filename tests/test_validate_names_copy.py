"""run_qupkake/validate_names.py is a copy of DRAP's stage4/validate_names.py
(DRAP is the source of truth). With DRAP_DIR set, the copy must match it apart
from the header comment under the shebang; a difference means re-copy."""
import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
COPY = ROOT / "run_qupkake" / "validate_names.py"
DRAP_DIR = os.environ.get("DRAP_DIR", "")


def without_header(text):
    """Drop the comment lines right after the shebang (DRAP's line 2 is the docstring)."""
    lines = text.splitlines(keepends=True)
    i = 1
    while i < len(lines) and lines[i].startswith("#"):
        i += 1
    return lines[0] + "".join(lines[i:])


@pytest.mark.skipif(not DRAP_DIR, reason="needs DRAP_DIR (a DRAP checkout)")
def test_copy_matches_drap():
    drap = Path(DRAP_DIR) / "stage4" / "validate_names.py"
    assert drap.is_file(), f"no {drap}"
    assert without_header(COPY.read_text()) == drap.read_text(), (
        f"{COPY} differs from {drap}: re-copy it and update the commit in its header")


def test_copy_names_its_source():
    head = COPY.read_text().splitlines()[1:6]
    assert any("Garn0123/drap" in x and "stage4/validate_names.py" in x for x in head)
    assert any("source of truth" in x for x in head)
