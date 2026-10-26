"""Shared setup for the tests of the kept tools (naming/, docking/, common/env.sh).

    pytest                       # from the repo root; mol_compiler/ has its own suite

The library-construction tests (pka_triage/, run_qupkake/, pipeline/) went with
that code; they are at tag llb-full-2026-10.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
for sub in ("docking", "naming"):
    sys.path.insert(0, str(ROOT / sub))
