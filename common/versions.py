#!/usr/bin/env python3
"""
versions.py -- print, as JSON, the versions that decide what a library run built.

    llb env PREP python $(llb root)/common/versions.py
    llb env PKA  python $(llb root)/common/versions.py

Run once per environment by `llb pipeline start` (each env has its own RDKit);
whatever cannot be imported in that env is simply left out.
"""
from __future__ import annotations

import json
import platform
import sys
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT / "run_qupkake"), str(ROOT / "pka_triage")]


def main() -> int:
    out: dict = {"python": platform.python_version()}
    for dist in ("rdkit", "torch", "pandas", "numpy", "qupkake", "dimorphite_dl"):
        try:
            out[dist] = metadata.version(dist)
        except metadata.PackageNotFoundError:
            try:
                out[dist] = __import__(dist).__version__
            except Exception:  # noqa: BLE001 -- not in this env
                pass
    try:
        import qupkake_protomers as q
        out.update(site_filter_version=str(q.SITE_FILTER_VERSION),
                   coupling_version=q.COUPLING_VERSION,
                   aromatic_n_version=q.AROMATIC_N_VERSION)
    except Exception:  # noqa: BLE001
        pass
    try:
        import molgpka
        out["molgpka_weights_sha256"] = molgpka.KNOWN_WEIGHTS
    except Exception:  # noqa: BLE001
        pass
    print(json.dumps(out, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
