#!/usr/bin/env python3
"""check_bcl.py -- can db2_converter find BCL::Conf? Exit 0 if it can.

Run inside the db2_converter env, after DB2C_SETUP (which sets BCL_ROOT and
PATH on our clusters). `-m conformator_bcl` with a missing BCL does not stop
build_ligand: every molecule just loses its BCL half, or fails as 1generate,
which a gate on build rate would only catch hours later. So submit_db2.sh and
each array task check up front.

db2_converter (63d6656) has no single documented BCL key, so this matches every
config entry whose level, name or value mentions BCL and requires each absolute
path among them to exist.
"""
from __future__ import annotations

import os
import sys


def bcl_entries(config) -> list[tuple[str, str, str]]:
    out = []
    for level, section in config.items():
        if not hasattr(section, "items"):
            continue
        for name, value in section.items():
            if "bcl" in f"{level} {name} {value}".lower():
                out.append((str(level), str(name), str(value)))
    return out


def main() -> int:
    try:
        from db2_converter.config import config
    except Exception as e:                       # noqa: BLE001 -- report anything
        print(f"check_bcl: cannot import db2_converter.config: {e}", file=sys.stderr)
        return 1
    entries = bcl_entries(config)
    if not entries:
        print("check_bcl: no BCL entry in db2_converter.config (levels: "
              f"{', '.join(map(str, config))})", file=sys.stderr)
        return 1
    bad = 0
    for level, name, value in entries:
        missing = value.startswith("/") and not os.path.exists(value)
        bad += missing
        print(f"  bcl  {level}.{name} = {value}{'   MISSING' if missing else ''}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
