#!/usr/bin/env python3
# Copied from DRAP (git@github.com:Garn0123/drap.git), stage4/validate_names.py at
# commit cb6bd9c. DRAP is the source of truth: change it there and re-copy; edit
# nothing below this header. tests/test_validate_names_copy.py checks the copy
# against a DRAP checkout when DRAP_DIR is set.
"""
Check a candidate name set against NAMING_CONTRACT.md, before any conformer is built.

The enumeration step runs in its own environment and cannot see what the matrix
build, the split and the reporting do with a name afterwards. This is that knowledge,
made executable: four invariants, each reported with the failure it prevents.

Run it on the UNION of everything that will be built, not on one file. Two runs that
are each internally unique can still collide with one another, and the collision is
silent -- db2 keeps the last 16 characters and one molecule simply becomes another.

    python3 stage4/validate_names.py --smi protomers.smi
    python3 stage4/validate_names.py --smi a.smi --smi b.smi
    python3 stage4/validate_names.py --names names.txt --col 0
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

NAME_LEN = 16
SENTINEL = "C"
BASE62 = set("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
SAFE_NAME = re.compile(r"^[A-Za-z0-9_-]+$")
DB2_TAIL = 16          # db2 writes name[-16:] (mol2db2/hierarchy.py)

FAILURES: list[str] = []


def fail(inv: str, msg: str, examples=()):
    FAILURES.append(inv)
    print(f"\n  [FAIL] {inv}\n         {msg}")
    for e in list(examples)[:5]:
        print(f"           {e}")


def ok(inv: str, msg: str):
    print(f"  [ok  ] {inv}  -- {msg}")


def read_names(paths, col):
    names, seen_files = [], []
    for p in paths:
        p = Path(p)
        if not p.exists():
            raise SystemExit(f"no such file: {p}")
        n = 0
        for line in p.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if col >= len(parts):
                continue
            names.append(parts[col])
            n += 1
        seen_files.append((p.name, n))
    return names, seen_files


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smi", action="append", default=[],
                    help="2-column .smi; the name is column 1. Repeat for the union.")
    ap.add_argument("--names", action="append", default=[],
                    help="plain list of names, one per line")
    ap.add_argument("--col", type=int, default=None,
                    help="0-based column holding the name (default: 1 for --smi, "
                         "0 for --names)")
    args = ap.parse_args()

    if not args.smi and not args.names:
        raise SystemExit("give at least one --smi or --names")

    names, files = [], []
    if args.smi:
        n, f = read_names(args.smi, 1 if args.col is None else args.col)
        names += n; files += f
    if args.names:
        n, f = read_names(args.names, 0 if args.col is None else args.col)
        names += n; files += f

    print(f"{len(names):,} names from {len(files)} file(s)")
    for nm, c in files:
        print(f"  {nm}: {c:,}")
    if not names:
        raise SystemExit("no names read -- wrong --col?")
    print()

    # ---- invariant 1: the LAST 16 characters are unique --------------------
    tails = Counter(n[-DB2_TAIL:] for n in names)
    dupes = {t: c for t, c in tails.items() if c > 1}
    if dupes:
        by_tail = defaultdict(list)
        for n in names:
            t = n[-DB2_TAIL:]
            if t in dupes:
                by_tail[t].append(n)
        ex = [f"{t!r} <- {sorted(set(v))[:3]}" for t, v in list(by_tail.items())[:5]]
        fail("1 last-16 unique",
             f"{len(dupes):,} truncated name(s) are shared by {sum(dupes.values()):,} "
             f"molecules. db2 keeps name[-16:], so these become ONE row -- the exact "
             f"collapse that cost +2.560 RMSE across 137,140 molecules.", ex)
    else:
        ok("1 last-16 unique", f"{len(tails):,} distinct tails, no collisions")

    # ---- format: length and sentinels -------------------------------------
    wrong_len = [n for n in names if len(n) != NAME_LEN]
    if wrong_len:
        fail("1 length == 16",
             f"{len(wrong_len):,} name(s) are not {NAME_LEN} characters. Anything "
             f"longer is truncated from the FRONT by db2; anything shorter loses the "
             f"sentinel check that detects truncation.",
             [f"{n!r} (len {len(n)})" for n in wrong_len])
    else:
        ok("1 length == 16", "every name is exactly 16 characters")

    bad_sent = [n for n in names
                if len(n) == NAME_LEN and not (n[0] == SENTINEL and n[-1] == SENTINEL)]
    if bad_sent:
        fail("1 sentinels", f"{len(bad_sent):,} name(s) do not start AND end with "
             f"{SENTINEL!r}. The sentinels are what make a damaged name detectable "
             f"instead of silently aliasing onto another molecule.", bad_sent)
    else:
        ok("1 sentinels", f"all names start and end with {SENTINEL!r}")

    body_bad = [n for n in names
                if len(n) == NAME_LEN and not set(n[1:15]) <= BASE62]
    if body_bad:
        fail("1 base62 body", f"{len(body_bad):,} name(s) have non-base62 characters "
             f"between the sentinels.", body_bad)
    else:
        ok("1 base62 body", "base id, protomer and stereo fields are all base62")

    # ---- invariant 2: filesystem- and shell-safe --------------------------
    unsafe = [n for n in names if not SAFE_NAME.match(n)]
    if unsafe:
        fail("2 shell-safe", f"{len(unsafe):,} name(s) fall outside [A-Za-z0-9_-]. "
             f"db2_converter interpolates the name unquoted into a path and into "
             f"`rm -r {{name}}`.", unsafe)
    else:
        ok("2 shell-safe", "all names match [A-Za-z0-9_-]+")

    # ---- invariant 3: normalize_id round-trips, base_id groups -------------
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "stage1"))
    try:
        from ids import normalize_id
    except ImportError as e:
        fail("3 normalize_id", f"could not import stage1/ids.py ({e})")
    else:
        changed = [n for n in names if normalize_id(n) != n]
        if changed:
            fail("3 normalize_id round-trip",
                 f"{len(changed):,} name(s) are rewritten by normalize_id. Each "
                 f"protomer must be its own row, so nothing should be stripped -- a "
                 f"rewrite here merges protomers back together.",
                 [f"{n!r} -> {normalize_id(n)!r}" for n in changed])
        else:
            ok("3 normalize_id round-trip", "every name maps to itself")

    good = [n for n in names if len(n) == NAME_LEN]
    groups = defaultdict(set)
    for n in good:
        groups[n[1:13]].add(n)
    sizes = Counter(len(v) for v in groups.values())
    ok("3 base_id grouping",
       f"{len(groups):,} distinct base ids; "
       f"{sizes.get(1,0):,} with one variant, "
       f"{sum(c for s,c in sizes.items() if s>1):,} with several "
       f"(max {max(sizes) if sizes else 0})")

    # ---- invariant 4 is a data check, not a name check ---------------------
    print(f"\n  [note] 4 charge recoverability is NOT checkable from names alone.")
    print(f"         The acceptance test runs after the matrix is built: every newly")
    print(f"         named molecule should come out ambiguity_class = 0, because")
    print(f"         protomers now have their own names and nothing collapses.")
    print(f"         Any class 1 or 2 in the new data means the naming did not work.")

    print()
    if FAILURES:
        raise SystemExit(f"{len(FAILURES)} invariant(s) FAILED: {', '.join(FAILURES)}\n"
                         f"Do not generate conformers from this name set.")
    print("all name invariants hold -- safe to build conformers")


if __name__ == "__main__":
    main()
