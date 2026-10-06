#!/usr/bin/env python3
"""
db2_to_mol2.py -- one DOCK-ready mol2 per molecule, taken FROM the db2.

For mol2 (flexible-ligand) docking with DOCK6, alongside or instead of db2
docking, without a second 3D tool such as CORINA. The mol2 is not a new
conformation: it is a conformer read back out of the db2 by db2tool tomol2
(dock6_claude, pinned 294d235), which uses DOCK's own db2 reader and mol2
writer. So the mol2 carries, per atom, exactly what the db2 carries:

    coordinates   one conformer set from the db2 (lowest internal energy)
    charges       the db2's AMSOL partial charges -- same as db2 docking uses
    atom types    SYBYL types from the db2
    title         the db2 M-line name, i.e. the 16-character contract name

What this wrapper adds to tomol2:

  * ONE conformer per molecule. db2_converter writes one db2 hierarchy per rigid
    fragment and concatenates them into NAME.db2.gz, each a separate record with
    the same name, so `tomol2 --best` emits one molecule per hierarchy. Here the
    lowest HDB_Internal_energy1 across all of a name's hierarchies is kept.
  * Never silently missing. A molecule whose every set is flagged broken is
    skipped by tomol2 --best; it is retried with --keep-broken and its lowest
    broken set kept, flagged in the report (as `db2tool audit --write-mol2`
    does). Mol2 docking builds its own conformers from this one, so a broken
    set is still a valid starting geometry. --no-broken-fallback turns this off.
  * The invariant-4 check, again. With --library, the mol2's summed partial
    charge must match the protomer's net formal charge in library.tsv, and every
    library name must come out. Exit 1 otherwise.

Usage (db2tool from config/hpc.env via common/with_env.sh DOCK):
    common/with_env.sh DOCK python3 run_qupkake/db2_to_mol2.py build/ \\
        --library library/library.tsv -o library/library.mol2

    inputs: .db2.gz / .db2 files, or directories searched for them
    writes: -o FILE (multi-molecule mol2) and FILE.tsv (one row per molecule)

stdlib only.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HEADER = re.compile(r"^##########\s+([A-Za-z0-9_]+):\s*(.*?)\s*$")
# db2tool: "<path>\tSKIPPED\t<name>: every set is flagged broken" -- the name stops
# at the colon (contract and ZINC names never contain one)
SKIPPED = re.compile(r"^(.*?)\tSKIPPED\t([^:\s]+)")


def find_db2(inputs: list[str]) -> list[Path]:
    out = []
    for s in inputs:
        p = Path(s)
        if p.is_dir():
            out += sorted(p.rglob("*.db2.gz")) + sorted(p.rglob("*.db2"))
        elif p.exists():
            out.append(p)
        else:
            raise SystemExit(f"no such file or directory: {p}")
    if not out:
        raise SystemExit("no .db2 / .db2.gz inputs found")
    return out


def run_tomol2(db2tool: str, files: list[Path], keep_broken: bool):
    """(records, skipped names). One record per db2 hierarchy emitted."""
    cmd = [db2tool, "tomol2", "--best", "-o", "-", "-q"]
    if keep_broken:
        cmd.append("--keep-broken")
    res = subprocess.run(cmd + [str(f) for f in files],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise SystemExit(f"{' '.join(cmd[:6])} ... exited {res.returncode}:\n{res.stderr[-2000:]}")
    skipped = {m.group(2) for m in map(SKIPPED.match, res.stderr.splitlines()) if m}
    return parse_records(res.stdout), skipped, res.stderr


def parse_records(text: str) -> list[dict]:
    """Split tomol2 output on its '##########' header blocks."""
    records, cur, in_header = [], None, False
    for line in text.splitlines(keepends=True):
        m = HEADER.match(line)
        if m:
            if not in_header:
                cur = {"header": {}, "lines": []}
                records.append(cur)
                in_header = True
            cur["header"][m.group(1)] = m.group(2)
            cur["lines"].append(line)
            continue
        in_header = False
        if cur is not None:
            cur["lines"].append(line)
    for r in records:
        h = r["header"]
        r["name"] = h.get("Name", "")
        r["e1"] = float(h.get("HDB_Internal_energy1", "inf"))
        r["e2"] = float(h.get("HDB_Internal_energy2", "nan"))
        r["conf"] = h.get("HDB_Conf_Number", "")
        r["broken"] = h.get("HDB_Set_Broken", "") == "yes"
        r["charge"], r["n_atoms"] = charge_sum(r["lines"])
    return records


def charge_sum(lines: list[str]) -> tuple[float, int]:
    """Sum of the partial-charge column of the ATOM section."""
    total, n, in_atoms = 0.0, 0, False
    for line in lines:
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.startswith("@<TRIPOS>ATOM")
            continue
        if in_atoms and line.strip():
            total += float(line.split()[-1])
            n += 1
    return total, n


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="db2 files or directories holding them")
    ap.add_argument("-o", "--output", required=True, help="multi-molecule mol2 to write")
    ap.add_argument("--library", default=None,
                    help="library.tsv from assign_names.py: charge + completeness check")
    ap.add_argument("--db2tool", default=os.environ.get("DB2TOOL_EXE") or "db2tool",
                    help="db2tool binary (default: $DB2TOOL_EXE from config/hpc.env, "
                         "else db2tool on PATH)")
    ap.add_argument("--no-broken-fallback", action="store_true",
                    help="leave out molecules whose every set is flagged broken")
    ap.add_argument("--charge-tol", type=float, default=0.05)
    args = ap.parse_args(argv)

    files = find_db2(args.inputs)
    records, skipped, _ = run_tomol2(args.db2tool, files, keep_broken=False)

    # every-set-broken molecules: one more pass, over just the files that hold them
    if skipped and not args.no_broken_fallback:
        stems = {f: f.name.split(".")[0] for f in files}
        retry = [f for f in files if stems[f] in skipped] or files
        more, _, _ = run_tomol2(args.db2tool, retry, keep_broken=True)
        records += [r for r in more if r["name"] in skipped]

    by_name: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        by_name[r["name"]].append(r)
    chosen = {n: min(rs, key=lambda r: (r["broken"], r["e1"])) for n, rs in by_name.items()}

    expected: dict[str, int] = {}
    if args.library:
        with open(args.library) as fh:
            expected = {r["name"]: int(r["net_charge"])
                        for r in csv.DictReader(fh, delimiter="\t")}

    out = Path(args.output)
    rows, bad = [], 0
    with open(out, "w") as fh:
        for name in sorted(chosen):
            r = chosen[name]
            fh.write("".join(r["lines"]).rstrip("\n") + "\n\n")
            want = expected.get(name)
            status = "OK"
            if want is not None and abs(r["charge"] - want) > args.charge_tol:
                status, bad = "CHARGE_MISMATCH", bad + 1
            elif args.library and want is None:
                status = "NOT_IN_LIBRARY"
            rows.append([name, len(by_name[name]), r["conf"], f"{r['e1']:.3f}",
                         f"{r['e2']:.3f}", int(r["broken"]), r["n_atoms"],
                         f"{r['charge']:+.4f}", "" if want is None else f"{want:+d}", status])
    missing = sorted(set(expected) - set(chosen))
    with open(str(out) + ".tsv", "w") as fh:
        fh.write("name\tdb2_hierarchies\tconf_set\tinternal_energy1\tinternal_energy2\t"
                 "broken_set\tatoms\tmol2_charge\texpected_charge\tstatus\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
        fh.writelines(f"{n}\t0\t\t\t\t\t\t\t{expected[n]:+d}\tMISSING\n" for n in missing)

    n_broken = sum(1 for r in chosen.values() if r["broken"])
    multi = sum(1 for rs in by_name.values() if len(rs) > 1)
    print(f"{len(files):,} db2 file(s) -> {len(chosen):,} molecule(s) in {out}")
    print(f"  {multi:,} had several db2 hierarchies (one kept: lowest internal energy)")
    print(f"  {n_broken:,} came from a broken set (every set broken; flagged in the .tsv)")
    if args.no_broken_fallback and skipped:
        print(f"  {len(skipped):,} left out: every set broken (--no-broken-fallback)")
    if args.library:
        print(f"  charge check vs library.tsv: {bad:,} mismatch(es), "
              f"{len(missing):,} library name(s) missing")
    print(f"  report -> {out}.tsv")
    return 1 if (bad or missing) else 0


if __name__ == "__main__":
    sys.exit(main())
