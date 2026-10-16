#!/usr/bin/env python3
"""
verify_db2_build.py
===================

Cross-checks a db2_converter run against the protomer set that was requested.

Why this is needed
------------------
db2_converter's own `chemistrycheck` compares each generated conformer to the
input SMILES by **InChI** (db2_converter/pipeline.py):

    if Chem.MolToInchi(mol) == Chem.MolToInchi(Chem.MolFromSmiles(canonical)):

Standard InChI normalizes mobile hydrogens and zwitterions, so a molecule that
comes back from the conformer generator / UNICON round trip in a *different*
protonation state can still pass -- glycine's zwitterion and neutral form share
an InChI, and the same is true for many amino-acid-like and heteroaromatic
tautomer pairs.  For a run whose entire point is enumerating protonation
states, that check is exactly the one that cannot be relied on.

This script instead reads the AMSOL total charge that mol2db2 writes into the
db2 itself (the second M line, `solvdata.totalCharge`) and compares it to the
net formal charge you asked for.

Usage
-----
    python verify_db2_build.py library/library.tsv --outputpath db2_run/out
    (or protomers.csv; --outputpath and its immediate subdirectories are
    searched, --deep for the whole tree, and every *.faillist found is read
    unless --faillist is given)

Exit 1 on a charge mismatch, a renamed or mangled name, or a molecule MISSING
with no faillist entry -- a silent loss. Molecules build_ligand reported as
failed are BUILD_FAILED: recorded attrition, not a failed run.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# NAMING_CONTRACT.md: 'C' + 12 base id + protomer + stereo + 'C'
CONTRACT = re.compile(r"^C[0-9A-Za-z]{14}C$")


def read_db2_header(path: Path) -> dict:
    """Pull name / total charge / SMILES off the leading M lines of a db2."""
    info: dict = {"file": path.name}
    m_lines = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("M "):
                m_lines.append(line)
                if len(m_lines) >= 3:
                    break
            elif m_lines:
                break
    if len(m_lines) >= 1:
        info["db2_name"] = m_lines[0][2:18].strip()
    if len(m_lines) >= 2:
        try:
            info["total_charge"] = float(m_lines[1].split()[1])
        except (IndexError, ValueError):
            info["total_charge"] = None
    if len(m_lines) >= 3:
        info["smiles"] = m_lines[2][2:].strip()
    return info


def scan(outdir: Path, deep: bool) -> tuple[list[Path], list[Path]]:
    """(*.db2.gz, *.faillist) under outdir.

    build_ligand writes both flat into --outputpath, and a submit_db2.sh run's
    out/ holds one such directory per shard (out/NNNN/). So by default only
    outdir and its immediate subdirectories are listed: ~1,000 listings for a
    1,000-shard run. A recursive walk also descends every molecule's working
    directory and on NFS took over two hours on 22k names (2026-10-16); it also
    finds db2_converter's working copies (NAME/db2/output.N.db2.gz), which are
    not results. --deep restores it.
    """
    if deep:
        return sorted(outdir.rglob("*.db2.gz")), sorted(outdir.rglob("*.faillist"))
    db2, fails = [], []
    dirs = [outdir]
    with os.scandir(outdir) as it:
        dirs += [Path(e.path) for e in it if e.is_dir(follow_symlinks=False)]
    for d in dirs:
        with os.scandir(d) as it:
            for e in it:
                if e.name.endswith(".db2.gz") and e.is_file():
                    db2.append(Path(e.path))
                elif e.name.endswith(".faillist") and e.is_file():
                    fails.append(Path(e.path))
    return sorted(db2), sorted(fails)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("protomers_csv",
                    help="protomers.csv (qupkake_protomers.py) or library.tsv (assign_names.py)")
    ap.add_argument("--outputpath", required=True, action="append",
                    help="repeatable: e.g. the main run's out/ and a fallback run's out/ "
                         "(a name that failed in one and built in another is OK). "
                         "db2_converter --outputpath, or a submit_db2.sh run's out/ "
                         "(searched recursively)")
    ap.add_argument("--faillist", action="append", default=[],
                    help=".faillist written by build_ligand (repeatable); default: "
                         "every *.faillist under --outputpath")
    ap.add_argument("--deep", action="store_true",
                    help="search --outputpath recursively (slow on NFS; default: it "
                         "and its immediate subdirectories, where build_ligand writes)")
    ap.add_argument("--threads", type=int, default=32,
                    help="parallel db2 reads (I/O-bound; NFS latency, not CPU)")
    ap.add_argument("--charge-tol", type=float, default=0.05,
                    help="tolerance when comparing AMSOL charge to formal charge")
    ap.add_argument("-o", "--report", default="db2_verification.csv")
    args = ap.parse_args(argv)

    outdirs = [Path(o) for o in args.outputpath]
    # protomers.csv (qupkake_protomers.py) or library.tsv (assign_names.py)
    with open(args.protomers_csv) as fh:
        tsv = args.protomers_csv.endswith(".tsv")
        expected = [{"protomer_name": r.get("protomer_name") or r["name"],
                     "parent_name": (r.get("parent_name") or r.get("input_name")
                                     or r.get("zinc_id", "")),
                     "net_charge": r["net_charge"]}
                    for r in csv.DictReader(fh, delimiter="\t" if tsv else ",")]

    fail_reasons: dict[str, str] = {}
    # One index over the whole tree: a submit_db2.sh run keeps each shard's
    # results in out/NNNN/, and globbing per name would be N directory scans.
    all_db2, found_faillists = [], []
    for outdir in outdirs:
        db2s, fls = scan(outdir, args.deep)
        all_db2 += db2s
        found_faillists += fls
        print(f"  {len(db2s):,} db2 file(s), {len(fls):,} faillist(s) under {outdir}",
              flush=True)
    exact_index: dict[str, list[Path]] = {}
    split_index: dict[str, list[Path]] = {}
    for f in all_db2:
        stem = f.name[:-len(".db2.gz")]
        exact_index.setdefault(stem, []).append(f)
        m = re.match(r"^(.*)\.\d+$", stem)        # --checkstereo split: NAME.0, NAME.1
        if m:
            split_index.setdefault(m.group(1), []).append(f)

    faillists = [Path(f) for f in args.faillist] or found_faillists
    for fl in faillists:
        if not fl.exists():
            continue
        for line in fl.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 4:
                fail_reasons[parts[1]] = parts[3]

    # Every db2's size and M lines, read up front with a thread pool: each read
    # is a few NFS round trips, and one at a time ran ~0.5 s per name on Turbo
    # (2026-10-16). The work is waiting, not computing, so threads hide it.
    def probe(f: Path):
        try:
            size = f.stat().st_size
            return f, size, (read_db2_header(f) if size else {})
        except (OSError, EOFError) as e:
            return f, None, {"error": f"{type(e).__name__}: {e}"}
    probed: dict[Path, tuple] = {}
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        for k, (f, size, hdr) in enumerate(pool.map(probe, all_db2), 1):
            probed[f] = (size, hdr)
            if k % 5000 == 0:
                print(f"  read {k:,} of {len(all_db2):,} db2 headers "
                      f"({time.time() - t0:.0f} s)", flush=True)
    print(f"  read {len(all_db2):,} db2 headers in {time.time() - t0:.0f} s "
          f"({args.threads} threads)", flush=True)

    rows, tally = [], Counter()
    for exp in expected:
        name = exp["protomer_name"]
        want_q = int(exp["net_charge"])

        # --checkstereo splits unspecified stereocentres into NAME.0, NAME.1 ...
        exact = exact_index.get(name, [])
        hits = exact or split_index.get(name, [])
        # A contract name (NAMING_CONTRACT.md) must reach the db2 unchanged: a
        # .N suffix makes it 18 characters and mol2db2 keeps only the last 16.
        contract = CONTRACT.match(name) is not None

        if not hits:
            # A failure build_ligand reported is attrition, recorded and
            # expected; a molecule gone with no reason anywhere is a silent
            # loss, and only that fails the run.
            reported = name in fail_reasons
            tally["build_failed" if reported else "missing"] += 1
            rows.append({"protomer_name": name, "parent": exp["parent_name"],
                         "expected_charge": want_q,
                         "status": "BUILD_FAILED" if reported else "MISSING",
                         "detail": fail_reasons.get(name, "no db2 and no faillist entry"),
                         "db2_file": "", "db2_name": "", "amsol_charge": ""})
            continue

        for h in hits:
            size, hdr = probed[h]
            if size is None:
                tally["unreadable"] += 1
                status, detail = "UNREADABLE", hdr["error"]
            elif size == 0:
                tally["empty"] += 1
                status, detail, hdr = "EMPTY", "zero-length db2", {}
            else:
                q = hdr.get("total_charge")
                if q is None:
                    tally["unreadable"] += 1
                    status, detail = "UNREADABLE", "no parsable M charge line"
                elif contract and not exact:
                    tally["renamed"] += 1
                    status = "RENAMED"
                    detail = (f"db2_converter split this into {h.name}; stereo was not "
                              f"closed upstream (run check_db2_stereo.py)")
                elif contract and hdr.get("db2_name") != name:
                    tally["name_mismatch"] += 1
                    status = "NAME_MISMATCH"
                    detail = f"M line holds {hdr.get('db2_name')!r}, requested {name!r}"
                elif abs(q - want_q) > args.charge_tol:
                    tally["charge_mismatch"] += 1
                    status = "CHARGE_MISMATCH"
                    detail = f"built {q:+.3f}, requested {want_q:+d}"
                else:
                    tally["ok"] += 1
                    status, detail = "OK", ""
            rows.append({"protomer_name": name, "parent": exp["parent_name"],
                         "expected_charge": want_q, "status": status,
                         "detail": detail, "db2_file": h.name,
                         "db2_name": hdr.get("db2_name", ""),
                         "amsol_charge": hdr.get("total_charge", "")})

    # db2 files present that nothing asked for
    wanted = {r["protomer_name"] for r in expected}
    for f in all_db2:
        stem = f.name[:-len(".db2.gz")]
        if stem not in wanted and stem.rsplit(".", 1)[0] not in wanted:
            tally["unexpected"] += 1
            rows.append({"protomer_name": "", "parent": "", "expected_charge": "",
                         "status": "UNEXPECTED", "detail": "no matching request",
                         "db2_file": f.name, "db2_name": "", "amsol_charge": ""})

    with open(args.report, "w", newline="") as fh:
        w = csv.DictWriter(fh, lineterminator="\n", fieldnames=["protomer_name", "parent",
                                           "expected_charge", "status", "detail",
                                           "db2_file", "db2_name", "amsol_charge"])
        w.writeheader()
        w.writerows(rows)

    # db2 truncates the M-line name to the last 16 characters
    trunc = Counter(r["db2_name"] for r in rows if r["db2_name"])
    collisions = {k: v for k, v in trunc.items() if v > 1}

    print("\n" + "=" * 56)
    print("DB2 BUILD VERIFICATION")
    print("=" * 56)
    print(f"  requested protomers          {len(expected)}")
    for k in ["ok", "build_failed", "missing", "empty", "charge_mismatch", "renamed",
              "name_mismatch", "unreadable", "unexpected"]:
        print(f"  {k:<28} {tally[k]}")
    if collisions:
        print(f"\n  !! db2 M-line names collide after 16-char truncation:")
        for k, v in collisions.items():
            print(f"     {k!r} x{v}")
    if fail_reasons:
        print("\n  build_ligand failure reasons:")
        for reason, n in Counter(fail_reasons.values()).items():
            print(f"    {reason:<20} {n}")
    if tally["missing"]:
        print(f"\n  !! {tally['missing']} molecule(s) have no db2 and no faillist entry: "
              f"a silent loss (killed task? wrong --outputpath?). See {args.report}.")
    print(f"\n  -> {args.report}\n")

    return 1 if (tally["missing"] or tally["charge_mismatch"] or tally["renamed"]
                 or tally["name_mismatch"]) else 0


if __name__ == "__main__":
    sys.exit(main())
