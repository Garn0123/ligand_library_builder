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
    python verify_db2_build.py protomers/protomers.csv --outputpath example/
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("protomers_csv",
                    help="protomers.csv (qupkake_protomers.py) or library.tsv (assign_names.py)")
    ap.add_argument("--outputpath", required=True,
                    help="db2_converter --outputpath directory")
    ap.add_argument("--faillist", default=None,
                    help="optional .faillist written by build_ligand")
    ap.add_argument("--charge-tol", type=float, default=0.05,
                    help="tolerance when comparing AMSOL charge to formal charge")
    ap.add_argument("-o", "--report", default="db2_verification.csv")
    args = ap.parse_args(argv)

    outdir = Path(args.outputpath)
    # protomers.csv (qupkake_protomers.py) or library.tsv (assign_names.py)
    with open(args.protomers_csv) as fh:
        tsv = args.protomers_csv.endswith(".tsv")
        expected = [{"protomer_name": r.get("protomer_name") or r["name"],
                     "parent_name": (r.get("parent_name") or r.get("input_name")
                                     or r.get("zinc_id", "")),
                     "net_charge": r["net_charge"]}
                    for r in csv.DictReader(fh, delimiter="\t" if tsv else ",")]

    fail_reasons: dict[str, str] = {}
    if args.faillist and Path(args.faillist).exists():
        for line in Path(args.faillist).read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 4:
                fail_reasons[parts[1]] = parts[3]

    rows, tally = [], Counter()
    for exp in expected:
        name = exp["protomer_name"]
        want_q = int(exp["net_charge"])

        # --checkstereo splits unspecified stereocentres into NAME.0, NAME.1 ...
        exact = sorted(outdir.glob(f"{name}.db2.gz"))
        hits = exact or sorted(outdir.glob(f"{name}.[0-9]*.db2.gz"))
        # A contract name (NAMING_CONTRACT.md) must reach the db2 unchanged: a
        # .N suffix makes it 18 characters and mol2db2 keeps only the last 16.
        contract = CONTRACT.match(name) is not None

        if not hits:
            tally["missing"] += 1
            rows.append({"protomer_name": name, "parent": exp["parent_name"],
                         "expected_charge": want_q, "status": "MISSING",
                         "detail": fail_reasons.get(name, "no db2 produced"),
                         "db2_file": "", "db2_name": "", "amsol_charge": ""})
            continue

        for h in hits:
            if h.stat().st_size == 0:
                tally["empty"] += 1
                status, detail, hdr = "EMPTY", "zero-length db2", {}
            else:
                hdr = read_db2_header(h)
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
    for f in outdir.glob("*.db2.gz"):
        stem = f.name[:-len(".db2.gz")]
        if stem not in wanted and stem.rsplit(".", 1)[0] not in wanted:
            tally["unexpected"] += 1
            rows.append({"protomer_name": "", "parent": "", "expected_charge": "",
                         "status": "UNEXPECTED", "detail": "no matching request",
                         "db2_file": f.name, "db2_name": "", "amsol_charge": ""})

    with open(args.report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["protomer_name", "parent",
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
    for k in ["ok", "missing", "empty", "charge_mismatch", "renamed",
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
    print(f"\n  -> {args.report}\n")

    return 1 if (tally["missing"] or tally["charge_mismatch"] or tally["renamed"]
                 or tally["name_mismatch"]) else 0


if __name__ == "__main__":
    sys.exit(main())
