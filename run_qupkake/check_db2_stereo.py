#!/usr/bin/env python3
"""
check_db2_stereo.py -- will db2_converter keep every name exactly as given?

Run IN THE db2_converter ENVIRONMENT, on library.smi, before build_ligand.

build_ligand --checkstereo calls rdk_enumerate_smi on each SMILES and, when it
returns more than one isomer, renames them NAME.0, NAME.1 ... An 18-character
name is then stored as its last 16 characters, which removes the leading
sentinel and part of the base id (NAMING_CONTRACT.md, invariant 1).

assign_names.py already closed the library under the same enumeration. This
repeats the check with db2_converter's OWN function, under db2_converter's own
RDKit, because stereo perception changes between RDKit releases -- the laptop
env that ran assign_names.py and the cluster env that runs build_ligand need
not agree. It imports the function rather than copying it for that reason.

Exit 0: every SMILES enumerates to exactly itself, every name is 16 characters.
Exit 1: lists the lines db2_converter would rename. Do not build.

    python check_db2_stereo.py library/library.smi
    python check_db2_stereo.py library/library.smi --db2c-src ~/src/db2_converter
"""
from __future__ import annotations

import argparse
import contextlib
import os
import re
import sys
from pathlib import Path

CONTRACT = re.compile(r"^C[0-9A-Za-z]{14}C$")


@contextlib.contextmanager
def quiet_stderr():
    """Silence C++ stderr (RDKit prints a full stack trace for every invariant
    violation, which buried the summary on the 22k-name drug-like library)."""
    sys.stderr.flush()
    saved = os.dup(2)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved, 2)
        os.close(devnull)
        os.close(saved)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("smi")
    ap.add_argument("--db2c-src", default=os.environ.get("DB2C_SRC") or None,
                    help="db2_converter checkout, if the package is not installed "
                         "(default: $DB2C_SRC, set by common/with_env.sh from config/hpc.env)")
    args = ap.parse_args(argv)

    if args.db2c_src:
        sys.path.insert(0, args.db2c_src)
    try:
        from db2_converter.utils.rdkit_gen import rdk_enumerate_smi
    except ImportError as e:
        raise SystemExit(f"cannot import db2_converter ({e}). Run this in the "
                         f"db2_converter env, or pass --db2c-src.")
    from rdkit import Chem, rdBase

    split, unstable, badname, crashed, n = [], [], [], [], 0
    for line in open(args.smi):
        parts = line.split()
        if not parts:
            continue
        smi, name = parts[0], parts[1]
        n += 1
        if not CONTRACT.match(name):
            badname.append(name)
        try:
            with quiet_stderr():
                isomers = rdk_enumerate_smi(smi)
        except Exception as e:  # noqa: BLE001
            # RDKit 2024.03.5 EnumerateStereoisomers can raise "Range Error idx2"
            # on a double bond with no controlling atom on one end. The build
            # catches it too (pipeline.write_enumerated_smifile), logs "cannot
            # be parsed" and leaves the molecule out: no db2, nothing else lost.
            crashed.append((name, smi, str(e).strip().splitlines()[0][:60]))
            continue
        # db2_converter's test is len(outsmis) == 1, and it does not dedupe
        if len(isomers) != 1:
            split.append((name, len(isomers), smi))
        elif Chem.CanonSmiles(isomers[0]) != Chem.CanonSmiles(smi):
            # One isomer, so the name survives (db2_converter keeps the line
            # as given) -- but enumeration moved the stereo labels. Seen on
            # adamantane-like cages, where two stereo indices of one protomer
            # each enumerate to the other: probably one molecule, two rows.
            unstable.append((name, smi, isomers[0]))

    print(f"{n:,} SMILES checked with db2_converter.rdk_enumerate_smi (RDKit {rdBase.rdkitVersion})")
    if badname:
        print(f"  [FAIL] {len(badname):,} name(s) not in contract form, e.g. {badname[:3]}")
    # every split name, not just the ten printed (llb pipeline excludes them all)
    split_path = Path(args.smi).with_name(Path(args.smi).stem + ".stereo_split.tsv")
    with open(split_path, "w") as fh:
        fh.write("name\tn_isomers\tsmiles\n")
        fh.writelines(f"{a}\t{k}\t{b}\n" for a, k, b in split)
    if split:
        print(f"  [FAIL] {len(split):,} SMILES would be renamed NAME.<i> by --checkstereo "
              f"(all in {split_path}):")
        for name, k, smi in split[:10]:
            print(f"           {name}  {k} isomer(s)  {smi}")
        print("         Re-run `llb names` with an RDKit of THIS version (PREP_ENV in "
              "the config), or investigate the perception difference, before building.")
    if unstable:
        print(f"  [warn] {len(unstable):,} SMILES enumerate to a different single SMILES. "
              f"Names are safe; the structure may be duplicated within its base id "
              f"(a wasted dock, not a split leak):")
        for name, a, b in unstable[:5]:
            print(f"           {name}  {a}\n                             -> {b}")
    out = Path(args.smi).with_name(Path(args.smi).stem + ".stereo_crashed.tsv")
    with open(out, "w") as fh:                   # written even when empty: a clean run
        fh.write("name\tsmiles\terror\n")
        fh.writelines(f"{a}\t{b}\t{c}\n" for a, b, c in crashed)
    if crashed:
        print(f"  [warn] {len(crashed):,} SMILES crash RDKit's stereo enumeration; "
              f"build_ligand will skip them (logged 'cannot be parsed', no db2). "
              f"All in {out}; first 10:")
        for name, smi, err in crashed[:10]:
            print(f"           {name}  {smi}  ({err})")
    if badname or split:
        return 1
    print("  [ok  ] every SMILES is stereo-closed; --checkstereo will leave all names alone")
    return 0


if __name__ == "__main__":
    sys.exit(main())
