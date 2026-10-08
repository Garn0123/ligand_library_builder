#!/usr/bin/env python3
"""
mol2_library.py -- the final multi-mol2 library: one conformer per molecule,
DOCK-format headers, cut into shards of N molecules.

    llb mol2 db2_run/out --library library/library.tsv -o library/library.mol2
    llb mol2-library library/library.mol2 --library library/library.tsv -o mol2_lib --shard-size 1000

Input is mol2 from `llb mol2` (db2_to_mol2.py): one record per molecule, the
lowest-energy conformer read back out of its db2, with db2tool's header block.
Several input files (or directories of .mol2) are fine; if a name appears more
than once, the record with the lowest HDB_Internal_energy1 is kept.

Each record is written with a header block in DOCK6's own format, the one its
library_file.cpp writes and Read_Mol2_retain reads back:

    ##########                                Name:    CZbLeIWZR3d9900C
     ^ "########## " + label right-aligned in 36 + value right-aligned in 20

DOCK labels, computed the way DOCK computes them, from the mol2 itself, so a
DOCK run that reads them back (dbfilter, Read_Mol2_retain) sees its own numbers:
    Name               the 16-character contract name (mol2 title)
    Molecular_Weight   sum of atomic masses over the mol2's atoms, H included
    Formal_Charge      sum of the partial charges, rounded to 0.1 (as DOCK prints)
    HBond_Acceptors    atoms whose SYBYL type starts O, N, S, F or Cl
                       (amber_typer.cpp flag_acceptor -- every such atom counts)
    HBond_Donors       H atoms bonded to one of those (flag_donator)
    Heavy_Atoms        non-H atoms
DOCK_Rotatable_Bonds is NOT written: DOCK's definition comes from its flex.defn
typing, and a different number under DOCK's label would be read back as DOCK's.

db2tool's own HDB_ lines (conformer set, internal energies, broken flag) are
kept. Then LLB_ lines from library.tsv and RDKit on the protomer SMILES, under
labels no DOCK reader parses:
    LLB_SMILES  LLB_Parent_ID  LLB_Input_Name  LLB_Base_ID  LLB_Protomer
    LLB_Stereo  LLB_Net_Charge  LLB_pH  LLB_Population  LLB_Site_Source
    LLB_Invertomer_Of  LLB_HBA_Lipinski  LLB_HBD_Lipinski  LLB_Rotatable_Bonds
    LLB_TPSA  LLB_cLogP
Every value is one whitespace-free token ("-" when empty): DOCK splits these
lines on whitespace and expects exactly three fields.

Writes OUT/<prefix>_NNNNN.mol2 (N molecules each, sorted by name, so a parent's
protomers and stereoisomers sit together), OUT/index.tsv (every header value
plus its shard), and OUT/MANIFEST.tsv (shard, count, first/last name, sha256).
Exit 1 if a library name is missing or a mol2's charge does not match its
protomer -- the same checks `llb mol2` makes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

# db2core::DELIMITER / STRING_WIDTH / FLOAT_WIDTH (dock6_claude db2_io.h)
DELIMITER, STRING_WIDTH, FLOAT_WIDTH = "########## ", 36, 20
HEADER = re.compile(r"^##########\s+([A-Za-z0-9_]+):\s*(.*?)\s*$")
HB_ELEMENTS = ("O", "N", "S", "F", "Cl", "CL")
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "P": 30.974,
        "S": 32.06, "Cl": 35.45, "Br": 79.904, "I": 126.904, "B": 10.81, "Si": 28.085,
        "Se": 78.971, "Na": 22.990, "K": 39.098, "Li": 6.94, "Mg": 24.305, "Ca": 40.078,
        "Zn": 65.38, "Fe": 55.845}
DOCK_ORDER = ["Name", "Molecular_Weight", "Formal_Charge", "HBond_Acceptors",
              "HBond_Donors", "Heavy_Atoms"]
LLB_ORDER = ["LLB_SMILES", "LLB_Parent_ID", "LLB_Input_Name", "LLB_Base_ID",
             "LLB_Protomer", "LLB_Stereo", "LLB_Net_Charge", "LLB_pH", "LLB_Population",
             "LLB_Site_Source", "LLB_Invertomer_Of", "LLB_HBA_Lipinski",
             "LLB_HBD_Lipinski", "LLB_Rotatable_Bonds", "LLB_TPSA", "LLB_cLogP"]


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------
def read_records(path: Path) -> list[dict]:
    """Records of a multi-mol2: header dict (## lines) + body from @<TRIPOS>MOLECULE."""
    records, header, body = [], {}, None
    with open(path) as fh:
        for line in fh:
            m = HEADER.match(line)
            if m:
                if body is not None:                 # a header starts a new record
                    records.append({"header": header, "body": body})
                    header, body = {}, None
                header[m.group(1)] = m.group(2)
                continue
            if line.startswith("@<TRIPOS>MOLECULE"):
                if body is not None:
                    records.append({"header": header, "body": body})
                    header = {}
                body = [line]
            elif body is not None:
                body.append(line)
    if body is not None:
        records.append({"header": header, "body": body})
    for r in records:
        r["name"] = r["body"][1].strip() if len(r["body"]) > 1 else ""
        r["e1"] = float(r["header"].get("HDB_Internal_energy1", "inf"))
    return records


def sections(body: list[str]) -> dict[str, list[str]]:
    out, cur = defaultdict(list), None
    for line in body:
        if line.startswith("@<TRIPOS>"):
            cur = line.strip()[9:]
            continue
        if cur and line.strip():
            out[cur].append(line)
    return out


def element(sybyl: str) -> str:
    """'C.ar' -> C, 'Cl' -> Cl, 'N.pl3' -> N."""
    e = sybyl.split(".")[0]
    return e if e in MASS or len(e) == 1 else e[0] + e[1:].lower()


def dock_descriptors(body: list[str]) -> dict:
    """DOCK's header numbers, computed from the mol2 the way DOCK does."""
    s = sections(body)
    types, charge, mw = {}, 0.0, 0.0
    for line in s["ATOM"]:
        f = line.split()
        types[f[0]] = f[5]
        charge += float(f[-1])
        mw += MASS.get(element(f[5]), 0.0)
    hb = {i for i, t in types.items() if t.startswith(HB_ELEMENTS)}
    donors = set()
    for line in s["BOND"]:
        _, a, b = line.split()[:3]
        for h, x in ((a, b), (b, a)):
            if types.get(h, "").startswith("H") and x in hb:     # DOCK: type[0] == 'H'
                donors.add(h)
    heavy = sum(1 for t in types.values() if element(t) != "H")
    return {"Molecular_Weight": mw, "Formal_Charge": charge, "HBond_Acceptors": len(hb),
            "HBond_Donors": len(donors), "Heavy_Atoms": heavy,
            "_unknown_elements": sorted({element(t) for t in types.values()} - set(MASS))}


def llb_fields(row: dict | None) -> dict:
    if row is None:
        return {}
    mol = Chem.MolFromSmiles(row["smiles"])
    d = {"LLB_SMILES": row["smiles"], "LLB_Parent_ID": row["parent_id"],
         "LLB_Input_Name": row.get("input_name", ""), "LLB_Base_ID": row["base_id"],
         "LLB_Protomer": row["protomer_index"], "LLB_Stereo": row["stereo_index"],
         "LLB_Net_Charge": row["net_charge"], "LLB_pH": row.get("ph_values", ""),
         "LLB_Population": row.get("population_estimate", ""),
         "LLB_Site_Source": row.get("site_source", ""),
         "LLB_Invertomer_Of": row.get("invertomer_of", "")}
    if mol is not None:
        d.update({"LLB_HBA_Lipinski": Lipinski.NumHAcceptors(mol),
                  "LLB_HBD_Lipinski": Lipinski.NumHDonors(mol),
                  "LLB_Rotatable_Bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
                  "LLB_TPSA": f"{rdMolDescriptors.CalcTPSA(mol):.2f}",
                  "LLB_cLogP": f"{Crippen.MolLogP(mol):.3f}"})
    return d


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------
def token(v) -> str:
    s = "-" if v is None or v == "" else str(v)
    return re.sub(r"\s+", "_", s)


def fmt(label: str, value) -> str:
    if label == "Molecular_Weight":
        value = f"{value:.6g}"                       # ostream default, as DOCK prints
    elif label == "Formal_Charge":
        value = f"{round(value * 10) / 10:g}"        # round(x*10)/10, as DOCK prints
    v = token(value)
    # DOCK's setw(20) lets a longer value run into the label ("Name:CZ..."); keep
    # one space so every line splits into exactly three whitespace fields
    v = f"{v:>{FLOAT_WIDTH}}" if len(v) < FLOAT_WIDTH else " " + v
    return f"{DELIMITER}{label + ':':>{STRING_WIDTH}}{v}\n"


def header_block(fields: dict, hdb: dict) -> str:
    lines = [fmt(k, fields[k]) for k in DOCK_ORDER]
    lines += [fmt(k, v) for k, v in hdb.items() if k.startswith("HDB_")]
    lines += [fmt(k, fields[k]) for k in LLB_ORDER if k in fields]
    return "\n" + "".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="mol2 files (from llb mol2) or directories of them")
    ap.add_argument("--library", required=True, help="library.tsv from llb names")
    ap.add_argument("-o", "--outdir", type=Path, required=True)
    ap.add_argument("--shard-size", type=int, default=1000, help="molecules per mol2 file")
    ap.add_argument("--prefix", default="library")
    ap.add_argument("--charge-tol", type=float, default=0.05)
    args = ap.parse_args(argv)
    if args.shard_size < 1:
        raise SystemExit("--shard-size must be >= 1")
    if args.outdir.exists() and any(args.outdir.glob("*.mol2")):
        raise SystemExit(f"{args.outdir} already holds mol2 shards; use a new -o")

    files = []
    for s in map(Path, args.inputs):
        files += sorted(s.rglob("*.mol2")) if s.is_dir() else [s]
    by_name: dict[str, list[dict]] = defaultdict(list)
    for f in files:
        for r in read_records(f):
            by_name[r["name"]].append(r)
    if not by_name:
        raise SystemExit("no mol2 records found")
    chosen = {n: min(rs, key=lambda r: r["e1"]) for n, rs in by_name.items()}

    with open(args.library) as fh:
        library = {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}

    args.outdir.mkdir(parents=True, exist_ok=True)
    names = sorted(chosen)
    index_rows, manifest, problems = [], [], Counter()
    for s, start in enumerate(range(0, len(names), args.shard_size)):
        part = names[start:start + args.shard_size]
        path = args.outdir / f"{args.prefix}_{s:05d}.mol2"
        with open(path, "w") as fh:
            for name in part:
                r = chosen[name]
                row = library.get(name)
                d = dock_descriptors(r["body"])
                fields = {"Name": name, **{k: d[k] for k in DOCK_ORDER[1:]}, **llb_fields(row)}
                status = "OK"
                if row is None:
                    status = "NOT_IN_LIBRARY"
                elif round(d["Formal_Charge"]) != int(row["net_charge"]):
                    status = "CHARGE_MISMATCH"       # a different integer charge
                if d["_unknown_elements"]:
                    status = f"UNKNOWN_ELEMENT:{','.join(d['_unknown_elements'])}"
                problems[status] += status != "OK"
                fh.write(header_block(fields, r["header"]))
                fh.write("".join(r["body"]).rstrip("\n") + "\n")
                index_rows.append({"name": name, "shard": path.name, "status": status,
                                   "copies_in_input": len(by_name[name]),
                                   **{k: fields.get(k, "") for k in DOCK_ORDER[1:]},
                                   **{k: r["header"].get(k, "") for k in
                                      ("HDB_Conf_Number", "HDB_Internal_energy1", "HDB_Set_Broken")},
                                   **{k: fields.get(k, "") for k in LLB_ORDER}})
        manifest.append((path.name, len(part), part[0], part[-1],
                         hashlib.sha256(path.read_bytes()).hexdigest()))

    missing = sorted(set(library) - set(chosen))
    cols = list(index_rows[0])
    with open(args.outdir / "index.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, lineterminator="\n", delimiter="\t", fieldnames=cols)
        w.writeheader()
        for r in index_rows:
            r = dict(r)
            r["Molecular_Weight"] = f"{r['Molecular_Weight']:.3f}"
            r["Formal_Charge"] = f"{r['Formal_Charge']:+.3f}"
            w.writerow(r)
        for n in missing:
            w.writerow({"name": n, "status": "MISSING"})
    with open(args.outdir / "MANIFEST.tsv", "w") as fh:
        fh.write("file\tmolecules\tfirst\tlast\tsha256\n")
        fh.writelines("\t".join(map(str, m)) + "\n" for m in manifest)

    dup = sum(1 for rs in by_name.values() if len(rs) > 1)
    print(f"{len(files):,} mol2 file(s), {len(chosen):,} molecule(s) -> {len(manifest):,} "
          f"shard(s) of <= {args.shard_size:,} in {args.outdir}")
    if dup:
        print(f"  {dup:,} name(s) appeared more than once; lowest HDB_Internal_energy1 kept")
    for k, v in sorted(problems.items()):
        if v:
            print(f"  ! {k}: {v:,}")
    if missing:
        print(f"  ! {len(missing):,} library name(s) not in any input mol2 (index.tsv: MISSING)")
    print(f"  index -> {args.outdir}/index.tsv, shards -> {args.outdir}/MANIFEST.tsv")
    return 1 if (missing or sum(problems.values())) else 0


if __name__ == "__main__":
    sys.exit(main())
