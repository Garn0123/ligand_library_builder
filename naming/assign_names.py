#!/usr/bin/env python3
"""
assign_names.py -- MolSanitizer's states and stereoisomers -> contract names.

The naming step of NAMING_CONTRACT.md (DRAP repo), for libraries whose chemistry
MolSanitizer decides. Runs ONCE over every MolSanitizer output of a library,
because uniqueness is a property of the whole library.

    llb names --inputs samples/H*.smi --ms ms_pass1/*_clean.smi -o library

    C  <12: base id>  <state>  <stereo>  C          exactly 16 characters

MolSanitizer names what it emits after its input name:

    ZINCaa000013yz56_2.3     state 2 (a tautomer x protonation state), stereoisomer 3
    ZINCaa000006kl12_1       several states, no stereo enumerated
    ZINCaa000011uv12.2       one state, several stereoisomers
    ZINCaa000009qr78         one state, one isomer

so its hierarchy is ours. Here each index is renumbered from 0 per parent (per
state for stereo), in MolSanitizer's order, and written as one base62 character.
The state character covers tautomers AND protomers: MolSanitizer numbers them
together (`--tautomers --protonation`), where LLB's old `protomer` character
counted protonation states only. library.tsv keeps the column name
protomer_index so docking/ and DRAP read it unchanged.

  base id   llb_ids.base_id_of(parent id): ZINC id [4:16] for a ZINC parent,
            'Z' + 11 base62 of sha256(name) for any other (see llb_ids.py).

Why rename before 3D: a db2 keeps only the last 16 characters of a name
(MolSanitizer's writer, run 2026-10-20: Ch300000BlLjl00C.1 -> 300000BlLjl00C.1),
and a ZINC22 id is already 16. So MolSanitizer runs twice: pass 1 on SMILES only
(tautomers, protonation, stereo) -> this step -> pass 2 on library.smi with
`--no-neutralize --no-stereoisomers --gen3d -f db2`, which keeps every name.

Parsing needs --inputs (what went into pass 1): an input name may itself
contain '_' or '.', so the suffix is what follows a KNOWN input name, never a
guess from the shape of the string. Inputs MolSanitizer did not return
(filtered, failed) are counted in stats.json and listed in dropped.tsv.

Enforced here, as in LLB's old names step (tag llb-full-2026-10):
* stereo closed: re-enumerating any output SMILES (RDKit, tryEmbedding=False)
  gives exactly one isomer, so a name means one structure in pass 2.
* one structure, one name: a SMILES reached from two parents keeps the first
  (sorted by parent id), and duplicates within a parent are dropped.
* one fragment per row (run MolSanitizer with --removesalts if inputs have salts).
* > 62 states for a parent, or > 62 isomers for a state: dropped and recorded.
* N-H+ invertomers: MolSanitizer enumerates both N configurations as stereo-
  isomers. library.tsv's invertomer_of names the row each one is an N-inversion
  partner of, so `llb collapse-invertomers` can keep the better score.
* invariants 1-3 (16 characters, unique names and db2 tails, one parent per
  base id) are asserted before anything is written.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

from llb_ids import BASE62, base_id_of, is_zinc, parent_id
from rdkit import Chem, RDLogger
from rdkit.Chem.EnumerateStereoisomers import (EnumerateStereoisomers,
                                                StereoEnumerationOptions)

RDLogger.DisableLog("rdApp.*")

NAME_LEN = 16
STEREO_OPTS = StereoEnumerationOptions(tryEmbedding=False)
COLUMNS = ["name", "base_id", "parent_id", "zinc_id", "input_name", "protomer_index",
           "stereo_index", "net_charge", "heavy_atoms", "ph_values", "population_estimate",
           "invertomer_of", "note", "smiles", "site_source", "ms_name"]


def read_smi(path: Path):
    """(smiles, name) per non-empty line of a 'SMILES NAME' file, plain or gzipped."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                yield parts[0], parts[1]


def split_ms_name(ms_name: str, inputs: set[str]) -> tuple[str, int, int] | None:
    """MolSanitizer's 'INPUT[_state][.stereo]' -> (input, state, stereo), 1-based; 0 = absent.

    Tries the longest known input name that the string starts with, so an input
    name containing '_' or '.' is never cut. None if nothing known matches.
    """
    cands = [ms_name] + [ms_name[:i] for i in range(len(ms_name) - 1, 0, -1)
                         if ms_name[i] in "_."]
    for head in cands:                          # longest first
        if head not in inputs:
            continue
        rest = ms_name[len(head):]
        state = stereo = 0
        if rest.startswith("_"):
            num, _, rest = rest[1:].partition(".")
            if not num.isdigit():
                continue
            state, rest = int(num), ("." + rest if rest else "")
        if rest.startswith("."):
            if not rest[1:].isdigit():
                continue
            stereo, rest = int(rest[1:]), ""
        if rest == "":
            return head, state, stereo
    return None


def stereo_count(mol: Chem.Mol) -> int:
    return len({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(mol, options=STEREO_OPTS)})


def invertomer_key(smi: str) -> str:
    """SMILES with the configuration of every N-H+ centre erased.

    Isomers sharing this key differ ONLY at protonated-amine centres, which
    interconvert in solution (deprotonate, invert, reprotonate). Quaternary N+
    (no H) cannot invert, so its isomers stay separate compounds.
    """
    m = Chem.MolFromSmiles(smi)
    for a in m.GetAtoms():
        if a.GetSymbol() == "N" and a.GetFormalCharge() == 1 and a.GetTotalNumHs() > 0:
            a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return Chem.MolToSmiles(m)


def contract_name(base: str, p: int, s: int) -> str:
    return f"C{base}{BASE62[p]}{BASE62[s]}C"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="+", type=Path, required=True,
                    help="the 'SMILES NAME' files given to MolSanitizer pass 1 (samples/H*.smi)")
    ap.add_argument("--ms", nargs="+", type=Path, required=True,
                    help="MolSanitizer pass 1 output (*_clean.smi)")
    ap.add_argument("-o", "--outdir", type=Path, required=True)
    ap.add_argument("--ph", default="",
                    help="the pH settings pass 1 used, recorded in ph_values (e.g. '7.4+-1')")
    ap.add_argument("--force", action="store_true", help="overwrite an existing library")
    args = ap.parse_args(argv)

    out = args.outdir
    if (out / "library.tsv").exists() and not args.force:
        raise SystemExit(f"{out}/library.tsv exists; --force to overwrite")

    # ---- parents: input name -> parent id -------------------------------------
    dropped: list[tuple[str, str, str, str]] = []     # input_name, ms_name, reason, detail
    pid_of: dict[str, str] = {}
    seen: set[str] = set()
    base_owner: dict[str, str] = {}
    for f in args.inputs:
        for _, name in read_smi(f):
            if name in seen:
                continue
            seen.add(name)
            pid, kind = parent_id(name)
            if pid is None:
                dropped.append((name, "", kind, str(f)))
                continue
            base = base_id_of(pid)
            if base_owner.setdefault(base, pid) != pid:
                dropped.append((name, "", "base_id_collision", base_owner[base]))
                continue
            pid_of[name] = pid
    inputs = set(pid_of)

    # ---- MolSanitizer rows, grouped parent -> state -> [(stereo, smiles, ms_name)]
    groups: dict[str, dict[int, list]] = defaultdict(lambda: defaultdict(list))
    n_ms = 0
    for f in args.ms:
        for smi, ms_name in read_smi(f):
            n_ms += 1
            hit = split_ms_name(ms_name, inputs)
            if hit is None:
                dropped.append(("", ms_name, "unknown_input_name", str(f)))
                continue
            inp, state, stereo = hit
            groups[inp][state].append((stereo, smi, ms_name))

    # ---- name -----------------------------------------------------------------
    library: list[dict] = []
    owner: dict[str, str] = {}                  # canonical SMILES -> name
    counts = Counter()
    for inp in sorted(groups, key=lambda n: pid_of[n]):
        pid = pid_of[inp]
        base = base_id_of(pid)
        p_new = 0
        for state in sorted(groups[inp]):
            kept = []
            for stereo, smi, ms_name in sorted(groups[inp][state]):
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    dropped.append((inp, ms_name, "unparsable_smiles", smi)); continue
                if len(Chem.GetMolFrags(mol)) > 1:
                    dropped.append((inp, ms_name, "multiple_fragments", smi)); continue
                if stereo_count(mol) != 1:
                    dropped.append((inp, ms_name, "stereo_not_closed", smi)); continue
                can = Chem.MolToSmiles(mol)
                if can in owner:
                    dropped.append((inp, ms_name, "duplicate_structure",
                                    f"{can} already named {owner[can]}")); continue
                if any(can == k[0] for k in kept):
                    dropped.append((inp, ms_name, "duplicate_structure", can)); continue
                kept.append((can, mol, ms_name))
            if not kept:
                continue
            if p_new >= len(BASE62):
                dropped.extend((inp, k[2], "too_many_states", "> 62") for k in kept)
                continue
            if len(kept) > len(BASE62):
                dropped.extend((inp, k[2], "too_many_stereoisomers", "> 62") for k in kept[62:])
                kept = kept[:62]
            first_of: dict[str, str] = {}
            for s_idx, (can, mol, ms_name) in enumerate(kept):
                name = contract_name(base, p_new, s_idx)
                inv = first_of.setdefault(invertomer_key(can), name)
                owner[can] = name
                library.append({
                    "name": name, "base_id": base, "parent_id": pid,
                    "zinc_id": pid if is_zinc(pid) else "", "input_name": inp,
                    "protomer_index": p_new, "stereo_index": s_idx,
                    "net_charge": Chem.GetFormalCharge(mol),
                    "heavy_atoms": mol.GetNumHeavyAtoms(),
                    "ph_values": args.ph, "population_estimate": "",
                    "invertomer_of": inv if inv != name else "", "note": "",
                    "smiles": can, "site_source": "molsanitizer", "ms_name": ms_name})
            p_new += 1
        counts["parents_named"] += p_new > 0

    # ---- invariants 1-3, asserted before anything final is written -------------
    names = [r["name"] for r in library]
    problems = []
    if any(len(n) != NAME_LEN or not (n[0] == n[-1] == "C") for n in names):
        problems.append("a name is not C + 14 base62 + C")
    for label, keys in (("names", names), ("db2 tails name[-16:]", [n[-16:] for n in names])):
        dup = [k for k, c in Counter(keys).items() if c > 1]
        if dup:
            problems.append(f"{len(dup)} duplicate {label}, e.g. {dup[:3]}")
    owners = defaultdict(set)
    for r in library:
        owners[r["base_id"]].add(r["parent_id"])
    if any(len(v) > 1 for v in owners.values()):
        problems.append("a base id maps to more than one parent")
    if problems:
        raise SystemExit("contract violated -- nothing written:\n  " + "\n  ".join(problems))

    # ---- write -----------------------------------------------------------------
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "library.smi", "w") as fh:
        fh.writelines(f"{r['smiles']} {r['name']}\n" for r in library)
    with open(out / "library.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, lineterminator="\n", delimiter="\t", fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(library)
    missing = sorted(inputs - set(groups), key=lambda n: pid_of[n])
    with open(out / "dropped.tsv", "w") as fh:
        fh.write("input_name\tms_name\treason\tdetail\n")
        fh.writelines("\t".join(d) + "\n" for d in dropped)
        fh.writelines(f"{n}\t\tnot_returned_by_molsanitizer\t\n" for n in missing)
    reasons = Counter(d[2] for d in dropped)
    reasons["not_returned_by_molsanitizer"] += len(missing)
    stats = {
        "inputs": len(seen),
        "molsanitizer_rows": n_ms,
        "parents_named": counts["parents_named"],
        "names": len(library),
        "names_per_parent": round(len(library) / max(counts["parents_named"], 1), 3),
        "n_h_invertomers": sum(1 for r in library if r["invertomer_of"]),
        "dropped": {k: v for k, v in reasons.items() if v},
        "ph_values": args.ph,
        "inputs_files": [str(p) for p in args.inputs],
        "ms_files": [str(p) for p in args.ms],
        "library_smi_sha256": hashlib.sha256((out / "library.smi").read_bytes()).hexdigest(),
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=2))

    print(f"{stats['names']:,} names for {stats['parents_named']:,} parents "
          f"({stats['names_per_parent']} per parent) -> {out}/library.smi")
    print(f"  of which {stats['n_h_invertomers']:,} are N-H+ invertomers of another row "
          f"(library.tsv: invertomer_of)")
    for k, v in sorted(stats["dropped"].items()):
        print(f"  dropped {k:<30} {v:,}")
    print(f"\nnext:  msani -i {out}/library.smi --no-neutralize --no-stereoisomers --gen3d -f db2")
    return 0


if __name__ == "__main__":
    sys.exit(main())
