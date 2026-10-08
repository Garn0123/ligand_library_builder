#!/usr/bin/env python3
"""
assign_names.py -- every shard's protomers -> one library under contract names.

The naming step of NAMING_CONTRACT.md (DRAP repo). Runs ONCE over the union of
all QupKake shards, because uniqueness is a property of the whole library and no
shard can see the others.

    C  <12: base id>  <protomer>  <stereo>  C          exactly 16 characters

  base id   llb_ids.base_id_of(parent id). For a ZINC parent, ZINC id [4:16]:
            ZINC20 (ZINC000012345678 -> 000012345678) and ZINC22
            (ZINCh10000007sNS -> h10000007sNS) both give 12 base62 characters.
            For any other name, 'Z' + 11 base62 of its sha256 -- a namespace no
            ZINC id can enter (see llb_ids.py). library.tsv maps back either way.
  protomer  base62, renumbered 0.. per parent in QupKake population order.
  stereo    base62, 0.. per (parent, protomer), over the stereo CLOSURE below.

Stereo closure -- why this step enumerates stereo at all
--------------------------------------------------------
db2_converter's --checkstereo re-enumerates every input SMILES and, if it finds
more than one isomer, renames them NAME.0, NAME.1 (pipeline.py,
write_enumerated_smifile). An 18-character name is then cut to its LAST 16 by
mol2db2 -- the leading 'C' and first base-id character are lost, and the
contract is broken inside the db2 where nothing upstream can see it.

Protonation creates stereocentres the neutral parent did not have. Measured with
RDKit's default enumeration (exactly what db2_converter calls):

    CCN(C)Cc1ccccc1          1 isomer
    CC[NH+](C)Cc1ccccc1      2 isomers   <- protonated tertiary amine
    C[C@H]1CC[NH+](C)CC1     2 isomers   <- ring cis/trans through N+

So stereo enumerated on the parent (qupkake_protomers.py --stereo enumerate) is
not enough. Here every protomer is expanded with the same options db2_converter
uses, each isomer gets its own stereo index, and the result is checked to be a
fixed point: re-enumerating any output SMILES gives exactly one isomer, so
--checkstereo leaves every name alone. check_db2_stereo.py re-runs that check
with db2_converter's own function, in its own environment, before any build.

N-H+ invertomers are genuine stereoisomers and are kept as rows, but
library.tsv's invertomer_of column names the row each one is an N-inversion
partner of, so the matrix can collapse them (they interconvert in solution).

Also enforced here
------------------
* one structure, one name: a SMILES reached from two parents keeps the first
  (sorted by ZINC id) and the other is dropped -- two rows with different base
  ids for one structure is a train/val leak the identity split cannot group.
* > 62 protomers, or > --max-stereo isomers for one protomer: dropped and
  recorded (db2_converter itself skips > 32 isomers as "0enumerate_many").
* invariants 1-3 are asserted on the result before anything is written as
  final. stage4/validate_names.py (DRAP side) should still be run on the output.

Usage:
    python assign_names.py --shards-dir protomers --parents parents -o library
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from llb_ids import base_id_of, is_zinc
from rdkit import Chem, RDLogger
from rdkit.Chem.EnumerateStereoisomers import (EnumerateStereoisomers,
                                                StereoEnumerationOptions)

RDLogger.DisableLog("rdApp.*")

BASE62 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
CONTRACT = re.compile(r"^C[0-9A-Za-z]{14}C$")
NAME_LEN = 16
# db2_converter/utils/rdkit_gen.py: rdk_enumerate_smi -> StereoEnumerationOptions(
# tryEmbedding=False), everything else default. Must stay identical to it.
DB2C_OPTS = StereoEnumerationOptions(tryEmbedding=False)


def closure(smi: str) -> list[str]:
    """Every stereoisomer db2_converter would split `smi` into (sorted)."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return []
    return sorted({Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m, options=DB2C_OPTS)})


def invertomer_key(smi: str) -> str:
    """SMILES with the configuration of every N-H+ centre erased.

    A protonated amine with three different substituents is a genuine
    stereocentre, so its two configurations are separate rows here. But they
    interconvert in solution (deprotonate, invert, reprotonate), so downstream
    should treat a pair as one species -- e.g. keep the better score. Isomers
    sharing this key differ ONLY at such centres. Quaternary N+ (no H) is left
    alone: it cannot invert, so its isomers are separate compounds.
    """
    m = Chem.MolFromSmiles(smi)
    for a in m.GetAtoms():
        if a.GetSymbol() == "N" and a.GetFormalCharge() == 1 and a.GetTotalNumHs() > 0:
            a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return Chem.MolToSmiles(m)


def contract_name(base: str, p: int, s: int) -> str:
    return f"C{base}{BASE62[p]}{BASE62[s]}C"


def base_id(name: str) -> str:
    return name[1:13]


def read_manifest_done(shards_dir: Path, parents_dir: Path | None):
    """Shard dirs to read, and which expected shards are missing or unfinished."""
    found = sorted(d for d in shards_dir.glob("shard_*") if (d / "DONE").exists())
    missing = []
    if parents_dir is not None:
        with open(parents_dir / "shards.tsv") as fh:
            want = [f"shard_{int(r['index']):05d}" for r in csv.DictReader(fh, delimiter="\t")]
        have = {d.name for d in found}
        missing = [w for w in want if w not in have]
    return found, missing


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shards-dir", required=True,
                    help="OUT_DIR of the QupKake array, or of llb merge (holds shard_NNNNN/)")
    ap.add_argument("--parents", default=None,
                    help="prepare_parents.py --outdir: completeness check + per-bin attrition")
    ap.add_argument("-o", "--outdir", default="library")
    ap.add_argument("--max-stereo", type=int, default=32,
                    help="isomers per protomer; db2_converter skips > 32")
    ap.add_argument("--allow-partial", action="store_true",
                    help="name what is finished even if some shards are not")
    args = ap.parse_args(argv)

    shards_dir, outdir = Path(args.shards_dir), Path(args.outdir)
    parents_dir = Path(args.parents) if args.parents else None
    found, missing = read_manifest_done(shards_dir, parents_dir)
    if missing and not args.allow_partial:
        raise SystemExit(f"{len(missing)} shard(s) not DONE, e.g. {missing[:5]}. Names "
                         f"must be assigned over the whole library; finish them or pass "
                         f"--allow-partial (and never merge two partial namings).")
    if not found:
        raise SystemExit(f"no finished shards under {shards_dir}")
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- gather ------------------------------------------------------------
    rows_by_parent: dict[str, list[dict]] = defaultdict(list)
    qk_failed: list[tuple[str, str]] = []
    for d in found:
        with open(d / "protomers.csv") as fh:
            for r in csv.DictReader(fh):
                rows_by_parent[r["parent_name"]].append(r)
        f = d / "qupkake_failed.tsv"
        if f.exists():
            with open(f) as fh:
                qk_failed += [(r["name"], r["reason"]) for r in csv.DictReader(fh, delimiter="\t")]

    dropped: list[tuple[str, str, str, str]] = []     # parent_id, protomer, reason, detail
    bases: dict[str, str] = {}
    for pid in rows_by_parent:
        try:
            bases[pid] = base_id_of(pid)
        except ValueError:
            raise SystemExit(f"parent {pid!r} is neither a ZINC id nor an LLB<base id>. "
                             f"prepare_parents.py assigns these; was it skipped?")

    # original names (prepare_parents.py --outdir); a ZINC id is its own name
    input_name: dict[str, str] = {}
    parents_rows: list[dict] = []
    if parents_dir is not None:
        with open(parents_dir / "parents.tsv") as fh:
            parents_rows = list(csv.DictReader(fh, delimiter="\t"))
        for r in parents_rows:
            pid = r.get("parent_id") or r["zinc_id"]     # older parents.tsv: zinc_id only
            input_name[pid] = r.get("input_name") or pid

    # ---- closure + naming ----------------------------------------------------
    library: list[dict] = []
    owner: dict[str, str] = {}                        # final SMILES -> name that holds it
    closure_added = 0
    for pid in sorted(rows_by_parent):
        by_state: dict[int, list[dict]] = defaultdict(list)
        for r in rows_by_parent[pid]:
            by_state[int(r["protomer_index"])].append(r)
        p_new = 0
        for p_old in sorted(by_state):
            rows = by_state[p_old]
            isomers: list[str] = []
            for r in sorted(rows, key=lambda r: int(r["stereo_index"])):
                for smi in closure(r["smiles"]):
                    if smi not in isomers:
                        isomers.append(smi)
            closure_added += len(isomers) - len(rows)
            head = rows[0]
            if not isomers:
                dropped.append((pid, str(p_old), "unparsable_protomer", head["smiles"]))
                continue
            if len(isomers) > args.max_stereo:
                dropped.append((pid, str(p_old), "too_many_stereoisomers",
                                f"{len(isomers)} > {args.max_stereo}"))
                continue
            not_closed = [s for s in isomers if len(closure(s)) != 1]
            if not_closed:
                dropped.append((pid, str(p_old), "stereo_not_closed", not_closed[0]))
                continue
            taken = [s for s in isomers if s in owner]
            for s in taken:
                dropped.append((pid, str(p_old), "duplicate_structure",
                                f"{s} already named {owner[s]}"))
            isomers = [s for s in isomers if s not in owner]
            if not isomers:
                continue
            if p_new >= len(BASE62):
                dropped.append((pid, str(p_old), "too_many_protomers", "> 62"))
                continue
            first_of: dict[str, str] = {}
            for s_idx, smi in enumerate(isomers):
                name = contract_name(bases[pid], p_new, s_idx)
                inv = first_of.setdefault(invertomer_key(smi), name)
                owner[smi] = name
                mol = Chem.MolFromSmiles(smi)
                library.append({
                    "name": name, "base_id": base_id(name), "parent_id": pid,
                    "zinc_id": pid if is_zinc(pid) else "",
                    "input_name": input_name.get(pid, pid),
                    "protomer_index": p_new, "stereo_index": s_idx,
                    "net_charge": Chem.GetFormalCharge(mol),
                    "heavy_atoms": mol.GetNumHeavyAtoms(),
                    "ph_values": head["ph_values"],
                    "population_estimate": head["population_estimate"],
                    "invertomer_of": inv if inv != name else "",
                    "note": head["note"], "smiles": smi,
                    # merged shards (pka_triage/merge_protomers.py) say which
                    # predictor decided the parent; plain QupKake shards don't
                    "site_source": head.get("site_source") or "qupkake"})
            p_new += 1

    # ---- invariants 1-3, asserted before anything final is written -----------
    names = [r["name"] for r in library]
    problems = []
    if any(len(n) != NAME_LEN or not CONTRACT.match(n) for n in names):
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
        raise SystemExit("contract violated -- nothing final written:\n  " + "\n  ".join(problems))

    # ---- write ---------------------------------------------------------------
    with open(outdir / "library.smi", "w") as fh:
        fh.writelines(f"{r['smiles']} {r['name']}\n" for r in library)
    cols = ["name", "base_id", "parent_id", "zinc_id", "input_name", "protomer_index",
            "stereo_index", "net_charge",
            "heavy_atoms", "ph_values", "population_estimate", "invertomer_of", "note",
            "smiles", "site_source"]
    with open(outdir / "library.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
        w.writeheader()
        w.writerows(library)
    with open(outdir / "dropped.tsv", "w") as fh:
        fh.write("parent_id\tprotomer_index\treason\tdetail\n")
        fh.writelines("\t".join(d) + "\n" for d in dropped)
        fh.writelines(f"{n}\t\t{r}\t\n" for n, r in qk_failed)

    # per-bin attrition: parents in -> parents with >= 1 named protomer
    attrition = {}
    if parents_rows:
        named = {r["parent_id"] for r in library}
        by_bin = defaultdict(lambda: Counter())
        for p in parents_rows:
            b = p["heavy_atoms"]
            by_bin[b]["parents"] += 1
            by_bin[b]["named"] += (p.get("parent_id") or p["zinc_id"]) in named
        attrition = {b: dict(c) for b, c in sorted(by_bin.items(), key=lambda x: int(x[0]))}

    reasons = Counter(d[2] for d in dropped) + Counter(r for _, r in qk_failed)
    stats = {
        "shards_read": len(found), "shards_missing": missing,
        "parents_named": len({r["parent_id"] for r in library}),
        "parents_hashed": len({r["parent_id"] for r in library if not r["zinc_id"]}),
        "names": len(library),
        "protomers": len({(r["parent_id"], r["protomer_index"]) for r in library}),
        "stereo_closure_added": closure_added,
        "n_h_invertomers": sum(1 for r in library if r["invertomer_of"]),
        "dropped": dict(reasons),
        "net_charge": dict(sorted(Counter(r["net_charge"] for r in library).items())),
        "attrition_by_heavy_atoms": attrition,
        "library_smi_sha256": hashlib.sha256((outdir / "library.smi").read_bytes()).hexdigest(),
    }
    (outdir / "stats.json").write_text(json.dumps(stats, indent=2, default=str))

    print(f"{stats['shards_read']} shard(s), {stats['parents_named']:,} parents -> "
          f"{stats['protomers']:,} protomers -> {stats['names']:,} names")
    print(f"  stereo closure added {closure_added:,} isomer(s) that --checkstereo would "
          f"otherwise have renamed")
    print(f"  of which {stats['n_h_invertomers']:,} are N-H+ invertomers of another row "
          f"(library.tsv: invertomer_of) -- one species in solution; collapse downstream")
    for k, v in sorted(reasons.items()):
        print(f"  dropped {k:<26} {v:,}")
    if missing:
        print(f"  ! PARTIAL: {len(missing)} shard(s) missing -- see stats.json")
    for b, c in attrition.items():
        print(f"  H{int(b):02d}  {c['named']:,} / {c['parents']:,} parents named")
    print(f"invariants 1-3 hold -> {outdir}/library.smi")
    print(f"\nnext:  llb validate-names --smi {outdir}/library.smi   "
          f"(DRAP's stage4/validate_names.py, found via DRAP_DIR in the config)")
    print(f"       llb check-stereo {outdir}/library.smi             (this repo, db2_converter env)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
