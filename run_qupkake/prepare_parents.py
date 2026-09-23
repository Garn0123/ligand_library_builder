#!/usr/bin/env python3
"""
prepare_parents.py -- ZINC SMILES -> validated, deduplicated, sharded parents.

The step between the download (mol_download/sample_2d.py or run.sh) and QupKake.
Everything QupKake is given should already be a molecule we intend to keep, under
a name that can become a contract name (NAMING_CONTRACT.md), because QupKake is
the expensive step and the naming step downstream refuses anything else.

What it does, in order:

  1. ID check.   Only ^ZINC[0-9A-Za-z]{12}$ passes. The contract's 12-character
                 base id IS name[4:16] of the ZINC id -- ZINC20 (12 digits) and
                 ZINC22 (12 base62) both fit -- so an id of any other shape has no
                 contract name and is rejected here, not three stages later.
  2. Standardize with qupkake_protomers.standardize (imported, not copied, so
                 the dedupe key is exactly the structure QupKake will see):
                 largest fragment, cleanup, neutralize.
  3. Dedupe.     Same id twice -> keep the first; if the two SMILES disagree that
                 is recorded as an id_conflict. Same neutral structure under two
                 ids -> keep the smallest id. Two ids for one structure would
                 become two matrix rows with DIFFERENT base ids, which the
                 identity split cannot group -- a train/val leak by construction.
  4. Bin check.  Heavy atoms of the neutral parent vs the H-bin in the input
                 filename (H17.smi -> 17). Mismatches are kept and flagged, not
                 dropped: the ladder bins on RDKit's count, and a shift is a
                 finding to report.
  5. Shard.      Shuffled with a fixed seed so every shard mixes sizes (xtb cost
                 grows with size; see mol_compiler/README.md on load balance),
                 then cut into --shard-size pieces. shards.tsv lists them with a
                 sha256 each; shards.tsv.meta carries the manifest's own checksum,
                 which the array job verifies before touching an index.

Does not do fingerprint/identity-group dedupe (SIZE_LADDER_SPEC.md 3); that is
stage4/make_split.py's grouping and belongs on the DRAP side.

Usage:
    python prepare_parents.py samples/H*.smi -o parents --shard-size 250
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem

from qupkake_protomers import standardize

ZINC_ID = re.compile(r"^ZINC[0-9A-Za-z]{12}$")
BIN_FROM_FILE = re.compile(r"H(\d\d)")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="'SMILES ZINCID' files, e.g. samples/H17.smi")
    ap.add_argument("-o", "--outdir", default="parents")
    ap.add_argument("--shard-size", type=int, default=250,
                    help="parents per QupKake array task")
    ap.add_argument("--seed", type=int, default=20260923)
    ap.add_argument("--preflight", type=int, default=20,
                    help="parents, spread across bins, written to preflight.smi")
    args = ap.parse_args(argv)

    outdir = Path(args.outdir)
    if (outdir / "shards.tsv").exists():
        raise SystemExit(f"{outdir}/shards.tsv exists. Regenerating it changes what "
                         f"every array index means; use a new --outdir.")
    (outdir / "shards").mkdir(parents=True, exist_ok=True)

    rejected: list[tuple[str, str, str]] = []
    by_id: dict[str, dict] = {}
    for path in map(Path, args.inputs):
        m = BIN_FROM_FILE.search(path.name)
        hbin = int(m.group(1)) if m else None
        for lineno, line in enumerate(open(path), 1):
            parts = line.split()
            if not parts or parts[0].startswith("#"):
                continue
            smi, zid = parts[0], parts[1] if len(parts) > 1 else ""
            where = f"{path.name}:{lineno}"
            if not ZINC_ID.match(zid):
                rejected.append((zid or where, "bad_id", f"{where} {smi}"))
                continue
            mol = standardize(smi)
            if mol is None:
                rejected.append((zid, "unparsable", smi))
                continue
            std = Chem.MolToSmiles(mol)
            if zid in by_id:
                prev = by_id[zid]
                reason = "duplicate_id" if prev["smiles"] == std else "id_conflict"
                rejected.append((zid, reason, f"{where} {std} (kept {prev['source']} "
                                              f"{prev['smiles']})"))
                continue
            ha = mol.GetNumHeavyAtoms()
            by_id[zid] = {"zinc_id": zid, "smiles": std, "input_smiles": smi,
                          "heavy_atoms": ha, "file_bin": hbin, "source": where,
                          "bin_match": "" if hbin is None else int(ha == hbin)}

    by_structure = defaultdict(list)
    for zid, rec in by_id.items():
        by_structure[rec["smiles"]].append(zid)
    for smi, ids in by_structure.items():
        if len(ids) > 1:
            keep, *drop = sorted(ids)
            for zid in drop:
                rejected.append((zid, "duplicate_structure", f"same as {keep}: {smi}"))
                del by_id[zid]

    parents = sorted(by_id.values(), key=lambda r: r["zinc_id"])
    with open(outdir / "parents.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(parents[0]) if parents
                           else ["zinc_id"])
        w.writeheader()
        w.writerows(parents)
    with open(outdir / "rejected.tsv", "w") as fh:
        fh.write("zinc_id\treason\tdetail\n")
        fh.writelines(f"{a}\t{b}\t{c}\n" for a, b, c in rejected)

    order = list(parents)
    random.Random(args.seed).shuffle(order)
    shards = [order[i:i + args.shard_size] for i in range(0, len(order), args.shard_size)]
    rows = []
    for i, shard in enumerate(shards):
        p = outdir / "shards" / f"shard_{i:05d}.smi"
        p.write_text("".join(f"{r['smiles']} {r['zinc_id']}\n" for r in shard))
        rows.append((i, p.relative_to(outdir), len(shard), sha256(p)))
    with open(outdir / "shards.tsv", "w") as fh:
        fh.write("index\tpath\tn_parents\tsha256\n")
        fh.writelines(f"{i}\t{p}\t{n}\t{h}\n" for i, p, n, h in rows)
    (outdir / "shards.tsv.meta").write_text(
        f"checksum={sha256(outdir / 'shards.tsv')}\nn_shards={len(rows)}\n"
        f"shard_size={args.shard_size}\nseed={args.seed}\nn_parents={len(parents)}\n")

    # Preflight set: round-robin across bins so the largest molecules -- the
    # slowest for xtb and the likeliest to fail -- are in the timing run.
    per_bin = defaultdict(list)
    for r in order:
        per_bin[r["heavy_atoms"]].append(r)
    pre, bins = [], sorted(per_bin, reverse=True)
    while len(pre) < min(args.preflight, len(order)):
        for b in bins:
            if per_bin[b] and len(pre) < args.preflight:
                pre.append(per_bin[b].pop())
    (outdir / "preflight.smi").write_text(
        "".join(f"{r['smiles']} {r['zinc_id']}\n" for r in pre))

    # ---- report -----------------------------------------------------------
    reasons = Counter(r for _, r, _ in rejected)
    print(f"parents kept        {len(parents):,}")
    for k in ("bad_id", "unparsable", "duplicate_id", "id_conflict", "duplicate_structure"):
        print(f"  rejected {k:<20} {reasons.get(k, 0):,}")
    bins = Counter((r["file_bin"], r["heavy_atoms"]) for r in parents)
    print("\n  file bin   heavy atoms   n")
    for (fb, ha), n in sorted(bins.items(), key=lambda x: (x[0][0] or 0, x[0][1])):
        flag = "" if fb is None or fb == ha else "   <- bin mismatch"
        print(f"  {fb!s:>8}   {ha:>11}   {n:,}{flag}")
    print(f"\n{len(rows)} shard(s) of <= {args.shard_size} -> {outdir}/shards.tsv")
    print(f"preflight set ({len(pre)}) -> {outdir}/preflight.smi")
    if reasons.get("id_conflict"):
        print(f"\n  ! {reasons['id_conflict']} id(s) arrived with two different "
              f"structures; see rejected.tsv", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
