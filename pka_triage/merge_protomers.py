#!/usr/bin/env python3
"""
merge_protomers.py -- MolGpKa-decided + QupKake-routed protomers -> one shard set.

    llb merge --parents parents --triage triage --qupkake protomers_qk -o merged
    llb names --shards-dir merged --parents parents -o library

assign_names.py groups protomers.csv rows by parent across every shard it reads,
so a parent must come from exactly ONE source. This writes merged/shard_NNNNN/
for every shard of the ORIGINAL parents directory (so `llb names --parents
parents` still checks completeness against the real manifest), each holding:

  protomers.csv       per parent, one source:
                        route molgpka            -> triage rows      site_source=molgpka
                        route qupkake, succeeded -> QupKake rows     site_source=qupkake
                        route qupkake, failed    -> triage rows      site_source=molgpka_fallback
                                                    (--on-qupkake-fail molgpka, default)
                                                 or dropped, listed in qupkake_failed.tsv
                                                    (--on-qupkake-fail drop)
  qupkake_failed.tsv  dropped parents (assign_names.py reports them as dropped)
  DONE                which triage / QupKake shards it came from, and their sha256

Refuses, rather than guesses, when:
  * a triage shard is not DONE, or a routed parent has no finished QupKake shard
    (--allow-partial writes only the complete original shards; names then
    needs --allow-partial too, and the missing ones are listed);
  * QupKake ran at different pH values, margin, min population or max states
    than the triage (protomer sets would mean different things for different
    parents);
  * a QupKake shard's DONE does not match the routed manifest's sha256, or
    QupKake output contains a parent that was not routed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "run_qupkake"))

from qupkake_protomers import PROTOMER_COLUMNS  # noqa: E402

OUT_COLUMNS = PROTOMER_COLUMNS + ["site_source", "route_reason"]


def read_done(path: Path) -> dict:
    return dict(line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)


def read_manifest(d: Path) -> list[dict]:
    with open(d / "shards.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def shard_dir(root: Path, idx) -> Path:
    return root / f"shard_{int(idx):05d}"


# DONE key written by qupkake_array.sbatch -> triage settings key
SAME_SETTINGS = {"margin": "margin", "min_population": "min_population",
                 "max_states": "max_states"}


def load_qupkake(args, routed_ids: set[str], settings: dict):
    """-> (rows by parent, failures by parent, routed parents with no finished shard,
    provenance lines)."""
    rows: dict[str, list[dict]] = defaultdict(list)
    failed: dict[str, str] = {}
    unfinished: set[str] = set()
    prov: list[str] = []
    rdir = args.triage / "routed"
    if not routed_ids:
        return rows, failed, unfinished, prov
    if args.qupkake is None:
        raise SystemExit(f"{len(routed_ids)} parent(s) were routed to QupKake; pass "
                         f"--qupkake (the OUT dir of llb submit-qupkake {rdir} ...)")
    for r in read_manifest(rdir):
        d = shard_dir(args.qupkake, r["index"])
        shard_ids = [l.split()[1] for l in (rdir / r["path"]).read_text().splitlines()
                     if l.strip()]
        if not (d / "DONE").exists():
            unfinished.update(shard_ids)
            continue
        done = read_done(d / "DONE")
        if done.get("shard_sha256") != r["sha256"]:
            raise SystemExit(f"{d}/DONE was made from a different shard than "
                             f"{rdir}/{r['path']} (sha256 differs). Was --qupkake run on "
                             f"a different routed/ directory?")
        triage_ph = settings["ph"]
        qk_ph = sorted(float(x) for x in done.get("ph", "").split())
        if qk_ph != sorted(triage_ph):
            raise SystemExit(f"{d} ran at pH {qk_ph}, the triage at {sorted(triage_ph)}. "
                             f"Rerun QupKake with PH=\"{' '.join(map(str, triage_ph))}\".")
        for key, skey in SAME_SETTINGS.items():
            if key in done and float(done[key]) != float(settings[skey]):
                raise SystemExit(f"{d} ran with {key}={done[key]}, the triage with "
                                 f"{settings[skey]}; the two halves would build protomers "
                                 f"differently. Rerun one so they agree.")
        prov.append(f"qupkake_shard_{int(r['index']):05d}={r['sha256']}")
        with open(d / "protomers.csv") as fh:
            for row in csv.DictReader(fh):
                rows[row["parent_name"]].append(row)
        f = d / "qupkake_failed.tsv"
        if f.exists():
            with open(f) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    failed[row["name"]] = row["reason"]
    stray = (set(rows) | set(failed)) - routed_ids
    if stray:
        raise SystemExit(f"QupKake output has {len(stray)} parent(s) that were not routed, "
                         f"e.g. {sorted(stray)[:3]}")
    # a routed parent absent from a finished shard's output: QupKake rejected it
    # at input (qupkake_protomers rejected_input.smi), so it is a failure too
    for pid in routed_ids - unfinished - set(rows) - set(failed):
        failed[pid] = "no_qupkake_output"
    return rows, failed, unfinished, prov


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parents", type=Path, required=True, help="prepare_parents.py --outdir")
    ap.add_argument("--triage", type=Path, required=True, help="triage.py -o")
    ap.add_argument("--qupkake", type=Path, default=None,
                    help="OUT dir of the QupKake array run on TRIAGE/routed")
    ap.add_argument("-o", "--outdir", type=Path, default=Path("merged"))
    ap.add_argument("--on-qupkake-fail", choices=["molgpka", "drop"], default="molgpka")
    ap.add_argument("--allow-partial", action="store_true")
    args = ap.parse_args(argv)

    settings = json.loads((args.triage / "settings.json").read_text())
    manifest = read_manifest(args.parents)
    if not (args.triage / "routed" / "shards.tsv").exists() and \
            not (args.triage / "stats.json").exists():
        raise SystemExit(f"{args.triage} has not been finished: llb triage "
                         f"{args.parents} -o {args.triage} --finish")
    if (args.outdir / "merge_stats.json").exists():
        raise SystemExit(f"{args.outdir} already holds a merge; use a new -o")

    # ---- triage side --------------------------------------------------------
    decisions: dict[str, dict] = {}
    with open(args.triage / "decisions.tsv") as fh:
        for d in csv.DictReader(fh, delimiter="\t"):
            decisions[d["parent_id"]] = d
    routed_ids = {p for p, d in decisions.items() if d["route"] == "qupkake"}
    qk_rows, qk_failed, unfinished, qk_prov = load_qupkake(args, routed_ids, settings)

    shard_lines = {r["index"]: [l.split()[:2] for l in
                                (args.parents / r["path"]).read_text().splitlines() if l.strip()]
                   for r in manifest}
    incomplete = [(r["index"], n) for r in manifest
                  if (n := sum(1 for _, p in shard_lines[r["index"]] if p in unfinished))]
    if incomplete and not args.allow_partial:
        n = sum(k for _, k in incomplete)
        raise SystemExit(f"{n} routed parent(s) in {len(incomplete)} shard(s) have no finished "
                         f"QupKake shard yet (e.g. original shard {incomplete[0][0]}). Finish "
                         f"the QupKake array, or pass --allow-partial.")
    waiting_shards = {i for i, _ in incomplete}

    counts = Counter()
    written = 0
    for r in manifest:
        if r["index"] in waiting_shards:
            continue
        tdir = shard_dir(args.triage / "shards", r["index"])
        tdone = tdir / "DONE"
        if not tdone.exists() or read_done(tdone).get("shard_sha256") != r["sha256"]:
            raise SystemExit(f"{tdir} is missing or was made from a different shard; "
                             f"finish the triage first")
        by_parent: dict[str, list[dict]] = defaultdict(list)
        with open(tdir / "protomers.csv") as fh:
            for row in csv.DictReader(fh):
                by_parent[row["parent_name"]].append(row)
        smiles_of = {p: smi for smi, p in shard_lines[r["index"]]}
        shard_ids = list(smiles_of)

        out_rows, dropped = [], []
        for pid in shard_ids:
            d = decisions[pid]
            if d["route"] == "molgpka":
                src, rows = "molgpka", by_parent.get(pid, [])
            elif pid in qk_rows:
                src, rows = "qupkake", qk_rows[pid]
            elif args.on_qupkake_fail == "molgpka" and by_parent.get(pid):
                src, rows = "molgpka_fallback", by_parent[pid]
                counts[f"fallback:{qk_failed.get(pid, 'unknown')}"] += 1
            else:
                dropped.append((pid, qk_failed.get(pid, "qupkake_failed"), smiles_of[pid]))
                counts["dropped_qupkake_failed"] += 1
                continue
            counts[src] += 1
            for row in rows:
                out = {c: row.get(c, "") for c in PROTOMER_COLUMNS}
                if src == "molgpka_fallback":
                    out["note"] = ";".join(x for x in (row.get("note"), "qupkake_failed") if x)
                out["site_source"] = src
                out["route_reason"] = d["reason"]
                out_rows.append(out)

        odir = shard_dir(args.outdir, r["index"])
        odir.mkdir(parents=True, exist_ok=True)
        with open(odir / "protomers.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS)
            w.writeheader()
            w.writerows(out_rows)
        with open(odir / "qupkake_failed.tsv", "w") as fh:
            fh.write("name\treason\tsmiles\n")
            fh.writelines(f"{a}\t{b}\t{c}\n" for a, b, c in dropped)
        tsha = hashlib.sha256((tdir / "protomers.csv").read_bytes()).hexdigest()
        (odir / "DONE").write_text(
            f"shard_index={r['index']}\nshard_sha256={r['sha256']}\n"
            f"triage_protomers_sha256={tsha}\nqupkake_dir={args.qupkake or ''}\n"
            f"on_qupkake_fail={args.on_qupkake_fail}\n"
            f"finished={time.strftime('%Y-%m-%dT%H:%M:%S')}\n")
        written += 1

    stats = {"shards_written": written, "shards_waiting_on_qupkake": len(incomplete),
             "parents_by_source": dict(counts), "settings": settings,
             "qupkake_shards": len(qk_prov), "on_qupkake_fail": args.on_qupkake_fail}
    (args.outdir / "merge_stats.json").write_text(json.dumps(stats, indent=2))

    print(f"{written} shard(s) -> {args.outdir}  (of {len(manifest)})")
    for k, v in sorted(counts.items()):
        print(f"  {k:<32} {v:,}")
    if incomplete:
        print(f"  ! PARTIAL: {len(incomplete)} shard(s) still waiting on QupKake; "
              f"names needs --allow-partial", file=sys.stderr)
    print(f"\nnext: llb names --shards-dir {args.outdir} --parents {args.parents} -o library")
    return 0


if __name__ == "__main__":
    sys.exit(main())
