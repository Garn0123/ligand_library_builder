#!/usr/bin/env python3
"""
triage.py -- MolGpKa on every parent; QupKake only for the ones it can't settle.

    llb triage parents -o triage                      # whole library, one process
    llb triage parents -o triage --part 3/20          # one slice (array task)
    llb triage parents -o triage --finish             # after every slice is DONE

Reads a prepare_parents.py directory (parents.tsv, shards.tsv, shards/) and, for
each parent shard, writes  triage/shards/shard_NNNNN/:

  sites.tsv       every MolGpKa site: parent_id, idx, kind, pka, kept (for the
                  harness; kept=0 = basic site on an amide N, see NON_BASIC_N)
  decisions.tsv   per parent: route (molgpka | qupkake), reason, min_gap, ...
  protomers.csv   MolGpKa protomers for EVERY parent, routed ones included (the
                  merge falls back to them if QupKake fails, and the harness
                  compares them with QupKake's). Same columns as
                  qupkake_protomers.py writes, plus site_source and route.
  DONE            shard sha256, settings hash, MolGpKa weight hashes

and, once every shard is DONE (automatic without --part):

  routed/         the parents sent to QupKake, as a prepare_parents.py directory
                  (parents.tsv, shards/, shards.tsv + .meta, preflight.smi):
                  llb submit-qupkake triage/routed protomers_qk --time=...
  decisions.tsv   all shards' decisions in one file
  stats.json      counts by reason, and the fraction routed at other windows

Routing (per parent, over every target pH):
  near_window   some site has |pKa - pH| <= --route-window
  coupled       two sites of the same kind, both on the charged side at some pH,
                within --coupling-bonds bonds (piperazine: both N predicted ~10.5,
                so the state model would build the dication). See TRIAGE.md:
                QupKake as we run it is ALSO independent-site, so routing these
                is a flag for review rather than a fix. --coupled ignore turns
                it off.
  molgpka_failed  MolGpKa raised on the molecule
Everything else is decided by MolGpKa; molecules with no MolGpKa site are kept
unchanged (note no_sites_predicted) and counted separately in stats.json.

Protomers are built by the same code as the QupKake path
(qupkake_protomers.build_protomers: Henderson-Hasselbalch per site, --margin,
populations, --max-states), so the two paths differ only in where the site
pKas come from. --route-window decides which MOLECULES go to QupKake; --margin
decides which SITES are enumerated both ways. They are separate on purpose.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "run_qupkake"))
sys.path.insert(0, str(HERE))

from rdkit import Chem  # noqa: E402

from molgpka import MolGpKa  # noqa: E402
from prepare_parents import write_shards  # noqa: E402
from qupkake_protomers import (PROTOMER_COLUMNS, Site, build_protomers,  # noqa: E402
                               make_name, protomer_row, unchanged_protomers)

TRIAGE_COLUMNS = ["site_source", "route"]

# MolGpKa's base SMARTS also match the N of amides (anilides, benzamides) and
# predicts pKa 3.4-6 there; N-protonation of an amide is about -1, so those
# sites are artifacts that build C(=O)[NH2+] protomers. Dropped by default
# (--keep-amide-bases to keep), and kept in sites.tsv with kept=0.
NON_BASIC_N = Chem.MolFromSmarts("[#7;$([#7]-[#6]=[#8,#16]),$([#7]-[#16](=[#8])=[#8])]")
DECISION_COLUMNS = ["parent_id", "route", "reason", "n_sites", "min_gap",
                    "nearest_site", "coupled_pairs", "n_protomers", "dropped_sites", "note"]
CURVE_WINDOWS = [0.25 * k for k in range(0, 17)]          # 0 .. 4 pKa units


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_manifest(parents: Path) -> list[dict]:
    with open(parents / "shards.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def settings_of(args) -> dict:
    return {"ph": args.ph, "route_window": args.route_window, "margin": args.margin,
            "coupling_bonds": args.coupling_bonds, "coupled": args.coupled,
            "keep_amide_bases": args.keep_amide_bases,
            "min_population": args.min_population, "max_states": args.max_states,
            "max_ambiguous": args.max_ambiguous}


def coupled_pairs(sites: list[Site], dist, ph_values, max_bonds: int) -> list[str]:
    """Same-kind sites both on the charged side at some pH, within max_bonds."""
    if max_bonds <= 0:
        return []
    out = set()
    for ph in ph_values:
        charged = [s for s in sites if s.changed_fraction(ph) >= 0.5]
        for i, a in enumerate(charged):
            for b in charged[i + 1:]:
                if a.kind == b.kind and a.idx != b.idx:
                    d = int(dist[a.idx][b.idx])
                    if d <= max_bonds:
                        lo, hi = sorted((a.idx, b.idx))
                        out.add(f"{a.kind[0]}{lo}-{hi}:{d}")
    return sorted(out)


def triage_shard(rows: list[tuple[str, str]], model: MolGpKa, args):
    """-> (decisions, site rows, protomers) for one shard's (smiles, parent_id)."""
    decisions, site_rows, protomers = [], [], []
    for smi, pid in rows:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            decisions.append({"parent_id": pid, "route": "qupkake",
                              "reason": "unparsable", "n_sites": 0})
            continue
        try:
            ref, found = model.predict(mol)
        except Exception as e:                       # noqa: BLE001 -- record, route on
            decisions.append({"parent_id": pid, "route": "qupkake",
                              "reason": "molgpka_failed", "n_sites": 0,
                              "note": f"{type(e).__name__}: {e}"[:200]})
            continue
        non_basic = set() if args.keep_amide_bases else {
            m[0] for m in ref.GetSubstructMatches(NON_BASIC_N)}
        sites, dropped = [], []
        for s in found:
            keep = not (s.kind == "basic" and s.idx in non_basic)
            (sites if keep else dropped).append(Site(idx=s.idx, kind=s.kind, pka=s.pka))
            site_rows.append((pid, s.idx, s.kind, f"{s.pka:.3f}", int(keep)))
        isomers = [mol]                              # --stereo carry, as the array runs

        if not sites:
            emitted = unchanged_protomers(pid, isomers, args.ph)
            gap, nearest, pairs, used_isomer = math.inf, "", [], True
        else:
            emitted, _, _, used_isomer = build_protomers(
                pid, isomers, ref, sites, args.ph, args.margin, args.min_population,
                args.max_states, args.max_ambiguous, "carry",
                fallback_note="built_on_molgpka_parent")
            gap, near_site = min((abs(s.pka - ph), s) for s in sites for ph in args.ph)
            nearest = f"{near_site.kind[0]}{near_site.idx}:{near_site.pka:.2f}"
            pairs = coupled_pairs(sites, Chem.GetDistanceMatrix(ref), args.ph,
                                  args.coupling_bonds)

        reasons = []
        if gap <= args.route_window:
            reasons.append("near_window")
        if pairs and args.coupled == "route":
            reasons.append("coupled")
        if sites and not emitted:
            reasons.append("no_state_built")
        route = "qupkake" if reasons else "molgpka"
        for p in emitted:
            p.name = make_name(pid, p.protomer_index, p.stereo_index, "index", p, 1)
            protomers.append((p, route))
        decisions.append({
            "parent_id": pid, "route": route, "reason": "+".join(reasons),
            "n_sites": len(sites),
            "min_gap": "" if gap == math.inf else f"{gap:.3f}",
            "nearest_site": nearest, "coupled_pairs": ";".join(pairs),
            "n_protomers": len(emitted),
            "dropped_sites": ";".join(f"{d.kind[0]}{d.idx}:{d.pka:.2f}" for d in dropped),
            "note": ("no_sites_predicted" if not sites else
                     "" if used_isomer else "built_on_molgpka_parent")})
    return decisions, site_rows, protomers


def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols, restval="")
        w.writeheader()
        w.writerows(rows)


def run_part(args, manifest, settings_sha) -> int:
    k, n = args.part
    mine = [r for r in manifest if int(r["index"]) % n == k]
    model = MolGpKa(args.molgpka_dir)
    t0, done_here, mols = time.time(), 0, 0
    for r in mine:
        idx = int(r["index"])
        out = args.outdir / "shards" / f"shard_{idx:05d}"
        done = out / "DONE"
        if done.exists():
            have = dict(l.split("=", 1) for l in done.read_text().split("\n") if "=" in l)
            if have.get("shard_sha256") == r["sha256"] and have.get("settings_sha256") == settings_sha:
                continue
            raise SystemExit(f"{done} is from a different shard or different settings; "
                             f"use a new -o, or delete {out} to redo it")
        shard = args.parents / r["path"]
        if sha256_file(shard) != r["sha256"]:
            raise SystemExit(f"{shard} does not match its manifest sha256")
        rows = [tuple(l.split()[:2]) for l in shard.read_text().splitlines() if l.strip()]
        decisions, site_rows, protomers = triage_shard(rows, model, args)

        out.mkdir(parents=True, exist_ok=True)
        write_tsv(out / "decisions.tsv", DECISION_COLUMNS, decisions)
        with open(out / "sites.tsv", "w") as fh:
            fh.write("parent_id\tidx\tkind\tpka\tkept\n")
            fh.writelines("\t".join(map(str, s)) + "\n" for s in site_rows)
        with open(out / "protomers.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(PROTOMER_COLUMNS + TRIAGE_COLUMNS)
            w.writerows(protomer_row(p) + ["molgpka", route] for p, route in protomers)
        (out / "DONE").write_text(
            f"shard_index={idx}\nshard_sha256={r['sha256']}\n"
            f"settings_sha256={settings_sha}\n"
            + "".join(f"{f}={h}\n" for f, h in sorted(model.weights_sha.items()))
            + f"molgpka_dir={model.root}\nfinished={time.strftime('%Y-%m-%dT%H:%M:%S')}\n")
        done_here += 1
        mols += len(rows)
    dt = time.time() - t0
    rate = f", {1000 * dt / mols:.1f} ms/parent" if mols else ""
    print(f"part {k}/{n}: {done_here} shard(s) triaged now, {len(mine) - done_here} "
          f"already DONE; {mols:,} parents in {dt:.0f} s{rate}")
    return 0


def finish(args, manifest) -> int:
    missing = [r["index"] for r in manifest
               if not (args.outdir / "shards" / f"shard_{int(r['index']):05d}" / "DONE").exists()]
    if missing:
        raise SystemExit(f"{len(missing)} shard(s) not triaged yet, e.g. {missing[:5]}; "
                         f"run the remaining --part slices first")
    decisions = []
    for r in manifest:
        with open(args.outdir / "shards" / f"shard_{int(r['index']):05d}" / "decisions.tsv") as fh:
            decisions += list(csv.DictReader(fh, delimiter="\t"))
    write_tsv(args.outdir / "decisions.tsv", DECISION_COLUMNS, decisions)

    with open(args.parents / "parents.tsv") as fh:
        parents = list(csv.DictReader(fh, delimiter="\t"))
    routed_ids = {d["parent_id"] for d in decisions if d["route"] == "qupkake"}
    routed = [p for p in parents if p["parent_id"] in routed_ids]
    rdir = args.outdir / "routed"
    if (rdir / "shards.tsv").exists():
        raise SystemExit(f"{rdir}/shards.tsv exists; QupKake array indices may already "
                         f"refer to it. Use a new -o to re-route.")
    rdir.mkdir(parents=True, exist_ok=True)
    if routed:
        with open(rdir / "parents.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(parents[0]))
            w.writeheader()
            w.writerows(routed)
        rows, pre = write_shards(rdir, routed, args.qupkake_shard_size, args.seed,
                                 args.preflight)
    else:
        rows, pre = [], []

    # stats, including what other windows would have routed
    n = len(decisions)
    reasons = Counter(d["reason"] or "decided" for d in decisions)
    gaps = [float(d["min_gap"]) for d in decisions if d["min_gap"]]
    coupled = sum(1 for d in decisions if d["coupled_pairs"])
    hard = sum(1 for d in decisions if d["reason"] in ("molgpka_failed", "unparsable",
                                                       "no_state_built"))
    curve = {}
    for w in CURVE_WINDOWS:
        near = {d["parent_id"] for d in decisions if d["min_gap"] and float(d["min_gap"]) <= w}
        with_c = near | {d["parent_id"] for d in decisions if d["coupled_pairs"]}
        curve[f"{w:g}"] = {"near_window": round(len(near) / n, 4) if n else 0,
                           "near_window_or_coupled": round(len(with_c) / n, 4) if n else 0}
    stats = {
        "parents": n, "routed_to_qupkake": len(routed),
        "routed_fraction": round(len(routed) / n, 4) if n else 0,
        "by_reason": dict(reasons),
        "no_sites_predicted": sum(1 for d in decisions if d["note"] == "no_sites_predicted"),
        "parents_with_coupled_pairs": coupled,
        "failed_or_unbuildable": hard,
        "min_gap_median": round(sorted(gaps)[len(gaps) // 2], 3) if gaps else None,
        "routed_fraction_by_window": curve,
        "settings": settings_of(args),
        "routed_shards": len(rows), "qupkake_shard_size": args.qupkake_shard_size,
    }
    (args.outdir / "stats.json").write_text(json.dumps(stats, indent=2))

    print(f"{n:,} parents: {n - len(routed):,} decided by MolGpKa, "
          f"{len(routed):,} routed to QupKake ({stats['routed_fraction']:.1%})")
    for k_, v in sorted(reasons.items(), key=lambda x: -x[1]):
        print(f"  {k_:<26} {v:,}")
    print(f"  ({stats['no_sites_predicted']:,} decided parents had no MolGpKa site)")
    print("  routed fraction if --route-window were:  " + "  ".join(
        f"{w}:{c['near_window_or_coupled' if args.coupled == 'route' else 'near_window']:.0%}"
        for w, c in curve.items() if float(w) in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)))
    if routed:
        print(f"\nnext: llb submit-qupkake {rdir} protomers_qk --time=...   "
              f"({len(rows)} shard(s) of <= {args.qupkake_shard_size}; time from the "
              f"preflight, {rdir}/preflight.smi)")
        print(f"then: llb merge --parents {args.parents} --triage {args.outdir} "
              f"--qupkake protomers_qk -o merged")
    else:
        print(f"\nnothing routed. next: llb merge --parents {args.parents} "
              f"--triage {args.outdir} -o merged")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("parents", type=Path, help="prepare_parents.py --outdir")
    ap.add_argument("-o", "--outdir", type=Path, default=Path("triage"))
    ap.add_argument("--ph", type=float, nargs="+", default=[6.4, 7.4, 8.4])
    ap.add_argument("--route-window", type=float, default=1.5,
                    help="route a parent when any site is within this of any pH")
    ap.add_argument("--margin", type=float, default=1.0,
                    help="enumerate a site both ways within this of the pH "
                         "(same meaning and default as the QupKake path)")
    ap.add_argument("--coupling-bonds", type=int, default=3,
                    help="same-kind charged sites this close count as coupled (0 = off)")
    ap.add_argument("--coupled", choices=["route", "ignore"], default="route")
    ap.add_argument("--keep-amide-bases", action="store_true",
                    help="keep MolGpKa's basic sites on amide/sulfonamide N")
    ap.add_argument("--min-population", type=float, default=0.01)
    ap.add_argument("--max-states", type=int, default=8)
    ap.add_argument("--max-ambiguous", type=int, default=6)
    ap.add_argument("--part", default=None, metavar="K/N",
                    help="only parent shards with index %% N == K; then --finish")
    ap.add_argument("--finish", action="store_true",
                    help="only collect: decisions, routed/, stats.json")
    ap.add_argument("--qupkake-shard-size", type=int, default=10,
                    help="parents per QupKake task in routed/ (one core each)")
    ap.add_argument("--preflight", type=int, default=10,
                    help="routed parents written to routed/preflight.smi")
    ap.add_argument("--seed", type=int, default=20260923)
    ap.add_argument("--molgpka-dir", default=os.environ.get("MOLGPKA_DIR"))
    args = ap.parse_args(argv)

    manifest = read_manifest(args.parents)
    settings = settings_of(args)
    settings_sha = hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()
    args.outdir.mkdir(parents=True, exist_ok=True)
    sfile = args.outdir / "settings.json"
    if sfile.exists():
        prev = json.loads(sfile.read_text())
        if prev != settings:
            raise SystemExit(f"{args.outdir} was triaged with {prev};\nthis run asks for "
                             f"{settings}. Use a new -o for new settings.")
    else:
        tmp = sfile.with_suffix(f".{os.getpid()}.tmp")
        tmp.write_text(json.dumps(settings, indent=2))
        os.replace(tmp, sfile)

    if args.finish:
        return finish(args, manifest)
    if args.part:
        k, n = (int(x) for x in args.part.split("/"))
        if not 0 <= k < n:
            raise SystemExit("--part K/N needs 0 <= K < N")
        args.part = (k, n)
        return run_part(args, manifest, settings_sha)
    args.part = (0, 1)
    run_part(args, manifest, settings_sha)
    return finish(args, manifest)


if __name__ == "__main__":
    sys.exit(main())
