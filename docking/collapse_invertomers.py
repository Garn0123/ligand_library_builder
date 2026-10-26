#!/usr/bin/env python3
"""
collapse_invertomers.py -- one score row per species: the best of each N-H+ invertomer group.

    llb collapse-invertomers scores.parquet --library library/library.tsv -o scores_collapsed.parquet
    llb collapse-invertomers scores.tsv --library library/library.tsv -o out.tsv --long name,target,grid_score

A protonated amine with three different substituents is a stereocentre in a
db2: conformers never invert it, so `llb names` names both configurations and
both are built and docked (a piperidinium's substituent axial or equatorial,
the N-H donor pointing one way or the other). In water they interconvert, so
they are one species, and a binding site picks whichever fits. This keeps the
docking of both and collapses afterwards: per species and per target, the best
score of the group (lower is better, DOCK's Grid_Score; --higher-is-better
flips it). library.tsv's invertomer_of links each group to its first member,
whose name the collapsed row carries.

Input, either shape:
  wide  one row per name, one numeric column per target; the name column is
        --id-col (default: molecule_id, name or the first column)
  long  --long NAME,TARGET,SCORE: one row per (name, target)
Format by extension: .parquet (needs pyarrow), .csv, anything else tab-separated.

Writes OUT (same shape as the input) and OUT's stem + .invertomers.tsv: for
every group of two or more names, per target, the winner and the margin over
the group's first member. The summary line answers whether docking both
configurations is worth it: how often a non-first invertomer wins, and by how
much. A name not in library.tsv passes through as its own species (counted).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="," if path.suffix == ".csv" else "\t")


def write_table(df: pd.DataFrame, path: Path) -> None:
    if path.suffix == ".parquet":
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, sep="," if path.suffix == ".csv" else "\t", index=False,
                  lineterminator="\n")


def species_map(library: Path) -> dict[str, str]:
    lib = pd.read_csv(library, sep="\t", dtype=str, keep_default_na=False)
    if "invertomer_of" not in lib.columns:
        raise SystemExit(f"{library} has no invertomer_of column (made before llb names "
                         f"marked invertomers?)")
    return {n: (inv or n) for n, inv in zip(lib["name"], lib["invertomer_of"])}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scores", type=Path)
    ap.add_argument("--library", type=Path, required=True, help="library.tsv from llb names")
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--long", default=None, metavar="NAME,TARGET,SCORE",
                    help="long input: the name, target and score columns")
    ap.add_argument("--id-col", default=None, help="wide input: the name column")
    ap.add_argument("--higher-is-better", action="store_true")
    args = ap.parse_args(argv)

    species = species_map(args.library)
    df = read_table(args.scores)

    if args.long:
        name_c, tgt_c, score_c = args.long.split(",")
        long = df[[name_c, tgt_c, score_c]].rename(
            columns={name_c: "name", tgt_c: "target", score_c: "score"})
    else:
        id_c = args.id_col or next((c for c in ("molecule_id", "name") if c in df.columns),
                                   df.columns[0])
        targets = [c for c in df.columns if c != id_c and pd.api.types.is_numeric_dtype(df[c])]
        long = df.melt(id_vars=[id_c], value_vars=targets, var_name="target",
                       value_name="score").rename(columns={id_c: "name"})
    long["name"] = long["name"].astype(str)
    df_names = long["name"].nunique()
    unknown = sorted(set(long["name"]) - set(species))
    long["species"] = long["name"].map(species).fillna(long["name"])
    all_species = list(dict.fromkeys(long["species"]))   # rows with no score stay, all-NaN
    long = long.dropna(subset=["score"])

    # best per (species, target); ties keep the group's first member
    sign = -1.0 if args.higher_is_better else 1.0
    long["_key"] = sign * long["score"]
    long["_first"] = (long["name"] != long["species"]).astype(int)
    best = (long.sort_values(["species", "target", "_key", "_first"])
                .drop_duplicates(["species", "target"]))

    # provenance + the number that says whether docking both is worth it
    sizes = long.groupby("species")["name"].nunique()
    multi = set(sizes[sizes > 1].index)
    first = long[long["_first"] == 0].set_index(["species", "target"])["score"]
    prov = best[best["species"].isin(multi)][["species", "target", "name", "score"]].copy()
    prov["first_member_score"] = first.reindex(
        list(zip(prov["species"], prov["target"]))).to_numpy()
    prov["gain"] = sign * (prov["first_member_score"] - prov["score"])
    prov = prov.rename(columns={"name": "winner"})
    prov_path = args.out.with_name(args.out.name.split(".")[0] + ".invertomers.tsv")
    prov.to_csv(prov_path, sep="\t", index=False, lineterminator="\n")

    if args.long:
        out = best[["species", "target", "score"]].rename(
            columns={"species": name_c, "target": tgt_c, "score": score_c})
    else:
        out = (best.pivot(index="species", columns="target", values="score")
                   .reindex(index=all_species, columns=targets).rename_axis(None, axis=1)
                   .reset_index().rename(columns={"species": id_c}))
    write_table(out, args.out)

    n_in, n_out = df_names, len(all_species)
    other = prov[prov["winner"] != prov["species"]]
    scored_pairs = prov.dropna(subset=["first_member_score"])
    print(f"{n_in:,} names -> {n_out:,} species ({len(multi):,} invertomer groups collapsed)"
          f" -> {args.out}")
    if len(prov):
        g = other["gain"].dropna()
        print(f"  a non-first invertomer won {len(other):,} of {len(prov):,} (group, target) "
              f"cells ({100 * len(other) / len(prov):.1f}%)"
              + (f"; gain over the first: median {g.median():.2f}, p90 {g.quantile(0.9):.2f}, "
                 f"max {g.max():.2f}" if len(g) else "")
              + (f"; {len(prov) - len(scored_pairs):,} cells had no score for the first"
                 if len(prov) != len(scored_pairs) else ""))
        print(f"  per (group, target) winner and margin -> {prov_path}")
    if unknown:
        print(f"  ! {len(unknown):,} name(s) not in {args.library} passed through as their "
              f"own species, e.g. {unknown[:3]}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
