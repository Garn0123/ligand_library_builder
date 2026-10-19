#!/usr/bin/env python3
"""
index.py -- find the exact ZINC22 chunk(s) holding each wanted id.

Two passes over small downloads, each followed by a fetch with run.sh:

  groups  Decide which chunks' sidecars to fetch. Writes
            plan/sidecar_urls.txt  a sidecar (.pdbqt.tgz, else .solv.tgz) per candidate chunk
            plan/groups.tsv        zinc_id gen tranche charge     (one row per hit)
            plan/no3d.txt          wanted ids in no list: ZINC has no 3D build of them
            plan/groups.mode       how the candidates were chosen (read by `chunks`)
          Candidates are, by default, every chunk of every (generation, tranche, charge)
          group whose id list (plan/idlist_urls.txt, fetched under --root) holds a
          wanted id. --no-idlists skips the lists: every chunk of every tranche a
          wanted id falls in. That is the right call for a large uniform draw (50k LLB
          parents hit nearly every group anyway, so the lists, ~40 GB for H18-22,
          narrow nothing); groups.tsv and no3d.txt are then not written.

          --fraction F keeps a seeded random F of the candidate chunks (hash of URL +
          --seed, so raising F later keeps every sidecar already fetched). Sidecars
          and db2s then both cost ~F of the full bill and the parents found are an
          unbiased ~F sample of those ZINC built: the way to fit a budget.

  chunks  List the members of the fetched sidecars (they hold exactly the db2's
          members, at ~0.5% of the bytes). Writes
            plan/index.tsv         zinc_id protomer charge gen tranche chunk member db2_url
                                   one row per db2 member (protomer) ZINC built for a
                                   wanted id; a member holds 1-3+ db2 records
            plan/db2_urls.txt      the db2.tgz chunks to fetch
            plan/coverage.tsv      per bin: wanted, in a list, in a chunk, chunks, GB
            plan/listed_not_built.txt  in an id list but in no chunk (~5%: build failures)
          With --budget-gb N, chunks are taken in seeded random order (same hash as
          --fraction) until N GB, keeping the covered ids an unbiased sample; with
          --greedy, by wanted ids per byte instead (more ids per GB, but it favours
          small chunks: charged and sparse groups). Ids left out: plan/over_budget.txt.

          --all-members also writes plan/members.tsv.gz: every member of every listed
          sidecar, wanted or not (a full id -> chunk index of what was fetched).

Usage:
    python3 index.py groups --plan plan --root zinc_mirror
    python3 index.py groups --plan plan --no-idlists --fraction 0.05
    python3 index.py chunks --plan plan --root zinc_mirror [--budget-gb 500]
"""
from __future__ import annotations

import argparse
import gzip
import heapq
import sys
import tarfile
import zlib
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from zid import MEMBER, local_path, rank, read_tsv, tranche_matches  # noqa: E402


def load_wanted(plan: Path) -> dict[str, str]:
    return {r["zinc_id"]: r["bin"] for r in read_tsv(plan / "wanted.tsv")}


def groups(a) -> None:
    wanted = load_wanted(a.plan)
    cat = read_tsv(a.catalog / "chunks.tsv")
    if a.no_idlists:
        by_bin = defaultdict(list)
        for z, b in wanted.items():
            by_bin[b].append(z)
        tranches = {c["tranche"] for c in cat}
        touched = {t for t in tranches
                   if any(tranche_matches(z, t) for z in by_bin.get(t[:3], ()))}
        cand = [c for c in cat if c["tranche"] in touched]
        what = f"every chunk of the {len(touched):,} tranches the ids fall in"
        (a.plan / "groups.tsv").write_text("zinc_id\tgen\ttranche\tcharge\n")
        (a.plan / "no3d.txt").unlink(missing_ok=True)
    else:
        urls = (a.plan / "idlist_urls.txt").read_text().split()
        hits, missing_lists = [], []
        for url in urls:
            p = local_path(url, a.root)
            if not p.exists():
                missing_lists.append(url)
                continue
            # .../H20P200/H20P200-N.f.txt
            tranche, rest = p.name.split("-", 1)
            charge, gen = rest.split(".")[0], rest.split(".")[1]
            with open(p) as fh:
                for line in fh:
                    z = line.strip()
                    if z in wanted:
                        hits.append((z, gen, tranche, charge))
        if missing_lists:
            (a.plan / "idlists_missing.txt").write_text("".join(f"{u}\n" for u in missing_lists))
            print(f"WARNING: {len(missing_lists):,} of {len(urls):,} id lists not on disk under "
                  f"{a.root} (plan/idlists_missing.txt). Fetch them first, or the groups are "
                  f"incomplete.", file=sys.stderr)
        with open(a.plan / "groups.tsv", "w") as fh:
            fh.write("zinc_id\tgen\ttranche\tcharge\n")
            fh.writelines("\t".join(h) + "\n" for h in sorted(hits))
        found = {h[0] for h in hits}
        (a.plan / "no3d.txt").write_text("".join(f"{z}\n" for z in sorted(set(wanted) - found)))
        print(f"{len(found):,} of {len(wanted):,} wanted ids are in a ZINC 3D id list "
              f"({len(wanted) - len(found):,} -> no3d.txt)", file=sys.stderr)
        keys = {(g, t, c) for _, g, t, c in hits}
        cand = [c for c in cat if (c["gen"], c["tranche"], c["charge"]) in keys]
        what = f"every chunk of the {len(keys):,} groups holding a wanted id"

    n_all = len(cand)
    db2_all = sum(int(c["db2_bytes"]) for c in cand)
    if a.fraction < 1:
        cand = [c for c in cand if rank(c["db2_url"], a.seed) < a.fraction]
    side, side_bytes = [], 0
    for c in cand:
        kind = "pdbqt" if int(c["pdbqt_bytes"]) else "solv"
        side.append(c["db2_url"][:-len("db2.tgz")] + f"{kind}.tgz")
        side_bytes += int(c[f"{kind}_bytes"])
    (a.plan / "sidecar_urls.txt").write_text("".join(f"{u}\n" for u in sorted(side)))
    (a.plan / "groups.mode").write_text(
        f"idlists={'no' if a.no_idlists else 'yes'}\nfraction={a.fraction}\nseed={a.seed}\n")
    print(f"candidates: {what}: {n_all:,} chunks, {db2_all / 1e12:.2f} TB of db2", file=sys.stderr)
    print(f"kept {len(cand):,} (fraction {a.fraction}, seed {a.seed!r}): {len(side):,} sidecars, "
          f"{side_bytes / 1e9:.1f} GB -> {a.plan}/sidecar_urls.txt", file=sys.stderr)


def chunks(a) -> None:
    wanted = load_wanted(a.plan)
    listed = {r["zinc_id"] for r in read_tsv(a.plan / "groups.tsv")}
    mode = dict(line.split("=", 1) for line in
                (a.plan / "groups.mode").read_text().split()) if (a.plan / "groups.mode").exists() else {}
    used_lists = mode.get("idlists", "yes") == "yes"
    frac = float(mode.get("fraction", 1))
    cat = {c["db2_url"][:-len(".db2.tgz")]: c for c in read_tsv(a.catalog / "chunks.tsv")}
    rows, bad, missing = [], [], 0
    all_fh = gzip.open(a.plan / "members.tsv.gz", "wt") if a.all_members else None
    if all_fh:
        all_fh.write("zinc_id\tprotomer\tcharge\tdb2_url\tmember\n")
    for url in (a.plan / "sidecar_urls.txt").read_text().split():
        p = local_path(url, a.root)
        stem = url.rsplit(".", 2)[0]                  # drop .pdbqt.tgz / .solv.tgz
        c = cat.get(stem)
        if c is None or not p.exists():
            missing += 1
            continue
        try:
            with tarfile.open(p, "r|gz") as tf:
                for m in tf:
                    mm = MEMBER.search(m.name)
                    if not mm:
                        continue
                    z, prot, charge = mm[1], mm[2], mm[3]
                    member = m.name.rsplit(".", 1)[0] + ".db2"
                    if all_fh:
                        all_fh.write(f"{z}\t{prot}\t{charge}\t{c['db2_url']}\t{member}\n")
                    if z in wanted:
                        rows.append((z, prot, charge, c["gen"], c["tranche"], c["chunk"],
                                     member, c["db2_url"]))
        except (tarfile.TarError, OSError, EOFError, zlib.error) as e:
            bad.append(f"{url}\t{e}")
    if all_fh:
        all_fh.close()
    if missing:
        print(f"WARNING: {missing:,} sidecars not on disk (or not in the catalog); "
              f"run status.sh on sidecar_urls.txt and fetch them first.", file=sys.stderr)
    if bad:
        (a.plan / "sidecars_unreadable.tsv").write_text("".join(f"{b}\n" for b in bad))
        print(f"WARNING: {len(bad):,} unreadable sidecars -> sidecars_unreadable.tsv; "
              f"verify.sh --purge and refetch.", file=sys.stderr)

    rows.sort()
    with open(a.plan / "index.tsv", "w") as fh:
        fh.write("zinc_id\tprotomer\tcharge\tgen\ttranche\tchunk\tmember\tdb2_url\n")
        fh.writelines("\t".join(r) + "\n" for r in rows)
    built = {r[0] for r in rows}
    if used_lists and frac == 1:
        (a.plan / "listed_not_built.txt").write_text(
            "".join(f"{z}\n" for z in sorted(listed - built)))

    # chunk -> wanted ids it holds; optionally greedy under a byte budget
    ids_in = defaultdict(set)
    for r in rows:
        ids_in[r[7]].add(r[0])
    size = {c["db2_url"]: int(c["db2_bytes"]) for c in cat.values()}
    chosen, covered = [], set()
    if a.budget_gb is None:
        chosen = sorted(ids_in)
        covered = set(built)
    else:
        budget, spent = a.budget_gb * 1e9, 0
        if a.greedy:
            # Lazy greedy: a chunk's gain only shrinks as others are chosen, so a popped
            # chunk whose recomputed gain still beats the next entry is the true best.
            heap = [(-len(s) / max(size[u], 1), u) for u, s in ids_in.items()]
        else:
            heap = [(rank(u, a.seed), u) for u in ids_in]      # seeded random order
        heapq.heapify(heap)
        while heap:
            _, u = heapq.heappop(heap)
            gain = ids_in[u] - covered
            if not gain:
                continue
            if a.greedy:
                score = -len(gain) / max(size[u], 1)
                if heap and score > heap[0][0]:
                    heapq.heappush(heap, (score, u))     # stale: re-queue at its real score
                    continue
            if spent + size[u] > budget:
                continue
            chosen.append(u)
            covered |= gain
            spent += size[u]
        chosen.sort()
        (a.plan / "over_budget.txt").write_text(
            "".join(f"{z}\n" for z in sorted(built - covered)))
    (a.plan / "db2_urls.txt").write_text("".join(f"{u}\n" for u in chosen))

    cov = defaultdict(Counter)
    for z, b in wanted.items():
        cov[b]["wanted"] += 1
        cov[b]["listed"] += z in listed
        cov[b]["built"] += z in built
        cov[b]["fetched"] += z in covered
    for u in chosen:
        b = cat[u[:-len(".db2.tgz")]]["bin"]
        cov[b]["chunks"] += 1
        cov[b]["bytes"] += size[u]
    with open(a.plan / "coverage.tsv", "w") as fh:
        fh.write("bin\twanted\tin_id_list\tbuilt_3d\tcovered_by_fetch\tchunks\tdb2_GB\tMB_per_id\n")
        tot = Counter()
        for b in sorted(cov) + ["all"]:
            c = tot if b == "all" else cov[b]
            if b != "all":
                tot.update(c)
            fh.write(f"{b}\t{c['wanted']}\t{c['listed']}\t{c['built']}\t{c['fetched']}\t"
                     f"{c['chunks']}\t{c['bytes'] / 1e9:.1f}\t"
                     f"{c['bytes'] / max(c['fetched'], 1) / 1e6:.0f}\n")
    print((a.plan / "coverage.tsv").read_text(), file=sys.stderr, end="")
    if not used_lists:
        print("(--no-idlists: in_id_list is not measured)", file=sys.stderr)
    if frac < 1:
        print(f"sidecars were a {frac} sample of the candidates: all of them would find "
              f"~{len(built) / frac:,.0f} built ids (~{len(built) / frac / max(len(wanted), 1):.0%} "
              f"of wanted) for ~{sum(size[u] for u in ids_in) / frac / 1e12:.1f} TB of db2",
              file=sys.stderr)
    print(f"{len(rows):,} db2 members for {len(built):,} ids -> {a.plan}/index.tsv; "
          f"{len(chosen):,} chunks -> {a.plan}/db2_urls.txt", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("groups", "chunks"):
        s = sub.add_parser(name)
        s.add_argument("--plan", type=Path, default=Path("plan"))
        s.add_argument("--catalog", type=Path, default=Path("catalog"))
        s.add_argument("--root", type=Path, default=Path("."),
                       help="the landing dir run.sh fetched into (holds zinc22/...)")
    for name in ("groups", "chunks"):
        sub.choices[name].add_argument("--seed", default="llb", help="for --fraction / --budget-gb")
    sub.choices["groups"].add_argument("--no-idlists", action="store_true")
    sub.choices["groups"].add_argument("--fraction", type=float, default=1.0)
    sub.choices["chunks"].add_argument("--budget-gb", type=float, default=None)
    sub.choices["chunks"].add_argument("--greedy", action="store_true")
    sub.choices["chunks"].add_argument("--all-members", action="store_true")
    a = ap.parse_args()
    {"groups": groups, "chunks": chunks}[a.cmd](a)


if __name__ == "__main__":
    main()
