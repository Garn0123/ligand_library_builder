#!/usr/bin/env python3
"""
plan.py -- which ZINC22 3D files a set of ZINC ids needs, and what they will cost.

Step 1 of the targeted fetch (RUNBOOK.md). Input: catalog/ from catalog.py and a
list of ZINC ids (an LLB parents.tsv or library.tsv -- the zinc_id column -- or
one id per line). Output, in --out (default plan/):

    wanted.tsv        zinc_id  bin  tranche          the ids this plan is for
    idlist_urls.txt   per-(tranche, charge, generation) id lists of every tranche
                      the ids fall in -- the cheap first download (~17 bytes/id)
    estimate.tsv      per bin: ids, tranches, id-list GB, and the db2 bill three ways:
                        all     every chunk in a touched tranche (upper bound)
                        expect  expected bytes of the chunks that actually hold an id
                        per_id  expect / expected ids with a 3D build
    rejected.txt      names that are not ZINC22 ids (hashed 'Z...' names, SMILES...)

The expectation needs the share of ids with a ZINC 3D build at all (LLB samples
the 2D set, which is larger). --p3d auto estimates it per tranche from the 2D
tranche sizes (one listing per bin of zinc22/2d/); --p3d 0.4 fixes it. It is only
an estimate: index.py groups measures it exactly.

Usage:
    python3 plan.py --catalog catalog --ids parents/parents.tsv --out plan
    python3 plan.py --catalog catalog --ids ids.txt --out plan --p3d 0.5
"""
from __future__ import annotations

import argparse
import base64
import re
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from zid import bin_of, read_ids, read_tsv, tranche_matches, tranche_of  # noqa: E402

BYTES_PER_2D_LINE = 11.7   # zinc22/2d .smi.gz, measured on H28 (sample_2d.py: 304 GB, ~26B lines)
TWOD_ROW = re.compile(r'href="(H\d\d[MP]\d{3})\.smi\.gz".*?<td align="right">[^<]*</td>\s*'
                      r'<td align="right">\s*([\d.]+)([KMGT]?)\s*</td>')


def twod_ids(bins: set[str]) -> dict[str, float]:
    """Approximate 2D ids per tranche from zinc22/2d/<bin>/ listing sizes."""
    out = {}
    for b in sorted(bins):
        req = urllib.request.Request(f"https://files.docking.org/zinc22/2d/{b}/")
        req.add_header("Authorization", "Basic " + base64.b64encode(b"gpcr:xtal").decode())
        req.add_header("User-Agent", "llb zinc3d/plan.py (one listing per bin)")
        page = urllib.request.urlopen(req, timeout=120).read().decode("utf-8", "replace")
        mult = {"": 1, "K": 1024, "M": 1024 ** 2, "G": 1024 ** 3, "T": 1024 ** 4}
        for t, num, unit in TWOD_ROW.findall(page):
            out[t] = float(num) * mult[unit] / BYTES_PER_2D_LINE
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--ids", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("plan"))
    ap.add_argument("--p3d", default="auto",
                    help="share of ids with a 3D build: 'auto' (from 2D sizes) or a number")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    ids, rejected = read_ids(a.ids)
    (a.out / "rejected.txt").write_text("".join(f"{x}\n" for x in rejected))
    chunks = read_tsv(a.catalog / "chunks.tsv")
    idlists = read_tsv(a.catalog / "idlists.tsv")
    cat_bins = {c["bin"] for c in chunks}
    tranches_by_bin = defaultdict(set)
    for c in chunks + idlists:
        tranches_by_bin[c["bin"]].add(c["tranche"])

    # each id -> its tranche(s); ids outside the catalogued bins are reported, not planned
    wanted, outside = [], Counter()
    per_tranche = Counter()
    for z in ids:
        b = bin_of(z)
        if b not in cat_bins:
            outside[b] += 1
            continue
        t = tranche_of(z)
        ts = [t] if t else sorted(x for x in tranches_by_bin[b] if tranche_matches(z, x))
        wanted.append((z, b, t or "P>500"))
        for x in ts:
            per_tranche[x] += 1 / len(ts)
    with open(a.out / "wanted.tsv", "w") as fh:
        fh.write("zinc_id\tbin\ttranche\n")
        fh.writelines(f"{z}\t{b}\t{t}\n" for z, b, t in wanted)

    touched = set(per_tranche)
    urls = sorted(r["url"] for r in idlists if r["tranche"] in touched)
    (a.out / "idlist_urls.txt").write_text("".join(f"{u}\n" for u in urls))

    ids3d = Counter()              # tranche -> molecules in ZINC 3D (id lists, ~5% high)
    for r in idlists:
        ids3d[r["tranche"]] += int(r["approx_ids"])
    if a.p3d == "auto":
        twod = twod_ids({t[:3] for t in touched})
        p3d = {t: min(1.0, ids3d[t] / twod[t]) if twod.get(t) else 0.0 for t in touched}
    else:
        p3d = {t: float(a.p3d) for t in touched}

    by_tranche = defaultdict(list)
    for c in chunks:
        if c["tranche"] in touched:
            side = int(c["pdbqt_bytes"]) or int(c["solv_bytes"]) or 1
            by_tranche[c["tranche"]].append((int(c["db2_bytes"]), side))

    est = defaultdict(lambda: Counter())
    for r in idlists:
        if r["tranche"] in touched:
            est[r["bin"]]["idlist_bytes"] += int(r["bytes"])
    for t, k in per_tranche.items():
        b = t[:3]
        e = est[b]
        e["ids"] += k
        e["tranches"] += 1
        e["ids_3d"] += k * p3d[t]
        cs = by_tranche.get(t, [])
        tot_side = sum(s for _, s in cs) or 1
        k3 = k * p3d[t]
        for db2, side in cs:
            e["all"] += db2
            e["sidecar_all"] += side
            # P(chunk holds >= 1 of the k3 ids), each landing in proportion to molecules
            e["expect"] += db2 * (1 - (1 - side / tot_side) ** k3)
    with open(a.out / "estimate.tsv", "w") as fh:
        fh.write("bin\tids\ttranches\tids_with_3d\tidlist_GB\tsidecar_all_GB\t"
                 "db2_all_TB\tdb2_expect_TB\tdb2_MB_per_id\n")
        tot = Counter()
        for b in sorted(est):
            e = est[b]
            tot.update(e)
            fh.write(f"{b}\t{e['ids']:.0f}\t{e['tranches']}\t{e['ids_3d']:.0f}\t"
                     f"{e['idlist_bytes'] / 1e9:.2f}\t{e['sidecar_all'] / 1e9:.1f}\t"
                     f"{e['all'] / 1e12:.2f}\t{e['expect'] / 1e12:.2f}\t"
                     f"{e['expect'] / max(e['ids_3d'], 1) / 1e6:.0f}\n")
        fh.write(f"all\t{tot['ids']:.0f}\t{tot['tranches']}\t{tot['ids_3d']:.0f}\t"
                 f"{tot['idlist_bytes'] / 1e9:.2f}\t{tot['sidecar_all'] / 1e9:.1f}\t"
                 f"{tot['all'] / 1e12:.2f}\t{tot['expect'] / 1e12:.2f}\t"
                 f"{tot['expect'] / max(tot['ids_3d'], 1) / 1e6:.0f}\n")

    print(f"{len(ids):,} ZINC ids ({len(rejected):,} other names -> rejected.txt); "
          f"{len(wanted):,} in catalogued bins", file=sys.stderr)
    for b, n in sorted(outside.items()):
        print(f"  {n:,} ids in {b}, which the catalog does not cover", file=sys.stderr)
    print(f"{len(touched)} tranches touched; {len(urls):,} id lists -> {a.out}/idlist_urls.txt",
          file=sys.stderr)
    print((a.out / "estimate.tsv").read_text(), file=sys.stderr, end="")


if __name__ == "__main__":
    main()
