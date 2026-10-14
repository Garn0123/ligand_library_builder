#!/usr/bin/env python3
"""
clash_report.py -- why are a molecule's db2 conformer sets flagged broken?

    llb env DOCK sh -c '"${DB2TOOL_EXE:-db2tool}" tomol2 --all --keep-broken -o sets.mol2 db2_run/out/0003/NAME.db2.gz'
    llb clash-report sets.mol2 -o clash.tsv

db2_converter (63d6656, mol2db2/clash.py) flags a conformer set broken by one
rule: a ROTATABLE hydrogen within 1.70 A of another hydrogen more than two
bonds away. Rotatable hydrogens are set by mol2db2/hydrogens.py: the H of an
O-H or S-H on carbon, and some N-H. This reads every set `db2tool tomol2 --all`
wrote and, for each, finds the closest such H...H contact, so you can see
whether "broken" means a near miss (1.6 A) or a collapsed conformer (1.0 A),
and which groups are involved.

Here "rotatable" is approximated as: H bonded to O or S, and H bonded to N
(N-H is reported separately, since db2_converter rotates only some of them).

Writes one row per set (name, conf, broken flag, closest contact, the two
groups) and prints a per-molecule summary. stdlib only.
"""
from __future__ import annotations

import argparse
import re
import statistics
import sys
from collections import defaultdict, deque

HEADER = re.compile(r"^##########\s+([A-Za-z0-9_]+):\s*(.*?)\s*$")
CUTOFF = 1.70


def records(path):
    """Yield (header dict, atoms, bonds) per mol2 record."""
    header, atoms, bonds, section = {}, [], [], None
    with open(path) as fh:
        for line in fh:
            m = HEADER.match(line)
            if m:
                if atoms:
                    yield header, atoms, bonds
                    header, atoms, bonds = {}, [], []
                header[m.group(1)] = m.group(2)
                continue
            if line.startswith("@<TRIPOS>"):
                tag = line.strip()[9:]
                if tag == "MOLECULE" and atoms:
                    yield header, atoms, bonds
                    header, atoms, bonds = {}, [], []
                section = tag
                continue
            f = line.split()
            if not f:
                continue
            if section == "ATOM":
                atoms.append((f[1], float(f[2]), float(f[3]), float(f[4]), f[5]))
            elif section == "BOND":
                bonds.append((int(f[1]) - 1, int(f[2]) - 1))
    if atoms:
        yield header, atoms, bonds


def element(t):
    return t.split(".")[0]


def closest_contact(atoms, bonds):
    """Closest H...H over rotatable-type H, > 2 bonds apart (bond path >= 3)."""
    n = len(atoms)
    nbr = defaultdict(list)
    for a, b in bonds:
        nbr[a].append(b)
        nbr[b].append(a)
    heavy_of = {i: nbr[i][0] for i in range(n) if element(atoms[i][4]) == "H" and nbr[i]}
    rot = [i for i, h in heavy_of.items() if element(atoms[h][4]) in ("O", "S", "N")]
    hs = list(heavy_of)
    best = None
    for i in rot:
        # bond distances from i, up to 2
        dist = {i: 0}
        q = deque([i])
        while q:
            u = q.popleft()
            if dist[u] >= 2:
                continue
            for v in nbr[u]:
                if v not in dist:
                    dist[v] = dist[u] + 1
                    q.append(v)
        xi = atoms[i][1:4]
        for j in hs:
            if j == i or j in dist:
                continue
            xj = atoms[j][1:4]
            d = sum((a - b) ** 2 for a, b in zip(xi, xj)) ** 0.5
            if best is None or d < best[0]:
                best = (d, i, j)
    if best is None:
        return None
    d, i, j = best
    group = lambda h: f"{element(atoms[heavy_of[h]][4])}-H"   # noqa: E731
    return d, group(i), group(j), atoms[i][0], atoms[j][0]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mol2", nargs="+", help="output of db2tool tomol2 --all --keep-broken")
    ap.add_argument("-o", "--output", default="clash.tsv")
    args = ap.parse_args(argv)

    per_mol = defaultdict(list)
    with open(args.output, "w") as out:
        out.write("name\tconf\tbroken\tclosest_HH\trotatable_H_on\tother_H_on\tatoms\n")
        for path in args.mol2:
            for h, atoms, bonds in records(path):
                name = h.get("Name", "?")
                c = closest_contact(atoms, bonds)
                broken = h.get("HDB_Set_Broken", "") == "yes"
                if c is None:
                    out.write(f"{name}\t{h.get('HDB_Conf_Number','')}\t{int(broken)}\t\t\t\t\n")
                    per_mol[name].append((None, broken, ""))
                    continue
                d, gi, gj, ai, aj = c
                out.write(f"{name}\t{h.get('HDB_Conf_Number','')}\t{int(broken)}\t{d:.3f}"
                          f"\t{gi}\t{gj}\t{ai}-{aj}\n")
                per_mol[name].append((d, broken, f"{gi}...{gj}"))

    print(f"{'name':<18} {'sets':>5} {'broken':>6} {'<1.70':>6} {'min':>6} {'median':>7}  closest pair (most common)")
    for name, rows in sorted(per_mol.items()):
        ds = [d for d, _, _ in rows if d is not None]
        pairs = defaultdict(int)
        for d, _, p in rows:
            if d is not None and d < CUTOFF:
                pairs[p] += 1
        top = max(pairs, key=pairs.get) if pairs else "-"
        print(f"{name:<18} {len(rows):>5} {sum(b for _, b, _ in rows):>6} "
              f"{sum(1 for d in ds if d < CUTOFF):>6} "
              f"{(min(ds) if ds else float('nan')):>6.2f} "
              f"{(statistics.median(ds) if ds else float('nan')):>7.2f}  {top}")
    print(f"-> {args.output}  (closest contact per set; rotatable H = H on O, S or N)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
