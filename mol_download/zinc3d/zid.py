"""zid.py -- ZINC22 id and path helpers shared by plan.py, index.py and extract.py.

A ZINC22 id is 'ZINC' + 12 base62 characters (0-9, a-z, A-Z):

    ZINCk5000003UMKO
        ^            heavy-atom count:  k = 20  ->  bin H20
         ^           logP bucket:       0-5 = M500..M000, 6 + n = P(10n) up to U = P500
          ^^^^^^^^^^ serial within the tranche

So an id names its tranche (H20M000) but not its generation (zinc-22a..z), its
charge or its chunk -- those come from the per-tranche id lists and the sidecars.
Checked 2026-10 against zinc-22a H20: M100=4, M000=5, P000=6, P010=7, P100=g,
P200=q, P300=A, P400=K, P500=U. Codes past U (P600, P900 directories exist) are
not decoded; ids carrying them match every tranche of their bin above P500.

Archive members are  H20/H20M000/0L/mS/ZINCk50000040LmS.0.L.db2  (id, protomer
index, charge letter), identical across the .db2.tgz, .pdbqt.tgz and .solv.tgz of
one chunk. The db2 M line holds only the LAST 16 characters of 'ZINC...id.n'
(here 'NCk50000040LmS.0'), so docking output names are truncated: join on the
member name (records.tsv), never on the M-line name.
"""
from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path

B62 = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
ZINC_ID = re.compile(r"^ZINC[0-9A-Za-z]{12}$")
MEMBER = re.compile(r"(ZINC[0-9A-Za-z]{12})\.(\d+)\.([A-Z])\.(db2|pdbqt|solv)(?:\.gz)?$")
BASE = "https://files.docking.org/"


def rank(key: str, seed: str = "llb") -> float:
    """Stable pseudo-random number in [0, 1) per key (a URL). Used for every
    --fraction and for --budget-gb's random order: the same items come first on
    every rerun, and raising a fraction keeps everything already selected."""
    h = hashlib.sha1(f"{seed}:{key}".encode()).digest()
    return int.from_bytes(h[:8], "big") / 2 ** 64


def bin_of(zid: str) -> str:
    return f"H{B62.index(zid[4]):02d}"


def tranche_of(zid: str) -> str | None:
    """'H20M000' for ZINCk5..., or None when the logP code is past P500."""
    n = B62.index(zid[5])
    if n <= 5:
        return f"{bin_of(zid)}M{(5 - n) * 100:03d}"
    if n <= 56:
        return f"{bin_of(zid)}P{(n - 6) * 10:03d}"
    return None


def tranche_matches(zid: str, tranche: str) -> bool:
    t = tranche_of(zid)
    if t is not None:
        return t == tranche
    return tranche.startswith(bin_of(zid) + "P") and int(tranche[4:]) > 500


def local_path(url: str, root: Path | str = ".") -> Path:
    """Where run.sh/fetch.sh (wget -x -nH) put a URL, relative to the landing dir."""
    return Path(root) / url[len(BASE):] if url.startswith(BASE) else Path(root) / url


def read_ids(path: Path) -> tuple[list[str], list[str]]:
    """ZINC ids out of an LLB parents.tsv / library.tsv (zinc_id or parent_id column),
    or a plain list (first whitespace field per line). Returns (zinc_ids, rejected)."""
    ids, rejected = [], []
    with open(path) as fh:
        head = fh.readline()
        fh.seek(0)
        cols = head.rstrip("\n").split("\t")
        col = next((c for c in ("zinc_id", "parent_id") if c in cols), None)
        if col:
            rows = (r.get(col) or "" for r in csv.DictReader(fh, delimiter="\t"))
        else:
            rows = (line.split()[0] if line.split() else "" for line in fh)
        for x in rows:
            x = x.strip()
            if not x:
                continue
            (ids if ZINC_ID.match(x) else rejected).append(x)
    return sorted(set(ids)), rejected


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
