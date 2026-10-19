#!/usr/bin/env python3
"""
catalog.py -- list every ZINC22 3D chunk on files.docking.org for some heavy-atom bins.

Reads Apache directory listings only; never downloads a data file. One row per
chunk (one tranche x charge x chunk id), with the server's sizes for the .db2.tgz
and its small sidecars, plus one row per per-tranche id list (.txt).

Layout (2026-10, varies by generation, so the walk does not assume one):

    zinc22/zinc-22<g>/H20/H20P200/H20P200-N.<g>.txt         ids of (tranche, charge, generation)
    zinc22/zinc-22<g>/H20/H20P200/<s>/H20P200-N-<g>aaaba.db2.tgz   one chunk
                                     .../...-N-<g>aaaba.pdbqt.tgz  same members, ~1% the bytes
                                     .../...-N-<g>aaaba.solv.tgz   same members

  <g>  generation a..z ("zinc-22f").  <s>  a chunk-group directory (a/, b/, ...).
  N    net charge: J=-4 K=-3 L=-2 M=-1 N=0 O=+1 P=+2 Q=+3.

Every member of a chunk is H20/H20P200/xx/yy/ZINCk...<id>.<n>.<charge>.<ext>, and
the db2, pdbqt and solv tarballs of one chunk hold the same members. That is what
makes index.py cheap: list the sidecar, not the db2.

Sizes are Apache's rounded human sizes ("276M"): good to ~1%, not exact.

Outputs (in --out, default catalog/):
    chunks.tsv   gen bin tranche charge chunk db2_url db2_bytes pdbqt_bytes solv_bytes
    idlists.tsv  gen bin tranche charge url bytes approx_ids   (ids = bytes / 17)
    dirs.done    listings already parsed (resume: re-running skips them)
    summary.tsv  per gen x bin: chunks, db2 bytes, sidecar bytes, approx ids

Usage:
    python3 catalog.py --heavy 18-30 --out catalog            # ~10-20k listings, ~1-2 h
    python3 catalog.py --heavy 18-30 --out catalog --summary  # re-summarise only

stdlib only, python3 >= 3.8. Two listing requests in flight at most.
"""
from __future__ import annotations

import argparse
import base64
import re
import sys
import threading
import time
import urllib.error
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BASE = "https://files.docking.org/zinc22/"
AUTH = "gpcr:xtal"          # published on ZINC's own index page (fetch.sh uses it too)
MAX_ATTEMPTS = 4
GENERATIONS = "abcdefghijklmnopqrstuvwxyz"

# <td><a href="NAME">...</a>   </td><td align="right">2024-08-24 23:41  </td><td align="right">276M</td>
ROW = re.compile(r'<a href="([^"?/][^"]*)">[^<]*</a>\s*</td>\s*<td align="right">[^<]*</td>\s*'
                 r'<td align="right">\s*([^<]*?)\s*</td>')
CHUNK = re.compile(r"^(H(\d\d)[MP]\d{3})-([A-Z])-([a-z]{6})\.(db2|pdbqt|solv)\.tgz$")
IDLIST = re.compile(r"^(H(\d\d)[MP]\d{3})-([A-Z])\.([a-z])\.txt$")
TRANCHE_DIR = re.compile(r"^H\d\d[MP]\d{3}/$")
GROUP_DIR = re.compile(r"^[a-z]/$")


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def human_bytes(s: str) -> int:
    s = s.strip()
    if s in ("", "-"):
        return 0
    mult = {"K": 1024, "M": 1024 ** 2, "G": 1024 ** 3, "T": 1024 ** 4}.get(s[-1], 1)
    return int(float(s.rstrip("KMGT")) * mult)


def fetch(url: str, delay: float) -> str:
    req = urllib.request.Request(url)
    req.add_header("Authorization", "Basic " + base64.b64encode(AUTH.encode()).decode())
    # Python's default User-Agent gets 403 from files.docking.org.
    req.add_header("User-Agent", "llb zinc3d/catalog.py (directory listings; 2 connections)")
    err = ""
    for attempt in range(1, MAX_ATTEMPTS + 1):
        time.sleep(delay)
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code in (404, 410):
                log(f"  {e.code} {url} -- skipped")
                return ""
            err = f"HTTP {e.code}"
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            err = f"{type(e).__name__}: {e}"
        if attempt < MAX_ATTEMPTS:
            log(f"  {err} on {url}; retry {attempt} in {20 * attempt}s")
            time.sleep(20 * attempt)
    raise RuntimeError(f"gave up on {url} ({err})")


def listing(url: str, delay: float) -> list[tuple[str, int]]:
    return [(name, human_bytes(size)) for name, size in ROW.findall(fetch(url, delay))]


def parse_bins(spec: str) -> list[str]:
    lo, _, hi = spec.partition("-")
    return [f"H{h:02d}" for h in range(int(lo), int(hi or lo) + 1)]


class Catalog:
    def __init__(self, out: Path, delay: float):
        self.out = out
        self.delay = delay
        self.lock = threading.Lock()
        out.mkdir(parents=True, exist_ok=True)
        self.done = set()
        if (out / "dirs.done").exists():
            self.done = set((out / "dirs.done").read_text().split())
        self.chunks_fh = open(out / "chunks.raw", "a")
        self.ids_fh = open(out / "idlists.raw", "a")
        self.done_fh = open(out / "dirs.done", "a")

    def walk_tranche(self, gen: str, hbin: str, tranche_url: str) -> None:
        """A tranche dir and its chunk-group subdirs. Recorded atomically per tranche."""
        if tranche_url in self.done:
            return
        chunk_rows, id_rows = [], []
        todo = [tranche_url]
        while todo:
            url = todo.pop()
            for name, size in listing(url, self.delay):
                if GROUP_DIR.match(name):
                    todo.append(url + name)
                elif m := CHUNK.match(name):
                    chunk_rows.append((gen, hbin, m[1], m[3], m[4], m[5], url + name, size))
                elif m := IDLIST.match(name):
                    id_rows.append((gen, hbin, m[1], m[3], url + name, size))
        with self.lock:
            for r in chunk_rows:
                self.chunks_fh.write("\t".join(map(str, r)) + "\n")
            for r in id_rows:
                self.ids_fh.write("\t".join(map(str, r)) + "\n")
            self.chunks_fh.flush(); self.ids_fh.flush()
            self.done_fh.write(tranche_url + "\n"); self.done_fh.flush()

    def run(self, bins: list[str], jobs: int) -> None:
        tranches = []
        for g in GENERATIONS:
            top = f"{BASE}zinc-22{g}/"
            have = {n.rstrip("/") for n, _ in listing(top, self.delay)}
            for hbin in bins:
                if hbin not in have:
                    continue
                url = f"{top}{hbin}/"
                ts = [url + n for n, _ in listing(url, self.delay) if TRANCHE_DIR.match(n)]
                log(f"zinc-22{g} {hbin}: {len(ts)} tranche dirs")
                tranches += [(g, hbin, t) for t in ts]
        todo = [t for t in tranches if t[2] not in self.done]
        log(f"{len(tranches)} tranche dirs, {len(todo)} still to list")
        n = 0
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            for _ in ex.map(lambda t: self.walk_tranche(*t), todo):
                n += 1
                if n % 200 == 0:
                    log(f"  {n}/{len(todo)} tranche dirs listed")


def finish(out: Path) -> None:
    """chunks.raw (one row per file) -> chunks.tsv (one row per chunk) + summary."""
    chunks: dict[tuple, dict] = {}
    for line in open(out / "chunks.raw"):
        gen, hbin, tranche, charge, chunk, kind, url, size = line.rstrip("\n").split("\t")
        key = (gen, hbin, tranche, charge, chunk)
        c = chunks.setdefault(key, {"db2_url": "", "db2": 0, "pdbqt": 0, "solv": 0})
        c[kind] = int(size)
        if kind == "db2":
            c["db2_url"] = url
    with open(out / "chunks.tsv", "w") as fh:
        fh.write("gen\tbin\ttranche\tcharge\tchunk\tdb2_url\tdb2_bytes\tpdbqt_bytes\tsolv_bytes\n")
        for key in sorted(chunks):
            c = chunks[key]
            if not c["db2_url"]:
                continue        # sidecar without a db2: nothing to fetch
            fh.write("\t".join(key) + f"\t{c['db2_url']}\t{c['db2']}\t{c['pdbqt']}\t{c['solv']}\n")

    ids: dict[str, tuple] = {}
    for line in open(out / "idlists.raw"):
        gen, hbin, tranche, charge, url, size = line.rstrip("\n").split("\t")
        ids[url] = (gen, hbin, tranche, charge, url, int(size))
    with open(out / "idlists.tsv", "w") as fh:
        fh.write("gen\tbin\ttranche\tcharge\turl\tbytes\tapprox_ids\n")
        for r in sorted(ids.values()):
            fh.write("\t".join(map(str, r)) + f"\t{r[5] // 17}\n")

    agg = defaultdict(lambda: [0, 0, 0, 0])
    for (gen, hbin, *_), c in chunks.items():
        if c["db2_url"]:
            a = agg[(gen, hbin)]
            a[0] += 1; a[1] += c["db2"]; a[2] += c["pdbqt"] or c["solv"]
    for gen, hbin, *_, size in ids.values():
        agg[(gen, hbin)][3] += size // 17
    with open(out / "summary.tsv", "w") as fh:
        fh.write("gen\tbin\tchunks\tdb2_GB\tsidecar_GB\tapprox_ids_2d\n")
        for (gen, hbin), (n, b, s, i) in sorted(agg.items()):
            fh.write(f"{gen}\t{hbin}\t{n}\t{b / 1e9:.1f}\t{s / 1e9:.2f}\t{i}\n")
    tot = [sum(a[k] for a in agg.values()) for k in range(4)]
    log(f"{tot[0]:,} chunks, {tot[1] / 1e12:.2f} TB db2, {tot[2] / 1e9:.1f} GB sidecars, "
        f"~{tot[3]:,} ids in the per-tranche lists -> {out}/summary.tsv")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--heavy", default="18-30", help="heavy-atom bins, e.g. 18-30 or 20")
    ap.add_argument("--out", type=Path, default=Path("catalog"))
    ap.add_argument("-j", "--jobs", type=int, default=2,
                    help="listing requests in flight (default 2; keep it small)")
    ap.add_argument("--delay", type=float, default=0.2, help="seconds before each request")
    ap.add_argument("--summary", action="store_true", help="only rebuild chunks.tsv/summary.tsv")
    a = ap.parse_args()
    if not a.summary:
        Catalog(a.out, a.delay).run(parse_bins(a.heavy), a.jobs)
    finish(a.out)


if __name__ == "__main__":
    main()
