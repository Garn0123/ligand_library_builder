#!/usr/bin/env python3
"""
sample_2d.py -- draw N SMILES per heavy-atom bin from ZINC22 2D.

Three modes
-----------
  --local ROOT   UNIFORM over the whole bin, from files already on disk (run.sh).
                 One reservoir pass across every tranche file in the bin, so
                 every molecule in the bin is equally likely and the bin's logP
                 mix comes out right without any allocation. THE production mode.
  (default)      STREAMING: reads the front of each remote tranche file and
                 stops. A few MB per bin. BIASED -- for plumbing runs only.
  --print-urls   the 2D URLs for the bins, for run.sh (the full download).

Why streaming is biased (measured 2026-09-23, H17M100, 3,600 molecules)
-----------------------------------------------------------------------
The file is a concatenation of registration batches, and its first few thousand
lines are one batch:

                           prefix   uniform   uniform (2nd draw)
    unique Murcko scaffolds  1,479     1,654     1,591
    mean NN Tanimoto         0.469     0.432     0.435
    rings / molecule         1.378     1.463     KS p 2e-5 vs 0.86 (null)
    amides / molecule        1.502     1.423     KS p 0.001 vs 1.0 (null)

Fewer scaffolds and more near-neighbours: a clustered, less diverse sample --
the wrong error for a chemical-space study. Simple string features (length,
stereo, charge) are flat across deciles and do NOT reveal it.

Cost of doing it properly (measured): one connection downloads ~14 MB/s, so
H14-H28 (~960 GB) is ~5 h at run.sh -j 4 if the server scales, ~19 h if not.
A local pass reads ~16M lines/s (gzip -dc), so H28 (~26B lines) is ~30 min.

Outputs, per bin, in --outdir:
    H17.smi            'SMILES ZINCID', the input to prepare_parents.py
    H17.manifest.tsv   per tranche file: lines, sampled (+ quota/window if streamed)
    H17.done           written last, records mode/seed/sha256; skipped on re-run

Usage:
    python3 sample_2d.py --heavy 14-28 --print-urls > urls_2d.txt      # then run.sh
    python3 sample_2d.py --heavy 14-28 --local zinc22 --urls urls_2d.txt \
        --outdir samples                                                # uniform
    python3 sample_2d.py --heavy 17 --per-bin 500 --outdir smoke       # pilot only

stdlib only: runs on a DTN, login or compute node with any python3 >= 3.8.
"""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import math
import shutil
import subprocess
import random
import re
import sys
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path

BASE = "https://files.docking.org/zinc22/"
# ZINC's own index page publishes this credential for the protected tranches
# (fetch.sh uses the same one).
AUTH = "gpcr:xtal"
ZINC_ID = re.compile(r"^ZINC[0-9A-Za-z]{12}$")
# One Apache autoindex row:
#   <td><a href="H17M000.smi.gz">...</a></td><td align="right">2025-11-04 15:56</td>
#   <td align="right"> 98M</td>
ROW = re.compile(r'href="(H\d\d[MP]\d{3}\.smi\.gz)".*?<td align="right">[^<]*</td>'
                 r'<td align="right">\s*(\d+(?:\.\d+)?)([KMG]?)\s*</td>')
TRANCHE_HREF = re.compile(r'href="H\d\d[MP]\d{3}\.smi\.gz"')
UNIT = {"": 1, "K": 1e3, "M": 1e6, "G": 1e9}
MAX_ATTEMPTS = 4


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def request(url: str, timeout: int = 120):
    req = urllib.request.Request(url)
    req.add_header("Authorization",
                   "Basic " + base64.b64encode(AUTH.encode()).decode())
    req.add_header("User-Agent", "sample_2d.py (size-ladder sampling; one connection)")
    return urllib.request.urlopen(req, timeout=timeout)


def with_retries(fn, what: str):
    """5xx / network errors back off and retry; 404/403/410 return None at once.

    Same policy as fetch.sh: a permanent error does not burn attempts, a
    transient one escalates the sleep, because bursty failures from this server
    mean it is at capacity (README, Gotcha 7).
    """
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            return fn()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404, 410):
                log(f"  {e.code} {what} -- permanent, skipped")
                return None
            err = f"HTTP {e.code}"
        except (urllib.error.URLError, TimeoutError, ConnectionError,
                EOFError, OSError) as e:
            err = f"{type(e).__name__}: {e}"
        if attempt < MAX_ATTEMPTS:
            log(f"  {err} on {what}; retry {attempt}/{MAX_ATTEMPTS - 1} "
                f"in {20 * attempt}s")
            time.sleep(20 * attempt)
    raise RuntimeError(f"gave up on {what} after {MAX_ATTEMPTS} attempts ({err})")


def list_bin(subset: str, hbin: str) -> list[tuple[str, float]]:
    """(filename, approx bytes) for every tranche file in one heavy-atom bin."""
    url = f"{BASE}{subset}/{hbin}/"
    body = with_retries(lambda: request(url).read().decode("utf-8", "replace"),
                        url)
    if body is None:
        return []
    out = [(m.group(1), float(m.group(2)) * UNIT[m.group(3)])
           for m in ROW.finditer(body)]
    # A row the regex misses is a tranche silently left out of the bin, skewing
    # its logP mix. Bytes-only sizes ('284 ', trailing space) did exactly that.
    n_href = len(TRANCHE_HREF.findall(body))
    if not out or len(out) != n_href:
        raise RuntimeError(f"parsed {len(out)} of {n_href} tranche rows from {url} "
                           f"-- listing format changed?")
    return out


def allocate(files: list[tuple[str, float]], n: int) -> dict[str, int]:
    """Quota per file proportional to size (largest remainder), summing to n."""
    total = sum(s for _, s in files)
    raw = {f: n * s / total for f, s in files}
    q = {f: int(v) for f, v in raw.items()}
    for f in sorted(raw, key=lambda f: raw[f] - q[f], reverse=True)[:n - sum(q.values())]:
        q[f] += 1
    return q


def stream_window(url: str, window: int) -> tuple[list[str], bool, int]:
    """First `window` valid lines of a remote .smi.gz; (lines, hit_eof, n_malformed)."""
    def go():
        lines, bad = [], 0
        with request(url, timeout=300) as resp, gzip.GzipFile(fileobj=resp) as gz:
            for raw in gz:
                parts = raw.decode("ascii", "replace").split()
                if len(parts) < 2 or not ZINC_ID.match(parts[1]):
                    bad += 1
                    continue
                lines.append(f"{parts[0]} {parts[1]}")
                if len(lines) >= window:
                    return lines, False, bad
        return lines, True, bad
    return with_retries(go, url)


def sample_bin(subset: str, hbin: str, per_bin: int, oversample: int,
               seed: int, outdir: Path) -> None:
    done = outdir / f"{hbin}.done"
    if done.exists():
        log(f"{hbin}: already done ({done}); --force to redo")
        return
    files = list_bin(subset, hbin)
    quota = allocate(files, per_bin)
    log(f"{hbin}: {len(files)} tranche files, {sum(s for _, s in files) / 1e9:.2f} GB "
        f"on the server, quota {per_bin}")

    # Smallest files first, so a file that runs out before its quota passes the
    # deficit on to larger ones, and the largest files absorb whatever remains.
    chosen, rows, carry = [], [], 0
    for fname, size in sorted(files, key=lambda x: x[1]):
        want = quota[fname] + carry
        if want == 0:
            rows.append((fname, size, 0, 0, "", 0, 0))
            continue
        url = f"{BASE}{subset}/{hbin}/{fname}"
        got = stream_window(url, want * oversample)
        if got is None:                          # permanent 404: pass the quota on
            carry = want
            rows.append((fname, size, want, 0, "404", 0, 0))
            continue
        lines, eof, bad = got
        rng = random.Random(f"{seed}:{subset}:{fname}")
        pick = rng.sample(lines, min(want, len(lines)))
        chosen += pick
        carry = want - len(pick)
        rows.append((fname, size, want, len(lines), "eof" if eof else "", len(pick), bad))
        time.sleep(0.3)                          # same courtesy delay as run.sh

    if carry:
        log(f"  ! {hbin}: short by {carry} -- the whole bin holds fewer molecules "
            f"than --per-bin x ... read, or files 404'd. Recorded, not hidden.")

    tmp = outdir / f"{hbin}.smi.tmp"
    tmp.write_text("".join(f"{l}\n" for l in chosen))
    tmp.rename(outdir / f"{hbin}.smi")
    with open(outdir / f"{hbin}.manifest.tsv", "w") as fh:
        fh.write("file\tserver_bytes\tquota\tlines_read\tstatus\tsampled\tmalformed\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in (r[0], int(r[1])) + r[2:]) + "\n")
    digest = hashlib.sha256((outdir / f"{hbin}.smi").read_bytes()).hexdigest()
    done.write_text(f"subset={subset}\nper_bin={per_bin}\nsampled={len(chosen)}\n"
                    f"oversample={oversample}\nseed={seed}\nsha256={digest}\n"
                    f"date={time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    log(f"  -> {hbin}.smi  {len(chosen)} molecules")


def sample_local(root: Path, subset: str, hbin: str, per_bin: int, seed: int,
                 outdir: Path, urls: Path | None) -> None:
    """Uniform sample of `per_bin` lines over every tranche file of a bin on disk.

    Reservoir sampling, Algorithm L (Li 1994): one pass, O(k) memory, and the
    random draws are only made at the (rare) replacement points, so the loop is
    a line counter. Lines are parsed only if selected.
    """
    done = outdir / f"{hbin}.done"
    if done.exists():
        log(f"{hbin}: already done ({done}); --force to redo")
        return
    bindir = root / subset / hbin
    files = sorted(bindir.glob(f"{hbin}[MP][0-9][0-9][0-9].smi.gz"))
    if not files:
        raise SystemExit(f"no {hbin} tranche files under {bindir} -- wrong --local root? "
                         f"run.sh writes zinc22/{subset}/{hbin}/...")
    # A missing tranche file is a hole in the bin's logP mix that no error
    # would ever report, so with --urls every expected file must be present.
    if urls is not None:
        want = {u.rsplit("/", 1)[1] for u in urls.read_text().split()
                if f"/{subset}/{hbin}/" in u}
        permanent = urls.parent / "permanent.tsv"
        gone = ({l.split()[1].rsplit("/", 1)[1] for l in permanent.read_text().splitlines()
                 if l.strip()} if permanent.exists() else set())
        absent = sorted(want - {f.name for f in files} - gone)
        if absent:
            raise SystemExit(f"{hbin}: {len(absent)} expected file(s) not on disk, e.g. "
                             f"{absent[:3]}. Finish run.sh / status.sh first.")
    else:
        log(f"  ! {hbin}: no --urls, so bin completeness is NOT checked")
    unzip = [shutil.which("pigz") or "gzip", "-dc"]
    rng = random.Random(f"{seed}:{subset}:{hbin}:local")
    k, i, nxt, w = per_bin, 0, -1, 0.0
    res: list[tuple[bytes, str]] = []
    rows = []
    for f in files:
        start = i
        proc = subprocess.Popen(unzip + [str(f)], stdout=subprocess.PIPE, bufsize=1 << 20)
        for line in proc.stdout:
            if i < k:
                res.append((line, f.name))
                if i == k - 1:
                    w = math.exp(math.log(rng.random()) / k)
                    nxt = i + int(math.log(rng.random()) / math.log(1 - w)) + 1
            elif i == nxt:
                res[rng.randrange(k)] = (line, f.name)
                w *= math.exp(math.log(rng.random()) / k)
                nxt += int(math.log(rng.random()) / math.log(1 - w)) + 1
            i += 1
        if proc.wait() != 0:
            raise SystemExit(f"{f}: decompression failed (exit {proc.returncode}). "
                             f"Run verify.sh --purge and re-fetch before sampling.")
        rows.append((f.name, i - start))
        log(f"  {f.name}: {i - start:,} lines")

    chosen, bad = [], 0
    for raw, _ in res:
        parts = raw.decode("ascii", "replace").split()
        if len(parts) < 2 or not ZINC_ID.match(parts[1]):
            bad += 1
            continue
        chosen.append(f"{parts[0]} {parts[1]}")
    per_file = Counter(fn for _, fn in res)

    tmp = outdir / f"{hbin}.smi.tmp"
    tmp.write_text("".join(f"{l}\n" for l in chosen))
    tmp.rename(outdir / f"{hbin}.smi")
    with open(outdir / f"{hbin}.manifest.tsv", "w") as fh:
        fh.write("file\tlines\tsampled\n")
        fh.writelines(f"{n}\t{c}\t{per_file.get(n, 0)}\n" for n, c in rows)
    digest = hashlib.sha256((outdir / f"{hbin}.smi").read_bytes()).hexdigest()
    done.write_text(f"mode=local\nroot={root.resolve()}\nsubset={subset}\n"
                    f"per_bin={per_bin}\nsampled={len(chosen)}\npopulation={i}\n"
                    f"malformed_selected={bad}\nn_files={len(files)}\nseed={seed}\n"
                    f"sha256={digest}\ndate={time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    log(f"  -> {hbin}.smi  {len(chosen):,} of {i:,} (uniform)")


def parse_heavy(spec: str) -> list[str]:
    out = []
    for part in spec.split(","):
        a, _, b = part.partition("-")
        out += [f"H{h:02d}" for h in range(int(a), int(b or a) + 1)]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--heavy", required=True, help="bins, e.g. '14-28' or '14,17,20'")
    ap.add_argument("--per-bin", type=int, default=20000)
    ap.add_argument("--oversample", type=int, default=10,
                    help="read this many times the quota from each file, sample from that")
    ap.add_argument("--seed", type=int, default=20260923)
    ap.add_argument("--subset", default="2d", help="'2d' (2d-all was identical for H17)")
    ap.add_argument("--outdir", default="samples")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--local", default=None, metavar="ROOT",
                    help="sample uniformly from files on disk under ROOT/<subset>/Hxx/ "
                         "(the zinc22/ directory run.sh writes). Production mode.")
    ap.add_argument("--urls", default=None,
                    help="with --local: the URL list run.sh fetched; every file in it "
                         "(minus permanent.tsv 404s beside it) must be on disk")
    ap.add_argument("--print-urls", action="store_true",
                    help="print every 2D tranche URL for the bins (for run.sh) and exit")
    args = ap.parse_args(argv)

    if args.print_urls:
        for hbin in parse_heavy(args.heavy):
            for fname, _ in list_bin(args.subset, hbin):
                print(f"{BASE}{args.subset}/{hbin}/{fname}")
        return 0

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for hbin in parse_heavy(args.heavy):
        if args.force:
            (outdir / f"{hbin}.done").unlink(missing_ok=True)
        if args.local:
            sample_local(Path(args.local), args.subset, hbin, args.per_bin, args.seed,
                         outdir, Path(args.urls) if args.urls else None)
        else:
            log(f"{hbin}: STREAMING mode -- a biased prefix sample, for plumbing runs "
                f"only (see docstring). Use --local for the real sample.")
            sample_bin(args.subset, hbin, args.per_bin, args.oversample, args.seed, outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
