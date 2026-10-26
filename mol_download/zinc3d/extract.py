#!/usr/bin/env python3
"""
extract.py -- pull the wanted ids' db2 records out of fetched ZINC22 chunks.

Streams each .db2.tgz in plan/db2_urls.txt (never untars) and keeps only the
members index.py put in plan/index.tsv. Per chunk it writes one gzipped db2 that
DOCK reads as-is, holding just those records:

    OUT/db2/<bin>/<chunk name>.db2.gz       e.g. H20/H20P200-N-faaaba.db2.gz
    OUT/records/<chunk name>.tsv            one row per archive member:
        zinc_id protomer charge m_name n_records out_file offset bytes archive member sha1

offset/bytes locate the member inside the UNcompressed out_file, so one member is
    gzip.open(out_file).seek(offset); read(bytes)

A member (one protomer of one id) can hold several db2 records -- 2-3 is common
-- each docked separately by DOCK under the same name; n_records counts them
('E' lines). Score per id = min over its records, as for LLB's db2s.

m_name is the name DOCK will report: the db2 M line keeps only the last 16
characters of 'ZINC<id>.<protomer>' ('NCk50000040LmS.0'), so join docking
output back to zinc_id through this column, not by string prefix.

Sharded for a SLURM array (--shard i --nshards n, archives kept whole). Each
chunk is written to a .tmp and renamed, then marked done, so a re-run skips it.
--delete-archive removes a chunk's .db2.tgz once its records are written and
every member index.tsv expected was found (the download is ~20-100x the
extract; keep it only if you want the other molecules too).

Afterwards, once every shard is done:
    python3 extract.py --plan plan --root zinc_mirror --out zinc_db2 --collect
merges records/ into OUT/records.tsv and checks every indexed record was found.

Usage:
    python3 extract.py --plan plan --root zinc_mirror --out zinc_db2 --shard 0 --nshards 1
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import sys
import tarfile
import zlib
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from zid import local_path, read_tsv  # noqa: E402

COLS = ["zinc_id", "protomer", "charge", "m_name", "n_records", "out_file", "offset",
        "bytes", "archive", "member", "sha1"]


def m_name(data: bytes) -> str:
    for line in data.split(b"\n", 4)[:4]:
        f = line.split()
        if len(f) >= 2 and f[0] == b"M" and not f[1].lstrip(b"+-").replace(b".", b"").isdigit():
            return f[1].decode()
    return ""


def extract_one(url: str, members: dict[str, tuple], root: Path, out: Path,
                delete: bool) -> tuple[int, int]:
    arch = local_path(url, root)
    name = url.rsplit("/", 1)[1][:-len(".db2.tgz")]
    hbin = name[:3]
    done = out / "records" / f"{name}.done"
    if done.exists():
        return 0, 0
    dest = out / "db2" / hbin / f"{name}.db2.gz"
    dest.parent.mkdir(parents=True, exist_ok=True)
    (out / "records").mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".gz.tmp")
    rows, offset = [], 0
    with gzip.open(tmp, "wb", compresslevel=6) as gz, tarfile.open(arch, "r|gz") as tf:
        for m in tf:
            if m.name not in members:
                continue
            data = tf.extractfile(m).read()
            if not data.endswith(b"\n"):
                data += b"\n"
            gz.write(data)
            z, prot, charge = members[m.name]
            n_rec = sum(1 for line in data.split(b"\n") if line.strip() == b"E")
            rows.append([z, prot, charge, m_name(data), str(n_rec), str(dest.relative_to(out)),
                         str(offset), str(len(data)), str(arch), m.name,
                         hashlib.sha1(data).hexdigest()])
            offset += len(data)
    os.replace(tmp, dest)
    rec = out / "records" / f"{name}.tsv"
    rec.write_text("".join("\t".join(r) + "\n" for r in rows))
    done.write_text(f"{len(rows)} of {len(members)}\n")
    if delete and len(rows) == len(members):
        arch.unlink()
    return len(rows), len(members)


def collect(a) -> None:
    expected = {(r["db2_url"], r["member"]) for r in read_tsv(a.plan / "index.tsv")}
    urls = set((a.plan / "db2_urls.txt").read_text().split())
    expected = {e for e in expected if e[0] in urls}
    found = set()
    with open(a.out / "records.tsv", "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for rec in sorted((a.out / "records").glob("*.tsv")):
            for line in open(rec):
                fh.write(line)
                f = line.rstrip("\n").split("\t")
                found.add(f[COLS.index("member")])
    lost = sorted(m for _, m in expected if m not in found)
    (a.out / "records_missing.txt").write_text("".join(f"{m}\n" for m in lost))
    ids = {line.split("\t", 1)[0] for line in open(a.out / "records.tsv")} - {"zinc_id"}
    nrec = sum(int(line.split("\t")[4]) for line in open(a.out / "records.tsv")
               if not line.startswith("zinc_id"))
    print(f"{len(found):,} members ({nrec:,} db2 records) for {len(ids):,} ids -> "
          f"{a.out}/records.tsv; {len(lost):,} indexed members not found "
          f"(records_missing.txt)", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--plan", type=Path, default=Path("plan"))
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--out", type=Path, default=Path("zinc_db2"))
    ap.add_argument("--shard", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", 0)))
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--delete-archive", action="store_true")
    ap.add_argument("--collect", action="store_true")
    a = ap.parse_args()
    if a.collect:
        return collect(a)

    by_url = defaultdict(dict)
    for r in read_tsv(a.plan / "index.tsv"):
        by_url[r["db2_url"]][r["member"]] = (r["zinc_id"], r["protomer"], r["charge"])
    urls = sorted((a.plan / "db2_urls.txt").read_text().split())
    mine = [u for i, u in enumerate(urls) if i % a.nshards == a.shard]
    got = want = 0
    for url in mine:
        if not local_path(url, a.root).exists():
            print(f"  missing archive (not fetched?): {url}", file=sys.stderr)
            continue
        try:
            g, w = extract_one(url, by_url[url], a.root, a.out, a.delete_archive)
        except (tarfile.TarError, OSError, EOFError, zlib.error) as e:
            print(f"  unreadable {url}: {e} -- verify.sh --purge, refetch, rerun", file=sys.stderr)
            continue
        got += g; want += w
        if g < w:
            print(f"  {url}: {g} of {w} indexed members found", file=sys.stderr)
    print(f"shard {a.shard}/{a.nshards}: {len(mine)} archives, {got:,} of {want:,} members "
          f"written this run", file=sys.stderr)


if __name__ == "__main__":
    main()
