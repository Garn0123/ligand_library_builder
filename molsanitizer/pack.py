#!/usr/bin/env python3
"""
pack.py -- a collected run's per-name db2.gz files -> ~1,000-name chunks plus a byte index.

    python3 pack.py RUN [--per-chunk 1000] [--workers 8]   # pack + verify; loose files kept
    python3 pack.py RUN --verify                           # re-verify an existing pack
    python3 pack.py RUN --delete-loose                     # pack/verify, THEN delete the loose files
    python3 pack.py RUN --extract NAMES.txt -o DIR         # one DIR/<name>.db2.gz per listed name

Pass 2 writes one db2.gz per name (db2_run/out/NNNNN/<name>.db2.gz): fine for 10^5
names, too many files for 10^6+. This concatenates them, byte for byte, into
db2_packed/cNNNNN.db2.gz. Gzip members concatenate into one valid gzip stream, and
dock6 / db2tool read through them (zlib gzread; rna_pipeline checked 2026-10-22).

Chunks follow the OK names of db2_verification.csv in name order, up to --per-chunk
names each, and never split a parent: a parent's names (same 12-character base id,
so its invertomers and stereoisomers) are always in one chunk, so a docking batch of
whole chunks never cuts an invertomer group. A parent with more names than
--per-chunk gets a chunk of its own.

Writes into RUN:
  db2_packed/cNNNNN.db2.gz   the chunks
  db2_index.tsv              name, parent, chunk, offset, length, sha256: one row per OK name;
                             bytes [offset, offset+length) of the chunk are that name's
                             loose .db2.gz exactly (a valid gzip file on their own)
  db2_packed/manifest.tsv    chunk, n_names, bytes, sha256, first_name, last_name
  pack.json                  per_chunk, counts, checksums of db2_verification.csv and the index

Verify (always, after packing; or --verify alone): every chunk's size and sha256 match
the manifest, every member's sha256 matches the index, every member decompresses and
its db2 records all carry its own name, and the index covers exactly the OK names.
--delete-loose deletes db2_run/out/ files only after that passed. After that, never
resubmit pass 2 for this run (build_shard.py refuses once pack.json exists): it would
not find its kept db2s and would rebuild them. Packing again needs the loose files.

Resumable: a chunk already on disk whose bytes match its planned members is kept.
stdlib only. Run it as a job (it reads every db2 once or twice), e.g.
  srun -A maom -p maom-h200 -c 8 --mem 4G -t 1:00:00 bash msani.sh pack RUN --workers 8
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

INDEX_COLS = ["name", "parent", "chunk", "offset", "length", "sha256"]
MANIFEST_COLS = ["chunk", "n_names", "bytes", "sha256", "first_name", "last_name"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, cols: list[str], rows) -> None:
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def ok_rows(run: Path) -> list[dict]:
    with open(run / "db2_verification.csv", newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["status"] == "OK"]
    return sorted(rows, key=lambda r: r["protomer_name"])


def plan_chunks(rows: list[dict], per_chunk: int) -> list[list[dict]]:
    """Up to per_chunk names per chunk, in name order, a parent's names never split."""
    groups: list[list[dict]] = []
    for r in rows:
        if groups and groups[-1][0]["parent"] == r["parent"]:
            groups[-1].append(r)
        else:
            groups.append([r])
    chunks: list[list[dict]] = [[]]
    for g in groups:
        if chunks[-1] and len(chunks[-1]) + len(g) > per_chunk:
            chunks.append([])
        chunks[-1].extend(g)
    return [c for c in chunks if c]


def db2_names(data: bytes) -> set[str]:
    """Names on the db2 record headers (M lines whose 2nd token is not a number)."""
    names = set()
    for line in gzip.decompress(data).decode().splitlines():
        t = line.split()
        if len(t) >= 11 and t[0] == "M" and not t[1].lstrip("+-").replace(".", "").isdigit():
            names.add(t[1])
    return names


def build_chunk(run: Path, k: int, members: list[dict]) -> tuple[list[dict], dict]:
    """Write chunk k (or keep it if it already holds exactly these members)."""
    tag = f"c{k:05d}"
    path = run / "db2_packed" / f"{tag}.db2.gz"
    loose = [run / "db2_run" / "out" / m["db2_file"] for m in members]
    missing = [str(p) for p in loose if not p.exists()]
    if missing:
        if path.exists():
            raise SystemExit(f"{tag}: loose files gone ({missing[0]}, ...); "
                             "an existing pack can be verified (--verify), not rebuilt")
        raise SystemExit(f"{tag}: {len(missing)} OK name(s) have no db2.gz, e.g. {missing[0]}")
    index, offset, whole = [], 0, hashlib.sha256()
    blobs = []
    for m, p in zip(members, loose):
        data = p.read_bytes()
        blobs.append(data)
        whole.update(data)
        index.append({"name": m["protomer_name"], "parent": m["parent"], "chunk": f"{tag}.db2.gz",
                      "offset": offset, "length": len(data), "sha256": hashlib.sha256(data).hexdigest()})
        offset += len(data)
    digest = whole.hexdigest()
    if not (path.exists() and path.stat().st_size == offset and sha256_file(path) == digest):
        tmp = path.with_name(path.name + ".tmp")
        with open(tmp, "wb") as fh:
            for data in blobs:
                fh.write(data)
        tmp.replace(path)
    man = {"chunk": f"{tag}.db2.gz", "n_names": len(members), "bytes": offset, "sha256": digest,
           "first_name": members[0]["protomer_name"], "last_name": members[-1]["protomer_name"]}
    return index, man


def verify_chunk(run: Path, man: dict, members: list[dict]) -> list[str]:
    path = run / "db2_packed" / man["chunk"]
    if not path.exists():
        return [f"{man['chunk']}: missing"]
    data = path.read_bytes()
    bad = []
    if len(data) != int(man["bytes"]) or hashlib.sha256(data).hexdigest() != man["sha256"]:
        bad.append(f"{man['chunk']}: size or sha256 differs from the manifest")
    for m in members:
        blob = data[int(m["offset"]): int(m["offset"]) + int(m["length"])]
        if hashlib.sha256(blob).hexdigest() != m["sha256"]:
            bad.append(f"{m['name']}: bytes differ from the index")
            continue
        try:
            names = db2_names(blob)
        except (OSError, EOFError, UnicodeDecodeError) as e:
            bad.append(f"{m['name']}: does not decompress ({e})")
            continue
        if names != {m["name"]}:
            bad.append(f"{m['name']}: db2 records named {sorted(names)[:3]}")
    return bad


def verify(run: Path, workers: int) -> int:
    index = read_tsv(run / "db2_index.tsv")
    manifest = read_tsv(run / "db2_packed" / "manifest.tsv")
    by_chunk: dict[str, list[dict]] = {}
    for r in index:
        by_chunk.setdefault(r["chunk"], []).append(r)
    problems = []
    want = {r["protomer_name"] for r in ok_rows(run)}
    have = [r["name"] for r in index]
    if len(have) != len(set(have)):
        problems.append("index lists a name more than once")
    if set(have) != want:
        problems.append(f"index vs OK names: {len(want - set(have))} missing, {len(set(have) - want)} extra")
    if set(by_chunk) != {m["chunk"] for m in manifest}:
        problems.append("index and manifest name different chunks")
    for m in manifest:
        if int(m["n_names"]) != len(by_chunk.get(m["chunk"], [])):
            problems.append(f"{m['chunk']}: manifest says {m['n_names']} names, index {len(by_chunk.get(m['chunk'], []))}")
    with ThreadPoolExecutor(max(1, workers)) as ex:
        for bad in ex.map(lambda m: verify_chunk(run, m, by_chunk.get(m["chunk"], [])), manifest):
            problems += bad
    meta = json.loads((run / "pack.json").read_text())
    if meta["verification_sha256"] != sha256_file(run / "db2_verification.csv"):
        problems.append("db2_verification.csv changed since packing (collect rerun?): pack again")
    if problems:
        print(f"pack verify FAILED ({len(problems)} problem(s)):", file=sys.stderr)
        for p in problems[:30]:
            print("  " + p, file=sys.stderr)
        return 1
    print(f"pack verify OK: {len(index)} names in {len(manifest)} chunks, "
          f"{sum(int(m['bytes']) for m in manifest) / 1e9:.2f} GB")
    return 0


def delete_loose(run: Path) -> None:
    out = run / "db2_run" / "out"
    n = 0
    names = {r["name"] for r in read_tsv(run / "db2_index.tsv")}
    for shard in sorted(p for p in out.iterdir() if p.is_dir()):
        for f in shard.iterdir():
            if f.name.endswith(".db2.gz") and f.name[: -len(".db2.gz")] in names:
                f.unlink()
                n += 1
        try:
            shard.rmdir()
        except OSError:
            pass                                        # holds something not packed: keep it
    meta = json.loads((run / "pack.json").read_text())
    meta["loose_deleted"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    (run / "pack.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"deleted {n} loose db2.gz files under {out}")


def extract(run: Path, names_file: Path, dest: Path) -> int:
    want = [l.split()[0] for l in open(names_file) if l.strip()]
    index = {r["name"]: r for r in read_tsv(run / "db2_index.tsv")}
    lost = [n for n in want if n not in index]
    if lost:
        print(f"{len(lost)} name(s) not in the index, e.g. {lost[:3]}", file=sys.stderr)
        return 1
    dest.mkdir(parents=True, exist_ok=True)
    handles: dict[str, object] = {}
    try:
        for n in want:
            r = index[n]
            fh = handles.get(r["chunk"]) or handles.setdefault(r["chunk"], open(run / "db2_packed" / r["chunk"], "rb"))
            fh.seek(int(r["offset"]))
            data = fh.read(int(r["length"]))
            if hashlib.sha256(data).hexdigest() != r["sha256"]:
                print(f"{n}: bytes differ from the index", file=sys.stderr)
                return 1
            (dest / f"{n}.db2.gz").write_bytes(data)
    finally:
        for fh in handles.values():
            fh.close()
    print(f"wrote {len(want)} db2.gz files to {dest}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--per-chunk", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=8, help="threads (file-system latency, not CPU)")
    ap.add_argument("--verify", action="store_true", help="only verify an existing pack")
    ap.add_argument("--delete-loose", action="store_true", help="after a passing verify, delete db2_run/out/ files")
    ap.add_argument("--extract", type=Path, metavar="NAMES", help="file of names (first column)")
    ap.add_argument("-o", "--out", type=Path, help="with --extract: destination directory")
    a = ap.parse_args(argv)
    run = a.run.resolve()

    if a.extract:
        if not a.out:
            ap.error("--extract needs -o DIR")
        return extract(run, a.extract, a.out)
    if not (run / "db2_verification.csv").exists():
        raise SystemExit(f"{run}: no db2_verification.csv; run `msani.sh collect` first")
    if a.verify:
        return verify(run, a.workers)

    meta_path = run / "pack.json"
    if meta_path.exists():
        old = json.loads(meta_path.read_text())
        if old.get("loose_deleted"):
            print(f"already packed and loose files deleted ({old['loose_deleted']}); verifying")
            return verify(run, a.workers)
        if old["per_chunk"] != a.per_chunk:
            raise SystemExit(f"{run} was packed with --per-chunk {old['per_chunk']}; "
                             "delete db2_packed/, db2_index.tsv and pack.json to repack differently")

    rows = ok_rows(run)
    if not rows:
        raise SystemExit("no OK names in db2_verification.csv")
    chunks = plan_chunks(rows, a.per_chunk)
    (run / "db2_packed").mkdir(exist_ok=True)
    t0 = time.time()
    with ThreadPoolExecutor(max(1, a.workers)) as ex:
        built = list(ex.map(lambda kc: build_chunk(run, kc[0], kc[1]), enumerate(chunks, start=1)))
    planned = {f"c{k:05d}.db2.gz" for k in range(1, len(chunks) + 1)}
    for f in (run / "db2_packed").glob("c*.db2.gz"):
        if f.name not in planned:
            f.unlink()                              # from an older, longer plan
    write_tsv(run / "db2_index.tsv", INDEX_COLS, [r for idx, _ in built for r in idx])
    write_tsv(run / "db2_packed" / "manifest.tsv", MANIFEST_COLS, [m for _, m in built])
    sizes = [int(m["n_names"]) for _, m in built]
    meta = {"per_chunk": a.per_chunk, "n_names": len(rows), "n_chunks": len(chunks),
            "names_per_chunk_min_max": [min(sizes), max(sizes)],
            "bytes": sum(int(m["bytes"]) for _, m in built),
            "verification_sha256": sha256_file(run / "db2_verification.csv"),
            "index_sha256": sha256_file(run / "db2_index.tsv"),
            "packed": time.strftime("%Y-%m-%dT%H:%M:%S"), "seconds": round(time.time() - t0, 1),
            "host": os.uname().nodename}
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    print(f"packed {len(rows)} names into {len(chunks)} chunks ({min(sizes)}-{max(sizes)} names each) "
          f"in {meta['seconds']} s")

    rc = verify(run, a.workers)
    if rc == 0 and a.delete_loose:
        delete_loose(run)
    elif a.delete_loose:
        print("verify failed: loose files NOT deleted", file=sys.stderr)
    return rc


if __name__ == "__main__":
    sys.exit(main())
