# db2pipe

Tools for repacking large libraries of ZINC `.db2.gz` files into fixed-size,
load-balanced chunks with a global manifest, so they can be fed into a docking
pipeline as evenly sized units of work. Built to handle a ~500 GB / 8.5 M
molecule ZINC 3D library, but works on anything from a handful of files up.

There are two ways to run it:

- **Serial pipeline** (`serial/`, stages `01`–`04`): simple, single-process,
  good up to a few million molecules or for development and testing.
- **Parallel pipeline** (`parallel/`: `make_shards` + `p1`–`p3` + `submit.slurm`):
  a SLURM job-array version of the same thing that scales to the full library.

Both produce the **same output format** — a directory of chunk files plus a
`manifest.tsv` — and both feed the same labeling stages (`02`, `03`).

---

## The problem this solves

A ZINC 3D download is a tree of gzipped `.db2` files. Two things make naive
splitting wrong:

1. **`.db2` is record-oriented, not line-oriented.** Each molecule is a block
   that begins with an `M` header line and ends with an `E` terminator, with a
   variable number of atom / bond / coordinate / conformer lines in between.
   Splitting on line counts (`split -l`, `zcat | split`) cuts molecules in
   half. Every tool here splits **only on record boundaries**.

2. **Tranche files are sorted by molecular size.** ZINC tranche directories are
   named by heavy-atom count and logP, so walking the tree in sorted order
   walks small molecules → large. Chunk *N* by simple accumulation ends up
   internally homogeneous, and since docking cost scales with atom and
   conformer count, per-chunk runtime can vary ~6× even when every chunk holds
   the same molecule count. The pipeline balances against this.

`grep -c ZINC` is **not** a reliable molecule count for these files: a record
can carry the ZINC token on more than one `M` line, and a truncated file has a
dangling record with no `E`. The record count (number of `E` terminators) is
the ground truth, and the tools report both and warn when they disagree.

---

## db2 record model (the core assumption)

Everything depends on this shape:

```
M ZINC000000000001 ...        <- header line(s), start with 'M', carry the id
M smiles ...
A ...                         <- atoms
B ...                         <- bonds
X ...                         <- coordinates
C ...                         <- conformers (many, on big molecules)
E                             <- terminator, a line that is exactly "E"
```

- A record runs from an `M` line to the next `E` line inclusive.
- The **id** is the `ZINC…` token on the `M` (header) lines (old format), or —
  for ZINC22, where the id overflows the fixed-width header field — the truncated
  first token of the first `M` line (see below).
- Records are addressed by **position** (which chunk, which index within it),
  **never by id** — ids can legitimately repeat (protonation states,
  tautomers, stereoisomers), so id is not a key.

The fast parser (`iter_records_bytes` in `db2common.py`) splits on the byte
sequence `\nE\n`. **If a real db2 file can contain a line that starts with `E`
but is not a terminator, this assumption breaks** — see "Open questions" below.
The slow parser (`iter_records`) splits on any line whose first character is
`E`; the two are verified to produce byte-identical records on the test data.

---

## ZINC22 tarballs and the new id scheme

Newer data (ZINC22) ships differently, and the parallel pipeline handles both
old and new transparently:

- **Inputs are `.db2.tgz` archives**, each holding thousands of tiny `.db2`
  members under a nested prefix
  (`H05/H05M000/2t/tm/ZINC550000002ttm.0.O.db2`). They are **read in place,
  member by member** (`iter_sources`) — never extracted, because millions of
  tiny files would exhaust inodes on a shared filesystem. `make_shards` shards on
  the archive list; `p1_collect` streams members; a member ending in `.gz` is
  decompressed on the fly.
- **The manifest `source_file` becomes `<archive>::<member>`**, giving exact
  per-molecule provenance the old format lacked. (No pipeline stage reopens
  source paths from the manifest, so nothing downstream needs to split on `::`.)
- **The id is truncated in the header.** The db2 id field is fixed-width, so a
  full id + conformer like `ZINC550000002ttm.0` (18 chars) is *left*-truncated to
  `NC550000002ttm.0` (16) in the record. The pipeline stores the **full** id
  (parsed from the member filename, `id_from_name`) in the manifest, while
  `03_apply` verifies and relabels against the truncated in-record form
  tolerantly (`id_match` — the record id is a suffix of the full one). Old bare
  `.db2.gz` inputs are unaffected: their id still comes from the header `ZINC`
  token.
- **Corrupt archives fail loudly, per archive** (a bad `.db2.tgz` is caught at
  the container level and skips only that shard's archive; a malformed member is
  caught while reading). Purge and re-fetch corrupt archives rather than working
  around them.

The **serial pipeline (`01`–`04`) still reads only bare `.db2.gz`** — use the
parallel pipeline for `.db2.tgz` archives (which is what you'd use at ZINC22
scale anyway).

---

## Repository layout

```
db2pipe/
  db2common.py            shared library: parsers, tar/bare source reader
                          (iter_sources), id extraction, relabeling, inode
                          dedup, weight functions, progress reporting
  serial/                 standalone single-process pipeline
    01_chunk.py           walk tree -> fixed-size chunks + manifest
    02_label.py           manifest -> labels.tsv (decide new ids; touches no data)
    03_apply.py           apply labels.tsv to chunks (position-addressed rewrite)
    04_rebalance.py       post-hoc: redistribute existing chunks to balance load
    05_reindex.py         post-QC: reindex a manifest after molecules were
                          removed from the chunks (applies a db2tool map)
  parallel/               SLURM job-array pipeline (same output format)
    check_inputs.py       preflight: scan the tree for corrupt/truncated files
    make_shards.py        split input tree into shard lists; size BINS from a target
    p1_collect.py         [array] parse a shard, stride records into per-bin parts
    p2_assemble.py        [array] concatenate a bin's parts into a finished chunk
    p3_finalize.py        stitch per-shard manifests into one global manifest
    submit.slurm          driver: chains all four phases with job dependencies
  tests/                  pytest suite (see tests/README.md)
  README.md               this file
  HPC_EXECUTION_GUIDE.md  step-by-step SLURM execution guide
```

`serial/02_label.py` and `serial/03_apply.py` are shared by both pipelines and
run serially in both — the parallel pipeline calls them after assembly. The
shared library `db2common.py` sits at the repo root; every stage script imports
it from there.

`serial/05_reindex.py` is optional and runs *after* QC rather than as part of
compilation. Quality tooling
([db2tool](https://github.com/Garn0123/dock6_claude/tree/feat/db2tool)) can drop
molecules from finished chunks — corrupt records, or molecules whose every
conformer set is flagged clashing — and that renumbers every molecule after the
first deletion. The manifest's `(chunk, chunk_idx)` keys then point at the wrong
rows, silently. db2tool writes a map of what moved where; this stage applies it:

```bash
db2tool subset --drop-all-broken --map remap.tsv -o clean/c1.db2.gz chunks/c1.db2.gz
python3 serial/05_reindex.py --map remap.tsv -m chunks/manifest.tsv -o clean
```

Rows whose `new_idx` is `-1` are dropped; the rest take `chunk_idx = new_idx`.
Chunks absent from the map pass through untouched, so reindexing one chunk of a
set is safe. Like stage 3 it verifies ids before rewriting and aborts on a stale
map. **If QC only *repairs* molecules rather than removing them, no reindexing
is needed** — repair preserves count and order, so the manifest stays valid.

---

## Output contract

Both pipelines produce, in the output directory:

- `chunk_00001.db2.gz`, `chunk_00002.db2.gz`, … — real gzip, each a whole
  number of records, every chunk starting with `M` and ending with `E`.
- `manifest.tsv` with header:
  `chunk<TAB>chunk_idx<TAB>source_file<TAB>source_idx<TAB>original_id`
  one row per molecule, where `chunk_idx` is the 0-based position of that
  molecule inside its chunk. This is the addressing scheme `03_apply.py` relies
  on.

The serial `01_chunk.py` also writes `chunk_sources.tsv` (a human-readable
chunk→source summary).

---

## Serial pipeline

```bash
# 1. chunk + manifest
python3 serial/01_chunk.py -i /path/to/3D -o chunks -n 50000

# 2. decide new ids for duplicates (reads text only, fast)
python3 serial/02_label.py -m chunks/manifest.tsv -o labels.tsv

# 3. apply the labels (parallelizable per chunk with --only)
python3 serial/03_apply.py -c chunks -L labels.tsv -o chunks_labelled
```

### 01_chunk.py — chunk and manifest

Key options:

- `-n / --chunk-size` — molecules per chunk (hard ceiling; default 50000).
- `--order {path,shuffle}` — `path` is sorted (= sorted by size in a tranche
  tree); `shuffle` mixes tranches so each chunk is a representative sample.
  Shuffle alone cuts conformer spread from ~6.4× to ~1.9× at zero cost and with
  no cost model.
- `--weight {count,bytes,lines:X}` + `--target-weight N` — close a chunk when
  summed weight reaches `N`. `lines:C` weights by conformer count (best proxy
  for docking time). **`-n` still applies as a ceiling and will silently defeat
  the balancer if set too low** — in weight mode set `-n` well above any
  expected chunk size. Run with `--dry-run --weight lines:C` first; it prints
  the total weight and suggested `--target-weight` values.
- `--suffix EXT` (repeatable) — collect more than `.db2.gz` (e.g. bare `.db2`).
- `--keep-duplicate-paths` — by default symlinks/hardlinks to a file already in
  the set are skipped (deduped by inode); this disables that.
- `--progress-interval SEC` — time-based progress to stderr (default 30; `0`
  fully silent). Percent/ETA are byte-based, not file-count-based.

### 02_label.py — decide new ids

Reads the manifest, writes `labels.tsv`
(`chunk<TAB>chunk_idx<TAB>original_id<TAB>new_id`) containing **only the records
whose id changes**. Modes:

- `duplicates` (default) — suffix only ids that occur 2+ times; unique
  molecules keep their bare catalog id. Two-pass over the manifest; for very
  large libraries precompute the dup set on disk and pass `--dup-ids`:
  ```bash
  tail -n +2 chunks/manifest.tsv | cut -f5 | sort | uniq -d > dups.txt
  python3 serial/02_label.py -m chunks/manifest.tsv -o labels.tsv --dup-ids dups.txt
  ```
- `occurrence` — every id gets an nth-seen suffix (`_1`, `_2`, …). Holds every
  unique id in memory.
- `serial` — every id gets a global running index. O(1) memory.

`--keep-first` leaves the first occurrence unsuffixed. `--sep`, `--width` format
the suffix. `--max-id-len N` warns on over-long ids (see the column-shift
caveat below).

### 03_apply.py — apply labels

Rewrites chunks, locating each record by `(chunk, chunk_idx)`, not by id — so
duplicate ids are never ambiguous. **Before rewriting a record it verifies the
id at that position matches what `labels.tsv` expects**; a mismatch means the
labels are stale (e.g. `01` was re-run with a different chunk size after `02`)
and the whole run aborts with a nonzero exit rather than corrupting data. Writes
through a `.tmp` + `os.replace`, so each chunk is all-or-nothing.

- `--only chunk_00007.db2.gz` — process a single chunk (for a job array).
- `--copy-unedited` — also copy chunks with no edits (default: leave originals
  in place).

### 04_rebalance.py — post-hoc rebalancing

Redistributes **already-chunked** data at the record level, without re-reading
the source tree. Use this if you chunked before adding balancing.
(Re-running `01` on the chunks with `--order shuffle` does **not** work — it
just permutes whole 50k blocks.)

```bash
python3 serial/04_rebalance.py -c chunks -m chunks/manifest.tsv -o rebalanced -N 10 --mode stride
```

- `--mode stride` — record *i* → bin *i mod N*; systematic sampling across the
  size-sorted input. No cost model; gives ~1.0× spread and equal counts.
- `--mode greedy --weight lines:C` — least-loaded-bin assignment against a cost
  metric.
- Pass `-m` to carry original provenance into the new manifest. A mismatched
  manifest is detected and rejected. Opens N gzip writers at once, so it checks
  `RLIMIT_NOFILE` up front and raises the soft limit where possible.

**Rebalancing changes every record's position, so any existing `labels.tsv` is
stale afterward — re-run `02` against the new manifest.**

---

## Parallel pipeline (SLURM)

Same result as the serial pipeline, scaled out. Edit the config block at the top
of `parallel/submit.slurm` and run `bash parallel/submit.slurm`. It submits four
dependency-chained jobs:

| Phase | Script | Array | Work |
|------|--------|-------|------|
| plan | `make_shards.py` | 1 | walk tree, pack files into size-balanced shard lists |
| collect | `p1_collect.py` | S | parse a shard, stride its records into a *part* of every bin |
| assemble | `p2_assemble.py` | N | concatenate a bin's parts into a finished chunk |
| final | `p3_finalize.py` | 1 | stitch per-shard manifests into one global manifest |

**Why no shuffle step is needed.** Each collect task strides its own records
across all N bins, so every bin receives an even slice of every shard. The bins
come out globally balanced with no cross-task communication, even though each
shard covers only one region of the size distribution.

**Why assembly is cheap.** Concatenated gzip members form a valid gzip stream
(verified with both `zcat` and Python's `gzip`). So `p2_assemble.py` is a pure
byte copy — **no decompression, no recompression, no parsing** — and runs at
disk speed. Parts are joined in shard order, which is the order `p3` assumes
when computing global indices.

**The subtle part: global indices in `p3_finalize.py`.** A record's index inside
a finished chunk is *not* its index inside the shard that produced it, because a
chunk is `part(shard 0) ++ part(shard 1) ++ …`. The global index is the sum of
every earlier shard's contribution to that bin (from the `counts/` files) plus
the record's local index. This arithmetic is the thing most likely to break
under changes; the verification below checks every position against the actual
assembled bytes.

### Intermediate layout (under `WORK_DIR`)

```
work/
  shards/shard_00000.txt ...      one file list per collect task
  plan.tsv                        shard -> (files, bytes)
  parts/bin_00001/part_00007.db2.gz   shard 7's slice of bin 1
  manifests/shard_00007.tsv       bin, local_idx, source_file, source_idx, id
  counts/shard_00007.tsv          bin, records, weight   (offsets for p3)
  logs/                           per-task stdout/stderr
```

`p2_assemble.py` deletes each bin's parts after a successful join (unless
`--keep-parts`); `p3_finalize.py` removes `manifest_parts/` at the end.

### Sizing

- **`BINS` (number of finished chunks) is auto-derived — don't guess it.** The
  molecule count isn't known until collect runs, but the number of bins has to be
  fixed before it strides records. So `make_shards --target-per-bin 50000`
  estimates the count from a sample of files (molecules-per-compressed-byte ×
  total bytes), computes `BINS = ceil(estimate / target)`, and writes it to
  `work/bins.txt`; `p1_collect` reads that when `-N` is omitted, and
  `submit.slurm` uses it to size the assemble array. `submit.slurm` is driven by
  `TARGET_PER_BIN`, not a hardcoded `BINS`. The estimate is ~±10-20% (fine for a
  "~50k" target); pass `-N` to override it, or a leave `BINS=` set explicitly in
  `submit.slurm`. This is what prevents the classic "I guessed 8 M, it was 32 M,
  now every chunk is 180 k" mistake.
- `SHARDS` = collect array width. ~200 is a reasonable start. Intermediates are
  S×N part files, and N scales with the target (200×640 ≈ 128 000 for 32 M at
  50 k/chunk) — fine on Lustre, but be mindful before pushing S much higher.
- Scratch: phase 1 writes parts totaling ~1× the input size, freed incrementally
  by phase 2. Budget ~input-size on `$SCRATCH`; do not point `WORK_DIR` at home.

### Failure recovery

- Collect tasks are **idempotent** (re-running a shard yields a byte-identical
  shard manifest), so a failed task can just be resubmitted:
  `sbatch --array=3,17,102 ... p1_collect.py ...`, then rerun assemble + final.
- `--dependency=afterok` on an array waits for **every** task, so one failed
  shard blocks assembly rather than silently producing short chunks.

---

## Performance notes

- Observed baseline was ~250 mol/s serial. At ~59 KB/molecule that's ~15 MB/s —
  too slow to be cluster I/O, so the bottleneck is **Python per-line iteration**
  over records with many conformer lines. `p1_collect.py` uses the bytes-level
  `iter_records_bytes` splitter, which was **~3× faster** on the test records
  (the gap grows with conformer count) and is verified byte-identical to the
  line parser.
- The serial `01`–`04` still use the line-based `iter_records`. Switching them
  to the bytes parser is a natural, unclaimed optimization (see below).
- All progress goes to **stderr**; explicit `flush()` after each line. On
  Python < 3.9 redirected stderr is block-buffered — run with `python3 -u` on
  older cluster interpreters (the `submit.slurm` wrap already does).

---

## How to verify a run

The invariants worth re-checking after any change:

```bash
# 1. every chunk is whole records: starts with M, ends with E
for f in chunks/*.db2.gz; do
  printf "%s first=%s last=%s\n" "$f" \
    "$(zcat "$f" | head -1 | cut -c1)" "$(zcat "$f" | tail -1 | cut -c1)"
done

# 2. nothing lost or duplicated vs the source tree (id multiset identical)
find /path/to/3D -name '*.db2.gz' | sort | xargs zcat | grep -o 'ZINC[0-9]*' | sort > /tmp/src.txt
zcat chunks/*.db2.gz | grep -o 'ZINC[0-9]*' | sort > /tmp/out.txt
diff -q /tmp/src.txt /tmp/out.txt && echo "OK: identical id multiset"

# 3. non-M lines are untouched by chunking/labeling (byte-identical)
zcat chunks/*.db2.gz | grep -v '^M' | md5sum

# 4. manifest positions match the assembled bytes
#    (walk each chunk, compare extract_id at each index to the manifest row)
```

An automated **pytest suite lives in `tests/`** (`python3 -m pytest` from the
repo root; needs only `pip install -r requirements-dev.txt`). It covers the
invariants above — record-boundary integrity, id-multiset preservation,
manifest-position correctness, stale-label rejection, the parallel global-index
arithmetic, and the truncated/non-gzip/duplicate-path edge cases — using
synthetic `M…E` records generated on the fly (`tests/db2gen.py`) with
controllable atom/conformer counts and deliberately lopsided per-tranche
populations. The ad-hoc shell checks above remain a useful final sanity pass on
real data. See `tests/README.md` for the layout.

---

## Known caveats / gotchas

- **Column shift on relabel.** Appending a suffix lengthens the id token, which
  shifts every column to its right on that `M` line. Fine if the db2 reader
  splits `M` lines on whitespace; **silent corruption if anything reads fixed
  columns.** Always push one labeled chunk through the real docking pipeline
  before committing to a full pass. `--max-id-len` warns on length.
- **Labeling is not idempotent.** Running `03` over already-labeled output
  stacks a second suffix (`ZINC…_1_1`). Keep labeled and unlabeled trees
  separate.
- **`E`-line collision (unverified).** The fast parser assumes no non-terminator
  line begins with `E` (and that `\nE\n` only ever delimits records). Not yet
  confirmed against real ZINC db2 output — verify before trusting the parallel
  pipeline at scale.
- **`grep -c ZINC` overcounts / miscounts.** Use the record (`E`) count.
- **Duplicate *files*** (symlinks/hardlinks) are deduped by inode by default;
  duplicate **content** in two distinct files is not detected. Duplicate **ids**
  are expected and preserved, not an error.
- **Truncated / corrupt / non-gzip inputs** are detected: truncated records and
  unreadable files cause a nonzero exit with a list; files that carry `.gz` but
  aren't gzipped are read as plain text and reported (the `.gz` extension is not
  trusted — content is sniffed for the `1f 8b` magic bytes). A file with a valid
  gzip header but a **damaged deflate body** (a partial download / bad transfer)
  raises `zlib.error` mid-stream — which is *not* an `OSError` — so the readers
  catch it explicitly; a corrupt file is named and the stage fails cleanly
  instead of crashing with a traceback. To find every bad file in one pass
  before submitting, run the preflight scanner and split the tree:
  ```bash
  python3 parallel/check_inputs.py -i /path/to/3D \
      --good-list good.txt --bad-list bad.txt
  # repair/re-fetch the bad ones, or just run over the clean set:
  python3 parallel/make_shards.py --file-list good.txt -o work -S 200
  ```
  To finish an already-planned run over a library with a few known-bad files
  without re-planning, rerun just the affected collect shards with
  `p1_collect.py --skip-corrupt`: it keeps whatever reads cleanly, drops the
  damaged files (still listing them), and exits 0 so assembly can proceed.
- **File-descriptor limits.** Rebalance and collect open N/bins gzip writers at
  once; they raise `RLIMIT_NOFILE` where the hard limit allows and otherwise
  fail fast. For thousands of bins, raise `ulimit -n` or reduce N.

---

## Open questions / next steps

Unclaimed work, roughly in priority order:

1. **Confirm the `E`-line assumption** against real ZINC db2 files, and if it can
   be violated, make `iter_records_bytes` robust to it.
2. **Extend the test suite.** `tests/` now covers parser equivalence,
   id-multiset preservation, manifest-position correctness, stale-label
   rejection, the parallel global-index arithmetic, and the
   truncated/non-gzip/duplicate-path edge cases. Still uncovered:
   `04_rebalance.py`, and running the invariants against a real ZINC fixture
   rather than only synthetic records.
3. **Switch the serial `01`–`04` to the bytes parser** for the ~3× speedup, once
   (1) is settled.
4. **Verify the column-shift question** for the specific db2 reader in use;
   document the answer here.
5. Consider whether `02`/`03` should be foldable into a single labeled-output
   pass for the parallel pipeline (currently labeling is a separate serial step
   after assembly).
6. A top-level wrapper that runs the serial pipeline end-to-end (`01`→`02`→`03`)
   with one command, for the small-library case.

---

## Requirements

Python 3 standard library only — no third-party packages. Developed against
CPython 3.12; needs ≥ 3.9 for line-buffered stderr behavior (or use `python3 -u`).
SLURM is only needed for the parallel pipeline; the serial pipeline runs
anywhere.
