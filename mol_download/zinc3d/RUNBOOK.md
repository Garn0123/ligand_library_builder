# ZINC22's own 3D db2s for LLB's molecules: runbook

Goal: for the ZINC22 parents of an LLB library (H18-22 now, H18-H30 next), get
the db2 records ZINC itself built, so the same molecules can be docked from both
sources and compared.

Nothing here has run on Great Lakes yet. Duo means John runs every command.
The tools were tested end to end on a laptop against files.docking.org, on real
H30 chunks: catalog, plan, both index passes, extract and collect. curl stood in
for wget, which the laptop lacks. The extracted members matched the archives
byte for byte.

## The short version

- **Downloading the full set is impossible.** ZINC22 3D for H18-H30 is
  **~1.8 PB** of `.db2.tgz`, H18-22 alone **~390 TB**. These are estimates from
  2026-10-19, good to perhaps ±50% (table below). Every bin is spread over up to
  25 generations (`zinc-22a` … `y`) and ~1M+ chunks of 0.1-0.7 GB each.
- **An id names its tranche, not its file.** `ZINCk5000003UMKO` → H20 (`k`),
  logP bucket M000 (`5`). The generation, charge and chunk are not in the id.
  Each chunk has a `.pdbqt.tgz` sidecar with exactly the db2's member names
  (`H20/H20M000/0L/mS/ZINCk50000040LmS.0.L.pdbqt`) at **0.48%** of the bytes.
  Listing sidecars gives the exact id → chunk index.
- **The catch: uniform parents mean one chunk per parent.** LLB's parents are a
  uniform draw from ZINC22 2D, so almost every parent sits in its own chunk.
  Each covered parent costs a whole chunk, **~0.43 GB** for H18-22. Covering
  every parent of the 50k H18-22 library would take **~17 TB** of db2 plus
  1.8 TB of sidecars and ~1M directory listings. Most H18-22 parents do have a
  ZINC 3D build: ~100% for H18-21 and ~55% for H22, so ~80% of the library.
- **Recommendation: sample cells.** `catalog.py --fraction G` lists only a
  seeded random share G of the (generation, tranche) cells. Every later cost
  scales with G, and the parents that fall in those cells are an unbiased
  (cluster) sample of the library. **G = 0.03 → ~1,200 parents for ~0.5 TB of
  db2 + ~55 GB of sidecars**, about a day end to end. Raise G later with the
  same seed and nothing already done is repeated.

## Sizes (2026-10-19)

| bin | ids listed (M) | measured TB | total TB (est.) | ZINC 3D / ZINC 2D |
|---|---:|---:|---:|---:|
| H18 | 259 | 1.6 | 22 | ~1 |
| H19 | 446 | 12.1 | 49 | ~1 |
| H20 | 866 | 0.2 | 80 | ~1 |
| H21 | 912 | 0.4 | 100 | 0.88 |
| H22 | 1,067 | 0.5 | 135 | 0.56 |
| H23 | 1,165 | 0.9 | 181 | 0.35 |
| H24 | 1,903 | 1.5 | 416 | 0.34 |
| H25 | 1,405 | 1.0 | 277 | 0.14 |
| H26 | 1,480 | 3.4 | 327 | 0.10 |
| H27 | 13 | 0.1 | 3 | ~0 |
| H28 | 7 | 0.1 | 2 | ~0 |
| H29 | 449 | 0.1 | 192 | 0.01 |
| H30 | 20 | 6.6 | 7 | 0.09 |
| **all** | **9,994** | **28.5** | **~1,790** | |

How these were obtained (inputs in `estimate_2026-10-19/`):

- **Measured** = 33 (generation, bin) cells crawled completely: 75,756 chunks
  in zinc-22a, b, c, d and f (H18-19).
- **Total** adds every other cell as ids × the bin's bytes per id × a
  per-generation factor. The factor (0.7-2.8) is measured on one mid-logP
  tranche of that generation.
- **ids listed** = the generations' per-bin id lists (`H20.txt`). These run
  ~5% above what was built (zinc-22a H20P200 N: 10,163 listed, 9,670 built).
- **3D / 2D** compares those ids with the zinc22/2d tranche sizes at 11.7
  bytes per SMILES line. It reads above 1 for H18-20, where SMILES are shorter.
  Read it as "nearly all".

The last column matters for the H18-H30 library. ZINC has built only ~1/3 of
H23-24, ~1/10 of H25-26 and almost none of H27-H29. Most heavy LLB parents
have no ZINC db2 to compare against, whatever is downloaded.

## What a budget buys (50k H18-22 library)

Parents per bin follow `sample_2d.py --total` (druglike_10k: H18 5%, H19 9%,
H20 16%, H21 24%, H22 46%). That gives ~40,000 parents with a ZINC 3D build,
at a molecule-weighted chunk of 0.25 GB (H18) to 0.5 GB (H22).

| G (`--fraction`) | parents found | listing | sidecars | db2 | download time |
|---:|---:|---:|---:|---:|---|
| 0.01 | ~400 | ~1 h | ~18 GB | ~0.17 TB | ~2-4 h |
| **0.03** | **~1,200** | **~2-3 h** | **~55 GB** | **~0.5 TB** | **~5-12 h** |
| 0.06 | ~2,400 | ~5 h | ~110 GB | ~1.0 TB | ~10-24 h |
| 0.12 | ~4,800 | ~10 h | ~220 GB | ~2.0 TB | ~1-2 days |
| 1 | ~40,000 | ~3 days | ~1.8 TB | ~17 TB | weeks |

Assumptions behind the columns:

- **Listing** assumes `-j 2`, from the crawl's measured ~500 chunks/min at
  `-j 4`.
- **Download time** runs from 4 × 14 MB/s (one connection's measured speed,
  sample_2d.py 2026-09) to ~15 MB/s total (if the server does not scale), plus
  ~1 s per sidecar file.
- **Extraction** is gzip-bound, ~50-100 MB/s of archive per core, so ~2-5
  core-h per TB, spread over the array.
- **Cluster sampling makes small G noisy.** Cell sizes are very uneven, so
  the parent count at small G varies around these numbers. `plan.py` prints
  the real expectation for the cells actually listed (`estimate.tsv`:
  `findable`, `db2_expect_TB`).

## One-time setup

```bash
# config/hpc.env (or $LLB_CONFIG): the account is maom2 since 2026-10-16
grep -E '^(LLB_ACCOUNT|LLB_PARTITION|WGET|PARALLEL|PREP_ENV)=' $(llb root)/config/hpc.env
#   LLB_ACCOUNT=maom2    WGET=wget    PARALLEL=parallel
#   PREP_ENV needs python >= 3.8; the scripts are stdlib only

# compute nodes must reach files.docking.org (fetch.sbatch checks again)
srun -A maom2 -p standard -t 5 curl -sI https://files.docking.org/zinc22/ | head -1
```

The code is on branch `feat/zinc22-db2` until merged:

```bash
cd /nfs/turbo/umms-maom/johnbick/opt/ligand_library_builder
git fetch && git checkout feat/zinc22-db2
```

## Run

One landing directory on scratch holds the mirror (`zinc22/...`), catalog,
plan and logs. Every command runs from it.

```bash
R=/scratch/maom_root/maom99/johnbick/2026_10_zinc3d_h18_22
mkdir -p $R/logs && cd $R
P=/scratch/maom_root/maom99/johnbick/2026_10_08_build_small_mols/<library run>/parents/parents.tsv
```

### 0. Catalog a sample of cells (SLURM, listings only)

```bash
sbatch -A maom2 -p standard -c 2 --mem=4G -t 24:00:00 -o logs/catalog_%j.out \
  --wrap "llb zinc3d-catalog --heavy 18-22 --fraction 0.03 --out catalog -j 2"
# when it ends:
cat catalog/summary.tsv
```

It is resumable: resubmit the same command after a timeout. Do not go above
`-j 4`. A 4-connection crawl on 2026-10-19 ended in ~15 minutes of TLS resets
from the server. If `logs/catalog_*.out` shows repeated retries, cancel, wait,
and resume with `-j 1`.

### 1. Plan

```bash
llb zinc3d-plan --catalog catalog --ids $P --out plan
cat plan/estimate.tsv     # findable parents and the db2 bill, before any download
```

`--ids` accepts an LLB `parents.tsv` or `library.tsv` (`zinc_id` column), or a
plain id list. Hashed non-ZINC names (`Z...`) go to `plan/rejected.txt`.

### 2. Sidecars

```bash
llb zinc3d-index groups --plan plan --catalog catalog --no-idlists
llb zinc3d-fetch plan/sidecar_urls.txt                # one SLURM job; logs/fetch_<id>.out
llb fetch-status plan/sidecar_urls.txt                # done at "-> retry.txt 0"
```

`--no-idlists` takes every catalogued chunk of every tranche a parent falls
in. For a large uniform draw that is the right call: the per-group id lists
would narrow almost nothing. For a small or concentrated id set (a few
hundred ids, a single bin), fetch the id lists first. They are ~17 bytes per
id and cut the sidecars to the groups that hold an id:

```bash
llb zinc3d-fetch plan/idlist_urls.txt
llb zinc3d-index groups --plan plan --catalog catalog   # also writes groups.tsv, no3d.txt
```

### 3. Exact index, then the db2s

```bash
llb zinc3d-index chunks --plan plan --catalog catalog
cat plan/coverage.tsv     # per bin: wanted, built, covered, chunks, GB, MB per parent
```

`plan/index.tsv` lists every db2 member ZINC built for a found parent.
`plan/db2_urls.txt` lists the chunks that hold them.

To cap the download below what the cells found, add `--budget-gb N`.
Chunks are then taken in seeded random order, so the selection stays
unbiased. `--greedy` covers more parents per GB, but favours small
(charged, sparse) chunks. Then fix the list under its own name and fetch:

```bash
cp plan/db2_urls.txt plan/db2_urls.G003.txt
llb zinc3d-fetch plan/db2_urls.G003.txt 4 -- --time=7-00:00:00
```

Watch scratch (`df -h .`; `lfs quota -h -u $USER /scratch`). Step 4 with
`--delete-archive` can run behind the download, so only the part not yet
extracted has to fit.

### 4. Extract

```bash
llb zinc3d-extract plan zinc_db2 --delete-archive     # array; logs/extract_<job>_<task>.out
# resubmit the same command after any failure or timeout: finished chunks are skipped
llb zinc3d-collect --plan plan --out zinc_db2         # -> zinc_db2/records.tsv + a loss check
```

```
zinc_db2/db2/H20/H20P200-N-faaaba.db2.gz   only the found parents' members: DOCK reads it as is
zinc_db2/records.tsv                       zinc_id protomer charge m_name n_records out_file
                                           offset bytes archive member sha1
```

`submit_extract.sh` pins `plan/index.tsv` and `db2_urls.txt` by sha256.
After re-planning, use a new OUT.

### 5. Keep

```bash
K=/nfs/turbo/umms-maom/johnbick/zzz.libs/zinc22_3d/h18_22_G003
mkdir -p $K && rsync -a zinc_db2 plan catalog/summary.tsv catalog/catalog.mode $K/
```

The scratch mirror (sidecars, any archives not deleted) can go. A larger G
later needs the sidecars again, so keep them if that is likely.

### More parents later

Raise G in the same catalog directory with the same seed. Only the new cells
are listed, and fetching skips every file already on disk. Then rerun steps
1-4. Extraction skips finished chunks only within one OUT, so extract into a
new OUT, or keep the old one and accept re-extraction.

```bash
llb zinc3d-catalog --heavy 18-22 --fraction 0.06 --out catalog -j 2   # under sbatch, as above
```

## Looking up one LLB parent

`records.tsv` is the join table:

```bash
grep -P '^ZINCi500000JIrMs\t' zinc_db2/records.tsv | cut -f1-8
```

```python
import csv, gzip
def zinc_db2(zid, out="zinc_db2"):
    """ZINC's db2 text for one ZINC id: one string per member (protomer)."""
    for r in csv.DictReader(open(f"{out}/records.tsv"), delimiter="\t"):
        if r["zinc_id"] == zid:
            with gzip.open(f"{out}/{r['out_file']}") as g:
                g.seek(int(r["offset"]))
                yield r, g.read(int(r["bytes"])).decode()
```

To read one member from an archive that is still on disk:

```bash
tar xzOf zinc22/zinc-22f/H20/H20P200/a/H20P200-N-faaaba.db2.tgz H20/H20P200/0A/MA/ZINCkq00000f0AMA.0.N.db2
```

A parent missing from `records.tsv` was dropped at one of these points:

| Where it is | Why |
|---|---|
| `plan/rejected.txt` | Not a ZINC id. |
| Absent from `plan/index.tsv` | Its cell was not sampled. With id lists: in `plan/no3d.txt` (no ZINC 3D at all) or `plan/listed_not_built.txt` (ZINC's build failed). |
| `plan/over_budget.txt` | Found, but left out by `--budget-gb`. |

## Before comparing

- **Docking names are truncated.** The db2 M line keeps the last 16 characters
  of `ZINC<id>.<protomer>`. `ZINCk50000040LmS.0` is written as
  `NCk50000040LmS.0`, and that is what DOCK reports. Join scores back through
  `records.tsv:m_name`, never by prefix.
- **There are several records per member and several members per id.** A
  member is one protomer at one charge, holding 1-3+ db2 records. In the H30
  test, 6 members held 16 records. An id can span charge groups:
  `ZINCu500000005KS` has `.0.L` and `.1.M`. Score per id = min over all its
  records, as for LLB's db2s.
- **ZINC protonates differently from LLB.** ZINC's protomers come from its own
  pipeline. The charge letter is the protomer's net charge: J=-4 … N=0 … Q=+3.
  LLB's protomers come from QupKake/MolGpKa. That difference is part of what
  the comparison measures.
- **Check the found subset before trusting it.** Compare its mix with the
  library's: bin from `coverage.tsv`, logP via `index.tsv:tranche` against
  `plan/wanted.tsv`. Cluster sampling is unbiased in expectation, but at small
  G a few huge cells can dominate. If the mix is off, raise G.

## Alternative: chunk-first

If same-molecule pairs need not be the production library's own parents,
turn it around. Sample ZINC chunks, take their ids and SMILES (the tranche dir
has `<tranche>-<charge>.<gen>.smi.gz` beside each id list) and build those
SMILES with LLB. One ~0.3 GB chunk holds ~5,000 molecules, so 20 chunks per bin
is ~6 GB for ~10^5 molecules per bin. That is ~10^4 times less download per
molecule.

The cost moves to LLB builds and to docking the LLB side again: ~1.1 core-h
per molecule over 3,150 pockets at the pilot's 1.28 s. Subsample to the
comparison's size. On the ZINC side, `catalog.py` → `index.py chunks
--all-members` → `extract.py` already does the work. Extracting every member
of a chunk only needs an `index.tsv` built from `members.tsv.gz`.

## Files

| File | What |
|---|---|
| `catalog.py` | Lists files.docking.org: chunks (db2 + sidecar sizes) and per-tranche id lists, for all or a seeded sample (`--fraction`) of (generation, tranche) cells. |
| `plan.py` | Ids → tranches → what to fetch, with an expected bill scaled to the catalogued fraction. |
| `index.py groups` | Picks the sidecars to fetch: every chunk of the touched tranches (`--no-idlists`), or only the groups whose id lists hold an id. `--fraction` subsamples chunks instead. |
| `index.py chunks` | Sidecar members → exact chunk per id (`index.tsv`), `db2_urls.txt` (optionally `--budget-gb`), `coverage.tsv`. |
| `extract.py` | Streams chunks: wanted members → one `db2.gz` per chunk, plus `records.tsv` with offsets. Sharded and resumable. `--collect` merges and checks for losses. |
| `submit_fetch.sh` / `fetch.sbatch` | One SLURM job: `../run.sh` plus up to 3 retry passes. The list is checksummed; never an array. |
| `submit_extract.sh` / `extract.sbatch` | The extraction array. The plan is checksummed. |
| `zid.py` | Id → bin/tranche, member parsing, the seeded `rank()` behind every fraction. |
| `estimate_2026-10-19/` | Inputs to the size tables: per-(generation, bin) id counts, crawled cells, per-generation calibration. |

Downloads go through `../fetch.sh` / `../run.sh` / `../status.sh`, which also
handle `.txt` id lists: valid = non-empty, ends in a newline, not HTML. Every
item in `../README.md`'s Gotchas still applies. Layout quirks found while
writing this:

- Chunk ids are 6 letters in most generations (`faaaba`) but 3 in zinc-22q
  (`qaa`).
- zinc-22e has id lists but almost no 3D.
- Python's default User-Agent gets HTTP 403 from the server. The scripts set
  their own.
