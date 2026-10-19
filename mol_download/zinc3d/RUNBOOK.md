# ZINC22's own 3D db2s for LLB's molecules: runbook

Goal: for every ZINC22 parent in an LLB library (H18-22 now, H18-H30 next),
get the db2 record(s) ZINC itself built, so the same molecules can be docked
from both sources and compared. Nothing here has been run on Great Lakes yet.
Duo means John runs every command. The tools were tested end to end on a laptop
against files.docking.org, on real H30 chunks (curl standing in for wget).

## The short version

- **The full H18-H30 set is out of reach: @FULL_TB@ TB** in @FULL_CHUNKS@ chunks
  (catalog of 2026-10-19, table below). Every bin is spread over several
  generations (zinc-22a … y), and zinc-22c (H30), d, f and others are huge.
- **So fetch only the chunks that hold LLB's ids.** A ZINC22 id gives its tranche
  (bin + logP) but not its generation, charge or chunk. Two cheap downloads
  find the exact chunk:
  1. per-(tranche, charge, generation) **id lists** (`H20P200-N.f.txt`, ~17 bytes per id), then
  2. each candidate chunk's **`.pdbqt.tgz` sidecar**. It holds the db2's exact
     member list (`H20/H20P200/0A/MA/ZINCkq00000f0AMA.0.K.pdbqt`) at ~1% of
     the db2's bytes.
- **The catch:** LLB's parents are a uniform draw from ZINC22 2D, so nearly every
  parent sits in a different ~100-300 MB chunk. Each covered parent costs a whole
  chunk, about **@MB_PER_ID@ MB per parent** (estimate below). For the 50k H18-22
  library that is **~@EXPECT_50K_TB@ TB** for every parent that has a ZINC 3D build.
  Not every parent has one: only **~@P3D_PCT@%** do (the rest are 2D-only in ZINC).
- **Recommendation:** fetch with a byte budget. `index.py chunks --budget-gb`
  picks chunks greedily by parents covered per byte. **500 GB covers ~@COVER_500@
  parents** of H18-22, enough for a same-molecule comparison, and takes ~@HOURS_500@
  of download. Spend more only if the comparison needs it. A cheaper design
  ("chunk-first", at the end) needs new LLB builds.

## Sizes (catalog of 2026-10-19, Apache's rounded sizes)

@SIZE_TABLE@

`db2_TB` is what files.docking.org holds for the bin across all generations.
`sidecar_GB` is the `.pdbqt.tgz` total, which is what an exact index of
*everything* would cost. `ids_M` counts the per-tranche id lists in millions.
Those lists run ~5% above what was actually built (zinc-22a H20P200 N: 10,163
listed, 9,670 built).

## Estimate for LLB-sized id sets

From `plan.py` on the 9,773 druglike_10k parents
(`llb_validation_2026-10/druglike_10k`), H18-22 and drawn the same way as the
production library (`sample_2d.py --total`). The 50k column scales by 5 under
the same model (it is not linear: chunks start to be shared):

@ESTIMATE_TABLE@

`per_id` is MB of db2 per parent with a 3D build. `p3d` is estimated from
the 2D tranche sizes. `index.py groups` measures it exactly in step 3, so the
real numbers are known after ~@IDLIST_GB@ GB of id lists, before any big
download.

Time: one connection to files.docking.org did ~14 MB/s (sample_2d.py,
2026-09). At 4 connections that is 1 TB in ~5 h if the server scales and
~20 h if it does not. Extraction is gzip-bound, ~50-100 MB/s of archive per
core, so 1 TB is ~5 core-h spread over the array.

## One-time setup

```bash
# config/hpc.env (or $LLB_CONFIG): the account is maom2 since 2026-10-16
grep -E '^(LLB_ACCOUNT|LLB_PARTITION|WGET|PARALLEL|PREP_ENV|QUPKAKE_ENV)=' $(llb root)/config/hpc.env
#   LLB_ACCOUNT=maom2   WGET=wget   PARALLEL=parallel (module load or conda)
#   PREP_ENV (or QUPKAKE_ENV) must have python >= 3.8; these scripts are stdlib only

# compute nodes must reach files.docking.org (fetch.sbatch checks again)
srun -A maom2 -p standard -t 5 curl -sI https://files.docking.org/zinc22/ | head -1
```

The code is on branch `feat/zinc22-db2` until merged:

```bash
cd /nfs/turbo/umms-maom/johnbick/opt/ligand_library_builder
git fetch && git checkout feat/zinc22-db2
```

## Run

One landing directory on scratch holds the mirror (`zinc22/...`), the plan and
the logs. Every command below runs from it.

```bash
R=/scratch/maom_root/maom99/johnbick/2026_10_zinc3d
mkdir -p $R && cd $R
```

### 0. Catalog (login node, ~1-1.5 h, listings only)

```bash
llb zinc3d-catalog --heavy 18-30 --out catalog -j 4     # resumable: rerun to continue
cat catalog/summary.tsv
```

Or start from the 2026-10-19 snapshot committed beside this file and skip
the crawl:

```bash
mkdir -p catalog && for f in chunks idlists; do
  gzip -dc $(llb root)/mol_download/zinc3d/catalog_2026-10-19/$f.tsv.gz > catalog/$f.tsv; done
```

### 1. Plan

```bash
P=/scratch/maom_root/maom99/johnbick/2026_10_08_build_small_mols/<library run>/parents/parents.tsv
llb zinc3d-plan --catalog catalog --ids $P --out plan
cat plan/estimate.tsv        # bytes before anything is fetched
```

`--ids` accepts an LLB `parents.tsv` or `library.tsv` (`zinc_id` column) or a
plain id list. Hashed non-ZINC names (`Z...`) go to `plan/rejected.txt`.

### 2. Id lists (~@IDLIST_GB@ GB for H18-22)

```bash
llb zinc3d-fetch plan/idlist_urls.txt          # one SLURM job; logs/fetch_<jobid>.out
# when it ends:
llb fetch-status plan/idlist_urls.txt          # "-> retry.txt 0"; permanent 404s are fine
```

### 3. Groups, then sidecars

```bash
llb zinc3d-index groups --plan plan --catalog catalog --root .
#  -> "N of M wanted ids are in a ZINC 3D id list"  (the measured p3d)
#  -> plan/sidecar_urls.txt and its size
llb zinc3d-fetch plan/sidecar_urls.txt
```

### 4. Exact index, then choose what to download

```bash
llb zinc3d-index chunks --plan plan --catalog catalog --root .               # everything
llb zinc3d-index chunks --plan plan --catalog catalog --root . --budget-gb 500
cat plan/coverage.tsv        # per bin: wanted, listed, built, covered, chunks, GB, MB/id
```

Rerun `chunks` with a different `--budget-gb` as often as needed. It only
reads the sidecars. `plan/index.tsv` always lists every built record of every
wanted id. `plan/db2_urls.txt` is the budgeted selection, and
`plan/over_budget.txt` holds the ids left out.

Once happy, copy the selection under a name that will not change, so the fetch's
checksum guard means something:

```bash
cp plan/db2_urls.txt plan/db2_urls.500G.txt
llb zinc3d-fetch plan/db2_urls.500G.txt 4 -- --time=7-00:00:00
```

Watch scratch while it runs (`df -h .`; `lfs quota -h -u $USER /scratch`).
The extraction in step 5 can run with `--delete-archive` behind the download,
so only the not-yet-extracted part has to fit.

### 5. Extract the wanted records

```bash
llb zinc3d-extract plan zinc_db2 --delete-archive      # array; logs/extract_<job>_<task>.out
# resubmit the same command after any failure or timeout; done chunks are skipped
llb zinc3d-collect --plan plan --out zinc_db2          # -> zinc_db2/records.tsv
```

Result:

```
zinc_db2/db2/H20/H20P200-N-faaaba.db2.gz   DOCK-ready: only the wanted records of that chunk
zinc_db2/records.tsv                       zinc_id protomer charge m_name n_records out_file
                                           offset bytes archive member sha1
                                           (offset/bytes: the member inside the unzipped out_file)
```

Extraction runs after `index.py chunks` and reads that run's `db2_urls.txt`.
If the selection changes, use a new `OUT`: `submit_extract.sh` refuses a
changed plan.

### 6. Keep

```bash
K=/nfs/turbo/umms-maom/johnbick/zzz.libs/zinc22_3d/h18_22_$(date +%Y%m%d)
mkdir -p $K && rsync -a zinc_db2 plan catalog/summary.tsv $K/
```

The mirror's `.pdbqt.tgz` sidecars (a few GB) can stay on scratch, or be
deleted. Re-indexing needs them, but they are quick to refetch.

## Looking up one LLB parent

`records.tsv` is the join table between LLB parents and ZINC's records:

```bash
grep -P '^ZINCi500000JIrMs\t' zinc_db2/records.tsv | cut -f1-6
# zinc_id  protomer  charge  m_name  n_records  out_file
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

To read one member straight from a still-present archive without extracting
the rest:

```bash
tar xzOf zinc22/zinc-22f/H20/H20P200/a/H20P200-N-faaaba.db2.tgz H20/H20P200/0A/MA/ZINCkq00000f0AMA.0.N.db2
```

Not in `records.tsv`? Check where the id was dropped, in this order:
- `plan/rejected.txt`: not a ZINC id.
- `plan/no3d.txt`: in no 3D id list. ZINC has no 3D build of this molecule.
- `plan/listed_not_built.txt`: listed, but ZINC's build failed.
- `plan/over_budget.txt`: built, but left out by `--budget-gb`.

## Things to know before comparing

- **Docking names are truncated.** The db2 M line keeps the last 16 characters
  of `ZINC<id>.<protomer>`. `ZINCk50000040LmS.0` is written as
  `NCk50000040LmS.0`, which is what DOCK reports. Join scores back through
  `records.tsv:m_name`, never by prefix.
- **Several records per member, several members per id.** A member is one
  protomer (`.0`, `.1`, …) at one charge. It holds 1-3+ db2 records, and an id
  can have members in several charge groups (`ZINCu500000005KS` has `.0.L` and
  `.1.M`). Score per id = min over all its records, the same rule as for
  LLB's db2s.
- **ZINC protonates differently from LLB.** ZINC's protomers come from its own
  pipeline (charge letter = the protomer's net charge). LLB's come from
  QupKake/MolGpKa. That difference is part of what the comparison measures,
  so do not "fix" it.
- **One chunk is one tranche × charge × registration batch.** A budgeted
  selection is biased toward parents that share a chunk with others (greedy
  by ids per byte), which happens more in small tranches. Check that the
  covered parents' bin/logP mix still matches the library's
  (`coverage.tsv`, and `index.tsv`'s tranche column against
  `plan/wanted.tsv`). If it drifts, rerun with a bigger budget or plain
  random subsampling.

## Alternative: chunk-first (cheapest per molecule)

If thousands of same-molecule pairs are enough and need not be the
production library's own parents, turn the scheme around. Take a stratified
random sample of ZINC chunks per bin, extract their ids and SMILES (each
tranche dir has `<tranche>-<charge>.<gen>.smi.gz` beside the id list), and
run those SMILES through LLB. One ~250 MB chunk gives ~5,000 molecules, so
20 chunks per bin for H18-H30 is ~65 GB for ~10^6 molecules, against
~@MB_PER_ID@ MB per molecule above. The cost moves to LLB builds and extra
docking of the LLB side (~1.1 core-h per molecule over 3,150 pockets at the
pilot's 1.28 s), so subsample to the size the comparison needs.
`catalog/chunks.tsv` plus `extract.py`, with an `index.tsv` written for all
members (`index.py chunks --all-members` gives `plan/members.tsv.gz`),
covers the ZINC side.

## Files

| File | What |
|---|---|
| `catalog.py` | Lists files.docking.org: every chunk under the given bins, with db2 and sidecar sizes, and every per-tranche id list. |
| `plan.py` | Ids → tranches → id lists to fetch, plus a cost estimate. |
| `index.py groups` | Id lists → (generation, tranche, charge) groups holding each id → sidecar URLs. |
| `index.py chunks` | Sidecar members → exact chunk per id (`index.tsv`), budgeted `db2_urls.txt`, `coverage.tsv`. |
| `extract.py` | Streams chunks; wanted members → `db2.gz` per chunk + `records.tsv`. Array-sharded, resumable. |
| `submit_fetch.sh` / `fetch.sbatch` | One SLURM job: `run.sh` + up to 3 retry passes. The list is checksummed. |
| `submit_extract.sh` / `extract.sbatch` | Extraction array. The plan is checksummed. |
| `zid.py` | ZINC22 id → bin/tranche decoding, member-name parsing, id-file reading. |
| `catalog_2026-10-19/` | Snapshot of the H18-H30 catalog (gzipped TSVs + summary). |

Downloads go through the existing `../fetch.sh` / `../run.sh` /
`../status.sh`. `fetch.sh` now also validates `.txt` id lists: non-empty,
ends in a newline, not HTML. Everything in `../README.md`'s Gotchas still
applies.
