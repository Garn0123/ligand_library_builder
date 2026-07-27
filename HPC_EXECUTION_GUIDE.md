# db2pipe — HPC Execution Guide (SLURM)

Exact, copy-pasteable steps to compile a large tree of ZINC `.db2.gz` files into
fixed-size, load-balanced chunks on a SLURM cluster, then label them for docking.

This guide covers the **parallel (SLURM) pipeline**. For the single-process path
(`01`–`04`), see `README.md`. Read the "db2 record model" and "Known caveats"
sections of `README.md` once before your first production run — this guide
assumes them.

The whole job is five scripts run as four dependency-chained SLURM jobs, followed
by two serial labeling steps:

```
make_shards.py  →  p1_collect.py  →  p2_assemble.py  →  p3_finalize.py
   (plan, 1)        (collect, S)       (assemble, N)       (final, 1)
                                                                │
                                              chunks + manifest.tsv in OUT_DIR
                                                                │
                                          02_label.py  →  03_apply.py
                                          (decide ids)     (rewrite chunks)
```

Result: `OUT_DIR/chunk_00001.db2.gz … chunk_000NN.db2.gz` + `OUT_DIR/manifest.tsv`,
each chunk a whole number of records, all chunks roughly equal in docking cost.

---

## 0. Before you start — preflight checklist

Run these from a login node in the pipeline directory (the folder containing
`submit.slurm`).

**1. Python 3.9+ is on the compute nodes.** Stdlib only, no pip installs needed.
The scripts are launched with `python3 -u` by `submit.slurm` (unbuffered stderr
so progress shows up live).

```bash
python3 --version        # need >= 3.9
```

If the cluster needs a module for a modern Python (e.g. `module load python`),
load it in your shell **and** add the same `module load` line near the top of
`submit.slurm` (after the shebang) so the compute nodes get it too.

**2. `$SCRATCH` is set and has room.** Intermediates total roughly **1× the input
library size**, freed incrementally during assembly. For a ~500 GB library,
budget ~500 GB free on scratch. **Do not put `WORK_DIR` on your home directory** —
it will be too small and too slow, and home quotas will kill the run mid-array.

```bash
echo "$SCRATCH"          # must be non-empty
df -h "$SCRATCH"         # confirm free space >= input size
```

**3. Know your account, partition, and input path.**

```bash
sacctmgr -n show assoc user=$USER format=account%30,partition%20   # your accounts/partitions
ls -d /path/to/ZINC_sync/published/3D                              # confirm the input tree
```

**4. (Recommended) Dry-run the planner to validate discovery and size the job.**
`make_shards.py` just walks the tree and stats file sizes — it's cheap and tells
you the total library size, which you need to choose `BINS` (step 1). Run it into
a throwaway work dir:

```bash
python3 make_shards.py -i /path/to/ZINC_sync/published/3D -o "$SCRATCH/db2_probe" -S 200
```

Read its summary: `found N files`, `total input: X GB`, and the shard-balance
spread. If it finds 0 files, your `INPUT_DIR` or `--suffix` is wrong. You can
delete `$SCRATCH/db2_probe` afterward; `submit.slurm` re-runs the planner itself.

---

## 1. Choose your parameters

You set six numbers. Only two need real thought — `BINS` and `SHARDS`.

### BINS — the number of finished chunks (and the array width for assembly)

`BINS = ceil(estimated_total_molecules / target_chunk_size)`.

- Target chunk size of **50,000 molecules** is the default assumption.
- **You must estimate the molecule count up front.** The parallel pipeline
  produces *exactly* `BINS` chunks; chunk size falls out as
  `total_molecules / BINS`. An estimate that's off by 20% just makes chunks 20%
  bigger or smaller — harmless. It does **not** need to be exact.
- If you don't know the count, estimate it from library size. ZINC 3D db2 runs
  ~59 KB/molecule uncompressed; from the planner's `total input: X GB`
  (compressed) a rough live-count is available from the ZINC tranche catalog, or
  just target chunk *size* by bytes. For the reference 8.5 M-molecule library,
  `BINS = 8.5e6 / 50000 ≈ 170`.

### SHARDS — the collect array width (how many parallel readers)

- This is the parallelism of the expensive phase. **~200 is a good default.**
- Intermediates are `SHARDS × BINS` part files (200 × 170 = 34,000). That's fine
  on Lustre/GPFS. Be cautious before pushing `SHARDS` far past a few hundred —
  file count grows as `SHARDS × BINS`.
- More shards = shorter per-task walls but more scheduler churn and more small
  files. Match it to how many array tasks your partition will actually run at
  once.

### The other four

| Var | Default | Meaning |
|-----|---------|---------|
| `MODE` | `stride` | `stride` = systematic sampling, needs no cost model, ~1.0× count spread. `greedy` = least-loaded-bin against `WEIGHT`. **Use `stride` unless you specifically want cost-weighted balancing.** |
| `WEIGHT` | `count` | Cost metric for `greedy` only (ignored by `stride`). `count`, `bytes`, or `lines:C` (weight by conformer count — the best proxy for docking time, **but confirm your db2 uses `C` for conformer lines first**). |
| `COLLECT_TIME` | `04:00:00` | Wall limit per collect task. Generous for ~200 shards; each task handles `total/SHARDS` molecules at ≥250 mol/s. Lower it to improve queue priority once you've seen a real run. |
| `ASSEMBLE_TIME` | `00:30:00` | Wall limit per assemble task. Assembly is a pure byte copy at disk speed; 30 min is ample. |

**Why `stride` needs no shuffle step:** each collect task strides its own records
across all `BINS` bins (offset by shard id), so every bin gets an even slice of
every shard. The chunks come out globally balanced with zero cross-task
communication, even though each shard covers only one region of the size-sorted
tree.

---

## 2. Configure `submit.slurm`

Edit the config block at the top of `submit.slurm` (lines ~15–30). Everything you
must change:

```bash
INPUT_DIR=/path/to/ZINC_sync/published/3D   # the tree of .db2.gz files
WORK_DIR=$SCRATCH/db2work                    # intermediates — MUST be on scratch, ~1x input size
OUT_DIR=$SCRATCH/db2chunks                   # finished chunks + manifest.tsv

SHARDS=200                                   # collect array width  (from step 1)
BINS=170                                      # finished chunk count (from step 1)
MODE=stride                                   # stride | greedy
WEIGHT=count                                  # count | bytes | lines:C   (greedy only)

ACCOUNT=your_account                          # <-- set to your real account
PARTITION=standard                            # <-- set to your real partition
CPUS=1                                        # each task is single-threaded; leave at 1
MEM=4G                                        # per-task memory; 4G is plenty for stride
COLLECT_TIME=04:00:00
ASSEMBLE_TIME=00:30:00
```

`PIPE` auto-detects the script directory — leave it. Do not point `WORK_DIR` at
`$HOME`.

**Concurrency throttle (optional).** The collect submit line uses
`--array=0-$((SHARDS-1))%$SHARDS`, i.e. `%SHARDS` = no throttle, all shards may
run at once. If your partition caps concurrent array tasks (or you want to be a
good neighbor), lower the number after `%`, e.g. edit it to `%50` to run at most
50 collect tasks simultaneously. This changes wall time, not correctness.

---

## 3. Submit

```bash
bash submit.slurm
```

This submits four jobs, each waiting on the previous with `--dependency=afterok`,
and prints their job IDs:

```
plan     : 1234567
collect  : 1234568  (array 0-199)
assemble : 1234569  (array 0-169)
final    : 1234570
```

`afterok` on an array waits for **every** task in it. So if any single collect
task fails, assembly does **not** start (this is deliberate — it prevents
silently producing short chunks). See step 5 for recovery.

Nothing runs on your login node; `bash submit.slurm` only submits and exits.

---

## 4. Monitor

```bash
squeue -u $USER                              # all four jobs / array tasks
squeue -u $USER -t RUNNING                   # what's actually running now

# live progress of the expensive phase (mol/s, file counts):
tail -f $SCRATCH/db2work/logs/collect_0.err

# after it finishes, per-phase logs:
ls  $SCRATCH/db2work/logs/
cat $SCRATCH/db2work/logs/final.out          # global chunk-balance summary
```

Every script writes progress and summaries to **stderr** (`*.err`); real results
also go there. Each phase's logs:

- `logs/plan.{out,err}` — file count, total GB, shard-balance spread.
- `logs/collect_<task>.{out,err}` — per-shard molecule count and mol/s; truncated
  files reported here (a truncated file makes that task exit nonzero).
- `logs/assemble_<task>.{out,err}` — parts joined, MB written per chunk.
- `logs/final.{out,err}` — **the run summary**: total molecules and the
  `molecules min/max spread` across chunks. A spread near `1.00x` means balancing
  worked.

When `final` completes, your deliverables are in `OUT_DIR`:

```bash
ls   $SCRATCH/db2chunks/                      # chunk_00001.db2.gz ... + manifest.tsv
head $SCRATCH/db2chunks/manifest.tsv
```

---

## 5. If something fails

**Collect tasks are idempotent** — re-running a shard produces a byte-identical
shard manifest — so recovery is just "rerun the failed shards, then rerun
assemble and final." `submit.slurm` chains everything in one shot, so after a
partial failure you resubmit the tail by hand.

**1. Find which collect tasks failed:**

```bash
sacct -j <COLLECT_JOBID> --format=JobID,State,ExitCode,Elapsed | grep -v COMPLETED
```

Common causes: a truncated input file (the task lists it and exits 1), an
out-of-time task (raise `COLLECT_TIME`), or a node/scratch hiccup.

**2. Rerun only the failed shards** (say tasks 3, 17, 102). Reuse the exact same
`-w`, `-N`, `--mode`, `--weight` you configured:

```bash
sbatch -A your_account -p standard --cpus-per-task=1 --mem=4G --time=04:00:00 \
  --job-name=db2collect_rerun --array=3,17,102 \
  -o "$SCRATCH/db2work/logs/collect_%a.out" -e "$SCRATCH/db2work/logs/collect_%a.err" \
  --wrap "python3 -u $PWD/p1_collect.py -w '$SCRATCH/db2work' -s \$SLURM_ARRAY_TASK_ID \
          -N 170 --mode stride --weight count"
```

**3. Resubmit assemble and final** (chained on the rerun):

```bash
JOB_A=$(sbatch --parsable -A your_account -p standard --cpus-per-task=1 --mem=4G \
  --time=00:30:00 --job-name=db2assemble --array=0-169 \
  --dependency=afterok:<COLLECT_RERUN_JOBID> \
  -o "$SCRATCH/db2work/logs/assemble_%a.out" -e "$SCRATCH/db2work/logs/assemble_%a.err" \
  --wrap "python3 -u $PWD/p2_assemble.py -w '$SCRATCH/db2work' -b \$SLURM_ARRAY_TASK_ID \
          -o '$SCRATCH/db2chunks'")

sbatch -A your_account -p standard --cpus-per-task=1 --mem=8G --time=01:00:00 \
  --job-name=db2final --dependency=afterok:$JOB_A \
  -o "$SCRATCH/db2work/logs/final.out" -e "$SCRATCH/db2work/logs/final.err" \
  --wrap "python3 -u $PWD/p3_finalize.py -w '$SCRATCH/db2work' -o '$SCRATCH/db2chunks'"
```

Replace `170`/`169`, account, partition, and paths with your values. Assembly
deletes each bin's parts after a successful join, so **don't** rerun an assemble
task that already succeeded — its parts are gone. If you need to rerun assembly,
rerun the whole collect phase first, or pass `--keep-parts` to `p2_assemble.py`
on the first run.

---

## 6. Label the chunks (serial phases `02` and `03`)

Labeling disambiguates molecule ids that legitimately repeat (protonation states,
tautomers, stereoisomers). It's a two-step, position-addressed rewrite and is the
**same in both pipelines**. It runs after `OUT_DIR/manifest.tsv` exists.

### 6a. Decide new ids — `02_label.py`

For a large library, precompute the duplicate-id set on disk (cheap, streaming)
and pass it in, rather than making `02` hold every unique id in memory:

```bash
tail -n +2 $SCRATCH/db2chunks/manifest.tsv | cut -f5 | sort | uniq -d \
  > $SCRATCH/db2work/dups.txt
```

Then run `02` as its **own sbatch job** (don't run an 8.5 M-row pass on a login
node). Give the sort/label step real memory and time:

```bash
sbatch -A your_account -p standard --cpus-per-task=1 --mem=8G --time=02:00:00 \
  --job-name=db2label -o $SCRATCH/db2work/logs/label.out -e $SCRATCH/db2work/logs/label.err \
  --wrap "python3 -u $PWD/02_label.py -m $SCRATCH/db2chunks/manifest.tsv \
          -o $SCRATCH/db2work/labels.tsv --mode duplicates \
          --dup-ids $SCRATCH/db2work/dups.txt --max-id-len 40"
```

`labels.tsv` contains **only** the records whose id changes. `--max-id-len 40`
warns if any new id gets long enough to worry about the column-shift caveat
(below). Default mode `duplicates` leaves unique molecules with their bare
catalog id.

### 6b. Apply the labels — `03_apply.py`

`03` locates each record by `(chunk, chunk_idx)` and, before rewriting, verifies
the id at that position matches what `labels.tsv` expects. A mismatch (e.g. the
manifest is stale) aborts the run with a nonzero exit instead of corrupting data.
Writes are `.tmp` + `os.replace`, so each chunk is all-or-nothing.

Run it as a **job array, one chunk per task**, with `--copy-unedited` so the
output directory is a **complete** labeled tree (every chunk present, edited or
not) rather than a sparse overlay:

```bash
sbatch -A your_account -p standard --cpus-per-task=1 --mem=4G --time=01:00:00 \
  --job-name=db2apply --array=1-170 \
  -o $SCRATCH/db2work/logs/apply_%a.out -e $SCRATCH/db2work/logs/apply_%a.err \
  --wrap 'CH=$(printf "chunk_%05d.db2.gz" $SLURM_ARRAY_TASK_ID); \
          python3 -u '"$PWD"'/03_apply.py -c '"$SCRATCH"'/db2chunks \
          -L '"$SCRATCH"'/db2work/labels.tsv -o '"$SCRATCH"'/db2chunks_labelled \
          --only $CH --copy-unedited'
```

Set `--array=1-170` to `1-BINS`. Chunks are 1-based, 5-digit
(`chunk_00001.db2.gz`).

- **`--copy-unedited`** copies chunks with no duplicates through to the output dir
  too, giving one self-contained `db2chunks_labelled/` for docking. Omit it if you
  prefer a sparse overlay (only changed chunks written; read unchanged ones from
  the original dir) to save I/O and space.
- **Labeling is not idempotent** — never run `03` over already-labeled output, or
  you stack a second suffix (`ZINC…_1_1`). Keep labeled and unlabeled trees
  separate.

The finished, docking-ready library is `OUT_DIR_labelled/` (e.g.
`$SCRATCH/db2chunks_labelled/`) plus the `manifest.tsv` in `OUT_DIR`.

---

## 7. Verify a run

Do this once end-to-end before trusting a full production pass, and re-check
after any code change. **Critically: before committing to a full labeled pass,
push one labeled chunk through your real docking pipeline** — appending an id
suffix shifts every column to its right on the `M` line, which is fine for a
whitespace-splitting reader but silent corruption for a fixed-column one (see
`README.md` "Column shift on relabel").

```bash
CH=$SCRATCH/db2chunks

# 1. Every chunk is whole records: first line starts with M, last with E.
for f in $CH/chunk_*.db2.gz; do
  printf "%s first=%s last=%s\n" "$(basename "$f")" \
    "$(zcat "$f" | head -1 | cut -c1)" "$(zcat "$f" | tail -1 | cut -c1)"
done | grep -v 'first=M last=E' && echo "!! bad boundary above" || echo "OK: all chunks M..E"

# 2. Nothing lost or duplicated vs the source tree (identical id multiset).
find /path/to/3D -name '*.db2.gz' | sort | xargs zcat | grep -o 'ZINC[0-9]*' | sort > /tmp/src.txt
zcat $CH/chunk_*.db2.gz | grep -o 'ZINC[0-9]*' | sort > /tmp/out.txt
diff -q /tmp/src.txt /tmp/out.txt && echo "OK: identical id multiset"

# 3. Chunk balance — read the final job's summary; spread should be near 1.00x.
grep -A2 'chunk balance' $SCRATCH/db2work/logs/final.out
```

Use the **record (E-terminator) count** as ground truth for molecule counts;
`grep -c ZINC` overcounts (a record can carry the ZINC token on more than one `M`
line). The scripts report both and warn when they disagree.

---

## Quick reference — full run, start to finish

```bash
# 0. preflight
python3 --version                                   # >= 3.9
echo "$SCRATCH"; df -h "$SCRATCH"                    # scratch has >= input-size free
python3 make_shards.py -i /path/to/3D -o "$SCRATCH/db2_probe" -S 200   # confirm file count + total GB

# 1-3. edit submit.slurm config block, then:
bash submit.slurm

# 4. watch
squeue -u $USER
tail -f $SCRATCH/db2work/logs/collect_0.err
cat     $SCRATCH/db2work/logs/final.out             # when done

# 6. label
tail -n +2 $SCRATCH/db2chunks/manifest.tsv | cut -f5 | sort | uniq -d > $SCRATCH/db2work/dups.txt
sbatch ... 02_label.py  ... --dup-ids $SCRATCH/db2work/dups.txt        # see step 6a
sbatch ... 03_apply.py  ... --array=1-BINS --only $CH --copy-unedited  # see step 6b

# 7. verify (see step 7), and push ONE labeled chunk through docking first.
```

---

## Sizing & resource cheat sheet

| Resource | Rule of thumb | Notes |
|----------|--------------|-------|
| `BINS` | `ceil(total_mols / 50000)` | = number of finished chunks. 8.5 M → ~170. Estimate is fine. |
| `SHARDS` | ~200 | Collect parallelism. Intermediates = `SHARDS × BINS` files. |
| Scratch space | ~1× input size | Parts freed incrementally by assembly. Never use `$HOME`. |
| Per-task memory | 4 G (collect/assemble/apply), 8 G (final/label) | `stride` mode is lightweight; `greedy` and `02` hold more. |
| Open files | `BINS + 32` per collect task | 170 bins ≈ 202 FDs, under the usual 1024. Thousands of bins → raise `ulimit -n`. |
| Collect wall | ≥ `total/SHARDS ÷ 250 mol/s` | Default 4 h is generous for ~200 shards; tighten after a real run. |

---

## Gotchas that bite in production

- **`WORK_DIR` on `$HOME`** — fills the quota, kills the array mid-run. Scratch only.
- **One failed collect task blocks all of assembly** (`afterok`). That's by
  design; recover with step 5, don't work around it.
- **Don't rerun a succeeded `p2_assemble` task** — it deleted its input parts.
  Use `--keep-parts` on the first run if you expect to reassemble.
- **`03` is not idempotent** — labeling twice double-suffixes ids. Separate trees.
- **Column shift on relabel is unverified for your reader** — validate one labeled
  chunk in real docking before a full pass. `--max-id-len` warns on long ids.
- **The fast (`p1`) parser assumes no non-terminator line begins with `E`.** Not
  yet confirmed against real ZINC db2 output (`README.md` "Open questions"). Run
  verification check #2 (id-multiset diff) against the source tree on your first
  real library to confirm nothing was mis-split.
- **`--weight lines:C` assumes `C` marks conformer lines** — confirm the letter in
  your db2 format before relying on cost-weighted (`greedy`) balancing.
```
