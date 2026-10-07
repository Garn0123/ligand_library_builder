# ZINC SMILES → protomers → contract names → db2

The path from ZINC SMILES to chunked db2, in the order it runs. Each step writes a
file the next one checks, so a failure stops at the step that caused it rather
than surfacing as a wrong row in the matrix.

```
mol_download/run.sh + sample_2d  ZINC22 2D on disk, uniform N/bin     samples/H*.smi
run_qupkake/prepare_parents.py   validate ids, neutralize, dedupe      parents/shards/*.smi
run_qupkake/slurm/qupkake_array  QupKake + protomers, one task/shard   protomers/shard_*/
run_qupkake/assign_names.py      union -> stereo closure -> names      library/library.smi
stage4/validate_names.py (DRAP)  the naming contract, on the union      gate
run_qupkake/check_db2_stereo.py  db2_converter will not rename         gate (db2_converter env)
build_ligand                     conformers -> db2                     (existing)
run_qupkake/verify_db2_build.py  names + charges survived into db2     gate
mol_compiler/parallel            db2 -> balanced chunks + manifest     docking input
```

## 0. Point the repo at your software (once per machine)

```bash
cp config/hpc.env.example config/hpc.env      # gitignored; edit it
```

Every step reads tool locations from `config/hpc.env`: conda env name or
prefix, an optional executable path, and an optional setup hook run *before*
activation (`module purge`, `module load ...`, `source ~/env/bcl.sh`). No
script hard-codes a path. Run a step in a tool's environment with:

```bash
common/with_env.sh PREP  python run_qupkake/prepare_parents.py ...   # RDKit steps
common/with_env.sh DB2C  python run_qupkake/check_db2_stereo.py ...  # db2_converter env
```

All commands below run from the repo root.

### Code in one place, runs anywhere

Keep the checkout in `opt/` and never run inside it. `bin/llb` runs any step,
in the right environment, against your **current** directory:

```bash
export PATH=/nfs/turbo/umms-maom/johnbick/opt/ligand_library_builder/bin:$PATH   # ~/.bashrc
export LLB_CONFIG=~/.config/llb/hpc.env      # optional: config outside the checkout too
mkdir -p /scratch/.../ladder_h14_h28 && cd /scratch/.../ladder_h14_h28
llb help                                      # every step, in pipeline order
```

Every output path is relative to where you run it, and SLURM logs go to
`./logs`, so one run directory holds one run and the checkout stays read-only.
`llb root` prints which checkout you're using. Each db2 build records its
commit, uncommitted changes and config in `provenance.<jobid>.txt`.

| step | `llb` | script |
| --- | --- | --- |
| URL list | `llb print-urls` / `llb extract-urls` | `mol_download/sample_2d.py` / `extract_urls.sh` |
| download | `llb fetch`, `llb fetch-verify`, `llb fetch-status` | `mol_download/run.sh`, `verify.sh`, `status.sh` |
| sample | `llb sample` | `mol_download/sample_2d.py` |
| parents | `llb parents` | `run_qupkake/prepare_parents.py` |
| QupKake | `llb qupkake`, `llb submit-qupkake` | `qupkake_protomers.py`, `submit_qupkake.sh` |
| names | `llb names`, `llb check-stereo` | `assign_names.py`, `check_db2_stereo.py` |
| db2 build | `llb submit-db2`, `llb verify-db2` | `submit_db2.sh` + `slurm/db2_array.sbatch`, `verify_db2_build.py` |
| mol2 | `llb mol2` | `db2_to_mol2.py` |
| chunk | `llb chunk` | `mol_compiler/parallel/submit_chunks.sh` → `submit.slurm` |

The commands below show the scripts by path; `llb <step>` is the same thing.

## 1. Download and sample

**Production: download the bins, then sample uniformly from disk.** Streaming a
prefix of each remote file is biased. The first few thousand lines of a tranche
are one registration batch, with fewer scaffolds (1,479 vs 1,654 per 3,600) and
more near-neighbours (NN Tanimoto 0.469 vs 0.432) than a uniform draw. Numbers are
in `sample_2d.py`'s docstring.

```bash
python3 mol_download/sample_2d.py --heavy 14-28 --print-urls > urls_2d.txt  # 904 files, ~960 GB
# or, from a script CartBlanche22 gave you (curl/wget/PowerShell; don't run it as-is):
#   mol_download/extract_urls.sh zinc22-2D-download.curl > urls_2d.txt
cd $ZINC_DATA                                     # scratch; files land under $PWD/zinc22/
$LLB/mol_download/run.sh $LLB/urls_2d.txt 4       # DTN / tmux, ~5-19 h
$LLB/mol_download/verify.sh --purge
$LLB/mol_download/status.sh $LLB/urls_2d.txt      # re-run run.sh retry.txt 2 until retry = 0
cd $LLB
python3 mol_download/sample_2d.py --heavy 14-28 --local $ZINC_DATA/zinc22 \
    --urls urls_2d.txt --outdir samples
```

(`$LLB` = this repo. `run.sh` takes `WGET` and `PARALLEL` from the config.)

`--local` is one reservoir pass over every file in the bin: every molecule
equally likely, the logP mix right by construction, about 30 min for H28. With
`--urls` it refuses to sample a bin that is missing files.

**Pilot only:** without `--local`, the sampler streams a few MB per bin. Use it to
prove the plumbing end to end, never for results.

## 2. Parents

```bash
common/with_env.sh PREP python run_qupkake/prepare_parents.py samples/H*.smi \
    -o parents --shard-size 250
```

Neutralizes with the same `standardize` QupKake's wrapper uses, and dedupes by
name and by neutral structure. Also writes `preflight.smi`: 20 parents, largest
bins first.

**Your own SMILES, no sampling:** skip step 1 and pass the files straight in,
plain or gzipped, any names:

```bash
common/with_env.sh PREP python run_qupkake/prepare_parents.py my_ligands.smi.gz -o parents
```

Names are handled automatically by their shape (`run_qupkake/llb_ids.py`):

| Input name | Becomes |
| --- | --- |
| `ZINC` + 12 (`ZINCh10000007sNS`) | kept; base ID = its last 12 characters |
| `ZINC12345` (old, unpadded) | padded to `ZINC000000012345`, then as above |
| `ZINC…` in any other shape (`ZINC000012345678.0`) | **rejected** `malformed_zinc_id`: a damaged ID, not a new molecule |
| anything else (`CHEMBL25`, `my(ligand)/2`) | base ID `Z` + 11 characters of its sha256; travels as `LLB<base>` |

No ZINC base ID starts with `O`–`Z` (ZINC22's first character is the heavy-atom
count, at most H49 = `N`; ZINC20's is a digit), so hashed and ZINC base IDs can't
collide. Your original name is kept as `input_name` in `parents.tsv` and
`library.tsv`. When one structure appears under two names, the ZINC ID wins.

## 3. QupKake

Preflight first, interactively, on real hardware. QupKake's per-molecule time
isn't known yet, and xtb cost grows with size:

```bash
srun -A maom99 -p standard -c 8 --mem 16G -t 1:00:00 --pty bash
time common/with_env.sh QUPKAKE python run_qupkake/qupkake_protomers.py \
    parents/preflight.smi -o preflight_out --name-style index -mp 8
cat preflight_out/qupkake_failed.tsv            # crashed molecules, now caught
```

Scale `--shard-size` and `#SBATCH --time` from that timing, then:

```bash
run_qupkake/submit_qupkake.sh parents protomers                 # + any sbatch args, e.g. --time=04:00:00
```

It activates the configured QupKake env on the login node first and checks that
the CLI, the Python package and `XTBPATH` resolve. A wrong path fails there,
once, rather than in every array task. Account, partition and array
concurrency come from the config.

Resubmitting the same range resumes. A task refuses to run if `shards.tsv` has
changed since it was written, or if a shard's `DONE` came from different input.

**A QupKake crash is no longer a "no sites" molecule.** QupKake catches its own
per-molecule and per-site exceptions and drops the row. The wrapper used to
read "absent from the SDF" as "nothing ionizable" and carry the molecule
through **neutral**, which for an amine is a wrong protomer that looks like a
clean one. `qupkake_protomers.py` now audits QupKake's working tree
(`processed/*.pt`, `raw/` vs `output/` sites), excludes failures, and records
them in `qupkake_failed.tsv`.

## 4. Names

```bash
common/with_env.sh PREP python run_qupkake/assign_names.py \
    --shards-dir protomers --parents parents -o library
python3 $DRAP/stage4/validate_names.py --smi library/library.smi    # DRAP repo
```

Refuses to run until every shard is `DONE`. Uniqueness is over the whole library,
and a shard can't see the others. `library.tsv` holds the name → ZINC id,
charge and pH lookup. db2_converter writes `NO_LONG_NAME`, so this file is the
only way back.

`invertomer_of` marks N-H⁺ epimers. They're genuine stereoisomers, kept as rows,
but one species in solution, so collapse each pair downstream (e.g. to its best
score).

**Stereo closure.** db2_converter's `--checkstereo` renames any SMILES with more
than one stereoisomer to `NAME.0`, `NAME.1`. That makes 18 characters, and
mol2db2 keeps the last 16. Protonation *creates* stereocentres
(`CC[NH+](C)Cc1ccccc1` has 2 isomers, its neutral parent 1), so enumerating the
parent is not enough. `assign_names.py` expands every protomer with
db2_converter's exact options and gives each isomer its own stereo index.

## 5. Before and after the build

```bash
llb check-stereo library/library.smi

# build: creates RUNDIR (config.env, shards/, out/, logs/), preflights db2_converter
# on this node, submits slurm/db2_array.sbatch -- John's 2026-08-27 script, with its
# environment now taken from DB2C_SETUP / DB2C_ENV / BUILD_LIGAND_EXE
llb submit-db2 db2_run library/library.smi --ncpus 36 --nshards 360 -- --time=48:00:00

# resubmit the same command to resume; then
llb verify-db2 library/library.tsv --outputpath db2_run/out
```

Defaults follow the db2_converter notebook (2026-08-27): `-m conformator`,
`-n 600 --keep_max_conf --checkstereo --reseth --rotateh --dock38`, MMFF off,
round-robin shards. `--method`, `--nconf`, `--mmff on`, `--ncpus`, `--nshards`,
`--mem-per-cpu` and `--time` override them. Each array task runs `NCPUS` builds on
node-local scratch, rsyncs back every 15 minutes, and refuses to start if
singularity, build_ligand, RDKit or AMSOL's libraries are missing.

A run directory belongs to one library. Resubmitting resumes without re-splitting.
Pointing it at a different library is refused. `verify-db2` reads every shard's
faillist: molecules build_ligand reported failing are `BUILD_FAILED` (attrition,
recorded), and only a molecule missing with **no** reason fails the run.

Keep `--checkstereo`. It's required anyway, since without it every chiral
ligand fails `chemistrycheck`, and your notebook's corollary ("every
stereocenter must be specified") is exactly what stereo closure guarantees.
After closure it renames nothing. **Don't pass
`--mergeiso`**: with closure it's a no-op, and if a split ever slipped through
it would merge the isomers under one name silently rather than visibly.

`verify_db2_build.py` now fails on `RENAMED` (a `.N` file) and `NAME_MISMATCH`
(M-line name ≠ requested), alongside the AMSOL-vs-formal charge check that
backs invariant 4.

## 5b. Mol2 from the db2 (optional, for mol2 / flexible-ligand docking)

No second 3D tool such as CORINA. Each molecule's mol2 is a conformer read back
out of its db2 by `db2tool tomol2` (dock6_claude `294d235`, using DOCK's own db2
reader and mol2 writer), so it is consistent with the db2 by construction:

| | mol2 | db2 |
| --- | --- | --- |
| coordinates | one conformer set from the db2, lowest internal energy | all sets |
| partial charges | the db2's AMSOL charges | same |
| atom types | the db2's SYBYL types | same |
| title | the 16-character contract name | same |

```bash
common/with_env.sh DOCK python3 run_qupkake/db2_to_mol2.py $d \
    --library library/library.tsv -o library/library.mol2
```

`db2tool` comes from `DB2TOOL_EXE` / `DOCK_SETUP` in `config/hpc.env`. Over plain
`tomol2`, the wrapper adds three things:

- **One conformer per molecule.** db2_converter writes one db2 record per rigid
  fragment, all under the same name, and `tomol2 --best` would emit each.
- **Nothing silently missing.** A molecule whose every set is flagged broken is
  taken from its lowest broken set and flagged in the `.tsv`
  (`--no-broken-fallback` to drop them instead).
- **The charge check again.** With `--library`, each mol2's summed charge must
  match the protomer's net charge, and every library name must come out. It
  exits 1 otherwise.

It accepts any db2 files, including mol_compiler chunks: run it per chunk to get
mol2 chunks that line up with the db2 chunks.

## 6. Chunking for docking (mol_compiler)

```bash
llb chunk . db2_run/out --target-per-bin 50000 --shards 200
```

Writes `work/` (intermediates, logs) and `chunks/` (`chunk_NNNNN.db2.gz` +
`manifest.tsv`) into the run directory, with account and partition from
`config/hpc.env`. `chunks.env` records the input and the code commit, and a
different input into the same run directory is refused.
`parallel/submit.slurm` still works standalone. Every setting is now an
environment override (`INPUT_DIR=… RUNDIR=… bash parallel/submit.slurm`), so
the checkout is never edited.

The verified `.db2.gz` files go to `mol_compiler`'s parallel pipeline as bare
`.db2.gz` inputs, as before. It takes a record's id from the first `M` line when no
`ZINC` token is present, which is the case for every contract name (tested:
`Ch300000BlLjl00C.db2.gz` → `Ch300000BlLjl00C`). Records are addressed by
position, never by id, so the names pass through untouched. One edge: an id with
no digit at all reads as `NO_ID` (`_looks_like_id`). A contract name needs a
digit-free ZINC id *and* protomer and stereo indices ≥ 10 for that, which no
zero-padded ZINC22 id produces.

## What is not done here

- **Identity-group dedupe** (`stage4/make_split.py`, SIZE_LADDER_SPEC §3) runs
  on the DRAP side.
- db2_converter is pinned at `hnlab/db2_converter@63d6656` in `DEPENDENCIES.md`.
  After any update, re-run `check_db2_stereo.py`.
- **Acceptance test (invariant 4)**: after the matrix build every new molecule
  should be `ambiguity_class = 0`.
