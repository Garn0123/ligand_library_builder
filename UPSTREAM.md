# ZINC SMILES → samples → (MolSanitizer) → docking jobs

The path from ZINC22 to docking jobs, in the order it runs. Each step writes a
file the next one checks, so a failure stops at the step that caused it rather
than surfacing as a wrong row in the matrix.

```
llb fetch + llb sample     ZINC22 2D on disk, N per heavy-atom bin   samples/H*.smi
MolSanitizer               SMILES -> protonated 3D -> db2            see molsanitizer/
llb chunk                  db2 -> balanced chunks + manifest         docking input
llb mol2 / mol2-library    mol2 read back out of the db2 (optional)  mol2 docking input
  ... docking (rna_pipeline pipeline/library_docking/) ...
llb collapse-invertomers   best score per N-H+ invertomer group      score table
```

Building the ligands (protonation, tautomers, names, conformers, db2) is not
done here any more. LLB's own construction stages (MolGpKa triage, QupKake,
contract names, db2_converter, `llb pipeline`) are at tag `llb-full-2026-10`;
[DECISIONS.md](DECISIONS.md) records why each choice was made for `lib_h18_22`.

## 0. Point the repo at your software (once per machine)

```bash
# the site layout: one directory per cluster, the checkout inside it
cp $LLB_SITE/ligand_library_builder/config/hpc.env.example $LLB_SITE/hpc.env   # edit it
$EDITOR $LLB_SITE/users/$USER.env      # per-user: optional LLB_ACCOUNT, ZINC_DATA
```

`llb` looks for its config in this order: `LLB_CONFIG` if set, then
`../hpc.env` beside the checkout (`$LLB_SITE/hpc.env`), then `config/hpc.env`
inside the checkout (gitignored). Then it reads `users/$USER.env` from the same
directory, whose keys win, so the cluster settings are shared and only the
per-user ones are yours.

Every step reads tool locations from that file: conda env name or prefix, an
optional executable path, and an optional setup hook run *before* activation
(`module purge`, `module load ...`). No script hard-codes a path. There are two
tools: `PREP` (Python with RDKit and pandas, `environments/llb.yml`) and `DOCK`
(`db2tool`). `llb env TOOL CMD` runs any command in a tool's environment, e.g.
`llb env DOCK which db2tool`.

### Code in one place, runs anywhere

Keep the checkout in `opt/` and never run inside it. `bin/llb` runs any step,
in the right environment, against your **current** directory:

```bash
export PATH=/nfs/turbo/umms-maom/johnbick/opt/ligand_library_builder/bin:$PATH   # ~/.bashrc
mkdir -p /scratch/.../ladder_h14_h28 && cd /scratch/.../ladder_h14_h28
llb help                                      # every step, in order
```

Every output path is relative to where you run it, and SLURM logs go to
`./logs`, so one run directory holds one run and the checkout stays read-only.
`llb root` prints which checkout you're using.

| step | `llb` | script |
| --- | --- | --- |
| URL list | `llb print-urls` / `llb extract-urls` | `mol_download/sample_2d.py` / `extract_urls.sh` |
| download | `llb fetch`, `llb fetch-verify`, `llb fetch-status` | `mol_download/run.sh`, `verify.sh`, `status.sh` |
| sample | `llb sample` | `mol_download/sample_2d.py` |
| ZINC 3D db2s | `llb zinc3d-*` | `mol_download/zinc3d/` (its `RUNBOOK.md`) |
| chunk | `llb chunk` | `mol_compiler/parallel/submit_chunks.sh` → `submit.slurm` |
| mol2 | `llb mol2`, `llb mol2-library` | `docking/db2_to_mol2.py`, `docking/mol2_library.py` |
| after docking | `llb collapse-invertomers` | `docking/collapse_invertomers.py` |

The scripts are listed for reading; run them through `llb`, which picks the
environment and keeps outputs out of the checkout.

## 1. Download and sample

**Production: download the bins, then sample uniformly from disk.** Streaming a
prefix of each remote file is biased. The first few thousand lines of a tranche
are one registration batch, with fewer scaffolds (1,479 vs 1,654 per 3,600) and
more near-neighbours (NN Tanimoto 0.469 vs 0.432) than a uniform draw. Numbers are
in `sample_2d.py`'s docstring.

```bash
llb print-urls --heavy 14-28 > urls_2d.txt  # 904 files, ~960 GB
# or, from a script CartBlanche22 gave you (curl/wget/PowerShell; don't run it as-is):
#   llb extract-urls zinc22-2D-download.curl > urls_2d.txt
cd $ZINC_DATA                                     # scratch; files land under $PWD/zinc22/
llb fetch $OLDPWD/urls_2d.txt 4                   # DTN / tmux, ~5-19 h
llb fetch-verify --purge
llb fetch-status $OLDPWD/urls_2d.txt              # then llb fetch retry.txt 2 until retry = 0
llb tranche-pattern permanent.tsv                 # are the remaining 404s real absences?
cd -                                              # back to the run directory
llb sample --heavy 14-28 --local $ZINC_DATA/zinc22 --urls urls_2d.txt --outdir samples
```

(`llb fetch` takes `WGET` and `PARALLEL` from the config.)

**Proportional to bin size:** `--total N` in place of `--per-bin` draws N
molecules uniformly over all the bins together, so each bin's share matches its
share of ZINC's molecules (one pass; still one `Hxx.smi` per bin, plus
`POOLED.done` with each bin's population and count):

```bash
llb sample --heavy 18-22 --total 10000 --local $ZINC_DATA/zinc22 --urls urls_2d.txt --outdir samples
```

`--local` is one reservoir pass over every file in the bin: every molecule
equally likely, the logP mix right by construction, about 30 min for H28. With
`--urls` it refuses to sample a bin that is missing files.

**Pilot only:** without `--local`, the sampler streams a few MB per bin. Use it to
prove the plumbing end to end, never for results.


## 2. Build the ligands: MolSanitizer

`samples/H*.smi` (one `SMILES ZINC_ID` per line, one file per heavy-atom bin)
is MolSanitizer's input. How to run it on Great Lakes, and what it writes, is
in `molsanitizer/`. Its output is one or more `.db2.gz` files, which the next
steps take as they are.

## 3. Docking jobs

### Chunking (mol_compiler)

```bash
llb chunk . <db2_dir> --target-per-bin 50000 --shards 200
```

Writes `work/` (intermediates, logs) and `chunks/` (`chunk_NNNNN.db2.gz` +
`manifest.tsv`) into the run directory, with account and partition from
your `hpc.env`. `chunks.env` records the input and the code commit, and a
different input into the same run directory is refused.
`parallel/submit.slurm` still works standalone. Every setting is an
environment override (`INPUT_DIR=… RUNDIR=… bash parallel/submit.slurm`), so
the checkout is never edited.

`mol_compiler` takes a record's id from its `ZINC` token, or from the first
`M` line when there is none. Records are addressed by position, never by id, so
names pass through untouched. One edge: an id with no digit at all reads as
`NO_ID` (`_looks_like_id`).

`rna_pipeline`'s `pipeline/library_docking/` plans the docking run itself
(pockets × ligand batches) and has its own `make_db2_chunks.sh` for DRAP's
jobA. Its `prepare.py` still reads an `llb pipeline` run directory
(`db2_verification.csv`, `db2_run/out/`). Pointing it at MolSanitizer's output
is a change in `rna_pipeline`.

### Mol2 from the db2 (optional, for mol2 / flexible-ligand docking)

No second 3D tool such as CORINA. Each molecule's mol2 is a conformer read back
out of its db2 by `db2tool tomol2` (dock6_claude `294d235`, using DOCK's own db2
reader and mol2 writer), so it is consistent with the db2 by construction:
same coordinates (one conformer set, lowest internal energy), partial charges
and SYBYL atom types, and the db2's molecule name as the title.

```bash
llb mol2 <db2_dir> [--library library.tsv] -o library.mol2
llb mol2-library library.mol2 --library library.tsv -o mol2_lib --shard-size 1000
```

`db2tool` comes from `DB2TOOL_EXE` / `DOCK_SETUP` in your `hpc.env`. Details
(broken sets, one conformer per molecule, the charge check) are in
[docking/README.md](docking/README.md).

## 4. After docking: invertomer collapse

```bash
llb collapse-invertomers scores.parquet --library library.tsv -o scores_collapsed.parquet
```

Keeps the best score per N-H+ invertomer group. It needs a `library.tsv` with an
`invertomer_of` column, which LLB's names step wrote for libraries built
before 2026-10-26 (`lib_h18_22`). DRAP's `stage4/collapse_invertomers.py` is the
streaming port used at matrix scale. Its parity test loads this script by path
and still names the old `run_qupkake/` location, so until DRAP is updated to
`docking/` it skips that check.

## What is not done here

- **Ligand construction**: MolSanitizer (`molsanitizer/`).
- **Identity-group dedupe** (`stage4/make_split.py`, SIZE_LADDER_SPEC §3) runs
  on the DRAP side.
- **Docking**: `rna_pipeline` `pipeline/library_docking/` and DOCK6.
