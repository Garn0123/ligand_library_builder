# ligand_library_builder

The two ends of building a small-molecule docking library: choosing molecules
from ZINC22, and turning finished db2 files into docking jobs. The middle step,
building DOCK-ready db2s from SMILES, is MolSanitizer's job (Carlsson lab;
set up under `molsanitizer/`).

```
ZINC download  ->  sample bins  ->  MolSanitizer (see molsanitizer/)  ->  docking jobs
llb fetch          llb sample       SMILES -> protonated, 3D, db2         llb chunk / llb mol2
```

```
bin/llb         one command for every step; put bin/ on PATH and run from any run directory
mol_download/   ZINC22 tranches: fetch (run.sh), verify, and sample by heavy-atom bin;
                zinc3d/: ZINC22's own 3D db2s for a list of ZINC ids
mol_compiler/   db2 -> load-balanced chunks + manifest, the unit of work for docking
docking/        mol2 from a db2, DOCK-header mol2 shards, invertomer collapse after docking
common/         env.sh / with_env.sh: activate each tool the way hpc.env says
config/         hpc.env.example: the ONLY place tool paths live
environments/   llb.yml, the one conda env (PREP) this repo defines
```

Docking itself (DOCK6, `dock6_claude`; the run driver in `rna_pipeline`,
`pipeline/library_docking/`) and the models trained on its output live
elsewhere.

## Setup on a new machine

One site directory per cluster, on shared group storage (never `~`), with the
checkout inside it:

```bash
LLB_SITE=/scratch/group/<project>/llb        # or Turbo on Great Lakes
git clone git@github.com:Garn0123/ligand_library_builder.git $LLB_SITE/ligand_library_builder
cp $LLB_SITE/ligand_library_builder/config/hpc.env.example $LLB_SITE/hpc.env   # cluster settings; edit it
mkdir -p $LLB_SITE/users && $EDITOR $LLB_SITE/users/$USER.env   # yours: LLB_ACCOUNT, ZINC_DATA
export PATH=$LLB_SITE/ligand_library_builder/bin:$PATH          # ~/.bashrc; the only line there
llb help                                       # then: cd to a run directory and go
```

`llb` finds `hpc.env` as `../hpc.env` from its checkout, then reads
`users/$USER.env` beside it, so one `hpc.env` serves everyone on the cluster.
`export LLB_CONFIG=/path/to/hpc.env` overrides that, and `config/hpc.env`
inside the checkout still works.

Download, sampling and chunking need only the system `python3` and SLURM. The
`docking/` tools and the zinc3d steps run in `PREP`, a conda env made from
`environments/llb.yml`; `llb mol2` needs `db2tool` (`DOCK` in `hpc.env`).

## Running it

**[UPSTREAM.md](UPSTREAM.md)** is the runbook: every step in order, the file
each one writes, and the check that gates the next. Start there.

Per-component detail:
[mol_download/README.md](mol_download/README.md) (fetching from files.docking.org
politely, and the gotchas that actually happened),
[mol_download/zinc3d/RUNBOOK.md](mol_download/zinc3d/RUNBOOK.md) (ZINC22's 3D db2s),
[mol_compiler/README.md](mol_compiler/README.md) (the db2 record model and chunking),
[docking/README.md](docking/README.md) (mol2 and invertomer collapse).

## Tests

    python -m pytest                          # tests/: docking/ tools, hpc.env lookup
    (cd mol_compiler && python -m pytest)     # chunking

Run them in the `llb` env. GitHub Actions runs both on every push
(`.github/workflows/tests.yml`).

## History

This repo was `mol_compiler`, then briefly `small_molecule_prep`. Its history is
intact under `mol_compiler/` (`git log --follow mol_compiler/<file>`).
`mol_download` and `run_qupkake` joined it on 2026-09-23.

From September to October 2026 LLB also built the ligands: MolGpKa triage,
QupKake protonation, 16-character contract names, and db2 via db2_converter
(conformator + BCL), driven end to end by `llb pipeline`. `lib_h18_22` (50k
parents, the first production library) was built and docked with it. On
2026-10-26 that construction code was removed in favour of MolSanitizer
(compared in `docs/MOLSANITIZER_COMPARISON.md`).

**The full pipeline is at tag `llb-full-2026-10`** (branch `archive/llb-full`).
That tag merges every branch that was open then: the pipeline, the install design and
toolchain image, the MolSanitizer comparison, the decision record and the
validation records. `llb-full-2026-10-reindex` keeps the unmerged
`mol_compiler` `05_reindex` work.
`git checkout llb-full-2026-10` gives back a working `llb pipeline`.

Kept in this tree because they justify `lib_h18_22`:
[DECISIONS.md](DECISIONS.md) (why each construction choice was made) and
[validation/](validation/) (the runs that tested it). The code they cite is at
`llb-full-2026-10`.
