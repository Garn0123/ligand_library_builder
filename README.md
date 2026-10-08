# ligand_library_builder

Setting up a small-molecule library for docking: ZINC SMILES in, chunked
DOCK-ready db2 out, every molecule carried under a name that survives the db2
format (`NAMING_CONTRACT.md` in the DRAP project, enforced here by
`run_qupkake/assign_names.py`).

```
bin/llb         one command for every step; put bin/ on PATH and run from any run directory
mol_download/   ZINC22 tranches: fetch (run.sh), verify, and sample N per heavy-atom bin
run_qupkake/    parents -> QupKake micro-pKa -> protomers -> contract names -> db2 checks
pka_triage/     MolGpKa triage + merge: QupKake only where it matters (branch feat/pka-triage)
environments/   conda specs for envs this repo defines (llb_pka.yml)
mol_compiler/   db2 -> load-balanced chunks + manifest, the unit of work for docking
common/         env.sh / with_env.sh: activate each tool the way config/hpc.env says
config/         hpc.env.example -- the ONLY place tool paths live
```

Docking itself (DOCK6, `dock6_claude`) and the models trained on its output live
elsewhere. This repo stops at chunked db2.

## Setup on a new machine

```bash
git clone git@github.com:Garn0123/ligand_library_builder.git /nfs/turbo/.../opt/ligand_library_builder
mkdir -p ~/.config/llb
cp /nfs/turbo/.../opt/ligand_library_builder/config/hpc.env.example ~/.config/llb/hpc.env  # edit it
export PATH=/nfs/turbo/.../opt/ligand_library_builder/bin:$PATH   # ~/.bashrc
export LLB_CONFIG=~/.config/llb/hpc.env                            # ~/.bashrc
llb help                                       # then: cd to a run directory and go
```

Environments are separate on purpose: QupKake pins `xtb 6.4.1` and its own
torch stack, and db2_converter needs its own module set. `hpc.env` names one env
per tool plus a setup hook for each; `llb <step>` runs each step inside the right
one, and `llb env TOOL cmd` runs anything else there.

## Running it

**[UPSTREAM.md](UPSTREAM.md)** is the runbook: every step in order, the file each
one writes, and the check that gates the next. Start there.

Per-component detail:
[mol_download/README.md](mol_download/README.md) (fetching from files.docking.org
politely, and the gotchas that actually happened),
[run_qupkake/README.md](run_qupkake/README.md) (how protonation states are chosen),
[mol_compiler/README.md](mol_compiler/README.md) (the db2 record model and chunking).

## History

This repo was `mol_compiler`, then briefly `small_molecule_prep`. Its history is intact under `mol_compiler/`
(`git log --follow mol_compiler/<file>`). `mol_download` and `run_qupkake` joined
it on 2026-09-23, because the three run in sequence, share one naming contract,
and change together.
