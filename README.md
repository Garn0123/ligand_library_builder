# small_molecule_prep

Setting up a small-molecule library for docking: ZINC SMILES in, chunked
DOCK-ready db2 out, every molecule carried under a name that survives the db2
format (`NAMING_CONTRACT.md` in the DRAP project, enforced here by
`run_qupkake/assign_names.py`).

```
mol_download/   ZINC22 tranches: fetch (run.sh), verify, and sample N per heavy-atom bin
run_qupkake/    parents -> QupKake micro-pKa -> protomers -> contract names -> db2 checks
mol_compiler/   db2 -> load-balanced chunks + manifest, the unit of work for docking
common/         env.sh / with_env.sh: activate each tool the way config/hpc.env says
config/         hpc.env.example -- the ONLY place tool paths live
```

Docking itself (DOCK6, `dock6_claude`) and the models trained on its output live
elsewhere. This repo stops at chunked db2.

## Setup on a new machine

```bash
git clone git@github.com:Garn0123/small_molecule_prep.git
cd small_molecule_prep
cp config/hpc.env.example config/hpc.env      # then point it at QupKake, db2_converter, ...
```

Environments are separate on purpose: QupKake pins `xtb 6.4.1` and its own
torch stack, and db2_converter needs its own module set. `config/hpc.env` names
one env per tool plus a setup hook for each, and `common/with_env.sh TOOL cmd`
runs any step inside the right one.

## Running it

**[UPSTREAM.md](UPSTREAM.md)** is the runbook: every step in order, the file each
one writes, and the check that gates the next. Start there.

Per-component detail:
[mol_download/README.md](mol_download/README.md) (fetching from files.docking.org
politely, and the gotchas that actually happened),
[run_qupkake/README.md](run_qupkake/README.md) (how protonation states are chosen),
[mol_compiler/README.md](mol_compiler/README.md) (the db2 record model and chunking).

## History

This repo was `mol_compiler`. Its history is intact under `mol_compiler/`
(`git log --follow mol_compiler/<file>`). `mol_download` and `run_qupkake` joined
it on 2026-09-23, because the three run in sequence, share one naming contract,
and change together.
