# lib_h18_22: the first production library (2026-10-16)

> **Code at tag `llb-full-2026-10`.** The construction code this record cites
> (`pka_triage/`, `run_qupkake/`, `pipeline/`, `environments/`) was removed on
> 2026-10-26, when ligand construction moved to MolSanitizer. Paths and commits
> here refer to that tag (`git checkout llb-full-2026-10`). Kept in the tree
> because `lib_h18_22` was built and docked with that code.

`REPORT.md` is the run's own report, written by `llb pipeline` (commit e03b65a at
start; three fixes were pulled mid-run, before per-stage commits were recorded:
917e7b4, e935fc5, 3e0dcb7). 50,000 parents, 10,000 per heavy-atom bin H18-H22
(`llb sample --per-bin 10000 --local`), no `--validate`.

Headline: 48,936 parents built (97.9%) from 98,091 names; 1,064 lost, almost all at
db2 (1,341 `1generate`, 1,042 `3chemistrycheck`, 9 `5amsolfail` names); about 8 h
from start to finish; QupKake 147 CPU-h for 19,064 routed parents (38.1%), db2
1,246 CPU-h, 43 GB. One stop: shard_01233 (a difluorocyclopropane-fused molecule
that fails half of all embeds) crashed QupKake twice; fixed (e935fc5) and resumed.

Known gaps in this report, fixed for later runs (8c0cff3): the db2_converter env's
versions are missing (a module banner hid the JSON), and every stage is credited
to the start commit. The `names` step ran on RDKit 2025.9.2 (the QupKake env) while
db2_converter uses 2024.03.5; check-stereo excluded the 2 names this split.

Run directory (scratch, purged after 60 days unused):
`/scratch/maom_root/maom99/johnbick/2026_10_08_build_small_mols/lib_h18_22`.
