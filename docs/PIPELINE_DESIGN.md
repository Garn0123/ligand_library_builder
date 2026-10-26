# `llb pipeline`: design (for review, 2026-10-16)

Goal: point `llb` at SMILES files and get a checked library, top to bottom, with
every molecule's fate and every setting recorded. It orchestrates the existing
steps (`llb parents` ... `llb verify-db2`) unchanged; it does not reimplement them.

## Commands

| command | what it does |
| --- | --- |
| `llb pipeline start SMILES... -o RUN [--validate] [--set key=value ...]` | checks `hpc.env` (`LLB_ACCOUNT`), writes `RUN/pipeline.json`, submits the first control job; returns. Activates no env: on ACES the container runtime exists only on compute nodes |
| `llb pipeline status RUN` | stage table, job ids, counts so far, why it stopped (if it did) |
| `llb pipeline advance RUN` | one idempotent step forward (what the control jobs call; safe by hand) |
| `llb pipeline resume RUN` | after a stop: re-check that stage's gate and continue |
| `llb pipeline cancel RUN` | `scancel` every job the run submitted |

**How it runs (primary: chained).** Each stage ends in a small *control job*
(1 core, ~10 min) that runs `advance`: it checks the stage's outputs, applies the
gate, then submits the next stage plus the next control job with
`--dependency=afterany:<array>`. `start` is the only command you type; nothing has
to stay alive on a login node. The first control job starts with `preflight`
(`llb doctor` and every env's versions, on the compute node), so a broken
environment stops the run before any work is queued. **Manual mode** is the same `advance`, run by you.
A lock file stops two `advance` calls from acting at once. Stage sizes (routed
shards, names) are only known when the previous stage ends, which is why the
chain extends itself rather than being submitted all at once.

## Stages

| # | stage | runs as | writes | gate: stop if | recorded, not stopping |
| --- | --- | --- | --- | --- | --- |
| 0 | preflight | control job | `pipeline/doctor.json`, `pipeline/versions.json` | any `llb doctor --stage host,prep,pka,qupkake,db2` FAIL (the stop names the host and the FAIL lines; fix, then `resume`) | doctor warnings (e.g. a UNICON license under 30 days), host |
| 1 | parents | control job | `parents/` | >1% unparsable/conflicting ids | rejected rows, tautomer fixes |
| 2 | triage | control job | `triage/` | MolGpKa failed on >1% | route + reason per parent |
| 3 | QupKake | array, 1 core/task | `protomers_qk/` | shards still missing after one retry at 2x time | `qupkake_embed_failed` etc. |
| 4 | reenumerate + merge (+ compare if `--validate`) | control job | `protomers_qk_h/`, `merged/`, `compare/` | any merge refusal | aromatic-N swaps, fallbacks |
| 5 | names | 1 job, 16 cores | `library/` | contract violation; >5% of parents lost | dropped stereo, kept-whole protomers |
| 6 | checks | control job (db2 env) | `library/build.smi`, `check_stereo.log` | validate-names fails | check-stereo `[FAIL]` and crash names excluded |
| 7 | db2 | array, 8 cores/task, ~25 names/task | `db2_run/` | tasks missing after one retry at 2x time | build failures by reason |
| 8 | verify + report | control job | `db2_verification.csv`, `ledger.tsv`, `attrition.tsv`, `REPORT.md` | any charge or name mismatch; parent loss >5% | everything else |

`--validate` runs QupKake on every parent (stage 3) and adds `llb compare`
(stage 4), for a new chemical space or a changed setting.

Retries: a task that failed or timed out is resubmitted once at twice the time
limit. A second failure stops the stage and names the shards. A failed gate
(as opposed to a time limit) always stops; nothing retries it automatically.

Waves: an array that would push the account past `MaxSubmit` (5,000) goes out
in waves, each wave's control job submitting the next.

## Records

**`pipeline.json`** (written by `start`, never edited; a rerun with different
values refuses unless it is a new `-o`):
inputs (paths, sha256, line counts); every setting (pH list, route window,
margin, min population, max states, shard sizes, time limits, cores, MMFF mode,
nconf, method); cluster (account, partition, `LLB_ARRAY_PARALLEL`).

**`pipeline/versions.json`** (written by `preflight`, rewritten if `preflight`
is rerun after a fix): LLB commit, `SITE_FILTER_VERSION`, `COUPLING_VERSION`,
`AROMATIC_N_VERSION`, db2_converter commit, RDKit in each env, QupKake and
MolGpKa versions, MolGpKa weight sha256. Runs started before `preflight`
existed (2026-10-20) have these in `pipeline.json` instead; their state has no
`preflight` entry, which is treated as done, and `REPORT.md` reads
`pipeline.json`'s versions.

**`pipeline/doctor.json`**: every `llb doctor` check of the preflight
(stage, check, status, detail, fix) and the host it ran on.

**`pipeline/state.json`**: per stage: status, job ids, start/end, retry count,
gate result, counts.

**`ledger.tsv`**: one row per input molecule:

| column | from |
| --- | --- |
| `input_name`, `input_smiles`, `source_file`, `heavy_atoms` | input / parents |
| `parent_id`, `parent_status` (kept, or the rejection reason), `tautomer_fix` | parents |
| `route`, `route_reason`, `site_source`, `aromatic_n_swapped` | triage / merge |
| `n_protomers`, `n_names`, `n_stereo_dropped`, `n_kept_whole`, `n_excluded_checks` | names / checks |
| `n_built`, `n_build_failed`, `build_fail_reasons` | verify |
| `fate` (`built`, or the stage it was lost at), `loss_reason` | derived |

**`attrition.tsv`**: molecules in and out at each stage, by heavy-atom bin.

**`REPORT.md`**: preflight (host, doctor checks and warnings), attrition, cost
(`llb timing`), versions, gates passed, and
anything stopped and resumed: the page to attach to a methods section.

## Changes needed in existing steps

- `submit-qupkake` / `submit-db2`: print the job id(s) machine-readably
  (`--parsable`), take a task-list to retry only failed indices.
- `names`, `merge`, `verify-db2`: already write the counts the ledger needs;
  the ledger is assembled from their files in stage 8, not threaded through.
- A tiny `llb versions` that prints every version above as JSON (run once per
  env by `preflight`; `common/versions.py`).

## Settled in review (2026-10-16)

1. The pipeline stops at verified db2s. `REPORT.md` and the final `status` print the
   `llb mol2` / `llb mol2-library` commands for that run.
2. **Stop threshold** = the loss rate at a stage above which the pipeline halts
   instead of recording and continuing. Set at a bit over twice what the drug-like
   10k measured; per run, `--set stop.<stage>=<fraction>` overrides one.

   | stage | counts as lost | measured (10k) | stop above |
   | --- | --- | --- | --- |
   | parents | inputs rejected (unparsable, id conflicts; not duplicates) | 0% | 1% |
   | triage | parents MolGpKa failed on | 0% | 1% |
   | qupkake | routed parents QupKake failed on (they fall back to MolGpKa) | 2.0% | 10% |
   | names | parents with no name | ~0.04% | 2% |
   | checks | names excluded by check-stereo | 0.03% | 0.5% |
   | db2 | named parents with nothing built | 2.0% | 5% |
   | end to end | input parents with no db2 | 2.0% | 5% |

3. Email on stop and on finish to `LLB_MAIL` from `hpc.env` (none if empty): every
   control job carries `--mail-type=FAIL` (a stop exits non-zero), the last one
   `END,FAIL`.
4. One library per `start`; batching by heavy-atom range stays manual.

Not built yet: submitting in waves. `start` refuses an array over 4,500 tasks
(`MaxSubmit` is 5,000), which today's batches stay under.
