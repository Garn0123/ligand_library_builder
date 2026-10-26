# MolSanitizer as the ligand-construction step

2026-10-26. Replaces everything LLB did between "sampled SMILES" and "db2 files":
`parents` → `triage`/MolGpKa → QupKake → `merge` → `names` → db2_converter →
`verify-db2`. Input is the per-bin `SMILES NAME` files from `llb sample`. Output is
a run directory that `rna_pipeline` `pipeline/library_docking/prepare.py --library`
and DRAP stage 4 read unchanged.

MolSanitizer: Lam, Pach, …, Carlsson, Cabeza de Vaca, ChemRxiv 2026,
doi:10.26434/chemrxiv.15010189/v1; github.com/carlssonlab/MolSanitizer at
`de79803`, PyPI 0.7.4. For the paper summary and the comparison with LLB, see
`docs/MOLSANITIZER_COMPARISON.md` on branch `claude/molsanitizer-comparison-b31000`.

Great Lakes commands (install, selftest, smoke, DOCK6 preflight, production) are in
the runsheet: https://claude.ai/artifact/V7yUdTvymyZheaapV3LnKk

## 1. What MolSanitizer needs and produces

| | |
| --- | --- |
| **Install** | `pip install molsanitizer==0.7.4`. Prebuilt wheels (manylinux_2_28 x86_64, Python 3.10–3.14) carry the compiled conformer, stereo and AMSOL code; Great Lakes (RHEL 8, glibc 2.28) qualifies. RDKit comes from PyPI: 0.7.4 asks for ≥ 2025.3.1, its docs list 2025.9.1–2026.3.6 as tested, so we pin **2026.3.6**. Apache 2.0. |
| **External tools** | **None needed, none licensed.** Open Babel (conda-forge) is optional: it is the embedding fallback after a 2-min RDKit timeout, and we install it. CORINA and OpenEye (OEB output) are optional and unused. AMSOL 7.1 is built in as a C++ port. |
| **Input** | 2 columns `SMILES ID`, no header, whitespace or tab separated, plain or gz/bz2/xz. |
| **Steps** (order fixed) | metal complexes → salts → filters (all off unless asked) → neutralize → tautomers (2 layers) → PAINS/unwanted → protonation (SMARTS rules, integer pKa, pH ± range) → stereoisomers (RDKit, ≤ 8, 3D-embeddable only) → conformers (srETKDGv3 + ring library + TorLib sampling, MMFF94s, ≤ 2,000, 25 kcal/mol window, 0.5 Å RMSD) → AMSOL → db2. |
| **Outputs** | `*_clean` SMILES (one row per tautomer × protomer × stereoisomer), **db2 per molecule** (`-f db2`) or one `db2.tgz`, plus SDF/MOL2/PDBQT if asked; `msani_error.err` for per-molecule 3D failures; `msani_timing.csv` with `--timing`. |
| **db2** | **Written by MolSanitizer itself.** Nothing fills a gap. One file per state holds every ring conformation × rigid fragment as separate db2 records. The paper docked its db2s with DOCK 3.8 (DUDE-Z). Our docking runs **DOCK 6.13's `hdb_db2_input_file` reader** (`rna_pipeline` library_docking, see its notes), so the smoke test docks a few with it before anything scales (§5). |
| **States** | One state per pH value, both forms when pH equals a rule's pKa, and protomers enumerated for flagged groups (asymmetric piperazines). A protonated tertiary amine's N is sampled both ways *inside* one entry's conformers (paper §2.5). |
| **Parallelism** | `-j N` processes in one run, or `msani_batch` (its own SLURM array, 200 lines/task). We use neither: one process per 1-CPU array task, with our own sharding, retry and accounting (§4). |

**The two traps this design avoids:**

1. *Names.* MolSanitizer appends `_k` per state and `.j` per stereoisomer, and its
   db2 writer keeps `name[-16:]` (`db2/layout.py`). Any suffix on a 16-character
   ZINC22 id therefore loses the id's leading characters inside the db2. On
   2026-10-20 `Ch300000BlLjl00C.1` came out as `M 300000BlLjl00C.1`. That is
   NAMING_CONTRACT invariant 1, the collapse that cost +2.56 RMSE. See §2.
2. *Defaults files.* The `msani` CLI reads a personal
   `~/.config/msani/msani_configurations.yaml` if one exists. Every task points
   `MSANI_CONFIG` at the run's own copy and passes every option explicitly.

## 2. Names: MolSanitizer's states → the 16-character contract

`C` + 12-character base id + state + stereo + `C` (DRAP `NAMING_CONTRACT.md`,
checked by `stage4/validate_names.py`).

| field | from |
| --- | --- |
| base id | the parent id, the same rule as LLB (`naming.py`, a copy of `run_qupkake/llb_ids.py`): a ZINC id's last 12 characters verbatim, otherwise `Z` + 11 base62 characters of sha256(name). A damaged `ZINC…` id is rejected. |
| state (contract "protomer") | MolSanitizer's `_k`, renumbered `0..` per parent in its order. It numbers **tautomer × protomer** states: MolSanitizer gives one index to both. Nothing downstream reads meaning from it. |
| stereo | MolSanitizer's `.j`, renumbered `0..` per state. |

How: the 2D step sends MolSanitizer the **parent id** and *parses* its suffixes
(`naming.parse_state_id`; an unexpected format stops the parent with a reason, it
is not guessed at). The 3D step sends the finished **contract name** with
`--no-stereoisomers`, so no suffix is added and `name[-16:]` is the whole name.
Every db2's M-line name is then checked against the requested name.

- **Invertomers.** Paper §2.5 says MolSanitizer samples both N-H⁺
  configurations inside one entry. If its RDKit stereo step still emits them as
  separate isomers, `invertomer_of` links them with LLB's rule, so DRAP
  `collapse_invertomers.py` keeps the better score. The smoke test shows which
  case holds (`collect.json` `n_h_invertomers`).
- **One structure, one name.** A SMILES reached from two parents keeps the
  lowest-gidx parent's name. The other is `DUPLICATE` and is not docked, because
  two base ids for one structure would leak across the identity split.
- **Limits.** Up to 62 states and 62 isomers per state. MolSanitizer stops at 8
  isomers anyway; LLB allowed 32.

## 3. Install (Great Lakes)

```bash
llb msani install          # -> /nfs/turbo/umms-maom/johnbick/opt/molsanitizer/0.7.4
```

`environment/environment.yml` holds Python 3.12 and Open Babel 3.1.1 from
conda-forge. `requirements.txt` holds molsanitizer 0.7.4 and rdkit 2026.3.6 from
PyPI. The install writes `PREFIX/LOCK/conda-linux-64.lock.txt` (explicit, with
md5s) and `requirements.lock.txt`. Commit both into `environment/`, and
`install.sh PREFIX --from-lock` then rebuilds the same environment. The install
also checks that the compiled extensions load and that `obabel` is found.
Settings are environment variables, listed in `msani.sh` (`LLB_MSANI_PREFIX`,
`LLB_MSANI_ACCOUNT=maom2`, …).

## 4. Workflow

```bash
llb msani plan H18.smi H19.smi ... -o RUN --per-shard 40   # shard + run.json (stdlib)
llb msani submit RUN                                      # 1 CPU per shard, logs in RUN/logs/
llb msani status RUN                                      # done / missing, as array ranges
llb msani submit RUN --tasks missing                      # resume: done shards are skipped
llb msani collect RUN --drap-dir $DRAP_DIR                # verify + merge + validate_names
```

Each array task (`build_shard.py`, one shard):

1. **2D**: one parent at a time through the Python API, which takes a
   non-integer pH (the CLI does not). Each failure gets its own reason, and each
   parent its own seconds.
2. **Name** the states (§2), mark db2-incompatible ones (`EXOTIC`) and
   in-shard duplicates.
3. **3D**: one `msani … -3d -f db2 --no-stereoisomers --no-neutralize
   --metal off --timing` call for the shard. Molecules the call leaves without a
   db2 or an error, because it crashed or timed out, are retried **one at a
   time** under `mol_timeout`, so one bad molecule costs only itself.
4. **Check every db2 before keeping it**: each record's name is the requested
   name, AMSOL total charge equals the formal charge (± 0.05, as
   `llb verify-db2`), and conformations > 0. Then gzip it to
   `RUN/db2_run/out/NNNNN/<name>.db2.gz` (temp file + rename).
5. `DONE.json` is written last: versions, host, CPU, SLURM ids, seconds, and the
   input sha256. A shard whose `DONE.json` matches is skipped. An interrupted
   shard is redone, but it keeps every db2.gz whose structure is unchanged.

`collect.py` works over the **whole** run. It refuses if any shard is missing, if
any input line is unaccounted for or counted twice, if names or name tails are
not unique, if a base id spans two parents, or if shards ran different software
or options. It writes:

| file | for |
| --- | --- |
| `db2_verification.csv` | `prepare.py`: only `status == OK` is docked |
| `db2_run/out/NNNNN/*.db2.gz` | `prepare.py --db2-root` default |
| `library/library.tsv`, `library.smi` | DRAP `collapse_invertomers.py` / `join_library.py --llb-library`; `validate_names.py`. LLB's columns: `site_source=molsanitizer`, `population_estimate` empty, `ph_values` = the run's pH |
| `library/dropped.tsv` | `DUPLICATE` / `EXOTIC` states |
| `parents_ledger.tsv` | **every input line**: `built`, or `failed` + reason |
| `timing.tsv` | per name: MolSanitizer's embedding / AMSOL / sampling / db2 seconds, and ours |
| `collect.json` | counts, fail reasons, attrition per input file, s/name, CPU-h per 1,000 parents, db2 bytes |
| `provenance.json` | MolSanitizer / RDKit / Python / Open Babel versions, options, LLB commit (and whether the checkout was dirty), CPUs, SLURM jobs, time span |

Tests: `pytest molsanitizer/tests` or `python molsanitizer/tests/test_msani_build.py`.
They need RDKit and pandas, not MolSanitizer: `tests/fake_msani` stands in for its
API and CLI in 0.7.4's formats.

## 5. Smoke test, then production sizing

**Smoke** (~54 molecules, John runs it; the commands are in the runsheet):
3 parents from each of H18–H30 plus `probes/edge_cases.smi`. The probe file
covers:

- an N-H⁺ invertomer;
- piperazine, metformin, imidazole;
- a tautomer (2-hydroxypyridine);
- Meldrum's acid;
- a bicyclo[1.1.1]pentane amine;
- a selenium compound (expected `EXOTIC`);
- sodium acetate (removed whole by salt stripping);
- an invalid SMILES;
- a short ZINC id;
- a damaged ZINC id and a repeated id (plan-rejected);
- the same structure under two names (`DUPLICATE`).

The run passes when:

1. collect succeeds;
2. `validate_names.py` passes;
3. every ledger row has a status;
4. the edge cases land as expected;
5. DOCK 6.13 docks 5 OK names against 1 pocket with plausible scores, using
   `prepare.py` + `run_pairs.sbatch` from the library-docking README §1.

Read off: s/name, names/parent, `MaxRSS` (sacct), db2 size, and `n_h_invertomers`.

**Production H18–H30**: 13 bins × 10,000 = 130,000 parents. Until the smoke test
measures it, assume ~2 names/parent (LLB's median was ~1.9) and 10–60 s per name
at 2,000 conformers. That brackets 260k names at **~0.7k–4.3k CPU-h**: 2–11 h of
wall time at 400 CPUs. db2_converter cost ~50 s/name. Pick `--per-shard` so a task
takes ~1.5 h:

```
per_shard ≈ 5400 / (names_per_parent × seconds_per_name)      # e.g. 5400/(2×30) = 90
```

At 90 per shard there are ~1,450 tasks: under MaxSubmit 5,000, under
MaxArraySize, and in one submit. Run them `--parallel 400`, which leaves room
inside maom_root's 2,000-CPU cap. Set the time limit to 3× the smoke p90 task
time. Copy the run to Turbo when done (scratch purges after 60 days).

**Docking cost depends on these choices.** MolSanitizer's default is 2,000
conformers; LLB built ≤ 600. The H18–22 pilot docked at 1.28 s per
ligand-pocket. More conformers, or more states from `--pH-range 1`, multiply the
~10⁵ core-h docking bill, not just this step's.

## 6. Open choices (defaults are MolSanitizer's published settings)

| option | default | alternative | trade-off |
| --- | --- | --- | --- |
| `--pH` / `--pH-range` | 7 / 0 (benchmarked) | 7.4 / 1 (LLB's window) | ±1 adds the imidazole / 2-aminopyridine cations (rule pKa 6) and more names to dock |
| `--numconfs` | 2,000 | 600 (LLB) | docking time per name |
| `--max-stereoisomers` | 8 | up to 62 | LLB kept ≤ 32 |
| `--tautomers` | on | off | the reason to switch (TautoBase 70% vs LLB 5/13 probes) |
| filters (PAINS, unwanted) | off | on | a training set, not a hit list: off |
