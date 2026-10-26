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

Cluster commands (install, selftest, smoke, DOCK6 preflight, production) are in
the runsheet: https://claude.ai/artifact/V7yUdTvymyZheaapV3LnKk

## 1. What MolSanitizer needs and produces

| | |
| --- | --- |
| **Install** | `pip install molsanitizer==0.7.4`. Prebuilt wheels (manylinux_2_28 x86_64, Python 3.10–3.14) carry the compiled conformer, stereo and AMSOL code; Lighthouse and Great Lakes (RHEL 8, glibc 2.28) qualify. RDKit comes from PyPI: 0.7.4 asks for ≥ 2025.3.1, its docs list 2025.9.1–2026.3.6 as tested, so we pin **2026.3.6**. Apache 2.0. |
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

## 2. Names: two passes with `llb names` between them

John's decision (2026-10-26; `UPSTREAM.md` section 2). MolSanitizer runs twice
and LLB's naming step sits between the passes:

1. **Pass 1** does the chemistry only and writes SMILES. Each row is named
   `INPUT_state.stereo`, after its input name.
2. **`llb names`** (`naming/assign_names.py`, the only naming code) turns those
   rows into `C` + 12-character base id + state + stereo + `C` (DRAP
   `NAMING_CONTRACT.md`).
3. **Pass 2** builds 3D and db2 from the contract names, with
   `--no-stereoisomers --no-neutralize`. It adds no suffix, so `name[-16:]` is
   the whole name.

| field | from |
| --- | --- |
| base id | the parent id (`naming/llb_ids.py`): a ZINC id's last 12 characters verbatim, otherwise `Z` + 11 base62 characters of sha256(name). A damaged `ZINC…` id is rejected. |
| state (contract "protomer") | MolSanitizer's `_k`, renumbered `0..` per parent. It numbers **tautomer × protomer** states, because MolSanitizer gives one index to both. |
| stereo | MolSanitizer's `.j`, renumbered `0..` per state. |

`llb names` also enforces four things:

- every SMILES is stereo-closed, so pass 2 builds exactly the named structure;
- one structure gets one name across the library;
- one fragment per row;
- `invertomer_of` is set for N-H⁺ invertomers, which MolSanitizer 0.7.4 emits
  as separate stereoisomers (checked on John's laptop).

Every pass-2 db2's M-line name is then checked against the requested name.

## 3. Install (Turbo; used from Lighthouse)

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
`LLB_MSANI_ACCOUNT=maom`, `LLB_MSANI_PARTITION=maom-h200`,
`LLB_MSANI_PARALLEL=140`, hard cap `LLB_MSANI_MAX_CPUS=149`). Builds run on **Lighthouse**
(John, 2026-10-26), kept under 150 CPUs: `submit` refuses `--parallel` above the cap
and refuses while another msani array of yours is queued or running.

## 4. Workflow

```bash
llb msani plan H18.smi ... H30.smi -o RUN       # pass 1 shards (500 parents each) + run.json
llb msani submit RUN 1                          # pass 1 array: chemistry, SMILES out
llb msani status RUN                            # done / missing per pass, as array ranges
srun ... llb msani names RUN --per-shard 40     # llb names over ALL of pass 1, then pass 2 shards
llb msani submit RUN 2                          # pass 2 array: 3D + db2
llb msani submit RUN 2 --tasks missing          # resume (either pass): done shards are skipped
llb msani collect RUN --drap-dir $DRAP_DIR      # verify + account + validate_names
```

| step | what it does |
| --- | --- |
| **pass 1** (`build_shard.py --pass 1`) | One parent at a time through MolSanitizer's Python API, which takes pH 7.4 (the CLI takes only integers). The parent goes in under its **input name**, so `llb names` can match the suffixes. MolSanitizer's rows go to `pass1/out/NNNNN_clean.smi`, and a parent that returns nothing gets a reason. |
| **names** (`names_step.py`) | Refuses unless every pass-1 shard is done. Runs `naming/assign_names.py --inputs <the run's sample files> --ms pass1/out/*_clean.smi -o library --ph 7.4+-1`, then cuts `library/library.smi` into pass-2 shards. |
| **pass 2** (`build_shard.py --pass 2`) | Marks db2-incompatible states (`EXOTIC`). Then runs one `msani … --no-neutralize --no-stereoisomers --metal off -3d -f db2 --timing` call per shard. Molecules that call leaves without a db2 or an error are retried **one at a time** under `mol_timeout`. **Every db2 is checked**: names on the M lines, AMSOL charge = formal charge ± 0.05, conformations > 0. It is then gzipped to `db2_run/out/NNNNN/<name>.db2.gz`. |

Each task writes `DONE.json` last, recording versions, host, CPU, SLURM ids,
seconds and the input sha256. A shard with a matching `DONE.json` is skipped. An
interrupted pass-2 shard keeps, and re-checks, every db2.gz it already wrote.

`collect.py` refuses if any of these hold:

- a shard or the naming step is missing;
- an input line is unaccounted for;
- a library name was not built exactly once;
- shards ran different software or options.

It writes:

| file | for |
| --- | --- |
| `db2_verification.csv` | `prepare.py`: only `status == OK` is docked |
| `db2_run/out/NNNNN/*.db2.gz` | `prepare.py --db2-root` default |
| `library/` (from `llb names`) | `library.tsv` for DRAP `collapse_invertomers.py` / `join_library.py --llb-library`; `library.smi` for `validate_names.py`; `dropped.tsv`; `stats.json` |
| `parents_ledger.tsv` | **every input line**: `built`, or `failed` + the step that lost it (`plan:`, `pass1:`, `names:`, `pass2:`) and why |
| `timing.tsv` | per name: MolSanitizer's embedding / AMSOL / sampling / db2 seconds, and ours |
| `collect.json` | counts, fail reasons, attrition per input file, s/parent (pass 1), s/name (pass 2), CPU-h per pass, db2 bytes |
| `provenance.json` | MolSanitizer / RDKit / Python / Open Babel versions, options, LLB commit (and whether the checkout was dirty), CPUs, SLURM jobs, time span |

Tests: `pytest molsanitizer/tests` or `python molsanitizer/tests/test_msani_build.py`.
They run the real `naming/assign_names.py` between the passes and fake only
MolSanitizer (`tests/fake_msani`, which uses 0.7.4's formats). They need RDKit
and pandas.

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

Read off: s/parent (pass 1), s/name (pass 2), names/parent, `MaxRSS` (sacct), db2 size, and `n_h_invertomers`.

**Production H18–H30**: 13 bins × 10,000 = 130,000 parents. Until the smoke test
measures it, assume ~2 names/parent (LLB's median was ~1.9; pH range 1 adds some) and
10–60 s per name at 600 conformers. That brackets 260k names at **~0.7k–4.3k CPU-h**: 5–31 h of
wall time at 140 CPUs. db2_converter cost ~50 s/name. Pass 1 is ~seconds per
parent: 500 parents per task, ~260 tasks, under an hour. `llb names` over the
union is one job (srun, a few minutes). Pass 2 carries the cost. Pick its
`--per-shard` (names per task) so a task takes ~1.5 h:

```
per_shard ≈ 5400 / seconds_per_name      # e.g. 5400/30 = 180 names -> ~1,450 tasks for 260k names
```

`submit` checks the cluster's MaxArraySize and refuses more than 4,500 tasks in one
call. Run it with `--parallel 140` (the default; the cap is 149). Set
the time limit to 3× the smoke run's slowest pass-2 task. Copy the run to Turbo
when done (scratch purges after 60 days).

**Docking cost.** The settings (§6) match LLB's conformer budget (600) and pH
window, so the H18–22 pilot's 1.28 s per ligand-pocket is the right first guess;
the smoke run's names per parent is what moves the ~10⁵ core-h docking bill.

## 6. Chemistry settings (John, 2026-10-26)

| option | set | MolSanitizer default | why |
| --- | --- | --- | --- |
| `--pH` / `--pH-range` | **7.4 / 1** (states at 6.4, 7.4, 8.4, union kept) | 7 / 0 | LLB's window; keeps the imidazole / 2-aminopyridine cations (rule pKa 6) that pH 7 alone drops |
| `--numconfs` | **600** | 2,000 | LLB built ≤ 600; docking time per name scales with it |
| `--max-stereoisomers` | **8** | 8 | |
| `--tautomers` | on | on | the reason for the switch |
| filters (PAINS, unwanted) | off | off | a training set, not a hit list |

The 2D step takes the non-integer pH through MolSanitizer's Python API (its CLI
accepts integers only); `ph_values` in library.tsv reads `6.4-8.4`.
