# Installing LLB on a new cluster: design (for review, 2026-10-20)

Goal: getting `llb pipeline` running on a new SLURM cluster takes a handful of
commands, and every problem shows up in minutes, not after a 1–2 h smoke test.
Third-party tools may still need installing, but LLB's own glue must be simple.

Status: design agreed with John on 2026-10-20, in an interview; the decisions are
numbered Q1–Q29 below. Nothing here is implemented yet.

## Why

Each move (Great Lakes → Lighthouse → TAMU ACES, October 2026) cost about a day,
and the failures appeared one at a time:

| what broke | why |
| --- | --- |
| `DB2C_SETUP="source ~/env/bcl.sh"` | `bcl.sh` lives in `~`, which is not shared between clusters; its module names (`Chemistry`, `gcc/8.2.0`) don't exist on ACES (`GCC/11.3.0`, no Chemistry) |
| UNICON / conformator | wrappers in `~/bin` that `singularity exec` an `ubuntu2204.sif` with bind mounts, all in `~` |
| ZBH licenses | per-user, annual, kept in `~`; "license expired" on the new cluster |
| `PKA_ENV`, `PREP_ENV` | Lighthouse needed the absolute prefix, not the env name |
| AMSOL | needs `libgfortran.so.5` at runtime, which comes from a host module |
| 3 conda envs + xtb 6.4.1 + MolGpKa + DRAP + BCL | rebuilt by hand on each cluster; BCL's build is painful (OOM at link, bzip2 SONAME, static-init crash) |
| account, partition, PATH, LLB_CONFIG | in `~/.bashrc`, retyped |

Checks existed, but only inside the db2 stage (`submit_db2.sh`,
`db2_array.sbatch`), so a broken QupKake or MolGpKa env surfaced hours in.

## The shape

**Audience (Q1).** Anyone with a SLURM cluster, Apptainer *or* SingularityCE,
and shared storage, including people outside the lab. No lab paths in the docs or
defaults; lab-specific values live in cluster profiles (Q21).

**Packaging (Q2).** One public **toolchain image** holds everything third-party
that can be redistributed. LLB's own code stays a git checkout on the host and is
bind-mounted into the image, so changing LLB never needs a rebuild. The host
needs only SLURM, a container runtime and the image file: no module names, no
setup hooks, no libgfortran, no CONDA_BASE, no nested ZBH container.

**One site directory, nothing in `~` (Q3, Q9).**

```
$LLB_SITE/                       e.g. /scratch/group/<project>/llb (never ~)
  ligand_library_builder/        the checkout, a tagged release (Q20)
  llb-toolchain_<version>.sif    the image
  hpc.env                        cluster settings, shared by the site's users
  users/<user>.env               per-user: ZBH_LICENSE, LLB_MAIL, optional LLB_ACCOUNT
  zbh/                           UNICON + conformator binaries (each user brings their own, or the site shares them)
  licenses/<user>/               default home for each user's ZBH license file
```

`llb` reads `$LLB_SITE/hpc.env` (found as `../hpc.env` from its own checkout),
then `users/$USER.env`. `LLB_CONFIG` still overrides that. The only line anyone
puts in `~/.bashrc` is optional: `export PATH=$LLB_SITE/ligand_library_builder/bin:$PATH`.
Moving clusters means recreating one directory.

## The toolchain image (Q8, Q11, Q13, Q17, Q23, Q24)

`oras://ghcr.io/garn0123/llb-toolchain:YYYY.MM.N`. It is built from
`container/llb.def` in GitHub Actions, with Apptainer, whenever `container/` or
`environments/` change. GHCR stores container images for free. The fallback is a
local `apptainer build --fakeroot`, which ACES allows on compute nodes.

| inside | how it is pinned | licence |
| --- | --- | --- |
| Ubuntu 22.04 base + micromamba | base image digest | — |
| `qupkake` env (torch, PyG, Lightning; QupKake f7a294e via pip) | `environments/qupkake.linux-64.lock` | BSD-3 |
| xtb 6.4.1, built from source in `/opt/xtb`; `XTBPATH=/opt/xtb/bin/xtb` | release tag + sha256 | LGPL-3 |
| `db2_converter` env (RDKit 2024.03.5, openbabel, ambertools, posebusters 0.2.9; db2_converter 63d6656) | `environments/db2_converter.linux-64.lock` | — |
| AMSOL 7.1, DOCK-patched (`amsol7.1_patched` as shipped in db2_converter, SIZES.i raised) | db2_converter commit | Apache-2.0, plus NOTICE stating the patch |
| BCL 4.3.1 static (`bcl-apps-static.exe` + the data `BCL_ROOT` needs) in `/opt/bcl` | release asset + sha256 (below) | MIT |
| `llb_pka` env (RDKit 2024.03.5, CPU torch, dimorphite_dl) | `environments/llb_pka.linux-64.lock` | — |
| MolGpKa 4dc8352 (code, weights, SMARTS) in `/opt/MolGpKa` | commit; `molgpka.py` already checks the weights' sha256 | MIT |
| `conda` stub → micromamba | — | — |
| `/opt/llb/manifest.json` | every pin above, the lockfile hashes, the build date | — |
| `/opt/llb/NOTICE`, licences | — | — |

Not inside: **UNICON and conformator** (ZBH; the academic licence is annual,
per user, non-transferable and not redistributable), DOCK/db2tool (Q6), and
LLB's own code.

**Lockfiles.** `environments/*.yml` stay the human-edited specs. `conda-lock
--kind explicit -p linux-64` turns them into committed `*.linux-64.lock` files,
and the image installs from those with micromamba, which also honours the
`# pip` lines. QupKake's own `environment.yaml` (python 3.9, torch 2.1,
rdkit 2022.03) becomes our `environments/qupkake.yml`.

**QupKake (Q17).** At import, QupKake runs `conda list xtb` and uses `$XTBPATH`
unless conda's xtb 6.4.1 is installed. Without a `conda` binary that call raises
FileNotFoundError. The image therefore ships a `conda` stub that forwards to
micromamba and sets `XTBPATH`. QupKake itself stays unpatched.

**BCL (Q23, Q24).** BCL isn't built in CI. John's working static build is
packed once as `bcl-4.3.1-static-linux-x86_64.tar.gz`, attached to an LLB
GitHub release (`toolchain-assets`), and the image build downloads it and
checks its sha256. `container/bcl/BUILD.md` records how it was built, including
the OOM, bzip2 and static-init workarounds, so it can be rebuilt. The image
sets `BCL_ROOT=/opt/bcl`, puts the binary on PATH, and **never** adds BCL's
`lib/` to LD_LIBRARY_PATH: that is what shadowed the env's libstdc++/libgomp. A
`sed` in `%post` points db2_converter's `config.py` BCLBASE/BCL at `/opt/bcl`,
the same edit db2_converter's README tells every user to make. If the static
binary turns out to need `lib/`, a wrapper sets it for BCL's own process only.
*To measure during implementation:* the install's size. Above ~2 GB it gets its
own image layer or a separate download.

**Versions (Q13).** The checkout holds `container/TOOLCHAIN_VERSION`, the image
version this code expects. `doctor` and `pipeline start` fail on a mismatch
(`--allow-toolchain-mismatch` overrides). REPORT.md adds the image sha256 and
`manifest.json` to what it already records (per-env `versions.py`, LLB commit).
The ZBH versions come from the shims at run time.

## How a step runs (Q12, Q15, Q16)

With `LLB_IMAGE` set (the default after `llb setup`), `common/with_env.sh TOOL
cmd` becomes

```
$LLB_CONTAINER exec --bind <auto> $LLB_IMAGE micromamba run -n <env> cmd
```

- `LLB_CONTAINER` is whichever of `apptainer` or `singularity` is found, after
  an optional `LLB_CONTAINER_SETUP` (e.g. `module load singularity` on Great
  Lakes, which runs SingularityCE 4.4.1). Both read the same SIF.
  `container/llb.def` uses only sections both runtimes support (`%post`,
  `%files`, `%environment`, `%labels`). `llb` sets both `APPTAINER_CACHEDIR`
  and `SINGULARITY_CACHEDIR` to node scratch.
- `<auto>` binds the checkout, `$LLB_SITE`, `$PWD`, `ZINC_DATA`, node scratch
  (`SLURM_TMPDIR` / `TMPDIR`), the ZBH directory and the license, plus anything
  in `LLB_BINDS`.
- LLB's shims for `unicon` and `conformator` live in the repo
  (`container/shims/`), are put on PATH inside the image, and run the user's
  ZBH binaries with the user's license (`ZBH_LICENSE`, default
  `$LLB_SITE/licenses/$USER/`) (Q10). The ZBH binaries run in the Ubuntu 22.04
  image directly, so the `ubuntu2204.sif` wrappers go away.
- `<TOOL>_SETUP` / `_ENV` / `_EXE` stay as per-tool overrides, e.g. to point one
  tool at a native env to test a new QupKake.
- **Nothing touches the image on a login node.** ACES has no container runtime
  there. `pipeline start` only checks `hpc.env` and submits. A new first stage,
  **`preflight`** (a ~5 min job on `LLB_PARTITION`), runs `doctor` and
  `versions.py` and gates everything after it. `submit-db2` and
  `submit-qupkake` drop their login-node activation, and every array task still
  checks itself.
- The sbatch scripts and the pipeline's control jobs stay on the host and need
  only bash and python3.

## Commands

### `llb doctor` (Q4, Q25)

This is the only place checks live. It prints one line per check: pass/fail and,
on failure, the fix.

| flag | where | what |
| --- | --- | --- |
| *(none)* | here (login node OK) | host side: runtime found, image present + version, `hpc.env` keys, account/partition valid (`sacctmgr`, `sinfo`), site and run directory writable, `wget` / `parallel` (Q28); says "run `--node` for the rest" |
| `--node` | 2-min `srun` on `LLB_PARTITION` | everything, inside the image, on a compute node |
| `--stage db2` (`qupkake`, `pka`, `prep`) | here | just that stage's checks; what `db2_array.sbatch`, `submit-db2` and the `preflight` stage call |
| `--smoke` | ~3-min job | 2 molecules through `build_ligand -m conformator_bcl --PBfilter`, the real default |
| `--no-zbh` | CI | skip the ZBH checks |

Checks per stage, reusing what exists:

| stage | checks |
| --- | --- |
| prep / pka | RDKit 2024.03.5 importable; MolGpKa weights' sha256 (`molgpka.py`); `versions.py` runs |
| qupkake | `import qupkake` (exercises the `conda list` path); `$XTBPATH` runs `xtb --version` = 6.4.1 |
| db2 | `build_ligand` on PATH; `import rdkit, db2_converter`; `check_bcl.py`; `import posebusters` when the filter is on; AMSOL `ldd` has no "not found"; `conformator` and `unicon` shims run; UNICON license date via the `Your license is valid until:` line that db2_converter's `check_UNICON_license` parses (warn under 30 days); the ZBH download/renewal URL on failure |
| host | runtime + version recorded; image digest = the one `hpc.env` pins |

The preflight blocks in `submit_db2.sh` and `db2_array.sbatch` are replaced by
calls to `doctor --stage db2`.

### `llb setup` (Q4, Q5, Q15, Q21)

Re-runnable; it never overwrites an edited `hpc.env` and shows a diff instead.

1. Recognize the cluster (hostname / `sinfo` cluster name) and start from its
   profile, `config/sites/{greatlakes,lighthouse,aces}.env`: the container module,
   recommended storage, partition hints. Profiles hold no accounts or user paths.
2. Choose `$LLB_SITE`. Suggest group/project storage (`$SCRATCH`,
   `/scratch/group/...`, Turbo), never `~`.
3. Clone or locate the checkout at the latest release tag (`--ref BRANCH` for
   development) (Q20).
4. Pull the image (`$LLB_CONTAINER pull oras://...`), inside an `srun` when the
   login node has no runtime.
5. Ask for the ZBH binaries' directory and the license file; write
   `users/$USER.env`.
6. Write a minimal `hpc.env`: `LLB_ACCOUNT` from `sacctmgr show assoc
   user=$USER`, `LLB_PARTITION` from `sinfo` (cores, memory, time limit, allowed
   accounts), a suggested `LLB_ARRAY_PARALLEL`, `LLB_MAIL`, `ZINC_DATA`, plus
   `LLB_IMAGE`. It asks only when there is more than one choice. Each key gets
   a comment saying what to change and when.
7. Run `doctor --node`, then `doctor --smoke`.

`llb setup --native` (Q14) builds the same envs from the same lockfiles with
micromamba, and compiles xtb into `$LLB_SITE/native/`, for a cluster without a
container runtime. It is documented as best effort, and CI doesn't test it.

## Other decisions

- **Conformer methods (Q19′).** Supported: `conformator_bcl` with `--PBfilter`
  (the default since d84b045) and `conformator`. `rdkit` is allowed as the
  option needing no license. `ccdc` and `confgenx` are unsupported, though the
  overrides still let you wire one in.
- **DRAP's validator (Q26).** DRAP isn't public. `stage4/validate_names.py` is
  stdlib-only, so it is copied into LLB as `run_qupkake/validate_names.py`,
  with a header naming the DRAP commit. `DRAP_DIR` becomes an optional override,
  and a test compares the copy to DRAP whenever `DRAP_DIR` is set.
- **DOCK / db2tool (Q6).** Stays out of the image as the
  `DOCK_SETUP` / `DB2TOOL_EXE` override; `doctor` checks it only when it is set.
- **Downloads (Q28).** `wget` and GNU `parallel` stay on the host (login or
  transfer nodes, with no runtime on ACES). If `parallel` is missing, `llb
  fetch` falls back to `xargs -P`.
- **Old configs (Q27).** An `hpc.env` without `LLB_IMAGE` keeps today's
  behaviour: native envs, setup hooks, `LLB_CONFIG`. All libraries restart on
  v1.0 anyway, but nothing forces a cutover. Once the image path is proven,
  `hpc.env.example` is rewritten so the knobs read as advanced overrides.
- **Resuming on another cluster (Q7).** `pipeline resume` always uses the
  *current* site's config; absolute paths recorded at start are provenance only,
  and REPORT.md notes the move. No automatic data copying.
- **One library across several clusters (Q29).** Every site pins the same LLB
  release and image digest. `pipeline start` records both, and parts may be
  combined only when they match. *How* to split and merge (ledger, attrition,
  db2 set, names unique across parts) is a separate design and gets its own doc.

## CI (Q18)

| job | when | what |
| --- | --- | --- |
| `tests` (exists) | every push | pytest in llb_pka, unchanged |
| `image` | `container/`, `environments/` change; tags | lock check (specs ↔ lockfiles), build the SIF, push to GHCR on tags |
| `image-doctor` | after `image` | `llb doctor --no-zbh` inside the image under **Apptainer 1.3 and SingularityCE 4.4** |
| `smoke` | after `image` | 5 molecules: parents → triage → QupKake (real xtb) → merge → names → check-stereo, inside the image, without SLURM |

Not in CI: the db2 build. db2_converter needs UNICON even with `--method
rdkit`, and the ZBH license is non-transferable, so no license goes into
CI secrets. `llb doctor --smoke` on the cluster covers it, and `llb setup` runs it.

## Rollout (Q22′)

No library is running (H23–30 stopped), so there are no compatibility
constraints from a live run. Phases in order:

1. **`llb doctor` + the `preflight` stage**, against today's native setup;
   fold the db2 preflight into `doctor --stage db2`.
2. **Site layout**: `$LLB_SITE`, `users/$USER.env`, `ZBH_LICENSE`, ZBH shims,
   the copied `validate_names.py`, the cluster profiles.
3. **Image**: `container/llb.def`, lockfiles, the BCL release asset +
   `BUILD.md`, `with_env.sh` image mode, the CI jobs, `TOOLCHAIN_VERSION`, image
   digest in REPORT.md.
4. **`llb setup`** (+ `--native`), README "Install" replacing "Setup on a new
   machine", `docs/sites/`.

Then merge `feat/pipeline` to `main` and tag **v1.0 only after `llb setup` +
`llb doctor --smoke` pass on Great Lakes, Lighthouse and ACES**. Every library
restarts from v1.0.

## What a new user does, when this is done

```bash
git clone https://github.com/Garn0123/ligand_library_builder /tmp/llb-bootstrap
/tmp/llb-bootstrap/bin/llb setup            # asks for: site dir, account (if >1), ZBH dir + license
export PATH=$LLB_SITE/ligand_library_builder/bin:$PATH
llb pipeline start my.smi -o run1           # preflight stage runs doctor first
```

Manual steps that remain: downloading UNICON + conformator from ZBH,
activating a license, and having an account on the cluster.

## Open items

- BCL static install size (measure in phase 3).
- `container/shims/`: how ZBH reads a license from a path other than the
  default. John: "can be specified from anywhere"; confirm the exact flag or
  env var while writing the shim.
- Follow-up design: splitting one library across clusters and merging it back (Q29).
