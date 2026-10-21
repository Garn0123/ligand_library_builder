# The LLB toolchain image

`llb.def` builds one Apptainer/SingularityCE image (SIF) that holds every
third-party tool LLB runs and is allowed to redistribute. LLB's own code is not
in it: `llb` bind-mounts its checkout and runs each step as
`apptainer exec ... llb.sif micromamba run -n <env> <cmd>`, so changing LLB never
needs a rebuild. Design: `docs/INSTALL_DESIGN.md`, "The toolchain image".

| in the image | where | pinned by |
| --- | --- | --- |
| Ubuntu 22.04 + micromamba 2.8.1 | `/`, `/usr/local/bin/micromamba` | base image digest, micromamba sha256 |
| `qupkake` env: QupKake f7a294e, torch 2.1 CPU, PyG 2.4, Lightning 2.1, RDKit 2022.03.3 | `/opt/conda/envs/qupkake` | `environments/qupkake.linux-64.lock` |
| xtb 6.4.1, built from source; `XTBPATH=/opt/xtb/bin/xtb` | `/opt/xtb` | release tarball sha256 |
| `conda` stub → micromamba (QupKake's import-time `conda list xtb`) | `/usr/local/bin/conda` | — |
| `db2_converter` env: RDKit 2024.03.5, AmberTools, Open Babel, PoseBusters 0.2.9 | `/opt/conda/envs/db2_converter` | `environments/db2_converter.linux-64.lock` |
| db2_converter 63d6656 (editable install) + AMSOL 7.1 DOCK-patched | `/opt/db2_converter` | commit |
| BCL 4.3.1 static; `BCL_ROOT=/opt/bcl` | `/opt/bcl` | GHCR artifact + sha256 |
| `llb_pka` env: RDKit 2024.03.5, torch 2.4 CPU, PyG, dimorphite_dl | `/opt/conda/envs/llb_pka` | `environments/llb_pka.linux-64.lock` |
| MolGpKa 4dc8352; `MOLGPKA_DIR=/opt/MolGpKa` | `/opt/MolGpKa` | commit + weights' sha256 |
| manifest, NOTICE, licences, smoke test | `/opt/llb/` | — |

Not in the image: UNICON and conformator (ZBH; licensed per user, not
redistributable; LLB's shims run each user's own copy), DOCK/db2tool.

db2_converter's `config.py` is edited at build time to `BCLBASE=/opt/bcl/bin`,
`BCL="bcl-apps-static.exe molecule:ConformerGenerator"`; `UNICON_EXE` and
`CONF_EXE` stay the bare names `unicon` and `conformator`, found on PATH. BCL's
`lib/` is never on `LD_LIBRARY_PATH`.

## Check an image

```bash
apptainer exec llb.sif /opt/llb/smoke.sh          # one line per check, exit 1 on any FAIL
apptainer exec llb.sif cat /opt/llb/manifest.json  # every pin, lockfile sha256s, build date
apptainer inspect --labels llb.sif                 # org.llb.toolchain.version, org.llb.bcl.installed
```

`singularity` works the same way on the same SIF (CI tests both).

## Build it locally

On a compute node (ACES allows `--fakeroot` there; login nodes may have no
runtime), from the **repository root**, since `llb.def`'s `%files` paths are
relative to it:

```bash
oras pull ghcr.io/garn0123/llb-assets:bcl-4.3.1 -o container/assets   # the BCL tarball, ~2.7 GB
export APPTAINER_TMPDIR=$TMPDIR APPTAINER_CACHEDIR=$TMPDIR/apptainer-cache   # node scratch, not ~
apptainer build --fakeroot llb.sif container/llb.def
```

`oras` is a single static binary: https://github.com/oras-project/oras/releases
(CI uses v1.3.2). The pull needs no login once the package is public. The build
downloads packages from conda-forge, PyPI and GitHub, so the node needs outbound
network. Expect it to take the better part of an hour (an estimate), mostly the conda envs and xtb's compile.

Without the BCL tarball (or while `BCL_SHA256` in `llb.def` is still `TODO`) the
build stops at once with a `FATAL`. To build an image without BCL on purpose,
for testing everything else:

```bash
touch container/assets/ALLOW_MISSING_BCL
apptainer build --fakeroot llb.sif container/llb.def    # loud WARNING; manifest says bcl.installed=false
```

Such an image can't run `-m conformator_bcl` (`check_bcl.py` fails), and CI
never releases one.

**Size.** The three envs are ~1.4 GB of compressed packages (about 4–5 GB
installed: two CPU torch builds with MKL, AmberTools), and the BCL install is
~2.7 GB. Expect a SIF of roughly 3–5 GB (estimate; the first CI build prints the
real figure in its job summary), and give the build ~3× that in
`APPTAINER_TMPDIR`. If BCL dominates, the design's fallback is a separate
download for it (docs/INSTALL_DESIGN.md, "BCL").

## Versions

- `container/TOOLCHAIN_VERSION` (`YYYY.MM.N`) is the image version this checkout
  expects. The build writes it into `/opt/llb/manifest.json` and the
  `org.llb.toolchain.version` label; `llb doctor` and `pipeline start` are to
  compare it with the image they are given (Q13, not implemented yet).
- Bump it whenever `container/` or `environments/` change in a way that is
  released: a lockfile, a pin in `llb.def`, BCL.
- Releasing: tag the repo `v*`. CI (`.github/workflows/image.yml`) builds, runs
  the smoke checks under Apptainer 1.3 and SingularityCE 4.4, and pushes
  `oras://ghcr.io/garn0123/llb-toolchain:<TOOLCHAIN_VERSION>`. A tag build refuses
  to skip BCL. Pull it with
  `apptainer pull llb-toolchain_<version>.sif oras://ghcr.io/garn0123/llb-toolchain:<version>`.
- Pushes and PRs touching `container/` or `environments/` build and test the
  image but don't publish it.

## Files

| file | what |
| --- | --- |
| `llb.def` | the image definition; every pin is at the top of `%post` |
| `TOOLCHAIN_VERSION` | the image version |
| `smoke.sh` | the image's self-test (`/opt/llb/smoke.sh`), run at the end of the build and in CI |
| `NOTICE`, `licenses/` | third-party licences, copied to `/opt/llb/` |
| `assets/` | build-context drop for the BCL tarball (git-ignored contents) |
| `bcl/BUILD.md` | how the static BCL tarball was built and how to publish it |
| `../environments/` | env specs and lockfiles (`environments/README.md`) |
