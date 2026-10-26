# Environments

| spec (edit this) | lockfile (generated, committed) | used by |
| --- | --- | --- |
| `qupkake.yml` | `qupkake.linux-64.lock` | QupKake stage (`QUPKAKE`) |
| `db2_converter.yml` | `db2_converter.linux-64.lock` | db2 build, check-stereo, verify-db2 (`DB2C`) |
| `llb_pka.yml` | `llb_pka.linux-64.lock` | pKa triage, tests (`PKA`) |

The `*.yml` files are the specs people edit; each has a header saying what is
pinned and why. The `*.linux-64.lock` files are conda-lock "explicit" lockfiles:
every conda package as a URL plus md5, and the pip packages as `# pip name @ url`
lines with their sha256. The toolchain image (`container/llb.def`) installs only
from the lockfiles, so an image built today and one built next year hold the
same packages.

## Regenerate a lockfile after editing its spec

From the repository root, with conda-lock 4.0.3 (the version that made the
committed files) and network access:

```bash
conda-lock --kind explicit -p linux-64 -f environments/qupkake.yml       --filename-template 'environments/qupkake.{platform}.lock'
conda-lock --kind explicit -p linux-64 -f environments/db2_converter.yml --filename-template 'environments/db2_converter.{platform}.lock'
conda-lock --kind explicit -p linux-64 -f environments/llb_pka.yml       --filename-template 'environments/llb_pka.{platform}.lock'
```

This re-solves from scratch, so it can move every unpinned package to its newest
build; review the diff before committing. Commit the spec and its lockfile
together: CI (`.github/workflows/image.yml`, job `lock-check`) fails when a
lockfile's `# input_hash` no longer matches its spec. Comment-only edits to a spec
do not change the hash.

To check without re-rendering (what CI runs): add `--check-input-hash`. conda-lock
still solves, but leaves the lockfile untouched when the hash matches, so
`git diff --exit-code environments/` afterwards tells you whether it is current.

conda-lock itself can live in a throwaway env, e.g.
`micromamba create -p /tmp/cl -c conda-forge conda-lock=4.0.3` (it needs conda,
mamba or micromamba on PATH to solve; pass `--micromamba` to use micromamba).

## Install an env from a lockfile, without the image

What `container/llb.def` does, and what `llb setup --native` will do:

```bash
grep -v '^# pip ' environments/qupkake.linux-64.lock > /tmp/conda.lock
micromamba create -y -p $PREFIX -f /tmp/conda.lock
grep '^# pip ' environments/qupkake.linux-64.lock | sed 's/^# pip //' > /tmp/pip.txt
$PREFIX/bin/python -m pip install --no-deps -r /tmp/pip.txt
```

The pip lines are installed separately, with `--no-deps` (conda-lock has already
resolved them), so the result doesn't depend on whether a given micromamba
version handles `# pip` lines in explicit files. The `git+` line (QupKake) needs
`git`.

Not in any lockfile: db2_converter itself (a clone at 63d6656, installed with
`pip install --no-deps -e`, since its setup.py would leave AMSOL and its data
files out of a normal install), xtb (built from source), BCL, MolGpKa. See
`container/llb.def`.
