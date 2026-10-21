# Building and publishing the static BCL tarball

The toolchain image does not build BCL (docs/INSTALL_DESIGN.md, Q23/Q24): the
build is slow, memory-hungry and fragile. Instead, one working static build is
packed once, published to GHCR as an ORAS artifact, and `container/llb.def` copies
it in and checks its sha256.

| | |
| --- | --- |
| source | https://github.com/BCLCommons/bcl, tag `BCL-v4.3.1-linux` (MIT) |
| artifact | `ghcr.io/garn0123/llb-assets:bcl-4.3.1` |
| file | `bcl-4.3.1-static-linux-x86_64.tar.xz` (the static install is ~2.7 GB, over GitHub's 2 GiB release-asset limit, hence GHCR) |
| checksum | `BCL_SHA256` in `container/llb.def` (**TODO** until this is done once) |

This file is a template. No build notes for BCL were found in the project tree
(searched `*.md`, `*.txt`, `*.sh` under the research directory on 2026-10-21),
so every command marked **TODO(John)** must be replaced with the one actually used
for the working build.

## 1. Build

On a Linux x86_64 node with plenty of memory (the link is the problem, below):

```bash
git clone https://github.com/BCLCommons/bcl.git && cd bcl
git checkout BCL-v4.3.1-linux
git rev-parse HEAD                     # record it here: TODO(John)

# TODO(John): modules / compiler used, e.g. `module load gcc/...` and cmake version
# TODO(John): the configure command for the static apps binary, e.g.
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ... 
# TODO(John): the build command, with FEW jobs (see "OOM at link")
#   cmake --build build --target bcl-apps-static -j 2
```

### Known problems

- **OOM at link.** Linking the static apps binary takes more memory than a
  normal node share, and with many parallel jobs the linker gets OOM-killed
  (typically `ld terminated with signal 9`). Build with few jobs
  (`-j 2`, or compile with `-j N` and then link alone with `-j 1`), on a node or
  allocation with a large `--mem`. TODO(John): the job count and memory that worked.
- **bzip2 SONAME.** BCL links bzip2, whose SONAME differs between distributions
  (Debian/Ubuntu `libbz2.so.1.0`, RHEL-family `libbz2.so.1`), so a dynamically
  linked BCL built on one family can fail to load on the other. A fully static binary avoids it by linking
  `libbz2.a`. TODO(John): how the static libbz2 was supplied (system `bzip2-static`,
  or built from source and pointed at with `-DBZIP2_LIBRARIES=...`).
- **Static-init crash.** The static binary crashed at start-up, in static
  initialisation. (Our guess at the cause, unconfirmed: BCL's many global
  registries are constructed in a different order once everything is linked into
  one binary.) TODO(John): the symptom as seen, for searching, and the workaround
  (link order, a flag, a patch).

Check the result before packing:

```bash
file build/.../bcl-apps-static.exe     # "statically linked"
ldd  build/.../bcl-apps-static.exe     # "not a dynamic executable"
build/.../bcl-apps-static.exe -help    # runs, no crash (the static-init check)
build/.../bcl-apps-static.exe molecule:ConformerGenerator -help
```

## 2. Pack: the layout llb.def expects

Exactly one top-level directory, the binary at `bin/bcl-apps-static.exe`, and the
data BCL::Conf needs at run time beside it (`BCL_ROOT=/opt/bcl` points at the
top-level directory once unpacked):

```
bcl-4.3.1/
  bin/bcl-apps-static.exe       the only thing in bin/: no shared libraries
  ...                           data BCL needs (TODO(John): which directories, e.g.
                                histogram/, model/, rotamer_library/ for ConformerGenerator)
  LICENSE                       BCL's MIT licence
```

`bin/` must not hold `*.so*`: db2_converter prepends `BCLBASE` (`/opt/bcl/bin`) to
`LD_LIBRARY_PATH` for its BCL call, and BCL's own libraries on that path are what
shadowed the envs' libstdc++/libgomp before. llb.def refuses a tarball that has
them. If the binary turns out to need anything from a `lib/`, add a wrapper in
`bin/` that sets `LD_LIBRARY_PATH` for BCL's own process only.

```bash
mkdir -p stage/bcl-4.3.1/bin
cp build/.../bcl-apps-static.exe stage/bcl-4.3.1/bin/
# TODO(John): copy the data directories into stage/bcl-4.3.1/
cp LICENSE stage/bcl-4.3.1/
du -sh stage/bcl-4.3.1                                         # record it: TODO
XZ_OPT='-T0 -9' tar -C stage -cJf bcl-4.3.1-static-linux-x86_64.tar.xz bcl-4.3.1
sha256sum bcl-4.3.1-static-linux-x86_64.tar.xz                 # -> BCL_SHA256 in llb.def
```

## 3. Publish to GHCR with ORAS

Needs the `oras` CLI (CI pins v1.3.2) and a GitHub token with `write:packages`
(a classic personal access token; `GITHUB_TOKEN` in Actions works too):

```bash
echo "$GHCR_TOKEN" | oras login ghcr.io -u Garn0123 --password-stdin
oras push ghcr.io/garn0123/llb-assets:bcl-4.3.1 \
    bcl-4.3.1-static-linux-x86_64.tar.xz:application/x-xz
```

Then make the package public, so CI and every cluster can pull without a login:
GitHub → your profile → Packages → `llb-assets` → Package settings → Change
visibility → Public. (Optionally "Connect repository" to ligand_library_builder.)

Check the round trip from a clean directory:

```bash
oras pull ghcr.io/garn0123/llb-assets:bcl-4.3.1 -o /tmp/bcl-check
sha256sum /tmp/bcl-check/bcl-4.3.1-static-linux-x86_64.tar.xz   # = BCL_SHA256
```

## 4. Switch the image over

1. Put the sha256 in `container/llb.def` (`BCL_SHA256=...`, replacing `TODO`).
2. Bump `container/TOOLCHAIN_VERSION` if an image was already released.
3. Push: CI pulls the artifact, the build checks the sha256, and the smoke checks
   now require `/opt/bcl/bin/bcl-apps-static.exe` to run. From then on tag builds
   (which never skip BCL) can release.

A new BCL build gets a new tag (`bcl-4.3.1-2`, ...) rather than overwriting
`bcl-4.3.1`, so an old `llb.def` keeps pointing at the bytes its sha256 names.
