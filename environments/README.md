# Environments

| spec | used by |
| --- | --- |
| `llb.yml` | `PREP` in `hpc.env`: `llb mol2-library`, `llb collapse-invertomers`, the `zinc3d-*` steps, and the tests |

`llb mol2` runs in `DOCK`, which is not a conda env: it is `db2tool` from
dock6_claude, loaded by `DOCK_SETUP` (see `config/hpc.env.example`). Download,
sampling and chunking need only the system `python3`.

```bash
mamba env create -p $PREFIX -f environments/llb.yml
```

The QupKake, db2_converter and MolGpKa envs, their conda-lock lockfiles and the
toolchain image that installed them were for library construction. They are at
tag `llb-full-2026-10` (`environments/`, `container/`). Ligand construction is
now MolSanitizer's job (`molsanitizer/`).
