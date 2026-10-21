# Cluster profiles

One file per cluster, `<cluster>.env`: plain bash `KEY=value` like `hpc.env`,
holding only facts about the cluster (container module, storage, partition
hints). Never accounts, users or user paths: those go in `$LLB_SITE/hpc.env`
and `users/$USER.env`.

Nothing reads these yet. `llb setup` (docs/INSTALL_DESIGN.md, phase 4) will
start a new `hpc.env` from the matching profile. Lines marked `# check:` are
not yet verified on that cluster.
