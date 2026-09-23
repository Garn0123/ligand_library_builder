#!/bin/bash
# submit_qupkake.sh PARENTS_DIR OUT_DIR -- preflight on this node, then submit
# the QupKake array with the account, partition and environment from
# config/hpc.env.
#
#   run_qupkake/submit_qupkake.sh parents protomers
#   run_qupkake/submit_qupkake.sh parents protomers --time=04:00:00   # extra sbatch args
#
# The preflight activates the configured QupKake env HERE and checks that the
# CLI, the Python package and XTBPATH all resolve, so a wrong path fails once on
# the login node instead of once per array task. Resubmitting resumes: finished
# shards are skipped.

set -euo pipefail

[[ $# -ge 2 ]] || { sed -n '2,12p' "$0" >&2; exit 64; }
PARENTS_DIR="$(cd "$1" && pwd)"; OUT_DIR="$(mkdir -p "$2" && cd "$2" && pwd)"; shift 2

SMP_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export SMP_ROOT
# shellcheck source=common/env.sh
source "$SMP_ROOT/common/env.sh"
smp_load_config

META="$PARENTS_DIR/shards.tsv.meta"
[[ -f "$META" ]] || { echo "FATAL: $META missing -- run prepare_parents.py first" >&2; exit 2; }
N="$(awk -F= '$1=="n_shards"{print $2}' "$META")"
[[ "$N" =~ ^[0-9]+$ && "$N" -gt 0 ]] || {
    echo "FATAL: no n_shards in $META; keys present:" >&2; cut -d= -f1 "$META" | sed 's/^/  /' >&2; exit 2; }

echo "== preflight (config: $SMP_CONFIG)"
# Every check exits explicitly: `set -e` is SUSPENDED inside a subshell on the
# left of `||`, so a failing command there would otherwise fall through and
# the array would be submitted with a broken environment.
( smp_activate QUPKAKE || exit 2
  exe="$(smp_require_exe QUPKAKE_EXE qupkake)" || exit 2
  echo "  qupkake CLI   $exe"
  python -c 'import qupkake, rdkit; print("  import        qupkake, rdkit", rdkit.__version__)' \
      || { echo "FATAL: qupkake or rdkit not importable in ${QUPKAKE_ENV}" >&2; exit 2; }
  if [[ -n "${XTBPATH:-}" ]]; then
      [[ -e "$XTBPATH" ]] || { echo "FATAL: XTBPATH=$XTBPATH does not exist" >&2; exit 2; }
      echo "  XTBPATH       $XTBPATH"
  else
      echo "  XTBPATH       (unset: QupKake's bundled xtb)"
  fi
) || exit 2

mkdir -p "$SMP_ROOT/run_qupkake/logs"
cd "$SMP_ROOT/run_qupkake"          # #SBATCH --output=logs/... is relative to here
echo "== submitting $N shard(s), account=$SMP_ACCOUNT partition=$SMP_PARTITION"
sbatch -A "$SMP_ACCOUNT" -p "$SMP_PARTITION" \
    --array="0-$((N - 1))%${SMP_ARRAY_PARALLEL:-50}" \
    --export=ALL,SMP_ROOT="$SMP_ROOT",SMP_CONFIG="$SMP_CONFIG",PARENTS_DIR="$PARENTS_DIR",OUT_DIR="$OUT_DIR" \
    "$@" slurm/qupkake_array.sbatch
