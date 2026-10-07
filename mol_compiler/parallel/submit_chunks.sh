#!/bin/bash
# submit_chunks.sh RUNDIR INPUT_DIR [options] -- chunk a db2 library from a run directory.
#
#   llb chunk my_run my_run/db2_run/out
#   llb chunk my_run my_run/db2_run/out --target-per-bin 20000 --shards 50
#
# Runs parallel/submit.slurm with WORK_DIR=RUNDIR/work, OUT_DIR=RUNDIR/chunks and
# its logs in RUNDIR/work/logs, so the checkout is never edited or written to.
# Account and partition come from config/hpc.env (LLB_ACCOUNT, LLB_PARTITION).
#
# Options (submit.slurm's settings; defaults there):
#   --target-per-bin N   molecules per finished chunk        (50000)
#   --shards N           collect array width                 (200)
#   --bins N             pin the chunk count, skip the estimate
#   --mode stride|greedy (stride)     --weight count|bytes|lines:C (count)
#   --collect-parallel N %N throttle on the collect array    (LLB_ARRAY_PARALLEL)
#   --collect-time T  --assemble-time T  --plan-time T  --mem M
#
# RUNDIR/chunks.env records what was chunked. Re-running with the same INPUT_DIR
# re-plans and resubmits from scratch; a different INPUT_DIR is refused -- give a
# new RUNDIR, or the manifest would describe a mix of two libraries.

set -euo pipefail

[[ $# -ge 2 ]] || { sed -n '2,/^$/p' "$0" >&2; exit 64; }
RUN="$1"; IN="$2"; shift 2
while [[ $# -gt 0 ]]; do
    case "$1" in
        --target-per-bin)   export TARGET_PER_BIN="$2"; shift 2 ;;
        --shards)           export SHARDS="$2"; shift 2 ;;
        --bins)             export BINS="$2"; shift 2 ;;
        --mode)             export MODE="$2"; shift 2 ;;
        --weight)           export WEIGHT="$2"; shift 2 ;;
        --collect-parallel) export COLLECT_PARALLEL="$2"; shift 2 ;;
        --collect-time)     export COLLECT_TIME="$2"; shift 2 ;;
        --assemble-time)    export ASSEMBLE_TIME="$2"; shift 2 ;;
        --plan-time)        export PLAN_TIME="$2"; shift 2 ;;
        --mem)              export MEM="$2"; shift 2 ;;
        *) echo "unknown option $1" >&2; exit 64 ;;
    esac
done

[[ -d "$IN" ]] || { echo "FATAL: INPUT_DIR $IN is not a directory" >&2; exit 2; }
mkdir -p "$RUN"
RUN="$(cd "$RUN" && pwd)"; IN="$(cd "$IN" && pwd)"

PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LLB_ROOT="$(cd "$PIPE/../.." && pwd)"
if [[ -f "$LLB_ROOT/common/env.sh" ]]; then
    # shellcheck source=../../common/env.sh
    source "$LLB_ROOT/common/env.sh"
    llb_load_config
    export LLB_ACCOUNT LLB_PARTITION
    export COLLECT_PARALLEL="${COLLECT_PARALLEL:-${LLB_ARRAY_PARALLEL:-}}"
    [[ -n "$COLLECT_PARALLEL" ]] || unset COLLECT_PARALLEL
fi

rec="$RUN/chunks.env"
if [[ -f "$rec" ]]; then
    prev="$(sed -n 's/^INPUT_DIR=//p' "$rec")"
    if [[ "$prev" != "$IN" ]]; then
        echo "FATAL: $RUN already chunked $prev; refusing to mix in $IN. Use a new RUNDIR." >&2
        exit 2
    fi
    echo "== re-chunking $IN into $RUN (re-plans from scratch)"
fi
{
    echo "# written by submit_chunks.sh $(date '+%Y-%m-%d %H:%M:%S')"
    echo "INPUT_DIR=$IN"
    echo "LLB_ROOT=$LLB_ROOT"
    echo "commit=$(git -C "$LLB_ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
    for v in TARGET_PER_BIN SHARDS BINS MODE WEIGHT COLLECT_PARALLEL; do
        echo "$v=${!v:-default}"
    done
} > "$rec"

export RUNDIR="$RUN" INPUT_DIR="$IN"
exec bash "$PIPE/submit.slurm"
