#!/bin/bash
# submit_qupkake.sh PARENTS_DIR OUT_DIR -- preflight on this node, then submit
# the QupKake array with the account, partition and environment from
# config/hpc.env.
#
#   run_qupkake/submit_qupkake.sh parents protomers
#   run_qupkake/submit_qupkake.sh parents protomers --time=04:00:00   # extra sbatch args
#   run_qupkake/submit_qupkake.sh parents protomers --only-missing    # retry: only shards
#                                                                     # with no DONE
#
# SLURM logs go to ./logs of the directory you run this from (or $LLB_LOGDIR),
# so the repo stays read-only.
#
# Each task gets 1 CPU and 4G (LLB_QUPKAKE_MEM to change); concurrency per array
# is LLB_ARRAY_PARALLEL from config/hpc.env -- with 1-core tasks, raise it.
#
# The preflight (`llb doctor --stage qupkake`) checks the QupKake env, its
# imports and xtb HERE, so a wrong path fails once on the login node instead of
# once per array task. Resubmitting resumes: finished shards are skipped.

set -euo pipefail

[[ $# -ge 2 ]] || { sed -n '2,12p' "$0" >&2; exit 64; }
PARENTS_DIR="$(cd "$1" && pwd)"; OUT_DIR="$(mkdir -p "$2" && cd "$2" && pwd)"; shift 2
ONLY_MISSING=0
SBATCH_ARGS=()
for a in "$@"; do
    if [[ "$a" == --only-missing ]]; then ONLY_MISSING=1; else SBATCH_ARGS+=("$a"); fi
done

LLB_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export LLB_ROOT
# shellcheck source=common/env.sh
source "$LLB_ROOT/common/env.sh"
llb_load_config

META="$PARENTS_DIR/shards.tsv.meta"
[[ -f "$META" ]] || { echo "FATAL: $META missing -- run llb parents first (or use a triage/routed directory)" >&2; exit 2; }
N="$(awk -F= '$1=="n_shards"{print $2}' "$META")"
[[ "$N" =~ ^[0-9]+$ && "$N" -gt 0 ]] || {
    echo "FATAL: no n_shards in $META; keys present:" >&2; cut -d= -f1 "$META" | sed 's/^/  /' >&2; exit 2; }

echo "== preflight (config: $LLB_CONFIG)"
"$LLB_ROOT/bin/llb" doctor --stage qupkake \
    || { echo "FATAL: llb doctor failed (see above); nothing submitted" >&2; exit 2; }

# Logs go to the run directory you submit from, never into the repo.
LOGDIR="${LLB_LOGDIR:-$PWD/logs}"
mkdir -p "$LOGDIR"
# One core per task (see qupkake_array.sbatch), so a library is many small
# tasks -- possibly more than SLURM's MaxArraySize, which bounds array INDICES.
# Larger libraries go out as several arrays, each told its shard offset.
max_array="$(scontrol show config 2>/dev/null | awk -F= '/^MaxArraySize/{gsub(/ /,"",$2); print $2}')"
[[ "$max_array" =~ ^[0-9]+$ ]] || max_array=1001          # SLURM's default
# Task indices within one array (offset ... offset+n-1) as SLURM ranges: all of
# them, or with --only-missing only shards with no DONE ("3,17-19,1233").
array_spec() {   # OFFSET N
    if (( ! ONLY_MISSING )); then echo "0-$(( $2 - 1 ))"; return; fi
    local i
    for (( i = 0; i < $2; i++ )); do
        [[ -f "$(printf '%s/shard_%05d/DONE' "$OUT_DIR" $(( $1 + i )))" ]] || echo "$i"
    done | awk 'NR==1{s=p=$1; next} $1==p+1{p=$1; next} {o=o (s==p?s:s"-"p) ","; s=p=$1}
                END{if (NR) print o (s==p?s:s"-"p)}'
}
LLB_SHARD_OFFSET=0
n_arrays=$(( (N + max_array - 1) / max_array ))
if (( ONLY_MISSING )); then
    echo "== submitting only shards with no DONE in $OUT_DIR, 1 CPU each,"
else
    echo "== submitting $N shard(s) as $n_arrays array(s) (MaxArraySize $max_array), 1 CPU each,"
fi
echo "   account=$LLB_ACCOUNT partition=$LLB_PARTITION, <= ${LLB_ARRAY_PARALLEL:-50} running per array"
submitted=0
while (( LLB_SHARD_OFFSET < N )); do
    n=$(( N - LLB_SHARD_OFFSET )); (( n > max_array )) && n=$max_array
    spec="$(array_spec "$LLB_SHARD_OFFSET" "$n")"
    if [[ -z "$spec" ]]; then LLB_SHARD_OFFSET=$(( LLB_SHARD_OFFSET + n )); continue; fi
    (( ONLY_MISSING )) && echo "   shards $LLB_SHARD_OFFSET + [$spec]"
    submitted=1
    sbatch -A "$LLB_ACCOUNT" -p "$LLB_PARTITION" \
        --cpus-per-task=1 --mem="${LLB_QUPKAKE_MEM:-4G}" \
        --array="${spec}%${LLB_ARRAY_PARALLEL:-50}" \
        --export=ALL,LLB_ROOT="$LLB_ROOT",LLB_CONFIG="$LLB_CONFIG",PARENTS_DIR="$PARENTS_DIR",OUT_DIR="$OUT_DIR",LLB_SHARD_OFFSET="$LLB_SHARD_OFFSET" \
        --output="$LOGDIR/qupkake_%A_%a.out" \
        "${SBATCH_ARGS[@]}" "$LLB_ROOT/run_qupkake/slurm/qupkake_array.sbatch"
    LLB_SHARD_OFFSET=$(( LLB_SHARD_OFFSET + n ))
done
(( submitted )) || echo "== nothing to submit: every shard has a DONE"
