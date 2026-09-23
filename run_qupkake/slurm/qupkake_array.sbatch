#!/bin/bash
#SBATCH --job-name=qupkake
#SBATCH --account=maom99
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/qupkake_%A_%a.out
#
# One QupKake + protomer-enumeration task per parent shard.
#
# Submit with run_qupkake/submit_qupkake.sh, which preflights the QupKake env
# on the login node and passes account, partition, SMP_ROOT, SMP_CONFIG,
# PARENTS_DIR and OUT_DIR. Software locations come from config/hpc.env
# (QUPKAKE_SETUP / QUPKAKE_ENV / QUPKAKE_EXE / XTBPATH), not from this file.
#
# Resubmitting the same range resumes: finished shards are skipped (FORCE=1 to
# redo).
#
# Names written here are throwaway (--name-style index). Contract names are
# assigned once, over the union of every shard, by assign_names.py -- a shard
# cannot see the others, so it cannot promise global uniqueness.

set -euo pipefail

: "${PARENTS_DIR:?set PARENTS_DIR to prepare_parents.py --outdir}"
: "${OUT_DIR:?set OUT_DIR}"
PH="${PH:-6.4 7.4 8.4}"
MARGIN="${MARGIN:-1.0}"
MIN_POP="${MIN_POP:-0.01}"
MAX_STATES="${MAX_STATES:-8}"
FORCE="${FORCE:-0}"

# ---- repo + config ------------------------------------------------------------
# sbatch runs a spooled COPY of this script, so BASH_SOURCE cannot find the
# repo. submit_qupkake.sh exports SMP_ROOT; SLURM_SUBMIT_DIR is the fallback.
env_sh=""
for d in "${SMP_ROOT:-}" "${SLURM_SUBMIT_DIR:-}" "${SLURM_SUBMIT_DIR:-}/.."; do
    [[ -n "$d" && -f "$d/common/env.sh" ]] && { env_sh="$d/common/env.sh"; break; }
done
if [[ -z "$env_sh" ]]; then
    echo "FATAL: common/env.sh not found under SMP_ROOT='${SMP_ROOT:-}' or" \
         "SLURM_SUBMIT_DIR='${SLURM_SUBMIT_DIR:-}'. Submit with submit_qupkake.sh." >&2
    exit 2
fi
# shellcheck disable=SC1090
source "$env_sh"
smp_load_config
SCRIPT="$SMP_ROOT/run_qupkake/qupkake_protomers.py"
[[ -f "$SCRIPT" ]] || { echo "FATAL: $SCRIPT missing" >&2; exit 2; }

# ---- manifest integrity ------------------------------------------------------
MANIFEST="$PARENTS_DIR/shards.tsv"
META="$MANIFEST.meta"
[[ -f "$MANIFEST" && -f "$META" ]] || { echo "FATAL: $MANIFEST or $META missing" >&2; exit 2; }
want=""
while IFS='=' read -r k v; do
    [[ "$k" == *checksum ]] && want="$v"
done < "$META"
if [[ -z "$want" ]]; then
    echo "FATAL: no *checksum key in $META; keys present:" >&2
    cut -d= -f1 "$META" | sed 's/^/  /' >&2
    exit 2
fi
have="$(sha256sum "$MANIFEST" | awk '{print $1}')"
if [[ "$have" != "$want" ]]; then
    echo "FATAL: $MANIFEST changed since it was written (sha256 $have, meta says $want)." >&2
    echo "Array index $SLURM_ARRAY_TASK_ID no longer means what it meant at submit." >&2
    exit 2
fi

IDX="${SLURM_ARRAY_TASK_ID:?run as an array task, or set SLURM_ARRAY_TASK_ID}"
row="$(awk -F'\t' -v i="$IDX" 'NR>1 && $1==i' "$MANIFEST")"
[[ -n "$row" ]] || { echo "FATAL: index $IDX not in $MANIFEST" >&2; exit 2; }
shard="$PARENTS_DIR/$(cut -f2 <<<"$row")"
shard_sha="$(cut -f4 <<<"$row")"
[[ "$(sha256sum "$shard" | awk '{print $1}')" == "$shard_sha" ]] \
    || { echo "FATAL: $shard does not match its manifest sha256" >&2; exit 2; }

out="$OUT_DIR/$(printf 'shard_%05d' "$IDX")"
if [[ -f "$out/DONE" && "$FORCE" != "1" ]]; then
    if grep -q "^shard_sha256=$shard_sha$" "$out/DONE"; then
        echo "shard $IDX already done ($out/DONE); FORCE=1 to redo"
        exit 0
    fi
    echo "FATAL: $out/DONE is from a different shard input; refusing to mix." >&2
    exit 2
fi

# ---- environment ---------------------------------------------------------------
smp_activate QUPKAKE
QK_EXE="$(smp_require_exe QUPKAKE_EXE qupkake)"
python -c 'import qupkake, rdkit' \
    || { echo "FATAL: env ${QUPKAKE_ENV} lacks qupkake or rdkit" >&2; exit 2; }
export OMP_NUM_THREADS=1       # -mp supplies the parallelism; xtb must not also thread
export MKL_NUM_THREADS=1

echo "shard $IDX  $(wc -l < "$shard") parents  host $(hostname)  $(date -Is)"
rm -rf "$out"
mkdir -p "$out"
# shellcheck disable=SC2086  # PH is a deliberate word list
python "$SCRIPT" "$shard" -o "$out" \
    --ph $PH --margin "$MARGIN" --min-population "$MIN_POP" --max-states "$MAX_STATES" \
    --stereo carry --name-style index -mp "${SLURM_CPUS_PER_TASK:-1}" \
    --qupkake-exe "$QK_EXE"

[[ -s "$out/protomers.csv" ]] || { echo "FATAL: no protomers.csv written" >&2; exit 1; }
{
    echo "shard_index=$IDX"
    echo "shard_sha256=$shard_sha"
    echo "ph=$PH"; echo "margin=$MARGIN"; echo "min_population=$MIN_POP"
    echo "max_states=$MAX_STATES"
    echo "qupkake_exe=$QK_EXE"; echo "xtbpath=${XTBPATH:-bundled}"
    echo "config=$SMP_CONFIG"
    echo "finished=$(date -Is)"
} > "$out/DONE"
echo "done: $out"
