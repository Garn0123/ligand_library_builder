#!/bin/bash
# submit_extract.sh PLAN OUT [--nshards N] [--delete-archive] [-- extra sbatch args]
#
# Run extract.py as a SLURM array from the CURRENT directory (the landing dir that
# holds ./zinc22/...): N tasks (default: one per 50 chunks, at most 500), 1 CPU
# each. Each task streams its share of PLAN/db2_urls.txt and writes the wanted
# records to OUT/db2/<bin>/<chunk>.db2.gz plus OUT/records/<chunk>.tsv.
#
#   llb zinc3d-extract plan zinc_db2
#   llb zinc3d-extract plan zinc_db2 --delete-archive      # free scratch as it goes
#
# Then, once every task is done:
#   llb zinc3d-collect --plan plan --out zinc_db2
#
# Resumable: a chunk with OUT/records/<chunk>.done is skipped, so resubmit the
# same command after failures or a timeout. OUT/plan.sha256 pins index.tsv and
# db2_urls.txt: a changed plan is refused rather than mixed into the same OUT.

set -euo pipefail

[[ $# -ge 2 ]] || { sed -n '2,/^$/p' "$0" >&2; exit 64; }
PLAN="$(cd "$1" && pwd)"; OUT="$2"; shift 2
NSHARDS=""; EXTRA=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --nshards)        NSHARDS="$2"; shift 2 ;;
        --delete-archive) EXTRA+=(--delete-archive); shift ;;
        --)               shift; break ;;
        *) echo "unknown option $1" >&2; exit 64 ;;
    esac
done

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LLB_ROOT="${LLB_ROOT:-$(cd "$HERE/../.." && pwd)}"
export LLB_ROOT
# shellcheck source=common/env.sh
source "$LLB_ROOT/common/env.sh"
llb_load_config

for f in index.tsv db2_urls.txt; do
    [[ -s "$PLAN/$f" ]] || { echo "FATAL: $PLAN/$f missing; run index.py chunks first" >&2; exit 2; }
done
mkdir -p "$OUT" logs
OUT="$(cd "$OUT" && pwd)"
sha() { sha256sum "$@" 2>/dev/null || shasum -a 256 "$@"; }
now="$(cd "$PLAN" && sha index.tsv db2_urls.txt)"
if [[ -f "$OUT/plan.sha256" ]]; then
    if [[ "$(cat "$OUT/plan.sha256")" != "$now" ]]; then
        echo "FATAL: $PLAN changed since $OUT was started (see $OUT/plan.sha256). Use a new OUT." >&2
        exit 2
    fi
    NSHARDS="${NSHARDS:-$(cat "$OUT/nshards")}"
    [[ "$NSHARDS" == "$(cat "$OUT/nshards")" ]] \
        || { echo "FATAL: $OUT was started with --nshards $(cat "$OUT/nshards")" >&2; exit 2; }
else
    n=$(grep -c . "$PLAN/db2_urls.txt")
    NSHARDS="${NSHARDS:-$(( (n + 49) / 50 ))}"
    (( NSHARDS > 500 )) && NSHARDS=500
    echo "$now" > "$OUT/plan.sha256"
    echo "$NSHARDS" > "$OUT/nshards"
fi

# preflight: the env the tasks will use can run the script
( llb_activate PREP || exit 2
  python -c 'import sys; assert sys.version_info >= (3, 8), sys.version' \
      || { echo "FATAL: PREP env python is older than 3.8" >&2; exit 2; }
) || exit 2

echo "== $NSHARDS task(s) over $(grep -c . "$PLAN/db2_urls.txt") chunks -> $OUT; account=$LLB_ACCOUNT"
sbatch -A "$LLB_ACCOUNT" -p "$LLB_PARTITION" \
    --job-name=zinc3d-extract --array="0-$((NSHARDS - 1))%${LLB_ARRAY_PARALLEL:-50}" \
    --cpus-per-task=1 --mem=2G --time=12:00:00 \
    --output="$PWD/logs/extract_%A_%a.out" \
    --export=ALL,LLB_ROOT="$LLB_ROOT",PLAN="$PLAN",OUT="$OUT",NSHARDS="$NSHARDS",LANDING="$PWD",EXTRA="${EXTRA[*]:-}" \
    "$@" "$HERE/extract.sbatch"
