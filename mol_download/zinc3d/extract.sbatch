#!/bin/bash
# One array task of extract.py. Submit with submit_extract.sh (llb zinc3d-extract),
# which creates logs/, pins the plan's sha256 and passes PLAN, OUT, NSHARDS,
# LANDING, LLB_ROOT and EXTRA.

set -uo pipefail

cd "${LANDING:?LANDING not set}" || exit 78
# shellcheck disable=SC1091
source "${LLB_ROOT:?LLB_ROOT not set}/common/env.sh" || exit 78
llb_load_config || exit 78
llb_activate PREP || exit 78

now="$(cd "$PLAN" && sha256sum index.tsv db2_urls.txt)"
[[ "$now" == "$(cat "$OUT/plan.sha256")" ]] \
    || { echo "FATAL: $PLAN changed since submission"; exit 78; }

# shellcheck disable=SC2086
python "$LLB_ROOT/mol_download/zinc3d/extract.py" --plan "$PLAN" --root "$LANDING" \
    --out "$OUT" --shard "${SLURM_ARRAY_TASK_ID:?}" --nshards "$NSHARDS" ${EXTRA:-}
