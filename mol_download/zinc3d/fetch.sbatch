#!/bin/bash
# One job: fetch LIST into LANDING with run.sh, then retry passes until status.sh
# has nothing left to retry. Submit with submit_fetch.sh (llb zinc3d-fetch), which
# creates logs/ and passes LIST, LIST_SHA, JOBS, LANDING and LLB_ROOT.

set -uo pipefail

cd "${LANDING:?LANDING not set}" || exit 78
# shellcheck disable=SC1091
source "${LLB_ROOT:?LLB_ROOT not set}/common/env.sh" || exit 78
llb_load_config || exit 78
export WGET="${WGET:-wget}"
D="$LLB_ROOT/mol_download"

now_sha="$(sha256sum "$LIST" | awk '{print $1}')"
if [[ "$now_sha" != "${LIST_SHA:?}" ]]; then
    echo "FATAL: $LIST changed between submission and start ($LIST_SHA -> $now_sha)"; exit 78
fi
{
    echo "job      ${SLURM_JOB_ID:-local} on $(hostname) $(date '+%F %T')"
    echo "list     $LIST"
    echo "sha256   $LIST_SHA"
    echo "urls     $(grep -c . "$LIST")"
    echo "wget     $("$WGET" --version | head -1)"
    echo "llb      $(git -C "$LLB_ROOT" rev-parse --short HEAD 2>/dev/null || echo '?')"
} | tee "logs/fetch.${SLURM_JOB_ID:-local}.manifest"

curl -sfI --max-time 30 https://files.docking.org/zinc22/ >/dev/null \
    || { echo "FATAL: no HTTPS to files.docking.org from $(hostname)"; exit 78; }

export URLS="$LIST"                 # run.sh reports status against the FULL list
"$D/run.sh" "$LIST" "${JOBS:-4}"
for pass in 2 3 4; do
    "$D/status.sh" "$LIST" > /dev/null
    left=$(grep -c . retry.txt || true)
    (( left == 0 )) && break
    echo "== pass $pass: $left URL(s) to retry, after a pause"
    sleep 300
    cp retry.txt "retry.pass$pass.txt"
    "$D/run.sh" "retry.pass$pass.txt" 2
done

"$D/status.sh" "$LIST"
left=$(grep -c . retry.txt || true)
echo "== done $(date '+%F %T'): $left still to retry (retry.txt), $(cut -f2 permanent.tsv 2>/dev/null | sort -u | grep -c . || true) permanent 404s"
(( left == 0 )) || exit 1
