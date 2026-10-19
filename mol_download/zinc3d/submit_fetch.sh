#!/bin/bash
# submit_fetch.sh LIST [JOBS] [-- extra sbatch args]
#
# Fetch every URL in LIST into the CURRENT directory (the landing dir; files go
# to ./zinc22/... as wget -x -nH mirrors them) as ONE SLURM job: run.sh at JOBS
# connections (default 4), then up to three retry passes at 2, until status.sh
# reports nothing left to retry. Never an array: files.docking.org is one server.
#
#   cd /scratch/.../zinc3d_h18_30
#   llb zinc3d-fetch plan/idlist_urls.txt            # then sidecar_urls.txt, db2_urls.txt
#   llb zinc3d-fetch plan/db2_urls.txt 4 -- --time=7-00:00:00
#
# Resumable: fetch.sh skips every file already on disk and valid, so resubmitting
# the same LIST after a timeout costs local CPU and no requests. The job records
# LIST's sha256 in logs/fetch.<jobid>.manifest and refuses to start if LIST
# changed since the last submission from this directory (a new plan means a new
# list file name).
#
# Compute nodes need outbound HTTPS: check once with
#   srun -A maom2 -p standard -t 5 curl -sI https://files.docking.org/zinc22/ | head -1

set -euo pipefail

[[ $# -ge 1 ]] || { sed -n '2,/^$/p' "$0" >&2; exit 64; }
LIST="$1"; shift
JOBS=4
if [[ $# -gt 0 && "$1" != "--" ]]; then JOBS="$1"; shift; fi
[[ "${1:-}" == "--" ]] && shift

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LLB_ROOT="${LLB_ROOT:-$(cd "$HERE/../.." && pwd)}"
export LLB_ROOT
# shellcheck source=common/env.sh
source "$LLB_ROOT/common/env.sh"
llb_load_config

[[ -s "$LIST" ]] || { echo "FATAL: $LIST is empty or missing" >&2; exit 2; }
LIST="$(cd "$(dirname "$LIST")" && pwd)/$(basename "$LIST")"
sha() { sha256sum "$1" 2>/dev/null | awk '{print $1}' || shasum -a 256 "$1" | awk '{print $1}'; }
LIST_SHA="$(sha "$LIST")"

mkdir -p logs
stamp="logs/fetch.$(basename "$LIST").sha256"
if [[ -f "$stamp" && "$(cat "$stamp")" != "$LIST_SHA" ]]; then
    echo "FATAL: $LIST changed since it was last submitted here (sha256 in $stamp)." >&2
    echo "       Give a changed list a new file name, or delete $stamp if you mean it." >&2
    exit 2
fi
echo "$LIST_SHA" > "$stamp"

for t in "${WGET:-wget}" "${PARALLEL:-parallel}"; do
    command -v "$t" >/dev/null || { echo "FATAL: $t not found (set WGET/PARALLEL in $LLB_CONFIG)" >&2; exit 2; }
done
chmod +x "$LLB_ROOT"/mol_download/*.sh 2>/dev/null || true   # README Gotcha 1

n=$(grep -c . "$LIST")
echo "== $n URLs from $LIST -> $PWD/zinc22 at -j $JOBS; account=$LLB_ACCOUNT partition=$LLB_PARTITION"
sbatch -A "$LLB_ACCOUNT" -p "$LLB_PARTITION" \
    --job-name=zinc3d-fetch --cpus-per-task="$JOBS" --mem=4G --time=3-00:00:00 \
    --output="$PWD/logs/fetch_%j.out" \
    --export=ALL,LLB_ROOT="$LLB_ROOT",LIST="$LIST",LIST_SHA="$LIST_SHA",JOBS="$JOBS",LANDING="$PWD" \
    "$@" "$HERE/fetch.sbatch"
