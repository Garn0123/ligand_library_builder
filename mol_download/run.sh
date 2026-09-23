#!/bin/bash
# run.sh -- fetch a URL list through GNU parallel at polite concurrency.
#
# Usage:  ./run.sh urls.txt        first pass  (defaults to -j 4)
#         ./run.sh retry.txt 2     retry pass  (use -j 2, gentler)
#
# files.docking.org is one academic server at UCSF, not a CDN. Four concurrent
# connections is the ceiling for a first pass; drop to two on retries. Raising
# it does not make things faster -- bandwidth per connection is the bottleneck,
# and the TLS drops in wget.log are the server telling you it is at capacity.

set -u

# Tool paths (WGET, PARALLEL) and ZINC_DATA from config/hpc.env, if present.
# Optional here: on a DTN with wget and parallel on PATH, no config is needed.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # scripts; data lands in $PWD
SMP_ROOT="${SMP_ROOT:-$(dirname "$HERE")}"
cfg="${SMP_CONFIG:-$SMP_ROOT/config/hpc.env}"
if [[ -f "$cfg" ]]; then
    # shellcheck disable=SC1090
    source "$cfg"
fi
export WGET="${WGET:-wget}"
PARALLEL="${PARALLEL:-parallel}"

LIST="${1:-urls.txt}"
JOBS="${2:-4}"
DELAY="${DELAY:-0.3}"

[[ -r "$LIST" ]] || { echo "no such list: $LIST" >&2; exit 1; }

command -v "$PARALLEL" >/dev/null || { echo "GNU parallel not found ($PARALLEL); set PARALLEL in config/hpc.env" >&2; exit 1; }
[[ -x "$HERE/fetch.sh" ]] || { echo "$HERE/fetch.sh is not executable -- run: chmod +x $HERE/*.sh" >&2; exit 1; }

stamp="$(date +%Y%m%d-%H%M%S)"
joblog="joblog-${stamp}.tsv"

# Fresh failure ledgers per pass. permanent.tsv accumulates across passes on
# purpose -- once a tranche 404s it stays 404.
: > failed.tsv

# wget -x -nH mirrors the URL path under the working directory, so run from
# $ZINC_DATA when it is configured.
if [[ -n "${ZINC_DATA:-}" && "$PWD" != "$ZINC_DATA" ]]; then
    echo "note: ZINC_DATA=$ZINC_DATA but running in $PWD; files land under $PWD" >&2
fi
echo "list=$LIST  jobs=$JOBS  joblog=$joblog"
"$PARALLEL" -j "$JOBS" --delay "$DELAY" --joblog "$joblog" -a "$LIST" "$HERE/fetch.sh" {}

echo
echo "--- pass complete ---"
"$HERE/status.sh" "${URLS:-urls.txt}"   # the FULL list, not retry.txt
