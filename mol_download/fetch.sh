#!/bin/bash
# fetch.sh -- download one ZINC22 tranche file (3D .db2.tgz or 2D .smi.gz, or a
# per-tranche id list .txt), idempotently and safely.
#
# Usage:  ./fetch.sh https://files.docking.org/zinc22/zinc-22a/H05/.../X.db2.tgz
#         ./fetch.sh https://files.docking.org/zinc22/2d/H17/H17M100.smi.gz
#
# Exit codes:
#   0  file present and valid (either already had it, or fetched it): gzip -t
#      for archives; for a .txt id list, non-empty, ends in a newline, not HTML
#   8  permanent HTTP failure (404/403/410) -- tranche does not exist, logged
#      to permanent.tsv, do not retry
#   *  transient failure after MAX_ATTEMPTS -- logged to failed.tsv, retry later
#
# Written for wget 1.19, which lacks --retry-on-http-error. The retry loop
# below supplies that behaviour: 5xx and TLS drops are retried with escalating
# backoff, 404s bail out immediately without burning attempts.

set -u

url="${1:-}"
[[ -z "$url" ]] && { echo "usage: fetch.sh URL" >&2; exit 64; }

BASE="${ZINC_BASE:-https://files.docking.org/}"
MAX_ATTEMPTS="${MAX_ATTEMPTS:-4}"
WGET_LOG="${WGET_LOG:-wget.log}"
FAILED_TSV="${FAILED_TSV:-failed.tsv}"
PERMANENT_TSV="${PERMANENT_TSV:-permanent.tsv}"

out="${url#"$BASE"}"

# gzip -t for archives. A plain-text id list has no checksum: require it to be
# non-empty, end in a newline (wget stops mid-line when a connection drops) and
# not be an HTML error page.
valid() {
    case "$1" in
        *.txt) [[ -s "$1" && -z "$(tail -c1 "$1")" ]] && ! head -c 64 "$1" | grep -qi '<' ;;
        *)     [[ -s "$1" ]] && gzip -t "$1" 2>/dev/null ;;
    esac
}

# --- already complete? -------------------------------------------------------
# gzip -t is the real test. An existence check is not enough: a dropped
# connection leaves a short file that looks fine to [[ -s ]].
if valid "$out"; then
    exit 0
fi

# Anything left over from a previous failure is garbage. Remove it. This is
# load-bearing -- see "The -c append trap" in README.md.
rm -f "$out"

# --- fetch with backoff ------------------------------------------------------
rc=1
code=""
for (( attempt=1; attempt<=MAX_ATTEMPTS; attempt++ )); do

    err="$(mktemp)"
    "${WGET:-wget}" --user=gpcr --password=xtal --auth-no-challenge -nv -x -nH \
         --tries=3 --waitretry=15 --retry-connrefused --timeout=60 \
         "$url" 2>"$err"
    rc=$?

    cat "$err" >> "$WGET_LOG"
    code="$(grep -o 'ERROR [0-9]\+' "$err" | tail -1 | awk '{print $2}')"
    rm -f "$err"

    # Success means exit 0 AND a valid archive on disk.
    if (( rc == 0 )) && valid "$out"; then
        exit 0
    fi

    # Permanent failures: stop immediately, do not sleep, do not retry.
    case "$code" in
        404|403|410)
            printf '%s\t%s\n' "$code" "$url" >> "$PERMANENT_TSV"
            rm -f "$out"
            exit 8
            ;;
    esac

    # Transient (5xx, TLS drop, timeout, truncation). Purge and back off.
    rm -f "$out"
    if (( attempt < MAX_ATTEMPTS )); then
        sleep $(( 20 * attempt ))
    fi
done

printf '%s\t%s\n' "${code:-$rc}" "$url" >> "$FAILED_TSV"
exit "$rc"
