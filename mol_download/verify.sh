#!/bin/bash
# verify.sh -- integrity-check every downloaded archive; optionally purge bad ones.
#
# Usage:  ./verify.sh          report only
#         ./verify.sh --purge  delete corrupt archives so they get re-fetched
#
# A file can transfer with wget exit 0 and still be truncated. gzip -t is the
# only check that actually catches it. Run this before handing anything to the
# Python pipeline -- finding corruption here costs seconds, finding it three
# hours into a DOCK run does not.

set -u

ROOT="${ZINC_ROOT:-zinc22}"
PURGE=0
[[ "${1:-}" == "--purge" ]] && PURGE=1

total=0
bad=0
empty=0
: > corrupt.txt

while IFS= read -r f; do
    total=$(( total + 1 ))
    if [[ ! -s "$f" ]]; then
        empty=$(( empty + 1 ))
        echo "$f" >> corrupt.txt
        (( PURGE )) && rm -f "$f"
        continue
    fi
    if ! gzip -t "$f" 2>/dev/null; then
        bad=$(( bad + 1 ))
        echo "$f" >> corrupt.txt
        (( PURGE )) && rm -f "$f"
    fi
done < <(find "$ROOT" \( -name '*.tgz' -o -name '*.gz' \) -type f 2>/dev/null)

printf '%-24s %8d\n' \
    "archives checked" "$total" \
    "zero-byte"        "$empty" \
    "failed gzip -t"   "$bad"

if (( empty + bad > 0 )); then
    echo
    if (( PURGE )); then
        echo "Purged $(( empty + bad )) file(s); listed in corrupt.txt."
        echo "Re-fetch them: llb fetch-status <url list>, then llb fetch retry.txt 2."
    else
        echo "Listed in corrupt.txt. Re-run with --purge to delete them."
    fi
fi
