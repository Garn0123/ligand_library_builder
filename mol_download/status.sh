#!/bin/bash
# status.sh -- what do we actually have, and what still needs fetching?
#
# Usage:  ./status.sh [urls.txt]
#
# The filesystem is the source of truth, NOT the wget log. TLS-level failures
# ("GnuTLS: Error in the pull function") never name a URL in the log, so a
# log-derived retry list silently drops them. Diffing expected-vs-present
# catches every failure mode at once.
#
# Writes:
#   missing.txt  every expected file not currently on disk (full URLs)
#   retry.txt    missing.txt minus known-permanent 404s -- feed this to run.sh

set -u

URLS="${1:-urls.txt}"
BASE="${ZINC_BASE:-https://files.docking.org/}"
ROOT="${ZINC_ROOT:-zinc22}"
PERMANENT_TSV="${PERMANENT_TSV:-permanent.tsv}"

[[ -r "$URLS" ]] || { echo "no url list: $URLS" >&2; exit 1; }

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

sed "s|^${BASE}||" "$URLS" | sort -u > "$tmp/expected"
# .tgz (3D db2 archives) or .smi.gz (2D SMILES). Not -name '*.{tgz,gz}':
# find does no brace expansion, so that pattern matches nothing at all.
find "$ROOT" \( -name '*.tgz' -o -name '*.gz' \) -type f 2>/dev/null | sed 's|^\./||' | sort -u > "$tmp/have"

comm -23 "$tmp/expected" "$tmp/have" | sed "s|^|${BASE}|" | sort -u > missing.txt

# Exclude tranches already confirmed absent upstream.
if [[ -s "$PERMANENT_TSV" ]]; then
    cut -f2 "$PERMANENT_TSV" | sort -u > "$tmp/permanent"
else
    : > "$tmp/permanent"
fi
comm -23 missing.txt "$tmp/permanent" > retry.txt

expected=$(wc -l < "$tmp/expected")
have=$(wc -l < "$tmp/have")
perm=$(wc -l < "$tmp/permanent")
miss=$(wc -l < missing.txt)
retry=$(wc -l < retry.txt)

printf '%-28s %8d\n' \
    "expected (urls.txt)"   "$expected" \
    "present on disk"       "$have" \
    "missing"               "$miss" \
    "  of which permanent"  "$perm" \
    "  -> retry.txt"        "$retry"

if (( retry == 0 && miss > 0 )); then
    echo
    echo "All remaining gaps are confirmed-permanent. Run ./tranche_pattern.sh"
    echo "to sanity-check that they cluster the way real absences should."
fi
