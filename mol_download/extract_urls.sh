#!/bin/bash
# extract_urls.sh -- the file URLs out of a ZINC download script, for run.sh.
#
# Usage:  ./extract_urls.sh zinc22-2D-download.curl > urls.txt
#         ./extract_urls.sh a.curl b.wget > urls.txt      # several scripts: union
#         cat script.sh | ./extract_urls.sh > urls.txt
#
# CartBlanche22's tranche browser hands out a curl, wget or PowerShell script
# with one command per file, e.g.
#
#   curl --user 'gpcr:xtal' --retry 3 ... -o H17/H17M300.smi.gz https://files.docking.org/zinc22/2d/H17/H17M300.smi.gz
#   wget -nH -r -l7 -np -A '*-L-*db2.tgz' https://files.docking.org/zinc22/zinc-22a/H11/...
#
# Don't run it as-is: the wget flavour crawls recursively (README, step 1), and
# the curl flavour writes H17/H17M300.smi.gz rather than the zinc22/2d/H17/...
# tree that status.sh and sample_2d.py --local expect. Fetch with run.sh instead.
#
# Why not `grep -o 'https://[^ ]*'`: it keeps a closing quote when the URL is
# quoted, and a trailing \r when the script has Windows line endings. A \r makes
# every URL differ from the path on disk, so status.sh reports every file
# missing, on every pass, forever.
#
# URLs go to stdout, sorted and unique. A summary goes to stderr. Exits 1 if no
# URL was found, or if any URL is not a file on files.docking.org (a directory,
# an index page, a different host) -- those are listed rather than passed on.

set -euo pipefail

if [[ $# -eq 0 && -t 0 ]]; then sed -n '2,24p' "$0" >&2; exit 64; fi
for f in "$@"; do [[ -r "$f" ]] || { echo "cannot read: $f" >&2; exit 2; }; done

tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT

# Strip CRs first, then take every http(s) token up to whitespace or a quote,
# backslash, semicolon or closing bracket.
cat "$@" | tr -d '\r' \
    | grep -oE "https?://[^][:space:]\"'\\\\;<>()]+" \
    | sed 's|^http://|https://|' \
    | sort -u > "$tmp/all" || true

grep -E '^https://files\.docking\.org/.+\.(smi\.gz|db2\.tgz|db2\.gz|sdf\.gz|mol2\.gz)$' \
    "$tmp/all" > "$tmp/good" || true
comm -23 "$tmp/all" "$tmp/good" > "$tmp/bad"

cat "$tmp/good"

n=$(wc -l < "$tmp/good" | tr -d ' ')
{
    echo "extract_urls: $n file URL(s)"
    if (( n > 0 )); then
        sed -E 's|.*/||; s|^(H[0-9]{2}).*\.([a-z0-9]+\.[a-z]+)$|\1 \2|' "$tmp/good" \
            | sort | uniq -c | awk '{printf "  %-4s %-8s %6d\n", $2, $3, $1}'
    fi
} >&2

if [[ -s "$tmp/bad" ]]; then
    echo "extract_urls: $(wc -l < "$tmp/bad" | tr -d ' ') URL(s) are not files on files.docking.org; NOT written:" >&2
    sed 's/^/  /' "$tmp/bad" | head -10 >&2
    exit 1
fi
(( n > 0 )) || { echo "extract_urls: no file URLs found" >&2; exit 1; }
