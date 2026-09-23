#!/bin/bash
# tranche_pattern.sh -- are the 404s real absences, or a bug in the URL list?
#
# Usage:  ./tranche_pattern.sh [permanent.tsv]
#
# ZINC22 filenames encode chemistry:
#
#     H05M100-O-aaaaaa.db2.tgz
#     ^^^                        H<heavy atom count>
#        ^                       M = negative cLogP, P = positive cLogP
#         ^^^                    |cLogP| x 100
#             ^                  net charge: L=-2 M=-1 N=0 O=+1 P=+2
#               ^^^^^^           subset + chunk id
#
#     H17M100.smi.gz                  (2D: no charge letter -- 2D tranches
#                                      hold every charge in one file)
#
# The tranche browser emits the full combinatorial grid whether or not a given
# cell was ever populated, so 404s are expected. What matters is WHERE they
# fall. Absences concentrated in M tranches (hydrophilic) at low heavy-atom
# count are chemically real: charged, very polar, very small molecules are a
# genuinely sparse corner of make-on-demand space.
#
# A flat spread across all logP signs and charges, or every file under one
# directory, means the URL list is wrong -- not the data.

set -u
SRC="${1:-permanent.tsv}"
[[ -r "$SRC" ]] || { echo "no such file: $SRC" >&2; exit 1; }

urls="$(cut -f2 "$SRC" 2>/dev/null | grep -c . || true)"
echo "permanent failures: $urls"
echo

parse() {
    cut -f2 "$SRC" | sed 's|.*/||' \
        | sed -nE 's|^H([0-9]+)([MP])([0-9]+)(-([A-Z])-)?.*|\1\t\2\t\5|p' \
        | awk -F'\t' -v OFS='\t' '$3==""{$3="none (2D file)"} 1'
}

echo "== by logP sign =="
parse | cut -f2 | sort | uniq -c | sort -rn
echo
echo "== by net charge =="
parse | cut -f3 | sort | uniq -c | sort -rn
echo
echo "== by heavy atom count =="
parse | cut -f1 | sort | uniq -c | sort -k2,2
echo
echo "== by directory =="
cut -f2 "$SRC" | sed 's|/[^/]*$||' | sort | uniq -c | sort -rn | head -15
echo
echo "Expected shape: heavily skewed toward M, concentrated at low HAC,"
echo "spread across several directories. Anything else -- especially all"
echo "failures sharing one directory or one charge letter -- means check"
echo "the URL list before discarding these."
