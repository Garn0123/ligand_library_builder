#!/bin/bash
# install.sh -- create the MolSanitizer env at PREFIX and freeze it.
#
#   install.sh [PREFIX]              from environment.yml + requirements.txt
#   install.sh [PREFIX] --from-lock  from conda-linux-64.lock.txt + requirements.lock.txt
#                                    (byte-for-byte the same packages as the first install)
#
# PREFIX defaults to $LLB_MSANI_PREFIX, else
#   /nfs/turbo/umms-maom/johnbick/opt/molsanitizer/0.7.4
# Refuses an existing PREFIX: an env is replaced by installing a new one beside it.
#
# Writes PREFIX/LOCK/{conda-linux-64.lock.txt,requirements.lock.txt,versions.txt}.
# Copy the two lock files into molsanitizer/environment/ and commit them; that
# is what makes the next install reproducible.
#
# Then checks the compiled extensions load and Open Babel is on PATH. It does not
# build any molecule (that needs a compute node: msani.sh selftest).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${LLB_MSANI_PREFIX:-/nfs/turbo/umms-maom/johnbick/opt/molsanitizer/0.7.4}"
MODE=spec
for a in "$@"; do
    case "$a" in
        --from-lock) MODE=lock ;;
        -*) echo "install.sh: unknown option $a" >&2; exit 64 ;;
        *) PREFIX="$a" ;;
    esac
done
[[ -e "$PREFIX" ]] && { echo "install.sh: $PREFIX exists; remove it or pick another prefix" >&2; exit 1; }

CONDA_BASE="${LLB_MSANI_CONDA_BASE:-${CONDA_BASE:-$(conda info --base 2>/dev/null || true)}}"
[[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]] || {
    echo "install.sh: no conda at '$CONDA_BASE'; set LLB_MSANI_CONDA_BASE" >&2; exit 1; }
set +u; source "$CONDA_BASE/etc/profile.d/conda.sh"; set -u
mkdir -p "$(dirname "$PREFIX")"

if [[ $MODE == lock ]]; then
    for f in conda-linux-64.lock.txt requirements.lock.txt; do
        [[ -s "$HERE/$f" ]] || { echo "install.sh: $HERE/$f missing; install from spec first" >&2; exit 1; }
    done
    conda create -y -p "$PREFIX" --file "$HERE/conda-linux-64.lock.txt"
    "$PREFIX/bin/python" -m pip install --no-deps -r "$HERE/requirements.lock.txt"
else
    conda env create -p "$PREFIX" -f "$HERE/environment.yml"
fi

mkdir -p "$PREFIX/LOCK"
conda list -p "$PREFIX" --explicit --md5 > "$PREFIX/LOCK/conda-linux-64.lock.txt"
# pip-installed packages only (conda's are in the explicit list above)
"$PREFIX/bin/python" - > "$PREFIX/LOCK/requirements.lock.txt" <<'EOF'
import importlib.metadata as md
rows = []
for d in md.distributions():
    inst = (d.read_text("INSTALLER") or "").strip()
    if inst == "pip":
        rows.append(f"{d.metadata['Name']}=={d.version}")
print("\n".join(sorted(rows, key=str.lower)))
EOF

set +u; conda activate "$PREFIX"; set -u
{
    echo "installed $(date +%Y-%m-%dT%H:%M:%S%z) on $(hostname) by $(whoami), mode=$MODE"
    echo "glibc: $(ldd --version 2>&1 | head -1)"
    msani --version 2>&1 | sed -n '/version/p'
    obabel -V
    python - <<'EOF'
import amsolcpp, msani_confgen_cpp, msani_stereoisomers  # the compiled extensions
from msani.conformers.conformers import AMSOLCPP_AVAILABLE, OBABEL_AVAILABLE
assert AMSOLCPP_AVAILABLE, "AMSOLcpp did not import: no db2 output possible"
assert OBABEL_AVAILABLE, "obabel not found: the embedding fallback is off"
print("extensions: amsolcpp, msani_confgen_cpp, msani_stereoisomers load; obabel found")
EOF
} | tee "$PREFIX/LOCK/versions.txt"

echo
echo "OK: $PREFIX"
echo "lock files: $PREFIX/LOCK/conda-linux-64.lock.txt  $PREFIX/LOCK/requirements.lock.txt"
