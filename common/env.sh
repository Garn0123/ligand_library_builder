# common/env.sh -- sourced, never executed. Finds the repo, reads
# config/hpc.env, and activates a tool's environment the way this cluster needs.
#
#   source "$LLB_ROOT/common/env.sh"
#   llb_load_config
#   llb_activate QUPKAKE        # or DB2C, PREP
#
# Written for scripts running under `set -euo pipefail`: conda's activation
# scripts dereference unset variables, so nounset is lifted around them only.

# The repo root. A SLURM job runs a spooled COPY of its script, so BASH_SOURCE
# is useless there; the submit wrappers export LLB_ROOT, and SLURM_SUBMIT_DIR
# is searched upward as a fallback.
llb_find_root() {
    local tried=() d
    for d in "${LLB_ROOT:-}" "${SLURM_SUBMIT_DIR:-}" "$PWD" \
             "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)"; do
        [[ -z "$d" ]] && continue
        while [[ "$d" != "/" && -n "$d" ]]; do
            tried+=("$d/common/env.sh")
            if [[ -f "$d/common/env.sh" && -f "$d/config/hpc.env.example" ]]; then
                LLB_ROOT="$d"; export LLB_ROOT; return 0
            fi
            d="$(dirname "$d")"
        done
    done
    echo "FATAL: ligand_library_builder root not found. Looked for:" >&2
    printf '  %s\n' "${tried[@]}" >&2
    echo "Set LLB_ROOT=/path/to/ligand_library_builder." >&2
    return 2
}

llb_load_config() {
    llb_find_root || return 2
    local cfg="${LLB_CONFIG:-$LLB_ROOT/config/hpc.env}"
    if [[ ! -f "$cfg" ]]; then
        echo "FATAL: no config at $cfg" >&2
        echo "  cp $LLB_ROOT/config/hpc.env.example $LLB_ROOT/config/hpc.env  and edit it" >&2
        return 2
    fi
    set +u
    # shellcheck disable=SC1090
    source "$cfg"
    set -u
    LLB_CONFIG="$cfg"
    # Resolve conda now, before any <TOOL>_SETUP can `module purge` it away.
    if [[ -z "${CONDA_BASE:-}" ]]; then
        CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    fi
    export LLB_CONFIG CONDA_BASE
}

# llb_activate TOOL -- run ${TOOL}_SETUP, then activate ${TOOL}_ENV.
# PREP falls back to QUPKAKE when PREP_ENV is empty.
llb_activate() {
    local tool="$1" setup_var="${1}_SETUP" env_var="${1}_ENV"
    local setup="${!setup_var:-}" env="${!env_var:-}"
    if [[ "$tool" == PREP && -z "$env" ]]; then
        setup="${QUPKAKE_SETUP:-}"; env="${QUPKAKE_ENV:-}"
    fi
    if [[ -z "$env" ]]; then
        echo "FATAL: ${env_var} is empty in $LLB_CONFIG" >&2; return 2
    fi
    if [[ -z "${CONDA_BASE:-}" || ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
        echo "FATAL: conda not found (CONDA_BASE='${CONDA_BASE:-}'). Set CONDA_BASE in $LLB_CONFIG." >&2
        return 2
    fi
    local rc=0
    set +u
    if [[ -n "$setup" ]]; then
        eval "$setup" || rc=$?
    fi
    if (( rc == 0 )); then
        # shellcheck disable=SC1091
        source "$CONDA_BASE/etc/profile.d/conda.sh" && conda activate "$env" || rc=$?
    fi
    set -u
    if (( rc != 0 )); then
        echo "FATAL: activating $tool failed (setup: '${setup}', env: '$env', rc=$rc)" >&2
        return 2
    fi
    [[ -n "${XTBPATH:-}" ]] && export XTBPATH
    echo "[$tool] env=$env python=$(command -v python)"
}

# llb_require_exe VAR DEFAULT -- resolve an executable: explicit path from the
# config, else DEFAULT on PATH. Prints the path; fails with what was tried.
llb_require_exe() {
    local var="$1" default="$2" want="${!1:-}"
    local exe="${want:-$default}"
    if [[ "$exe" == */* ]]; then
        [[ -x "$exe" ]] && { echo "$exe"; return 0; }
    elif command -v "$exe" >/dev/null 2>&1; then
        command -v "$exe"; return 0
    fi
    if [[ -n "$want" ]]; then
        echo "FATAL: ${var}='$want' (in $LLB_CONFIG) is not an executable" >&2
    else
        echo "FATAL: '$default' is not on PATH; set ${var} in $LLB_CONFIG" >&2
    fi
    return 2
}
