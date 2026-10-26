# common/env.sh -- sourced, never executed. Finds the repo, reads
# hpc.env (+ users/$USER.env), and activates a tool's environment the way this cluster needs.
#
#   source "$LLB_ROOT/common/env.sh"
#   llb_load_config
#   llb_activate PREP           # or DOCK
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

# Config lookup:
#   1. LLB_CONFIG, if set (an override: nothing else is tried)
#   2. $LLB_ROOT/../hpc.env   the site layout: the checkout sits inside $LLB_SITE
#   3. $LLB_ROOT/config/hpc.env   the old place, still read (Q27)
# then <dir of hpc.env>/users/$USER.env if it exists: per-user keys (LLB_ACCOUNT,
# ZINC_DATA) over the site's. LLB_SITE is the site directory when the
# config is the site one (empty otherwise), set before sourcing so hpc.env can use it.
llb_load_config() {
    llb_find_root || return 2
    local cfg site_cfg tried=()
    site_cfg="$(dirname "$LLB_ROOT")/hpc.env"
    if [[ -n "${LLB_CONFIG:-}" ]]; then
        tried=("LLB_CONFIG=$LLB_CONFIG")
        cfg="$LLB_CONFIG"
    else
        tried=("$site_cfg" "$LLB_ROOT/config/hpc.env")
        cfg="$site_cfg"
        [[ -f "$cfg" ]] || cfg="$LLB_ROOT/config/hpc.env"
    fi
    if [[ ! -f "$cfg" ]]; then
        local looked; looked="$(printf '%s, ' "${tried[@]}")"
        echo "FATAL: no hpc.env; looked for: ${looked%, }" >&2
        echo "  cp $LLB_ROOT/config/hpc.env.example $site_cfg  (site layout)" >&2
        echo "  or cp it to $LLB_ROOT/config/hpc.env, and edit it" >&2
        return 2
    fi
    local cfg_dir; cfg_dir="$(cd "$(dirname "$cfg")" && pwd)"
    # By path, not by which branch found it: a job inherits LLB_CONFIG=$LLB_SITE/hpc.env
    # from the shell that loaded it, and must still see the same LLB_SITE.
    LLB_SITE=""
    if [[ "$(cd "$cfg_dir" && pwd -P)" == "$(cd "$LLB_ROOT/.." && pwd -P)" ]]; then
        LLB_SITE="$cfg_dir"
    fi
    LLB_USER_CONFIG="$cfg_dir/users/${USER:-$(id -un)}.env"
    export LLB_SITE
    set +u
    # shellcheck disable=SC1090
    source "$cfg"
    if [[ -f "$LLB_USER_CONFIG" ]]; then
        # shellcheck disable=SC1090
        source "$LLB_USER_CONFIG"
    else
        LLB_USER_CONFIG=""
    fi
    set -u
    LLB_CONFIG="$cfg"
    # Resolve conda now, before any <TOOL>_SETUP can `module purge` it away.
    if [[ -z "${CONDA_BASE:-}" ]]; then
        CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    fi
    export LLB_CONFIG LLB_USER_CONFIG CONDA_BASE
    return 0
}

# llb_activate TOOL -- run ${TOOL}_SETUP, then activate ${TOOL}_ENV.
# DOCK (db2tool, a compiled binary) needs no conda env: with DOCK_ENV empty only
# DOCK_SETUP runs.
llb_activate() {
    local tool="$1" setup_var="${1}_SETUP" env_var="${1}_ENV"
    local setup="${!setup_var:-}" env="${!env_var:-}"
    if [[ -z "$env" && "$tool" == DOCK ]]; then
        local rc=0
        set +u
        if [[ -n "$setup" ]]; then eval "$setup" || rc=$?; fi
        set -u
        (( rc == 0 )) || { echo "FATAL: DOCK_SETUP failed (rc=$rc): '$setup'" >&2; return 2; }
        echo "[DOCK] no conda env (DOCK_ENV empty); setup: '${setup:-none}'"
        return 0
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
    echo "[$tool] env=$env python=$(command -v python)"
}
