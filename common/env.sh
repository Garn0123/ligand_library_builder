# common/env.sh -- sourced, never executed. Finds the repo, reads
# hpc.env (+ users/$USER.env), and activates a tool's environment the way this cluster needs.
#
#   source "$LLB_ROOT/common/env.sh"
#   llb_load_config
#   llb_activate QUPKAKE        # or DB2C, PREP, DOCK, PKA
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

# Config lookup (docs/INSTALL_DESIGN.md, "One site directory"):
#   1. LLB_CONFIG, if set (an override: nothing else is tried)
#   2. $LLB_ROOT/../hpc.env   the site layout: the checkout sits inside $LLB_SITE
#   3. $LLB_ROOT/config/hpc.env   the old place, still read (Q27)
# then <dir of hpc.env>/users/$USER.env if it exists: per-user keys (ZBH_LICENSE,
# LLB_MAIL, LLB_ACCOUNT) over the site's. LLB_SITE is the site directory when the
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
    if [[ -z "${ZBH_LICENSE:-}" && -n "$LLB_SITE" ]]; then
        ZBH_LICENSE="$LLB_SITE/licenses/${USER:-$(id -un)}/ZBH"
    fi
    # Resolve conda now, before any <TOOL>_SETUP can `module purge` it away.
    if [[ -z "${CONDA_BASE:-}" ]]; then
        CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    fi
    export LLB_CONFIG LLB_USER_CONFIG CONDA_BASE
    [[ -n "${ZBH_DIR:-}" ]] && export ZBH_DIR
    [[ -n "${ZBH_LICENSE:-}" ]] && export ZBH_LICENSE
    return 0
}

# llb_activate TOOL -- run ${TOOL}_SETUP, then activate ${TOOL}_ENV.
# PREP falls back to QUPKAKE when PREP_ENV is empty. DOCK (db2tool, a compiled
# binary) needs no conda env: with DOCK_ENV empty only DOCK_SETUP runs.
llb_activate() {
    local tool="$1" setup_var="${1}_SETUP" env_var="${1}_ENV"
    local setup="${!setup_var:-}" env="${!env_var:-}"
    # inside the image, with_env.sh has already entered the env (micromamba run)
    if [[ -n "${LLB_IN_IMAGE:-}" ]] && llb_image_env "$tool" >/dev/null; then
        echo "[$tool] env=${CONDA_DEFAULT_ENV:-?} (image) python=$(command -v python)"
        return 0
    fi
    if [[ "$tool" == PREP && -z "$env" ]]; then
        setup="${QUPKAKE_SETUP:-}"; env="${QUPKAKE_ENV:-}"
    fi
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
    [[ -n "${XTBPATH:-}" ]] && export XTBPATH
    # ZBH_DIR set: unicon/conformator resolve to the shims, which add the license
    [[ "$tool" == DB2C && -n "${ZBH_DIR:-}" ]] && PATH="$LLB_ROOT/container/shims:$PATH"
    echo "[$tool] env=$env python=$(command -v python)"
}

# ---- the toolchain image (docs/INSTALL_DESIGN.md) ------------------------------
# A tool runs inside LLB_IMAGE when LLB_IMAGE is set and its <TOOL>_ENV is empty;
# setting <TOOL>_ENV keeps that one tool native. DOCK is never in the image.

# llb_image_env TOOL -- the env inside the image for TOOL; fails when TOOL is native.
llb_image_env() {
    local env_var="${1}_ENV"
    [[ -n "${LLB_IMAGE:-}" && -z "${!env_var:-}" ]] || return 1
    case "$1" in
        QUPKAKE)  echo qupkake ;;
        DB2C)     echo db2_converter ;;
        PKA|PREP) echo llb_pka ;;
        *)        return 1 ;;
    esac
}

# llb_container -- apptainer or singularity (LLB_CONTAINER picks), after
# LLB_CONTAINER_SETUP (e.g. `module load singularity`).
llb_container() {
    local c rc=0
    set +u
    if [[ -n "${LLB_CONTAINER_SETUP:-}" ]]; then eval "$LLB_CONTAINER_SETUP" >/dev/null 2>&1 || rc=$?; fi
    set -u
    for c in ${LLB_CONTAINER:-apptainer singularity}; do
        command -v "$c" >/dev/null 2>&1 && { echo "$c"; return 0; }
    done
    echo "FATAL: no apptainer or singularity here (LLB_CONTAINER_SETUP='${LLB_CONTAINER_SETUP:-}', rc=$rc)." \
         "A login node may have none (TAMU ACES): run this in a job, or check with \`llb doctor --node\`." >&2
    return 2
}

# llb_image_exec ENV CMD... -- exec CMD in conda env ENV inside LLB_IMAGE, with
# the checkout, the site, the run directory, data, node scratch and ZBH bound.
llb_image_exec() {
    local env="$1" rt d; shift
    rt="$(llb_container)" || return 2
    [[ -f "$LLB_IMAGE" ]] || { echo "FATAL: LLB_IMAGE=$LLB_IMAGE does not exist (in $LLB_CONFIG)" >&2; return 2; }
    local binds=() extra=()
    IFS=, read -ra extra <<<"${LLB_BINDS:-}"
    for d in "$LLB_ROOT" "${LLB_SITE:-}" "$PWD" "${ZINC_DATA:-}" "${SLURM_TMPDIR:-}" "${TMPDIR:-}" \
             "${ZBH_DIR:-}" "${ZBH_LICENSE:-}" ${extra[@]+"${extra[@]}"}; do
        # /tmp is bound by default; binding it twice only adds a warning
        [[ -n "$d" && "$d" != /tmp && -e "${d%%:*}" ]] && binds+=(--bind "$d")
    done
    # the runtime's cache on node scratch, never ~ (ACES: 10 GB home)
    local cache="${TMPDIR:-/tmp}/llb-container-cache"
    export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$cache}" SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$cache}"
    # the ZBH shims (in the bound checkout) first on PATH inside the image
    export APPTAINERENV_PREPEND_PATH="$LLB_ROOT/container/shims"
    export SINGULARITYENV_PREPEND_PATH="$LLB_ROOT/container/shims"
    export LLB_IN_IMAGE=1
    # ${a[@]+...}: an empty array is "unbound" to bash < 4.4 under set -u
    exec "$rt" exec ${binds[@]+"${binds[@]}"} "$LLB_IMAGE" micromamba run -n "$env" "$@"
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
