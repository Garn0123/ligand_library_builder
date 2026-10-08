#!/bin/bash
# One array task = NCPUS concurrent build_ligand processes, one per shard.
# Work happens on node-local scratch and is rsynced back periodically so a
# wall-clock kill doesn't lose everything.
#
# Submit with run_qupkake/submit_db2.sh (or `llb submit-db2`), which creates
# the run directory, shards the library and passes RUNDIR. Everything this job
# reads or writes is under RUNDIR; the repo is only read. Run settings (METHOD,
# NCONF, MMFF, NCPUS, NSHARDS) live in RUNDIR/config.env; where the software
# lives (DB2C_SETUP, DB2C_ENV, BUILD_LIGAND_EXE) lives in config/hpc.env.
#
# Adapted from John's db2_array.sbatch (2026-08-27): same preflight, provenance,
# scratch + checkpointing and resume logic; the hard-coded environment lines
# (~/env/bcl.sh, the miniforge mamba.sh, `mamba activate db2_converter`) now
# come from config/hpc.env via common/env.sh.

set -uo pipefail

RUN="${RUNDIR:?RUNDIR not set}"
# shellcheck disable=SC1091
source "$RUN/config.env"

# --- environment -------------------------------------------------------
# DB2C_SETUP runs BEFORE activation: bcl.sh runs `module purge`, and a purge
# after activation resets PATH/LD_LIBRARY_PATH out from under the conda env. It
# also loads gcc/8.2.0 (libgfortran.so.5 for AMSOL) and singularity (the
# conformator/UNICON wrappers).
# shellcheck disable=SC1091
source "${LLB_ROOT:?config.env has no LLB_ROOT}/common/env.sh" \
    || { echo "cannot source $LLB_ROOT/common/env.sh"; exit 78; }
llb_load_config || exit 78
llb_activate DB2C || exit 78
BUILD_LIGAND="$(llb_require_exe BUILD_LIGAND_EXE build_ligand)" || exit 78
export OMP_NUM_THREADS=1                      # one thread per build, not per node

# --- preflight ---------------------------------------------------------
# Fail the task immediately rather than burning the wall clock producing
# nothing. A silently-missing module is the single most expensive failure
# mode at this level of parallelism.
preflight_fail=0
check() {
    if ! eval "$2" >/dev/null 2>&1; then
        echo "PREFLIGHT FAIL: $1"
        preflight_fail=1
    fi
}
check "singularity not on PATH (module purge without reload?)" "command -v singularity"
check "build_ligand not executable ($BUILD_LIGAND)"             "[[ -x '$BUILD_LIGAND' ]]"
check "python cannot import rdkit"                             "python -c 'import rdkit'"

AMSOL=$(python -c "from db2_converter.config import config; print(config['all']['AMSOLEXE'])" 2>/dev/null)
if [[ -n "$AMSOL" ]]; then
    missing=$(ldd "$AMSOL" 2>/dev/null | grep 'not found' || true)
    if [[ -n "$missing" ]]; then
        echo "PREFLIGHT FAIL: AMSOL has unresolved libraries:"
        echo "$missing"
        preflight_fail=1
    fi
else
    echo "PREFLIGHT FAIL: could not resolve AMSOLEXE from db2_converter.config"
    preflight_fail=1
fi

if (( preflight_fail )); then
    echo "aborting task $SLURM_ARRAY_TASK_ID on $(hostname) -- environment is broken"
    exit 78
fi
echo "preflight OK"

flags=""
[[ "$MMFF" == "on" ]] && flags="--MMFFopt"

# --- provenance --------------------------------------------------------
# bcl.sh runs `module purge` and reloads, so what the job ran with is not
# what was loaded at submit time. Record the post-setup state.
envfile="$RUN/logs/env.task$(printf '%03d' "$SLURM_ARRAY_TASK_ID").txt"
{
    echo "=== task ==="
    echo "date        $(date -Is 2>/dev/null || date)"
    echo "host        $(hostname)"
    echo "jobid       ${SLURM_JOB_ID:-}  array_task ${SLURM_ARRAY_TASK_ID}"
    echo "cpu         $(lscpu 2>/dev/null | awk -F: '/Model name/{gsub(/^ +/,"",$2); print $2; exit}')"
    echo "cpus_task   ${SLURM_CPUS_PER_TASK:-$NCPUS}"
    echo "scratch     ${SLURM_TMPDIR:-/tmp (SLURM_TMPDIR unset)}"
    echo
    echo "=== modules (after DB2C_SETUP purge+reload) ==="
    module list 2>&1
    echo
    echo "=== run settings ==="
    echo "method      $METHOD"
    echo "nconf       $NCONF"
    echo "mmffopt     $MMFF"
    echo "flags       $flags"
    echo "build       $BUILD_LIGAND"
} > "$envfile"
cat "$envfile"

# Full toolchain provenance, written once per *submission* by whichever task
# arrives first. Keyed on the array job id, so a resubmission after an
# environment change records the new state instead of reusing stale output.
prov="$RUN/provenance.${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}.txt"
if [[ ! -f "$prov" ]]; then
    {
        echo "=== recorded $(date -Is 2>/dev/null || date) by task ${SLURM_ARRAY_TASK_ID} on $(hostname) ==="
        echo
        echo "=== ligand_library_builder ==="
        echo "LLB_ROOT     $LLB_ROOT"
        git -C "$LLB_ROOT" rev-parse HEAD 2>/dev/null | sed 's/^/commit       /'
        git -C "$LLB_ROOT" status --porcelain 2>/dev/null | sed 's/^/dirty        /'
        echo "LLB_CONFIG   $LLB_CONFIG"
        grep -E '^(CONDA_BASE|DB2C_[A-Z]+|BUILD_LIGAND_EXE)=' "$LLB_CONFIG" | sed 's/^/  /'
        echo
        echo "=== conda env ==="
        echo "CONDA_PREFIX  ${CONDA_PREFIX:-unset}"
        python -V 2>&1
        python -c "import rdkit; print('rdkit        ', rdkit.__version__)" 2>&1
        echo
        echo "=== db2_converter ==="
        pkgdir=$(python -c "import db2_converter,os; print(os.path.dirname(db2_converter.__file__))" 2>/dev/null)
        echo "package_dir  ${pkgdir:-UNKNOWN}"
        if [[ -n "$pkgdir" ]]; then
            git -C "$pkgdir" rev-parse HEAD 2>/dev/null \
                | sed 's/^/commit       /' || echo "commit       (not a git checkout)"
            git -C "$pkgdir" status --porcelain 2>/dev/null \
                | sed 's/^/dirty        /'
        fi
        echo
        echo "=== resolved tool paths (db2_converter.config) ==="
        python -c "from db2_converter.config import config; import json; print(json.dumps(config, indent=2, default=str))" 2>&1
        echo
        echo "=== binaries ==="
        for exe in "$AMSOL" "$BUILD_LIGAND" "$(command -v conformator)" "$(command -v unicon)" \
                   "$(command -v antechamber)" "$(command -v obabel)" \
                   "$(command -v singularity)"; do
            [[ -n "$exe" && -e "$exe" ]] || continue
            printf '%s\n' "$exe"
            ls -lL "$exe" | sed 's/^/  /'
            md5sum "$exe" 2>/dev/null | sed 's/^/  md5 /'
        done
        echo
        echo "=== BCL ==="
        echo "BCL_ROOT     ${BCL_ROOT:-unset}"
        ls -l "${BCL_ROOT:-/nonexistent}/bin" 2>/dev/null | head -20 | sed 's/^/  /'
    } > "$prov" 2>&1
    echo "provenance written to $prov"
fi

# Node-local scratch, unique per job AND per array task. Without the job id
# two tasks landing on one node would share a directory.
scratch="${SLURM_TMPDIR:-/tmp}/db2.${SLURM_JOB_ID:-manual}.${SLURM_ARRAY_TASK_ID}"
mkdir -p "$scratch"

# --- checkpointing -----------------------------------------------------
sync_back() {
    for k in $(seq 0 $((NCPUS-1))); do
        s=$(( SLURM_ARRAY_TASK_ID * NCPUS + k ))
        (( s < NSHARDS )) || continue
        sid=$(printf "%04d" "$s")
        [[ -d "$scratch/$sid" ]] && rsync -a "$scratch/$sid/" "$RUN/out/$sid/"
    done
}
trap 'echo "[signal] syncing before exit"; kill "${syncer:-}" 2>/dev/null; sync_back; exit 1' USR1 TERM
trap 'kill "${syncer:-}" 2>/dev/null; sync_back' EXIT

( while true; do sleep 900; sync_back; done ) >/dev/null 2>&1 &
syncer=$!

# GNU time's -v resource report where available (Linux); plain run elsewhere.
timer=(/usr/bin/time -v)
/usr/bin/time -v true >/dev/null 2>&1 || timer=()

# --- launch one build per core ----------------------------------------
pids=()
for k in $(seq 0 $((NCPUS-1))); do
    s=$(( SLURM_ARRAY_TASK_ID * NCPUS + k ))
    (( s < NSHARDS )) || continue
    sid=$(printf "%04d" "$s")
    shard="$RUN/shards/s.$sid.smi"
    [[ -s "$shard" ]] || continue

    wd="$scratch/$sid"
    mkdir -p "$wd" "$RUN/out/$sid"

    # Pull back any partial results so build_ligand's resume logic can skip
    # molecules finished in an earlier attempt. Requires workingpath ==
    # outputpath, which is how they're set below.
    rsync -a "$RUN/out/$sid/" "$wd/" 2>/dev/null

    (
      ${timer[@]+"${timer[@]}"} "$BUILD_LIGAND" -i "$shard" \
        -n "$NCONF" --keep_max_conf \
        -m "$METHOD" --checkstereo --reseth --rotateh --dock38 \
        $flags \
        --workingpath "$wd" --outputpath "$wd"
    ) > "$RUN/logs/shard.$sid.out" 2>&1 &
    pids+=($!)
done

echo "launched ${#pids[@]} builds"
rc=0
for p in "${pids[@]}"; do wait "$p" || rc=1; done

kill "$syncer" 2>/dev/null
sync_back
# The periodic sync copies a molecule's working directory (NAME/) while it
# builds; build_ligand deletes it from scratch when done, but rsync without
# --delete leaves the copy in out/. Every build has ended here, so remove any
# out/ subdirectory scratch no longer has.
# Only directories named like a contract name (C + 14 base62 + C) are touched,
# and only when RUN and scratch are both set.
if [[ -n "${RUN:-}" && -n "${scratch:-}" && -d "$RUN/out" ]]; then
    for k in $(seq 0 $((NCPUS-1))); do
        s=$(( SLURM_ARRAY_TASK_ID * NCPUS + k ))
        (( s < NSHARDS )) || continue
        sid=$(printf "%04d" "$s")
        for d in "$RUN/out/$sid"/*/; do
            n="$(basename "$d")"
            [[ -d "$d" && "$n" =~ ^C[0-9A-Za-z]{14}C$ && ! -d "$scratch/$sid/$n" ]] || continue
            echo "removing stale working copy $RUN/out/$sid/$n"
            rm -rf -- "$RUN/out/$sid/$n"
        done
    done
fi
echo "task done rc=$rc"
exit $rc
