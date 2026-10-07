#!/bin/bash
# submit_db2.sh RUNDIR LIBRARY_SMI [options] [-- extra sbatch args]
#
# Set up a db2 build run directory and submit slurm/db2_array.sbatch over it.
#
#   run_qupkake/submit_db2.sh db2_run library/library.smi
#   run_qupkake/submit_db2.sh db2_run library/library.smi --ncpus 36 --nshards 360 -- --time=48:00:00
#
# RUNDIR gets everything the run produces; the repo is only read:
#   config.env   run settings (method, nconf, mmff, ncpus, nshards) + LLB_ROOT/LLB_CONFIG
#   shards/      s.NNNN.smi, round-robin split of LIBRARY_SMI
#   out/NNNN/    NAME.db2.gz, conformer.NAME.fixed.mol2, the shard's faillist
#   logs/        SLURM logs, per-shard build_ligand logs, env.taskNNN.txt
#   provenance.<jobid>.txt
#
# Options (defaults from the db2_converter notebook, 2026-08-27):
#   --method M    conformator      --nconf N   600
#   --mmff on|off off              --ncpus N   8 builds per array task
#   --nshards N   100              --mem-per-cpu 2G   --time 24:00:00
#
# Resubmitting with an existing RUNDIR resumes: shards are NOT re-split (a new
# split would scatter half-finished work across different shard numbers), and
# each task rsyncs its finished molecules back to scratch so build_ligand skips
# them. Never pass --rerun to build_ligand in an array: it deletes both paths.

set -euo pipefail

[[ $# -ge 2 ]] || { sed -n '2,/^$/p' "$0" >&2; exit 64; }
RUN="$1"; LIB="$2"; shift 2
METHOD=conformator; NCONF=600; MMFF=off; NCPUS=8; NSHARDS=100
MEM_PER_CPU=2G; TIME=24:00:00
while [[ $# -gt 0 ]]; do
    case "$1" in
        --method)      METHOD="$2"; shift 2 ;;
        --nconf)       NCONF="$2"; shift 2 ;;
        --mmff)        MMFF="$2"; shift 2 ;;
        --ncpus)       NCPUS="$2"; shift 2 ;;
        --nshards)     NSHARDS="$2"; shift 2 ;;
        --mem-per-cpu) MEM_PER_CPU="$2"; shift 2 ;;
        --time)        TIME="$2"; shift 2 ;;
        --)            shift; break ;;
        *) echo "unknown option $1" >&2; exit 64 ;;
    esac
done
[[ "$MMFF" == on || "$MMFF" == off ]] || { echo "--mmff must be on or off" >&2; exit 64; }

LLB_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export LLB_ROOT
# shellcheck source=common/env.sh
source "$LLB_ROOT/common/env.sh"
llb_load_config

mkdir -p "$RUN"/{shards,out,logs}
RUN="$(cd "$RUN" && pwd)"
[[ -s "$LIB" ]] || { echo "FATAL: $LIB is empty or missing" >&2; exit 2; }
LIB="$(cd "$(dirname "$LIB")" && pwd)/$(basename "$LIB")"

# ---- shards: once per RUNDIR -------------------------------------------------
if [[ -f "$RUN/config.env" ]]; then
    # shellcheck disable=SC1091
    prev_lib="$(sed -n 's/^LIBRARY_SMI=//p' "$RUN/config.env")"
    prev_sha="$(sed -n 's/^LIBRARY_SHA256=//p' "$RUN/config.env")"
    now_sha="$(sha256sum "$LIB" 2>/dev/null | awk '{print $1}' || shasum -a 256 "$LIB" | awk '{print $1}')"
    if [[ "$prev_sha" != "$now_sha" ]]; then
        echo "FATAL: $RUN was set up from $prev_lib (sha256 $prev_sha);" >&2
        echo "       $LIB now has sha256 $now_sha. Use a new RUNDIR for a new library." >&2
        exit 2
    fi
    NSHARDS="$(sed -n 's/^NSHARDS=//p' "$RUN/config.env")"
    echo "== resuming $RUN ($NSHARDS existing shards; settings below override method/nconf/mmff/ncpus)"
else
    n_mol=$(grep -c . "$LIB")
    (( NSHARDS > n_mol )) && NSHARDS=$n_mol
    # Round-robin, not contiguous: per-molecule cost varies widely and input
    # files cluster similar compounds (db2_converter notebook, section 4).
    awk -v n="$NSHARDS" -v d="$RUN/shards" \
        'NF { f = sprintf("%s/s.%04d.smi", d, (k++) % n); print > f }' "$LIB"
    echo "== split $n_mol molecules into $NSHARDS shards under $RUN/shards"
fi
LIB_SHA="$(sha256sum "$LIB" 2>/dev/null | awk '{print $1}' || shasum -a 256 "$LIB" | awk '{print $1}')"
cat > "$RUN/config.env" <<EOF
# written by submit_db2.sh $(date '+%Y-%m-%d %H:%M:%S')
LLB_ROOT=$LLB_ROOT
LLB_CONFIG=$LLB_CONFIG
LIBRARY_SMI=$LIB
LIBRARY_SHA256=$LIB_SHA
METHOD=$METHOD
NCONF=$NCONF
MMFF=$MMFF
NCPUS=$NCPUS
NSHARDS=$NSHARDS
EOF

# ---- preflight on this node: the same checks every task makes ----------------
# Every check exits explicitly: `set -e` is suspended inside a subshell on the
# left of `||`, so a failing command would otherwise fall through to sbatch.
echo "== preflight (config: $LLB_CONFIG)"
( llb_activate DB2C || exit 2
  bl="$(llb_require_exe BUILD_LIGAND_EXE build_ligand)" || exit 2
  echo "  build_ligand  $bl"
  python -c 'import rdkit, db2_converter' \
      || { echo "FATAL: rdkit or db2_converter not importable in ${DB2C_ENV}" >&2; exit 2; }
  command -v singularity >/dev/null \
      || echo "  WARNING: singularity not on PATH here; conformator/UNICON need it on the nodes"
) || exit 2

ntasks=$(( (NSHARDS + NCPUS - 1) / NCPUS ))
cd "$RUN"
echo "== submitting $ntasks task(s) x $NCPUS builds, account=$LLB_ACCOUNT partition=$LLB_PARTITION"
sbatch -A "$LLB_ACCOUNT" -p "$LLB_PARTITION" \
    --job-name=db2 --array="0-$((ntasks - 1))%${LLB_ARRAY_PARALLEL:-50}" \
    --cpus-per-task="$NCPUS" --mem-per-cpu="$MEM_PER_CPU" --time="$TIME" \
    --output="$RUN/logs/task_%A_%a.out" \
    --export=ALL,RUNDIR="$RUN" \
    "$@" "$LLB_ROOT/run_qupkake/slurm/db2_array.sbatch"
