#!/bin/bash
# One MolSanitizer build shard per array task. Submitted by msani.sh submit, which
# sets --account/--partition/--array/--output and passes the arguments:
#
#   build.sbatch RUN OFFSET MSB_DIR PREFIX CONDA_BASE
#
# shard = SLURM_ARRAY_TASK_ID + OFFSET (the offset keeps array indices under
# Great Lakes' MaxArraySize when a run has more shards than that).
# build_shard.py skips a shard that is already done, so resubmitting is safe.
#SBATCH --job-name=msani
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=08:00:00
set -euo pipefail

RUN="$1"; OFFSET="$2"; MSB_DIR="$3"; PREFIX="$4"; CONDA_BASE="$5"
SHARD=$(( SLURM_ARRAY_TASK_ID + OFFSET ))

set +u
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$PREFIX"
set -u

echo "shard $SHARD on $(hostname) at $(date +%Y-%m-%dT%H:%M:%S%z); job ${SLURM_ARRAY_JOB_ID:-}_${SLURM_ARRAY_TASK_ID:-}"
exec python "$MSB_DIR/build_shard.py" --run "$RUN" --shard "$SHARD"
