#!/bin/bash
# =====================================================================
# db2pipe -- parallel chunking of a large db2 library on SLURM
#
#   bash submit.slurm
#
# The plan step runs up front (via srun, blocking) so the finished-chunk
# count is known before the arrays are sized; then three jobs are chained:
#   plan     now      walk the tree, build shard lists, size BINS from a sample
#   collect  S tasks  parse shards, stride records into bin parts   <-- the work
#   assemble N tasks  concatenate parts into chunks (pure byte copy)
#   final    1 task   stitch the global manifest
#
# SETTINGS -- each is overridable from the environment, so the checkout never
# has to be edited:
#
#   RUNDIR=/scratch/.../run INPUT_DIR=/scratch/.../run/db2_run/out bash parallel/submit.slurm
#   llb chunk RUNDIR INPUT_DIR [opts]          # same thing, from ligand_library_builder
#
# With RUNDIR set, WORK_DIR and OUT_DIR default to RUNDIR/work and RUNDIR/chunks.
# ACCOUNT/PARTITION fall back to LLB_ACCOUNT/LLB_PARTITION (config/hpc.env, which
# `llb chunk` loads). Editing the defaults below in place still works too.
# =====================================================================
RUNDIR="${RUNDIR:-}"
INPUT_DIR="${INPUT_DIR:-/path/to/ZINC_sync/published/3D}"   # tree of .db2.gz / .db2.tgz
WORK_DIR="${WORK_DIR:-${RUNDIR:+$RUNDIR/work}}"             # intermediates; ~1x input size
WORK_DIR="${WORK_DIR:-${SCRATCH:-}/db2work}"
OUT_DIR="${OUT_DIR:-${RUNDIR:+$RUNDIR/chunks}}"             # finished chunks + manifest.tsv
OUT_DIR="${OUT_DIR:-${SCRATCH:-}/db2chunks}"
PIPE=$(cd "$(dirname "$0")" && pwd)         # parallel/ -- these scripts live here
SERIAL=$(cd "$PIPE/../serial" && pwd)       # serial/  -- 02_label, 03_apply live here
PYTHON="${PYTHON:-python3}"                 # stdlib only; no conda env needed

SHARDS="${SHARDS:-200}"                     # collect tasks (array width)
TARGET_PER_BIN="${TARGET_PER_BIN:-50000}"   # target molecules per finished chunk
BINS="${BINS:-}"                            # EMPTY = auto-size from TARGET_PER_BIN;
                                            # or a count, to skip the estimate
MODE="${MODE:-stride}"                      # stride | greedy
WEIGHT="${WEIGHT:-count}"                   # count | bytes | lines:C  (greedy only)

ACCOUNT="${ACCOUNT:-${LLB_ACCOUNT:-your_account}}"
PARTITION="${PARTITION:-${LLB_PARTITION:-standard}}"
CPUS="${CPUS:-1}"                           # each task is single-threaded
MEM="${MEM:-4G}"
PLAN_TIME="${PLAN_TIME:-01:00:00}"
COLLECT_TIME="${COLLECT_TIME:-04:00:00}"
ASSEMBLE_TIME="${ASSEMBLE_TIME:-00:30:00}"
COLLECT_PARALLEL="${COLLECT_PARALLEL:-$SHARDS}"   # %N throttle on the collect array

# =====================================================================
set -euo pipefail
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "FATAL: INPUT_DIR '$INPUT_DIR' is not a directory -- set INPUT_DIR (or edit the default)." >&2
  exit 2
fi
if [[ "$ACCOUNT" == your_account ]]; then
  echo "FATAL: no SLURM account -- set ACCOUNT, or LLB_ACCOUNT in config/hpc.env." >&2
  exit 2
fi
mkdir -p "$WORK_DIR" "$OUT_DIR" "$WORK_DIR/logs"
echo "input    : $INPUT_DIR"
echo "work/out : $WORK_DIR  ->  $OUT_DIR"

SB="sbatch --parsable -A $ACCOUNT -p $PARTITION --cpus-per-task=$CPUS --mem=$MEM"

# ---- plan (runs now; BINS must be known before the arrays are sized) -
# srun runs the walk (+ sample, when auto-sizing) on a compute node and
# blocks until it finishes. No srun on your cluster? Run make_shards
# yourself, set BINS above, and delete this block.
if [ -n "$BINS" ]; then
  srun -A "$ACCOUNT" -p "$PARTITION" --cpus-per-task="$CPUS" --mem=8G \
    --time="$PLAN_TIME" -o "$WORK_DIR/logs/plan.out" -e "$WORK_DIR/logs/plan.err" \
    "$PYTHON" -u "$PIPE/make_shards.py" -i "$INPUT_DIR" -o "$WORK_DIR" -S "$SHARDS"
else
  srun -A "$ACCOUNT" -p "$PARTITION" --cpus-per-task="$CPUS" --mem=8G \
    --time="$PLAN_TIME" -o "$WORK_DIR/logs/plan.out" -e "$WORK_DIR/logs/plan.err" \
    "$PYTHON" -u "$PIPE/make_shards.py" -i "$INPUT_DIR" -o "$WORK_DIR" -S "$SHARDS" \
    --target-per-bin "$TARGET_PER_BIN"
  BINS=$(cat "$WORK_DIR/bins.txt")
fi
echo "plan     : done, BINS=$BINS  (~$TARGET_PER_BIN mol/chunk)"

# ---- collect (the expensive phase) ----------------------------------
JOB_COLLECT=$($SB --job-name=db2collect --time=$COLLECT_TIME \
  --array=0-$((SHARDS-1))%$COLLECT_PARALLEL \
  -o "$WORK_DIR/logs/collect_%a.out" -e "$WORK_DIR/logs/collect_%a.err" \
  --wrap "$PYTHON -u $PIPE/p1_collect.py -w '$WORK_DIR' -s \$SLURM_ARRAY_TASK_ID \
          -N $BINS --mode $MODE --weight $WEIGHT")
echo "collect  : $JOB_COLLECT  (array 0-$((SHARDS-1)))"

# ---- assemble -------------------------------------------------------
# afterok on an array job waits for EVERY task to succeed.
JOB_ASSEMBLE=$($SB --job-name=db2assemble --time=$ASSEMBLE_TIME \
  --array=0-$((BINS-1)) --dependency=afterok:$JOB_COLLECT \
  -o "$WORK_DIR/logs/assemble_%a.out" -e "$WORK_DIR/logs/assemble_%a.err" \
  --wrap "$PYTHON -u $PIPE/p2_assemble.py -w '$WORK_DIR' -b \$SLURM_ARRAY_TASK_ID \
          -o '$OUT_DIR'")
echo "assemble : $JOB_ASSEMBLE  (array 0-$((BINS-1)))"

# ---- finalize -------------------------------------------------------
JOB_FINAL=$($SB --job-name=db2final --time=01:00:00 --mem=8G \
  --dependency=afterok:$JOB_ASSEMBLE \
  -o "$WORK_DIR/logs/final.out" -e "$WORK_DIR/logs/final.err" \
  --wrap "$PYTHON -u $PIPE/p3_finalize.py -w '$WORK_DIR' -o '$OUT_DIR'")
echo "final    : $JOB_FINAL"

cat <<EOF

Submitted. Watch with:
  squeue -u \$USER
  tail -f $WORK_DIR/logs/collect_0.err

When it finishes, chunks and manifest.tsv are in:
  $OUT_DIR

Then label as usual (these are unchanged and still serial):
  python3 $SERIAL/02_label.py -m $OUT_DIR/manifest.tsv -o labels.tsv
  python3 $SERIAL/03_apply.py -c $OUT_DIR -L labels.tsv -o ${OUT_DIR}_labelled

If a few collect tasks fail, rerun just those:
  sbatch --array=3,17,102 ... p1_collect.py ...
then resubmit assemble and final.
EOF
