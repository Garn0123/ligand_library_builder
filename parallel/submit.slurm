#!/bin/bash
# =====================================================================
# db2pipe -- parallel chunking of a large db2 library on SLURM
#
#   bash submit.slurm
#
# Four jobs, chained by dependency:
#   plan     1 task    walk the tree, build size-balanced shard lists
#   collect  S tasks   parse shards, stride records into bin parts   <-- the work
#   assemble N tasks   concatenate parts into chunks (pure byte copy)
#   final    1 task    stitch the global manifest
#
# EDIT THESE:
# =====================================================================
INPUT_DIR=/path/to/ZINC_sync/published/3D
WORK_DIR=$SCRATCH/db2work           # intermediates; needs ~1x input size
OUT_DIR=$SCRATCH/db2chunks          # finished chunks
PIPE=$(cd "$(dirname "$0")" && pwd)         # parallel/ -- these scripts live here
SERIAL=$(cd "$PIPE/../serial" && pwd)       # serial/  -- 02_label, 03_apply live here

SHARDS=200                          # collect tasks (array width)
BINS=170                            # finished chunks (~50k mols each for 8.5M)
MODE=stride                         # stride | greedy
WEIGHT=count                        # count | bytes | lines:C  (greedy only)

ACCOUNT=your_account
PARTITION=standard
CPUS=1                              # each task is single-threaded
MEM=4G
COLLECT_TIME=04:00:00
ASSEMBLE_TIME=00:30:00

# =====================================================================
set -euo pipefail
mkdir -p "$WORK_DIR" "$OUT_DIR" "$WORK_DIR/logs"

SB="sbatch --parsable -A $ACCOUNT -p $PARTITION --cpus-per-task=$CPUS --mem=$MEM"

# ---- plan -----------------------------------------------------------
JOB_PLAN=$($SB --job-name=db2plan --time=01:00:00 \
  -o "$WORK_DIR/logs/plan.out" -e "$WORK_DIR/logs/plan.err" \
  --wrap "python3 -u $PIPE/make_shards.py -i '$INPUT_DIR' -o '$WORK_DIR' -S $SHARDS")
echo "plan     : $JOB_PLAN"

# ---- collect (the expensive phase) ----------------------------------
JOB_COLLECT=$($SB --job-name=db2collect --time=$COLLECT_TIME \
  --array=0-$((SHARDS-1))%$SHARDS --dependency=afterok:$JOB_PLAN \
  -o "$WORK_DIR/logs/collect_%a.out" -e "$WORK_DIR/logs/collect_%a.err" \
  --wrap "python3 -u $PIPE/p1_collect.py -w '$WORK_DIR' -s \$SLURM_ARRAY_TASK_ID \
          -N $BINS --mode $MODE --weight $WEIGHT")
echo "collect  : $JOB_COLLECT  (array 0-$((SHARDS-1)))"

# ---- assemble -------------------------------------------------------
# afterok on an array job waits for EVERY task to succeed.
JOB_ASSEMBLE=$($SB --job-name=db2assemble --time=$ASSEMBLE_TIME \
  --array=0-$((BINS-1)) --dependency=afterok:$JOB_COLLECT \
  -o "$WORK_DIR/logs/assemble_%a.out" -e "$WORK_DIR/logs/assemble_%a.err" \
  --wrap "python3 -u $PIPE/p2_assemble.py -w '$WORK_DIR' -b \$SLURM_ARRAY_TASK_ID \
          -o '$OUT_DIR'")
echo "assemble : $JOB_ASSEMBLE  (array 0-$((BINS-1)))"

# ---- finalize -------------------------------------------------------
JOB_FINAL=$($SB --job-name=db2final --time=01:00:00 --mem=8G \
  --dependency=afterok:$JOB_ASSEMBLE \
  -o "$WORK_DIR/logs/final.out" -e "$WORK_DIR/logs/final.err" \
  --wrap "python3 -u $PIPE/p3_finalize.py -w '$WORK_DIR' -o '$OUT_DIR'")
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
