#!/bin/bash
# msani.sh -- MolSanitizer as LLB's ligand-construction step (`llb msani ...`).
#
# Two MolSanitizer passes with LLB's naming step between them (UPSTREAM.md section 2):
#
#   msani.sh install [PREFIX] [--from-lock]   create + freeze the env (environment/install.sh)
#   msani.sh plan INPUT.smi... -o RUN [opts]  shard `SMILES NAME` files for pass 1 (plan.py -h)
#   msani.sh submit RUN 1 [opts]              pass 1: chemistry, SMILES out (array; logs -> RUN/logs/)
#   msani.sh names RUN [--per-shard N]        `llb names` over all of pass 1, then pass 2's shards
#                                             (needs RDKit: runs in the msani env; srun it at scale)
#   msani.sh submit RUN 2 [opts]              pass 2: 3D + db2 on contract names (array)
#   msani.sh status RUN                       pass 1 / names / pass 2: done and missing, as array ranges
#   msani.sh collect RUN [--drap-dir DIR]     verify + account -> db2_verification.csv, ledger
#   msani.sh pack RUN [--delete-loose] [opts] per-name db2.gz -> 1,000-name chunks + db2_index.tsv
#                                             (pack.py -h; srun it: it reads every db2)
#   msani.sh shard RUN PASS K                 run one shard here (inside an srun/salloc)
#   msani.sh selftest [DIR]                   all of it on probes/selftest.smi, serially (srun it)
#   msani.sh abtest RUN [K] [N]               time N names of pass-2 shard K (default 1, 40) with and
#                                             without the forkserver preload (srun it; writes nothing in RUN)
#   msani.sh versions                         what the env holds
#
#   submit opts: [--tasks all|missing|LIST] [--parallel N] [--time T] [--mem M] [-- sbatch args]
#
# Settings (environment variables; defaults are Lighthouse, John's install, 2026-10-26):
#   LLB_MSANI_PREFIX      the env       /nfs/turbo/umms-maom/johnbick/opt/molsanitizer/0.7.4
#   LLB_MSANI_CONDA_BASE  conda         $CONDA_BASE, else `conda info --base`
#   LLB_MSANI_ACCOUNT     SLURM account maom
#   LLB_MSANI_PARTITION                 maom-h200
#   LLB_MSANI_PARALLEL    tasks running at once per array (each 1 CPU)   140
#   LLB_MSANI_MAX_CPUS    hard cap on --parallel (John: stay under 150 CPUs)  149
#   LLB_MSANI_MAX_SUBMIT  most tasks one submit may queue                4500
#
# The CPU cap holds across arrays too: submit refuses while another msani1/msani2
# job of yours is queued or running (pass --force-concurrent to override).
set -euo pipefail

MSB_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
PREFIX="${LLB_MSANI_PREFIX:-/nfs/turbo/umms-maom/johnbick/opt/molsanitizer/0.7.4}"
ACCOUNT="${LLB_MSANI_ACCOUNT:-maom}"
PARTITION="${LLB_MSANI_PARTITION:-maom-h200}"
PARALLEL="${LLB_MSANI_PARALLEL:-140}"
MAX_CPUS="${LLB_MSANI_MAX_CPUS:-149}"
MAX_SUBMIT="${LLB_MSANI_MAX_SUBMIT:-4500}"
# array indices must stay below the cluster's MaxArraySize
MAX_ARRAY_INDEX=$(( $(scontrol show config 2>/dev/null | awk '$1=="MaxArraySize"{print $3}' || true) + 0 ))
(( MAX_ARRAY_INDEX > 0 )) && MAX_ARRAY_INDEX=$(( MAX_ARRAY_INDEX - 1 )) || MAX_ARRAY_INDEX=4999

usage() { sed -n '2,/^set -euo/p' "$MSB_DIR/msani.sh" | sed '$d' | sed 's/^# \{0,1\}//'; }

conda_base() {
    local b="${LLB_MSANI_CONDA_BASE:-${CONDA_BASE:-}}"
    [[ -n "$b" ]] || b="$(conda info --base 2>/dev/null || true)"
    [[ -f "$b/etc/profile.d/conda.sh" ]] || {
        echo "msani.sh: no conda at '$b'; set LLB_MSANI_CONDA_BASE" >&2; exit 1; }
    echo "$b"
}

activate() {
    local b; b="$(conda_base)"
    [[ -x "$PREFIX/bin/python" ]] || { echo "msani.sh: no env at $PREFIX (msani.sh install)" >&2; exit 1; }
    set +u; source "$b/etc/profile.d/conda.sh"; conda activate "$PREFIX"; set -u
}

step="${1:-help}"; [[ $# -gt 0 ]] && shift
case "$step" in
    install)  exec "$MSB_DIR/environment/install.sh" "$@" ;;
    plan)     exec python3 "$MSB_DIR/plan.py" "$@" ;;
    status)   exec python3 "$MSB_DIR/collect.py" "$1" --status ;;
    collect)  exec python3 "$MSB_DIR/collect.py" "$@" ;;
    pack)     exec python3 "$MSB_DIR/pack.py" "$@" ;;
    shard)    activate; exec python "$MSB_DIR/build_shard.py" --run "$1" --pass "$2" --shard "$3" ;;
    names)    activate; exec python "$MSB_DIR/names_step.py" "$@" ;;
    versions)
        activate
        msani --version 2>&1 | sed -n '/version/p'; obabel -V
        [[ -f "$PREFIX/LOCK/versions.txt" ]] && { echo "--- $PREFIX/LOCK/versions.txt"; cat "$PREFIX/LOCK/versions.txt"; }
        ;;
    abtest)
        run="$(cd "${1:?msani.sh abtest RUN [SHARD] [N]}" && pwd)"; k="${2:-1}"; n="${3:-40}"
        activate
        base="${TMPDIR:-/tmp}/msani_abtest_$$"; mkdir -p "$base"
        head -n "$n" "$run/pass2/shards/$(printf %05d "$k").smi" > "$base/todo.smi"
        export MSANI_CONFIG="$run/msani_configurations.yaml" OMP_NUM_THREADS=1
        read -r -a flags < <(PYTHONPATH="$MSB_DIR" python -c \
            "import json, build_shard as b; print(' '.join(b.msani_flags(json.load(open('$run/run.json'))['options'])))")
        for mode in preload no-preload; do
            mkdir -p "$base/$mode"; cp "$base/todo.smi" "$base/$mode/"
            extra=(); [[ $mode == no-preload ]] && extra=(--no-preload)
            t0=$(date +%s)
            ( cd "$base/$mode" && python "$MSB_DIR/msani3d.py" --todo todo.smi --chunk "$n" "${extra[@]}" \
                -- "${flags[@]}" > out.txt 2> err.txt ) || echo "($mode exited non-zero; see $base/$mode/err.txt)"
            echo "$mode: $(( $(date +%s) - t0 )) s wall for $n names"
            python "$MSB_DIR/timing_summary.py" "$base/$mode/msani_timing.csv"
        done
        echo "kept in $base"
        ;;
    selftest)
        dir="${1:-${TMPDIR:-/tmp}/msani_selftest_$$}"
        [[ -e "$dir" ]] && { echo "msani.sh selftest: $dir exists" >&2; exit 1; }
        python3 "$MSB_DIR/plan.py" "$MSB_DIR/probes/selftest.smi" -o "$dir" --per-shard 100 --label selftest
        activate
        python "$MSB_DIR/build_shard.py" --run "$dir" --pass 1 --shard 1
        python "$MSB_DIR/names_step.py" "$dir" --per-shard 100
        python "$MSB_DIR/build_shard.py" --run "$dir" --pass 2 --shard 1
        python3 "$MSB_DIR/collect.py" "$dir"
        echo; echo "selftest run kept in $dir"
        ;;
    submit)
        run="$(cd "${1:?msani.sh submit RUN 1|2 ...}" && pwd)"; npass="${2:?msani.sh submit RUN 1|2 ...}"; shift 2
        [[ "$npass" == 1 || "$npass" == 2 ]] || { echo "msani.sh submit: pass is 1 or 2, not '$npass'" >&2; exit 64; }
        [[ "$npass" == 1 || -f "$run/pass2/plan.json" ]] || {
            echo "msani.sh submit: pass 2 has no shards yet; run msani.sh names $run first" >&2; exit 1; }
        tasks=all; time=""; mem=""; extra=(); concurrent=0
        while [[ $# -gt 0 ]]; do
            case "$1" in
                --tasks) tasks="$2"; shift 2 ;;
                --parallel) PARALLEL="$2"; shift 2 ;;
                --time) time="$2"; shift 2 ;;
                --mem) mem="$2"; shift 2 ;;
                --force-concurrent) concurrent=1; shift ;;
                --) shift; extra=("$@"); break ;;
                *) echo "msani.sh submit: unknown option $1" >&2; exit 64 ;;
            esac
        done
        (( PARALLEL >= 1 && PARALLEL <= MAX_CPUS )) || {
            echo "msani.sh submit: --parallel $PARALLEL is outside 1..$MAX_CPUS (LLB_MSANI_MAX_CPUS)" >&2; exit 1; }
        if (( ! concurrent )); then
            busy=$(squeue -h -u "$USER" -n msani1,msani2 -o '%i %j %T' 2>/dev/null | head -5 || true)
            [[ -z "$busy" ]] || {
                echo "msani.sh submit: msani jobs already queued/running (CPU cap $MAX_CPUS across arrays):" >&2
                echo "$busy" >&2; echo "wait for them, or --force-concurrent" >&2; exit 1; }
        fi
        [[ -f "$run/run.json" ]] || { echo "msani.sh submit: $run has no run.json (msani.sh plan)" >&2; exit 1; }
        # -> "OFFSET ARRAYSPEC COUNT"
        read -r offset spec count < <(PYTHONPATH="$MSB_DIR" python3 - "$run" "$tasks" "$MAX_ARRAY_INDEX" "$npass" <<'EOF'
import json, sys
from pathlib import Path
run, tasks, max_ix, npass = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
n = (json.loads((run / "run.json").read_text())["n_shards"] if npass == 1
     else json.loads((run / "pass2" / "plan.json").read_text())["n_shards"])
def expand(spec):
    out = []
    for part in spec.split(","):
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b or a) + 1))
    return out
if tasks == "all":
    want = list(range(1, n + 1))
elif tasks == "missing":
    from collect import shard_status
    want = shard_status(run, npass)[2]
else:
    want = expand(tasks)
bad = [k for k in want if not 1 <= k <= n]
if bad:
    sys.exit(f"shards out of range 1..{n}: {bad[:5]}")
if not want:
    print("0 - 0"); sys.exit()
offset = min(want) - 1
if max(want) - offset > max_ix:
    sys.exit(f"shards {min(want)}..{max(want)} span more than {max_ix} array indices; "
             f"submit them in pieces with --tasks A-B")
from collect import ranges
print(offset, ranges([k - offset for k in want]), len(want))
EOF
)
        [[ "$count" == 0 ]] && { echo "nothing to submit: every requested shard is done"; exit 0; }
        (( count > MAX_SUBMIT )) && {
            echo "msani.sh submit: $count tasks > LLB_MSANI_MAX_SUBMIT=$MAX_SUBMIT (MaxSubmit is per user," \
                 "and other jobs count too); submit with --tasks in pieces" >&2; exit 1; }
        mkdir -p "$run/logs"
        args=(--account="$ACCOUNT" --partition="$PARTITION" --array="${spec}%${PARALLEL}"
              --job-name="msani$npass" --output="$run/logs/msani${npass}_%A_%a.out")
        [[ -n "$time" ]] && args+=(--time="$time")
        [[ -n "$mem" ]] && args+=(--mem="$mem")
        jid=$(sbatch --parsable "${args[@]}" "${extra[@]}" "$MSB_DIR/slurm/build.sbatch" \
              "$run" "$npass" "$offset" "$MSB_DIR" "$PREFIX" "$(conda_base)")
        printf '%s\tpass%s\t%s\t%s\t%s\t%s\n' "$(date +%Y-%m-%dT%H:%M:%S%z)" "$npass" "$jid" "$offset" "$spec" "$count" >> "$run/submissions.tsv"
        echo "submitted pass $npass as $jid: $count task(s) (shard = array index + $offset; array $spec), $PARALLEL at once; logs in $run/logs/"
        ;;
    help|-h|--help) usage ;;
    *) echo "msani.sh: unknown step '$step'" >&2; usage >&2; exit 64 ;;
esac
