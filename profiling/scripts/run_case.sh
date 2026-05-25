#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

CASE_NAME="${CASE_NAME:-case}"
RESULTS_ROOT="${RESULTS_ROOT:-$ROOT/profiling/results/manual}"
CASE_DIR="$RESULTS_ROOT/$CASE_NAME"
mkdir -p "$CASE_DIR"

RANKS="${RANKS:-1}"
NODES="${NODES:-1}"
GPUS_PER_TASK="${GPUS_PER_TASK:-1}"
CPUS_PER_TASK="${CPUS_PER_TASK:-8}"
TASKS_PER_NODE="${TASKS_PER_NODE:-$(( (RANKS + NODES - 1) / NODES ))}"
EXE="${EXE:-$ROOT/build_gpu/main}"
INPUT="${INPUT:-$CASE_DIR/input.ini}"
LAUNCHER="${LAUNCHER:-mpirun}"
NX="${NX:-200}"
NY="${NY:-$NX}"
NZ="${NZ:-$NX}"
DIMS="${DIMS:-0 0 0}"
NSTEPS="${NSTEPS:-20}"
NITER="${NITER:-3}"
IBM_ENABLED="${IBM_ENABLED:-true}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD:-MANDATORY}"
export NVCOMPILER_ACC_NOTIFY="${NVCOMPILER_ACC_NOTIFY:-0}"

if [ "${REGENERATE_INPUT:-1}" = "1" ] || [ ! -f "$INPUT" ]; then
    NX="$NX" NY="$NY" NZ="$NZ" DIMS="$DIMS" NSTEPS="$NSTEPS" NITER="$NITER" \
    IBM_ENABLED="$IBM_ENABLED" profiling/scripts/make_input.sh "$INPUT"
fi

{
    echo "date: $(date -Is)"
    echo "host: $(hostname)"
    echo "case: $CASE_NAME"
    echo "nodes: $NODES"
    echo "ranks: $RANKS"
    echo "tasks_per_node: $TASKS_PER_NODE"
    echo "cpus_per_task: $CPUS_PER_TASK"
    echo "exe: $EXE"
    echo "input: $INPUT"
    echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"
    echo "OMP_TARGET_OFFLOAD: $OMP_TARGET_OFFLOAD"
    echo "NVCOMPILER_ACC_NOTIFY: $NVCOMPILER_ACC_NOTIFY"
    echo "launcher: $LAUNCHER"
    echo "MPI_EXTRA_ARGS: ${MPI_EXTRA_ARGS:-}"
    if type module >/dev/null 2>&1; then module list 2>&1; fi
} > "$CASE_DIR/env.txt"

if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi > "$CASE_DIR/nvidia-smi.txt" || true
fi

if [ ! -x "$EXE" ]; then
    echo "executable not found or not executable: $EXE" >&2
    exit 1
fi

if [ "$LAUNCHER" = "srun" ]; then
    cmd=(srun --nodes "$NODES" --ntasks "$RANKS" --ntasks-per-node "$TASKS_PER_NODE" \
         --cpus-per-task "$CPUS_PER_TASK" --gpus-per-task "$GPUS_PER_TASK" \
         --gpu-bind="${GPU_BIND:-single:1}" "$EXE" "$INPUT")
else
    read -r -a mpi_extra_args <<< "${MPI_EXTRA_ARGS:-}"
    cmd=(mpirun -np "$RANKS" --map-by "ppr:${TASKS_PER_NODE}:node" --bind-to none \
         "${mpi_extra_args[@]}" "$EXE" "$INPUT")
fi

printf '%q ' "${cmd[@]}" > "$CASE_DIR/command.txt"
printf '\n' >> "$CASE_DIR/command.txt"

set +e
/usr/bin/time -f 'elapsed_seconds %e\nmaxrss_kb %M' -o "$CASE_DIR/time.txt" \
    "${cmd[@]}" > "$CASE_DIR/stdout.txt" 2> "$CASE_DIR/stderr.txt"
rc=$?
set -e

timing_line="$(grep 'timing:' "$CASE_DIR/stdout.txt" | tail -n 1 || true)"
nsteps="$(printf '%s\n' "$timing_line" | awk '{for(i=1;i<=NF;i++) if($i=="nsteps") print $(i+1)}')"
loop_seconds="$(printf '%s\n' "$timing_line" | awk '{for(i=1;i<=NF;i++) if($i=="loop_seconds") print $(i+1)}')"
seconds_per_step="$(printf '%s\n' "$timing_line" | awk '{for(i=1;i<=NF;i++) if($i=="seconds_per_step") print $(i+1)}')"

printf 'case,nodes,ranks,tasks_per_node,nx,ny,nz,dims,nsteps,loop_seconds,seconds_per_step,exit_code\n' > "$CASE_DIR/summary.csv"
printf '%s,%s,%s,%s,%s,%s,%s,"%s",%s,%s,%s,%s\n' \
    "$CASE_NAME" "$NODES" "$RANKS" "$TASKS_PER_NODE" \
    "$NX" "$NY" "$NZ" "$DIMS" \
    "${nsteps:-}" "${loop_seconds:-}" "${seconds_per_step:-}" "$rc" >> "$CASE_DIR/summary.csv"

cat "$CASE_DIR/summary.csv"
exit "$rc"
