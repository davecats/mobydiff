#!/bin/bash
#SBATCH --job-name=moby_merge
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=19
#SBATCH --gres=gpu:4
#SBATCH --time=00:55:00
#SBATCH --partition=dev_accelerated
#SBATCH --account=hk-project-exasim
#SBATCH --no-requeue
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=davide.gatti@kit.edu
#
# The GPU half of the 2026-09-25 branch-consolidation gate (the merges of
# origin/boundaryLayer and origin/scalar). The CPU half ran on the login node
# and is recorded in the two merge commit messages; this covers what a login
# node cannot: the GPU path, and the two cases (beltrami_slaby, les_ibm) whose
# CPU runs are too slow to gate there.
#
# THE MERGE HAS TWO PARENTS AND THEREFORE TWO REFERENCES. Which one a case is
# gated against is not a detail -- it is the whole content of the gate:
#
#   REF_A = f0a8fe0, this branch's pre-merge head.
#     Everything that neither branch's physics touches must be max_abs 0
#     against it: min_channel (1 and 4 ranks), beltrami_slaby, les_ibm, and
#     Pass G's two production cases.
#
#   REF_B = be1451e, the old `scalar` branch tip. The branch itself was deleted
#     on 2026-09-25 once it was fully merged; the commit stays reachable from
#     main, so `git worktree add --detach <dir> be1451e` still builds it.
#     turb180 / wf180_y30 / lam30t must be max_abs 0 against THIS one and
#     WILL DIFFER from REF_A by O(1e-2). That is correct: `scalar` carries the
#     2026-08-05 fix for the cold-started RANS initial condition (k a factor 4
#     low on the last plane of every block, because init_rans_transport read a
#     halo that had not been filled yet). Agreeing with the branch that owns
#     the fix, to the last bit, is what proves the fix was imported rather
#     than mangled. A max_abs 0 against REF_A on these three would mean the
#     fix was LOST.
#
# PASS G RUNS WITH THE TRIP OFF (`trip_amp = 0`), deliberately. The production
# configs are boundary-layer cases with a trip, and the trip MATHS CHANGED with
# the boundaryLayer merge (Schlatter & Orlu coefficients -> the exact
# CaNS/SIMSON phase form), so a trip-on Pass G against REF_A is guaranteed to
# differ and would say nothing. Zeroing the amplitude makes the force
# identically zero on both sides and gates everything else at production scale
# -- 138 M and 60 M cells, blocks, the 2:1 interface, 4 ranks. The trip's own
# resolution (their maths, our block-list + spanwise-table optimisation
# re-applied on top) was gated separately and exactly on the CPU: max_abs 0 at
# 1 and 4 ranks against a binary carrying origin/boundaryLayer's bodyforce.f90
# verbatim, with trip_ts small enough that the random walk redrew ten times
# inside the run.
#
# FLAGS: -Mnofma -gpu=nofma on every side. The merge MOVES expressions (the
# trip's spanwise sums out of the per-cell loop), so the production-flag
# comparison would fail at 1e-15 with nothing wrong -- CLAUDE.md's rule.
#
# RESUMABLE: both drivers skip any run whose run.log exists.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
REF_A="${REF_A:?pre-merge worktree at f0a8fe0}"
REF_B="${REF_B:?worktree at be1451e, the old scalar branch tip}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/merge_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC"/run_exchange.sh "$SRC"/run_mapgate.sh "$SRC"/collect_exchange.py "$STG/"

# Pass G configs with the trip amplitude zeroed -- see the header.
cp -a "$SRC/../configs" "$STG/configs_notrip"
sed -i 's/^trip_amp *=.*/trip_amp = 0.0/' "$STG"/configs_notrip/*.ini

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

# be1451e predates compile.sh's nofma modes and carries its own script.
build_nofma() {   # build_nofma <dir>
    local d="$1"
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD)) -- gpu nofma"
    if grep -q 'gpu_nofma' "$d/compile.sh"; then
        ( cd "$d" && ./compile.sh gpu_nofma )
    else
        ( cd "$d" && ./compile.sh gpu && ./validation/scalar/compile_nofma.sh gpu )
    fi
}
for d in "$REF_A" "$REF_B" "$CODE_DIR"; do build_nofma "$d" || exit 1; done

NEW="$CODE_DIR/build_gpu_nofma/moby_solve"
RA="$REF_A/build_gpu_nofma/moby_solve"
RB="$REF_B/build_gpu_nofma/moby_solve"
for b in "$NEW" "$RA" "$RB"; do [ -x "$b" ] || { echo "ERROR: missing $b" >&2; exit 1; }; done

RES="${RESDIR:-$RUN_DIR/results_merge}"; mkdir -p "$RES"
{
    echo "job    : ${SLURM_JOB_ID:-none}"
    echo "date   : $(date '+%F %T %Z')"
    echo "nodes  : ${SLURM_JOB_NUM_NODES:-?}  (${SLURM_JOB_NODELIST:-?})"
    echo "new    : $(git -C "$CODE_DIR" rev-parse HEAD)  dirty $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
    echo "ref A  : $(git -C "$REF_A" rev-parse HEAD)   (pre-merge head)"
    echo "ref B  : $(git -C "$REF_B" rev-parse HEAD)   (the old scalar branch tip)"
    echo "flags  : -Mnofma -gpu=nofma, every side"
    echo "gpu    : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

echo "############ 1: vs the PRE-MERGE head -- must be max_abs 0 ############"
MOBY_ROOT="$CODE_DIR" NSTEPS=200 CASES="
min_channel:tutorials/min_channel/input_gpu.ini:4
min_channel1:tutorials/min_channel/input_gpu.ini:1
beltrami_slaby:validation/beltrami/slab_y.ini:1
les_ibm:validation/channel_interface/les_ibm/channel_ibm.ini:1
" bash "$STG/run_mapgate.sh" "$NEW" "$RA" "$RES/vs_premerge"

echo "############ 2: vs ORIGIN/SCALAR -- the RANS cases, must be max_abs 0 ############"
MOBY_ROOT="$CODE_DIR" NSTEPS=200 CASES="
turb180:validation/rans_sst/turb180.ini:1
wf180_y30:validation/rans_sst/wf180_y30.ini:1
lam30t:validation/rans_sst/lam30t.ini:1
" bash "$STG/run_mapgate.sh" "$NEW" "$RB" "$RES/vs_scalar"

echo "############ 3: the SAME three vs the pre-merge head -- MUST DIFFER ############"
# Informational and load-bearing in the opposite direction: max_abs 0 here
# would mean the scalar branch's RANS initial-condition fix was lost.
MOBY_ROOT="$CODE_DIR" NSTEPS=200 CASES="
turb180:validation/rans_sst/turb180.ini:1
wf180_y30:validation/rans_sst/wf180_y30.ini:1
lam30t:validation/rans_sst/lam30t.ini:1
" bash "$STG/run_mapgate.sh" "$NEW" "$RA" "$RES/ransfix_expected_to_differ"

echo "############ 4: Pass G, production scale, trip off -- must be max_abs 0 ############"
REF="$RA" H5MAXDIFF="$CODE_DIR/tools/h5maxdiff" CFG_DIR="$STG/configs_notrip" \
    PASSES="G" NSTEPS_GATE=20 bash "$STG/run_exchange.sh" "$NEW" "$RES/passG"

echo "=== verdict ==="
echo "-- must all be OK:"
grep -h -E "^(OK|FAIL)" "$RES"/vs_premerge/*.txt "$RES"/vs_scalar/*.txt "$RES"/passG/*.txt 2>/dev/null | sort | uniq -c
echo "-- must all be FAIL (the RANS IC fix is present):"
grep -h -E "^(OK|FAIL)" "$RES"/ransfix_expected_to_differ/*.txt 2>/dev/null | sort | uniq -c
echo "=== merge gate job finished ==="
