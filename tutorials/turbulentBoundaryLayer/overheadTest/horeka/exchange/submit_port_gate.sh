#!/bin/bash
#SBATCH --job-name=moby_port
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
# The GPU gate for the jacobi-interface port (branch
# port/jacobi-interface-features). Everything in that branch was gated on CPU
# only; this is the GPU half, plus les_ibm, which is too slow to run on a login
# node and so was never gated at all.
#
# REF = 68e8f16, `main` before the port. See run_portgate.sh for why the
# comparison injects `[flow] convection = skew` on the ref side rather than
# running a bare A/B: the S3 lockdown is a physics change, and what must be
# bit-exact is the lockdown against the old binary running the configuration it
# makes unconditional.
#
# Pass 2 is the feature-positive half -- the checks that say each ported knob
# actually does something on a GPU, not merely that disabling it changes
# nothing. A gate that only ever shows "no difference" cannot distinguish a
# working feature from a dead one.
#
# nofma on BOTH sides: the port MOVES expressions (the skew corrections left
# their `if` block), so a production-flag comparison would fail at 1e-15 with
# nothing wrong -- CLAUDE.md's rule.
#
# RESUMABLE: run_portgate.sh skips any leg whose run.log exists.
set -uo pipefail

CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"; REF_DIR="${REF_DIR:?}"
SRC="$CODE_DIR/tutorials/turbulentBoundaryLayer/overheadTest/horeka/exchange"
STG="$RUN_DIR/port_staged"; rm -rf "$STG"; mkdir -p "$STG"
cp "$SRC"/run_portgate.sh "$STG/"

module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT="${HDF5_ROOT:-$HOME/hdf5}"
export LD_LIBRARY_PATH="$HDF5_ROOT/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1

for d in "$REF_DIR" "$CODE_DIR"; do
    echo "=== building $d ($(git -C "$d" rev-parse --short HEAD)) -- gpu_nofma"
    ( cd "$d" && ./compile.sh gpu_nofma ) || exit 1
done
NEW="$CODE_DIR/build_gpu_nofma/moby_solve"
REF="$REF_DIR/build_gpu_nofma/moby_solve"
for b in "$NEW" "$REF"; do [ -x "$b" ] || { echo "ERROR: missing $b" >&2; exit 1; }; done

RES="${RESDIR:-$RUN_DIR/results_port}"; mkdir -p "$RES"
{
    echo "job   : ${SLURM_JOB_ID:-none}"
    echo "date  : $(date '+%F %T %Z')"
    echo "node  : ${SLURM_JOB_NODELIST:-?}"
    echo "new   : $(git -C "$CODE_DIR" rev-parse HEAD)  dirty $(git -C "$CODE_DIR" status --porcelain -uno | wc -l)"
    echo "ref   : $(git -C "$REF_DIR" rev-parse HEAD)   (main before the port)"
    echo "flags : -Mnofma -gpu=nofma, both sides"
    echo "gpu   : $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
} | tee "$RES/provenance.txt"

echo "############ PASS 1: equivalence -- every leg must be max_abs 0 ############"
MOBY_ROOT="$CODE_DIR" H5MAXDIFF="$CODE_DIR/tools/h5maxdiff" NSTEPS="${NSTEPS:-200}" \
    bash "$STG/run_portgate.sh" "$NEW" "$REF" "$RES/equiv"
equiv_rc=$?

echo "############ PASS 2: the ported knobs must DO something on a GPU ############"
P2="$RES/positive"; mkdir -p "$P2"
run1() {  # run1 <tag> <ini> <extra-fragment>
    local tag="$1" ini="$2" extra="$3" d="$P2/$1"
    rm -rf "$d"; cp -a "$(dirname "$CODE_DIR/$ini")" "$d"
    sed -e 's/^nsteps *=.*/nsteps = 100/' -e 's/^t_final *=.*/t_final = 0.0/' \
        "$CODE_DIR/$ini" > "$d/cfg.ini"
    printf '%s' "$extra" >> "$d/cfg.ini"
    ( cd "$d" && mpirun -n 1 "$NEW" cfg.ini > run.log 2>&1 ) || echo "  $tag FAILED"
}
# kpin: pinning the whole domain must make the answer independent of tu --
# k is then exactly 0, so nothing downstream can depend on it.
run1 kpin_tu1  validation/rans_sst/turb180.ini $'\n[rans]\nkpin_box = -1e9 1e9 -1e9 1e9 -1e9 1e9\ntu = 1.0\n'
run1 kpin_tu10 validation/rans_sst/turb180.ini $'\n[rans]\nkpin_box = -1e9 1e9 -1e9 1e9 -1e9 1e9\ntu = 10.0\n'
run1 notu1     validation/rans_sst/turb180.ini $'\n[rans]\ntu = 1.0\n'
run1 notu10    validation/rans_sst/turb180.ini $'\n[rans]\ntu = 10.0\n'
echo "-- kpin, pinned everywhere, tu 1 vs 10 (expect max_abs 0 on u,k,nut):"
"$CODE_DIR/tools/h5maxdiff" "$P2/kpin_tu1/turb180_100.h5" "$P2/kpin_tu10/turb180_100.h5" un k nut | tail -1
echo "-- CONTROL, no pin, same tu change (expect a DIFFERENCE, or the test above is empty):"
"$CODE_DIR/tools/h5maxdiff" "$P2/notu1/turb180_100.h5" "$P2/notu10/turb180_100.h5" k | tail -1

# scalar convection: the three forms must be distinct, and f = 1/2 must land
# between f = 0 and f = 1.
for m in divergence skew advective; do
    run1 "sc_$m" validation/scalar/ibmwavy.ini "$(printf '\n[scalar]\nconvection = %s\n' "$m")"
done
echo "-- scalar convection, pairwise (expect all three DIFFERENT, skew at the midpoint):"
for pair in "divergence skew" "divergence advective" "skew advective"; do
    set -- $pair
    printf "   %-22s %s\n" "$1 vs $2" \
        "$("$CODE_DIR/tools/h5maxdiff" "$P2/sc_$1/ibmwavy_100.h5" "$P2/sc_$2/ibmwavy_100.h5" theta | tail -1)"
done

# CV forces: the empty domain must give identically zero.
run1 cv_empty validation/cylinder/empty.ini ''
echo "-- CV forces on an empty domain (expect 0 nonzero entries):"
awk 'NR>1 {for(i=3;i<=NF;i++) if ($i+0 != 0) bad++} END {print "   nonzero:", bad+0, "of", (NR-1)*3}' \
    "$P2/cv_empty/forces_empty.txt" 2>/dev/null || echo "   NO FORCE TRACE"

echo "=== verdict ==="
grep -hE "^  " "$RES/equiv/summary.txt" 2>/dev/null | sed 's/^/ /'
echo "pass 1 exit: $equiv_rc  (0 = every equivalence leg OK)"
echo "=== port gate job finished $(date '+%F %T') ==="
