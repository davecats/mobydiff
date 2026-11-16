#!/bin/bash
#SBATCH --job-name=sd_time
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:4
#SBATCH --time=00:30:00
#SBATCH --partition=dev-gpu-h100
#SBATCH --account=hk-project-p0024333
#SBATCH --no-requeue
#
# Step 3b, the H100 step time of the SD7003 case (52 M cells, c/1024) before
# any production run: 300 steps of sd_seed.ini on 1, 2 and 4 H100 of one
# node, Chebyshev 12, then 150 steps at niter 24 and under red-black 12 on
# 4, each with the solver's profile. The binary is built by build_hk2.sh
# beforehand (a queued GPU job must not die on a compile error).
#
#     W2=$(ws_find junction); mkdir -p $W2/sd7003_run
#     (ship sd.case.h5 and sd_seed.ini into $W2/sd7003_run)
#     export CODE_DIR=$PROJECT/moby-junction RUN_DIR=$W2/sd7003_run
#     cd $RUN_DIR && sbatch --export=ALL $CODE_DIR/tutorials/sailplane/junction/horeka/time_sd7003_hk2.sh
set -uo pipefail
CODE_DIR="${CODE_DIR:?}"; RUN_DIR="${RUN_DIR:?}"
module purge
module load NVHPC/25.9-CUDA-12.9.1 HDF5/1.14.6-NVHPC-25.9-CUDA-12.9.1
export LD_LIBRARY_PATH="$EBROOTHDF5/lib:${LD_LIBRARY_PATH:-}"
export UCX_MEMTYPE_CACHE=n OMP_NUM_THREADS=1
BIN="$CODE_DIR/build_gpu/moby_solve"
[ -x "$BIN" ] || { echo "no $BIN: run build_hk2.sh first"; exit 1; }
cd "$RUN_DIR"
echo "job ${SLURM_JOB_ID:-none} $(date '+%F %T') $(hostname) code $(git -C "$CODE_DIR" rev-parse --short HEAD) gpu $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)" | tee -a timing.txt
run() {   # ranks, label, nsteps, pressure lines
    local n=$1 tag=$2 steps=$3 pres=$4 ini=time_$2_$1.ini
    sed -e "s/^nsteps = .*/nsteps = $steps/" -e 's/^field_interval = .*/field_interval = 0/' \
        -e "s/^field_prefix = .*/field_prefix = time_$tag/" -e "s/^runtime_file = .*/runtime_file = forces_time_$tag.txt/" \
        -e '/^accel = chebyshev$/d' -e "s/^niter = 12$/$pres/" sd_seed.ini | sed 's/|/\n/g' > "$ini"
    printf 'profile = true\n' >> "$ini"
    mpirun -n "$n" --bind-to none "$BIN" "$ini" > "time_${tag}_$n.log" 2>&1
    echo "$tag ranks $n rc $? $(grep -o 'seconds_per_step .*' "time_${tag}_$n.log" | head -1) projection $(grep 'step_timing: projection' "time_${tag}_$n.log" | grep -o 'seconds_per_step *[0-9.E+-]*')" | tee -a timing.txt
}
# (--bind-to none: "--map-by numa --bind-to core" could not place 4 ranks on
# the node, job 133191; ranks 1 and 2 of that job are in timing.txt.)
for r in ${RUNS:-"1 c12" "2 c12" "4 c12" "1 c24" "1 rb12" "4 c24" "4 rb12"}; do
    set -- $r
    case $2 in c12) p='accel = chebyshev|niter = 12';; c24) p='accel = chebyshev|niter = 24';; rb12) p='solver = redblack|sor = 1.5|niter = 12';; esac
    run $1 $2 $([ $2 = c12 ] && echo 300 || echo 150) "$p"
done
