#!/bin/bash
#SBATCH --job-name=moby_interp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000
#SBATCH --time=03:00:00
#SBATCH --partition=cpuonly
#SBATCH --account=hk-project-p0024333
#
# The c/512 final field injected onto the c/2048 layout (interp_restart3d.py,
# python3 + h5py of the HoreKa image), 2026-10-29: about 1.5 M leaves, 42 GB
# written leaf by leaf (538 s for 372 k leaves on the workstation).
#     export RUN_DIR=$W/calib_run SRC=calib_a1.5_y3.5_34925.h5 SRC_CASE=calib_a1.5_y3.5.case.h5 \
#            TGT_CASE=calib_a1.5_y3.5_l7.case.h5 OUT=calib_a1.5_y3.5_l7_restart.h5
#     cd $RUN_DIR && sbatch --export=ALL submit_interp.sh
set -uo pipefail
cd "${RUN_DIR:?}"
echo "job ${SLURM_JOB_ID:-none} $(date '+%F %T') interp $SRC -> $OUT on $TGT_CASE" | tee -a provenance.txt
/usr/bin/time -f "interp wall %e s, max RSS %M kB" python3 interp_restart3d.py "$SRC" "$SRC_CASE" "$TGT_CASE" "$OUT" --dt 2e-5 > "${OUT%.h5}.log" 2>&1
rc=$?; tail -2 "${OUT%.h5}.log"
echo "exit $rc interp $OUT $(date '+%F %T')" | tee -a provenance.txt
