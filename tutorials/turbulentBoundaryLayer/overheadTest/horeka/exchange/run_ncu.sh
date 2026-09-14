#!/bin/bash
# Why does jacobi_apply cost 3x jacobi_compute_phi per cell?
#
#     run_ncu.sh <exe> <results_dir>
#
# The post-fix phase table (results_horeka_2026-09-14.md) puts `apply` at 33 % of
# the refined 16-rank step and `sweep` at 10 %, both called 18 times. The nsys
# sqlite already ruled out occupancy -- 88 / 94 registers, identical grid and
# block -- so the difference is memory traffic, and only ncu can say whether that
# is bytes the kernel genuinely needs or bytes it wastes on uncoalesced access.
#
# Three kernels: jacobi_apply's two (F1L470 `p += phi*idt`, F1L502 the face
# correction) and jacobi_compute_phi as the CONTROL -- the one kernel already
# known to run at its bandwidth limit, so the others are read against it and not
# against a hardware spec sheet.
#
# ncu REPLAYS each kernel several times, so it is run at 1 rank on a few launches
# only; nothing here is a timing, and no step time from this run may be quoted.
set -uo pipefail

EXE="${1:?usage: run_ncu.sh <exe> <results_dir>}"
RES="${2:?usage: run_ncu.sh <exe> <results_dir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="${CFG_DIR:-$HERE/../configs}"
NCU="${NCU:-/software/all/toolkit/nvidia_hpc_sdk/25.3/Linux_x86_64/25.3/profilers/Nsight_Compute/ncu}"
NSTEPS="${NSTEPS:-6}"
CONFIGS="${CONFIGS:-refined_yp82_rect_jacobi}"

# DRAM bytes actually moved, the sector/request ratio that exposes coalescing,
# and the two throughput percentages that say which limit is being hit.
METRICS="dram__bytes_read.sum,dram__bytes_write.sum,\
gpu__dram_throughput.avg.pct_of_peak_sustained_elapsed,\
sm__throughput.avg.pct_of_peak_sustained_elapsed,\
l1tex__average_t_sectors_per_request_pipe_lsu_mem_global_op_ld.ratio,\
l1tex__average_t_sectors_per_request_pipe_lsu_mem_global_op_st.ratio,\
lts__t_sector_hit_rate.pct,\
launch__occupancy_limit_registers,sm__warps_active.avg.pct_of_peak_sustained_active"

[ -x "$NCU" ] || { echo "ERROR: ncu not at $NCU" >&2; exit 1; }
"$NCU" --version | head -2
mkdir -p "$RES"

for cfg in $CONFIGS; do
    run="$RES/ncu_$cfg"
    [ -f "$run/ncu.csv" ] && { echo "--- skip $cfg (done)"; continue; }
    rm -rf "$run"; mkdir -p "$run"
    cp "$CFG/$cfg.ini" "$run/config.ini" || { echo "MISSING $cfg"; continue; }
    sed -i -e "s/^nsteps *=.*/nsteps = $NSTEPS/" \
           -e "s/^runtime_interval *=.*/runtime_interval = $NSTEPS/" "$run/config.ini"

    # --launch-skip past the first substage so the kernels profiled are steady
    # ones, --launch-count small because every launch is replayed.
    echo "=== ncu $cfg ($(date '+%F %T'))"
    ( cd "$run" && mpirun -n 1 "$NCU" --target-processes all --csv --page raw \
        --metrics "$METRICS" \
        --kernel-name 'regex:jacobi_apply|jacobi_compute_phi' \
        --launch-skip 20 --launch-count 12 \
        "$EXE" config.ini > ncu.csv 2> ncu.log )
    rc=$?
    if [ $rc -ne 0 ] || ! grep -q "jacobi" "$run/ncu.csv" 2>/dev/null; then
        echo "    FAILED (exit $rc) -- see $run/ncu.log"
        head -20 "$run/ncu.log" | sed 's/^/    /'
    else
        echo "    $(grep -c jacobi "$run/ncu.csv") metric rows"
    fi
    rm -f "$run"/overhead_*.h5
done
echo "=== ncu done $(date '+%F %T')"
