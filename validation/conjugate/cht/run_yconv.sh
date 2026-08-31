#!/usr/bin/env bash
# Wall-normal grid convergence for the conjugate channel: dy+ = 1.5, 1.0, 0.75
# at fixed dx+ = 14.1, dz+ = 5.1 and fixed geometry (interfaces at 0.6 / 2.6).
#
#   ./run_yconv.sh 1.5 [1.0 0.75 ...]
#
# Environment: BIN, PREP, TDEV (develop time, default 5), TSTAT (default 25).
#
# The point is the PEAK temperature variance, which is a fluid-side quantity:
# one scalar (k1 = Flageul's own conjugate case) is enough, and a single scalar
# roughly halves the step cost relative to the six-scalar campaign.
set -uo pipefail
# istmcorax has no modulefile: the toolchain must be put on PATH explicitly,
# exactly as run_corax.sh does. Without it the CPU moby_prepare cannot even
# MPI_Init (it is linked against nvhpc's MPI).
export PATH=/opt/Nvidia/nvhpc/Linux_x86_64/25.9/compilers/bin:/opt/Nvidia/nvhpc/Linux_x86_64/25.9/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH
cd "$(dirname "$0")"
ROOT=$(cd ../../.. && pwd)
BIN=${BIN:-$ROOT/build_gpu_corax/moby_solve}
PREP=${PREP:-$ROOT/build_cpu/moby_prepare}
PY=${PY:-python3}
TDEV=${TDEV:-5}
TSTAT=${TSTAT:-25}
YLO=0.6; YHI=2.6; NX=160; NZ=224; NB=32; LY=3.2
LOG=$PWD/yconv.log

for dyp in "$@"; do
    # dy = dyp/Re_tau; the fluid gap is 2 and each slab 0.6, so ny follows
    read NY DT <<<"$($PY -c "
dyp=$dyp; dy=dyp/180.0
ny=int(round(3.2/dy)); dt=1.44e-3*(dyp/1.5)**2
print(ny, '%.6e'%dt)")"
    tag="yc$(echo $dyp | tr . p)"
    ndev=$($PY -c "print(int(round($TDEV/$DT)))")
    nstat=$($PY -c "print(int(round($TSTAT/$DT)))")
    echo "[yconv] === dy+ = $dyp  ny = $NY  dt ~ $DT  develop $ndev  stats $nstat  $(date +%H:%M:%S)" >> $LOG

    gen() { sed -e "s|@NX@|$NX|" -e "s|@NY@|$NY|" -e "s|@NZ@|$NZ|" -e "s|@NB@|$NB|" \
                -e "s|@CASE@|$tag.h5|" -e "s|@PREFIX@|$1|" -e "s|@NSTEPS@|$2|" \
                -e "s|@WRITE@|$3|" -e "s|@RESTART@|$4|" -e "s|@SAMPLE@|$5|" \
                -e "s|@STATS@|$6|" -e "s|@STATSFILE@|$7|" yconv.ini > ".$1.ini"; }

    # geometry (the STLs are cht180's -- identical walls)
    gen "${tag}_p" 1 1 - 0 0 x
    sed -e '/^coeff_file/d' -e '/^\[restart\]/,$d' \
        -e 's|^\[ibm\]|[ibm]\nstl_file = wall_lo.stl\nstl_file = wall_hi.stl|' \
        ".${tag}_p.ini" > ".${tag}_prep.ini"
    [ -f "$tag.h5" ] || mpirun -n 8 "$PREP" ".${tag}_prep.ini" "$tag.h5" >> $LOG 2>&1 || {
        echo "[yconv] prepare FAILED" >> $LOG; continue; }

    # mint a template, map KMM180 into the gap, seed the scalar
    gen "${tag}_m" 1 1 - 0 0 x
    sed '/^\[restart\]/,$d' ".${tag}_m.ini" > ".${tag}_mint.ini"
    mpirun -n 1 "$BIN" ".${tag}_mint.ini" >> $LOG 2>&1
    [ -f "${tag}_m_1.h5" ] || { echo "[yconv] MINT FAILED (no ${tag}_m_1.h5)" >> $LOG; continue; }
    # Check the IC step's own exit status. Not doing so once cost a whole
    # study run: a broken make_cht_ic.py left no IC file and the failure was
    # reported three steps later as "develop FAILED".
    $PY ./make_cht_ic.py "${tag}_m_1.h5" "${tag}_ic.h5" \
        --interfaces $YLO $YHI --scalars k1 >> $LOG 2>&1 \
        || { echo "[yconv] IC FAILED" >> $LOG; continue; }

    gen "${tag}_d" "$ndev" "$ndev" "${tag}_ic.h5" 0 0 x
    mpirun -n 1 "$BIN" ".${tag}_d.ini" >> $LOG 2>&1 || { echo "[yconv] develop FAILED" >> $LOG; continue; }

    gen "${tag}_s" "$nstat" "$nstat" "${tag}_d_${ndev}.h5" 25 2000 "${tag}_stats.h5"
    # per-leg statistics files: the channel accumulator CONTINUES an existing
    # file and cross-checks its grid, so a shared name makes every leg whose
    # ny differs from the first one die at the stats phase.
    rm -f "${tag}_stats.h5" "${tag}_s_vel.h5"
    mpirun -n 1 "$BIN" ".${tag}_s.ini" >> $LOG 2>&1 || { echo "[yconv] stats FAILED" >> $LOG; continue; }
    grep -E "seconds_per_step" "${tag}_s.log" 2>/dev/null | tail -1 >> $LOG
    echo "[yconv] dy+ = $dyp done $(date +%H:%M:%S)" >> $LOG
    $PY ./check_cht.py thermal "${tag}_stats.h5" --interfaces $YLO $YHI >> $LOG 2>&1
done
echo "[yconv] ALL DONE $(date +%H:%M:%S)" >> $LOG
