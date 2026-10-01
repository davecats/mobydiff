#!/usr/bin/env bash
# Numerics review step 9 (F7, the exact penalization factor): the two
# cylinder measurements, each run with TWO binaries from the same state.
#
#   NEW=/abs/build_gpu/moby_solve REF=/abs/ref/build_gpu/moby_solve ./run_f7_gates.sh [steady|dt|all]
#
#   steady  Re 40 from the clean-p converged state (cvpz_20301.h5), 2000
#           steps at the production settings. Both penalization factors have
#           the SAME steady fixed point (lambda q = R), so the control-volume
#           C_D must agree between the binaries to the convergence tolerance.
#   dt      Re 100 shedding state, pressure cleaned first (README "Clean-p
#           protocol": zero pn, 300 steps at niter 60), then ONE time unit at
#           fixed dt = 5e-3 / 2^k, k = 0..3, at niter 60 so the projection
#           error does not mask the time error. dt_study.py reads the final
#           fields: error against the finest NEW run, in the cut cells and in
#           the fluid, and the observed order.
#
# Inputs that are NOT in git: the restarts, in a sibling checkout (SIB), and
# the case files from ./setup.sh. One GPU, one rank; ~25 min for `all`.
set -uo pipefail
cd "$(dirname "$0")"
NEW=${NEW:?absolute path of the binary under test}
REF=${REF:?absolute path of the reference binary}
SIB=${SIB:-$HOME/Codes/mobydiff.scalar/validation/cylinder}
what=${1:-all}
mkdir -p f7 && cd f7
for f in ibm_coeff_re40.h5 ibm_coeff_re100.h5; do
    [ -f "../$f" ] || { echo "missing ../$f -- run ./setup.sh"; exit 1; }
    ln -sf "../$f" .
done

# derive <out.ini> from ../<base.ini>: restart file, then key=value overrides
derive() {
    local base=$1 out=$2 restart=$3; shift 3
    cp "../$base" "$out"
    for kv in "$@"; do
        local k=${kv%%=*} v=${kv#*=}
        grep -qE "^$k *=" "$out" || { echo "derive: no key $k in $base"; exit 1; }
        sed -i -E "s|^$k *=.*|$k = $v|" "$out"
    done
    printf '\n[restart]\nfile = %s\n' "$restart" >> "$out"
}
run() {  # run <binary> <ini> <log>
    mpirun -n 1 "$1" "$2" > "$3" 2>&1 || { echo "RUN FAILED: $3"; tail -5 "$3"; exit 1; }
}

if [ "$what" = steady ] || [ "$what" = all ]; then
    echo "== steady: Re 40, 2000 steps from the clean-p state, both binaries"
    for side in ref new; do
        bin=$REF; [ $side = new ] && bin=$NEW
        derive cyl_re40.ini st_$side.ini "$SIB/cvpz_20301.h5" \
            nsteps=2000 t_final=0.0 field_prefix=st_$side runtime_file=forces_st_$side.txt
        run "$bin" st_$side.ini st_$side.log
        python3 - "$side" <<'PY'
import sys, numpy as np
d = np.loadtxt(f"forces_st_{sys.argv[1]}.txt", skiprows=1)
tail = d[len(d)//2:]
print(f"   {sys.argv[1]:3s}  C_D = {tail[:,3].mean():.8f} +- {tail[:,3].std():.1e}   C_L = {tail[:,2].mean():+.2e}   (last {len(tail)} samples)")
PY
    done
fi

if [ "$what" = dt ] || [ "$what" = all ]; then
    echo "== dt: clean-p Re 100 state (zero pn, 300 steps at niter 60, reference binary)"
    python3 - "$SIB/cyl_re100_40001.h5" <<'PY'
import sys, shutil, h5py
shutil.copy(sys.argv[1], "cp100_zero.h5")
with h5py.File("cp100_zero.h5", "r+") as f:
    f["pn"][...] = 0.0
PY
    derive cyl_re100.ini cp100.ini cp100_zero.h5 \
        nsteps=300 t_final=0.0 niter=60 field_interval=300 field_prefix=cp100 runtime_file=forces_cp100.txt
    run "$REF" cp100.ini cp100.log
    ic=$(ls -t cp100_*.h5 | grep -v zero | head -1)
    echo "   clean state: $ic"
    for side in ref new; do
        bin=$REF; [ $side = new ] && bin=$NEW
        for k in 0 1 2 3; do
            n=$(( 200 << k ))
            dt=$(python3 -c "print(5.0e-3/2**$k)")
            derive cyl_re100.ini dt_${side}_$k.ini "$ic" \
                dt=$dt dtmax=$dt nsteps=$n t_final=0.0 niter=60 field_interval=$n \
                field_prefix=dt_${side}_$k force_sample_interval=$(( 10 << k )) \
                runtime_file=forces_dt_${side}_$k.txt
            rm -f dt_${side}_${k}_*.h5
            run "$bin" dt_${side}_$k.ini dt_${side}_$k.log
            echo "   $side  dt = $dt  ($n steps)  done"
        done
    done
    python3 ../dt_study.py .
fi
