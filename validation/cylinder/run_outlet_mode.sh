#!/usr/bin/env bash
# The outlet pressure mode on the Re 100 cylinder (docs/next_session_after_step9.md
# item 3; README.md, last sections): from a clean pressure, ONE time unit at
# Chebyshev niter 60 takes the stored p rms from 0.30 to 0.9-1.6 and the
# control-volume lift to +-3.9. This driver repeats that leg with different
# PROJECTIONS, everything else fixed, to find what the mode depends on.
#
#   NEW=/abs/build_gpu/moby_solve [OUT=outlet] ./run_outlet_mode.sh [leg ...]
#
# A leg is <name>:<k>[:<T>], k the step size dt = 5e-3 / 2^k (200 << k steps
# per time unit), T the number of time units (default 1; a snapshot is
# written every time unit) and <name> one of the projection settings below.
# Default: every leg at k = 0 and k = 3. outlet_mode.py reads what is there.
#
#   cheb60   Chebyshev-Jacobi, niter 60      (the leg of the finding)
#   cheb240  Chebyshev-Jacobi, niter 240
#   jac60    plain damped Jacobi (sor 0.8), niter 60
#   jac240   plain damped Jacobi, niter 240
#   rb60     red-black SOR (sor 1.5), niter 60
#   rb240    red-black SOR, niter 240
#
# The clean state is the one of run_f7_gates.sh: the Re 100 shedding restart
# with pn zeroed, 300 steps at Chebyshev niter 60 (f7/cp100_*.h5; made here
# when missing). Inputs not in git: the restart (SIB) and ./setup.sh's case
# file. One GPU, one rank.
set -uo pipefail
cd "$(dirname "$0")"
NEW=${NEW:?absolute path of the solver binary}
SIB=${SIB:-$HOME/Codes/mobydiff.scalar/validation/cylinder}
[ -f ibm_coeff_re100.h5 ] || { echo "missing ibm_coeff_re100.h5 -- run ./setup.sh"; exit 1; }
OUT=${OUT:-outlet}
mkdir -p "$OUT" && cd "$OUT"
ln -sf ../ibm_coeff_re100.h5 .

# derive <out.ini> from ../cyl_re100.ini: restart file, key=value overrides
# (a key the ini does not have is added to [pressure])
derive() {
    local out=$1 restart=$2; shift 2
    cp ../cyl_re100.ini "$out"
    for kv in "$@"; do
        local k=${kv%%=*} v=${kv#*=}
        if grep -qE "^$k *=" "$out"; then
            sed -i -E "s|^$k *=.*|$k = $v|" "$out"
        else
            sed -i -E "s|^\[pressure\]|[pressure]\n$k = $v|" "$out"
        fi
    done
    printf '\n[restart]\nfile = %s\n' "$restart" >> "$out"
}
run() { mpirun -n 1 "$NEW" "$1" > "$2" 2>&1 || { echo "RUN FAILED: $2"; tail -5 "$2"; exit 1; }; }

ic=$(ls -t ../f7/cp100_[0-9]*.h5 2>/dev/null | head -1)
if [ -z "$ic" ]; then
    echo "== clean-p state: zero pn, 300 steps at Chebyshev niter 60"
    python3 - "$SIB/cyl_re100_40001.h5" <<'PY'
import sys, shutil, h5py
shutil.copy(sys.argv[1], "cp100_zero.h5")
with h5py.File("cp100_zero.h5", "r+") as f:
    f["pn"][...] = 0.0
PY
    derive cp100.ini cp100_zero.h5 nsteps=300 t_final=0.0 niter=60 field_interval=300 \
        field_prefix=cp100 runtime_file=forces_cp100.txt
    run cp100.ini cp100.log
    ic=$(ls -t cp100_[0-9]*.h5 | head -1)
fi
echo "clean state: $ic"

settings() {  # the [pressure] overrides of a leg name
    case $1 in
        cheb60)  echo "niter=60 accel=chebyshev" ;;
        cheb240) echo "niter=240 accel=chebyshev" ;;
        jac60)   echo "niter=60 accel=none sor=0.8" ;;
        jac240)  echo "niter=240 accel=none sor=0.8" ;;
        rb60)    echo "niter=60 accel=none solver=redblack sor=1.5" ;;
        rb240)   echo "niter=240 accel=none solver=redblack sor=1.5" ;;
        *) echo "unknown leg $1" >&2; exit 1 ;;
    esac
}

legs=("$@")
[ ${#legs[@]} -eq 0 ] && legs=(cheb60:0 jac60:0 rb60:0 cheb240:0 jac240:0 rb240:0 \
                               cheb60:3 jac60:3 rb60:3 cheb240:3 jac240:3 rb240:3)
for leg in "${legs[@]}"; do
    IFS=: read -r name k T <<< "$leg"; T=${T:-1}
    per=$(( 200 << k )); n=$(( per*T )); dt=$(python3 -c "print(5.0e-3/2**$k)")
    pfx=om_${name}_$k; [ "$T" != 1 ] && pfx=om_${name}T${T}_$k
    # shellcheck disable=SC2046
    derive $pfx.ini "$ic" dt=$dt dtmax=$dt nsteps=$n t_final=0.0 field_interval=$per \
        field_prefix=$pfx force_sample_interval=$(( 10 << k )) runtime_file=forces_$pfx.txt \
        $(settings "$name")
    rm -f ${pfx}_[0-9]*.h5
    t0=$(date +%s)
    run $pfx.ini $pfx.log
    echo "   $pfx  dt = $dt  ($n steps)  $(( $(date +%s) - t0 )) s"
done
python3 ../outlet_mode.py .
