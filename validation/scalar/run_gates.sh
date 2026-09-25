#!/usr/bin/env bash
# Passive-scalar S1 gates (docs/next_session_scalar.md, increment S1).
#
#   ./run_gates.sh [group]      groups: uniform conserve conduction wave pr det
#
# Environment: BIN   solver binary   (default ../../build_cpu/moby_solve)
#              GBIN  GPU binary      (default ../../build_gpu/moby_solve, det only)
#              RANKS ranks           (default 1)
#
# The `det` group's CPU == GPU comparison must use the NOFMA builds
# (BIN=../../build_cpu_nofma/moby_solve GBIN=../../build_gpu_nofma/moby_solve):
# default FMA contraction differs between host and device and shows up as
# 1-2 ulp differences in every field, scalars included.
#
# Every group prints its measurement; PASS/FAIL lines come from
# check_scalar.py where the gate is a hard criterion.
set -uo pipefail
cd "$(dirname "$0")"

BIN=${BIN:-../../build_cpu/moby_solve}
GBIN=${GBIN:-../../build_gpu/moby_solve}
RANKS=${RANKS:-1}
sel=${1:-all}
status=0

run() {  # ini ranks logname [binary]
    local ini=$1 ranks=$2 log=$3 bin=${4:-$BIN}
    echo "== $log (ranks $ranks)"
    if ! mpirun -n "$ranks" "$bin" "$ini" > "$log.log" 2>&1; then
        echo "   RUN FAILED -- see $log.log"; status=1; return 1
    fi
    grep -E "seconds_per_step" "$log.log" | tail -1 | sed 's/^/   /'
}
want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
newest() { ls -t "$1"_*.h5 2>/dev/null | head -1; }

# --- (a) uniform scalar through the 3-level refined layout -----------------
if want uniform; then
    rm -f uniform3_*.h5
    run uniform3.ini "$RANKS" uniform3 && {
        python3 check_scalar.py uniform "$(newest uniform3)" \
            --scalar theta=2.5 --scalar phi=-1.25 --levels 3 | sed 's/^/   /' || status=1
        # the underlying uniform-FLOW gate must still hold
        python3 ../multilevel_body/check_uniform.py "$(newest uniform3)" \
            --u0 0.9396926207859084 --v0 0.3420201433256687 --w0 0.2 --levels 3 \
            | sed 's/^/   /' || status=1
    }
fi

# --- (b) global conservation in a periodic box -----------------------------
if want conserve; then
    rm -f conserve_*.h5 cons_ic.h5
    sed -e 's/^nsteps = 200/nsteps = 1/' -e 's/^field_interval = 200/field_interval = 1/' \
        -e 's/^field_prefix = conserve/field_prefix = cons_seed/' conserve.ini > .cons_seed.ini
    rm -f cons_seed_*.h5
    run .cons_seed.ini "$RANKS" cons_seed && {
        python3 make_scalar_ic.py "$(newest cons_seed)" cons_ic.h5 --wave 1 2 1 --amp 1.0
        sed -e 's|^field_prefix = conserve|field_prefix = conserve|' conserve.ini > .cons.ini
        printf '\n[restart]\nfile = cons_ic.h5\n' >> .cons.ini
        run .cons.ini "$RANKS" conserve && \
            python3 check_scalar.py conserve cons_ic.h5 "$(newest conserve)" \
                --scalar s1 --tolerance 1e-12 | sed 's/^/   /' || status=1
    }
    rm -f .cons_seed.ini .cons.ini
fi

# --- (c) conduction: exact linear profile + order 2 on a stretched grid -----
if want conduction; then
    # (c1) source-free: the LINEAR profile is an exact fixed point of the
    # discrete operator, so start ON it (init_profile = linear_y) and check it
    # does not drift. (Relaxing to it from s = 0 instead needs t >> L^2/D:
    # measured 2.4e-5 at t = 4, exactly the analytic slowest-mode residual
    # exp(-D pi^2 t/L^2) = 5e-5 -- convergence, not discretisation.)
    sed -e 's/^source = 1.0/source = 0.0/' -e 's/^y_min_value = 0.0/y_min_value = 1.0/' \
        -e 's/^initial = 0.0/initial = 0.0\ninit_profile = linear_y/' \
        -e 's/^t_final = 4.0/t_final = 0.5/' \
        -e 's/^field_prefix = conduction/field_prefix = cond_lin/' conduction.ini > .cond_lin.ini
    rm -f cond_lin_*.h5
    run .cond_lin.ini "$RANKS" cond_lin && \
        python3 check_scalar.py parabola "$(newest cond_lin)" --scalar s1 \
            --source 0.0 --diffusivity 1.0 --s0 1.0 --s1 0.0 | sed 's/^/   /'
    rm -f .cond_lin.ini
    # (c2) constant source: order study on the stretched grid
    for ny in 16 32 64; do
        # t_final = 12 leaves an analytic transient residual of exp(-D pi^2
        # t/L^2) = 1.4e-13, far below the truncation error under test.
        sed -e "s/^ny = 16/ny = $ny/" -e 's/^t_final = 4.0/t_final = 12.0/' \
            -e "s/^field_prefix = conduction/field_prefix = cond_$ny/" conduction.ini > ".cond$ny.ini"
        rm -f "cond_${ny}_"*.h5
        run ".cond$ny.ini" "$RANKS" "cond_$ny" && \
            python3 check_scalar.py parabola "$(newest "cond_$ny")" --scalar s1 \
                --source 1.0 --diffusivity 1.0 | sed 's/^/   /'
        rm -f ".cond$ny.ini"
    done
fi

# --- (d) advection-diffusion MMS: order 2 in space -------------------------
if want wave; then
    for n in 16 32 64; do
        sed -e "s/^nx = 16/nx = $n/" -e "s/^ny = 16/ny = $n/" -e "s/^nz = 16/nz = $n/" \
            -e "s/^field_prefix = wave/field_prefix = wv${n}_seed/" \
            -e 's/^nsteps = 500/nsteps = 1/' -e 's/^field_interval = 500/field_interval = 1/' \
            wave.ini > ".wv${n}s.ini"
        rm -f "wv${n}_seed_"*.h5 "wave${n}_"*.h5 "wv${n}_ic.h5"
        run ".wv${n}s.ini" "$RANKS" "wv${n}_seed" && {
            python3 make_scalar_ic.py "$(newest "wv${n}_seed")" "wv${n}_ic.h5" \
                --wave 1 1 1 --amp 1.0 --velocity 1.0 0.5 0.25 --zero-pressure
            sed -e "s/^nx = 16/nx = $n/" -e "s/^ny = 16/ny = $n/" -e "s/^nz = 16/nz = $n/" \
                -e "s/^field_prefix = wave/field_prefix = wave$n/" wave.ini > ".wv$n.ini"
            printf '\n[restart]\nfile = wv%s_ic.h5\n' "$n" >> ".wv$n.ini"
            run ".wv$n.ini" "$RANKS" "wave$n" && \
                python3 check_scalar.py wave "$(newest "wave$n")" --scalar s1 \
                    --wave 1 1 1 --amp 1.0 --velocity 1.0 0.5 0.25 --diffusivity 0.01 \
                    --ic "wv${n}_ic.h5" | sed 's/^/   /'
            rm -f ".wv$n.ini"
        }
        rm -f ".wv${n}s.ini"
    done
fi

# --- (e) Pr sweep in a laminar channel -------------------------------------
if want pr; then
    for pr in 0.1 1.0 10.0; do
        tag=pr${pr/./p}
        sed -e "s/^pr = 1.0/pr = $pr/" -e "s/^field_prefix = prsweep/field_prefix = $tag/" \
            prsweep.ini > ".$tag.ini"
        rm -f "${tag}_"*.h5
        run ".$tag.ini" "$RANKS" "$tag" && \
            python3 check_scalar.py profile "$(newest "$tag")" --scalar s1 \
                --diffusivity "$(python3 -c "print(1.0/(10.0*$pr))")" --tolerance 2e-2 \
                | sed 's/^/   /' || status=1
        rm -f ".$tag.ini"
    done
fi

# --- (f) determinism: 1 == 4 ranks, CPU == GPU, on a refined transport case -
if want det; then
    rm -f det_*.h5 det_ic.h5
    sed -e 's/^nsteps = 20/nsteps = 1/' -e 's/^field_interval = 20/field_interval = 1/' \
        -e 's/^field_prefix = det/field_prefix = det_seed/' det.ini > .det_seed.ini
    rm -f det_seed_*.h5
    run .det_seed.ini 1 det_seed && {
        python3 make_scalar_ic.py "$(newest det_seed)" det_ic.h5 --wave 1 2 1 --amp 1.0
        for tag in r1 r4 gpu; do
            sed -e "s/^field_prefix = det/field_prefix = det_$tag/" det.ini > ".det_$tag.ini"
            printf '\n[restart]\nfile = det_ic.h5\n' >> ".det_$tag.ini"
            rm -f "det_${tag}_"*.h5
            case $tag in
                r1)  run ".det_$tag.ini" 1 "det_$tag" ;;
                r4)  run ".det_$tag.ini" 4 "det_$tag" ;;
                gpu) [ -x "$GBIN" ] && run ".det_$tag.ini" 1 "det_$tag" "$GBIN" \
                         || echo "   (no GPU binary, skipping)" ;;
            esac
            rm -f ".det_$tag.ini"
        done
        echo "   1 rank vs 4 ranks:"
        python3 ../../tools/compare_fields.py "$(newest det_r1)" "$(newest det_r4)" \
            un vn wn pn s1 --tolerance 0 | sed 's/^/     /' || status=1
        if [ -f "$(newest det_gpu)" ]; then
            echo "   CPU vs GPU:"
            python3 ../../tools/compare_fields.py "$(newest det_r1)" "$(newest det_gpu)" \
                un vn wn pn s1 --tolerance 0 | sed 's/^/     /' || status=1
        fi
        echo "   interface conservation (INFORMATIVE, not a pass/fail: the flux"
        echo "   form telescopes exactly only across MATCHING faces, so a 2:1"
        echo "   interface leaks at its transfer error):"
        python3 check_scalar.py conserve det_ic.h5 "$(newest det_r1)" --scalar s1 \
            --tolerance 1e30 | sed 's/^/     /'
    }
    rm -f .det_seed.ini
fi

# --- (g) S0 restart round-trip: present vs absent scalar dataset -----------
if want restart; then
    rm -f smoke_*.h5 rst_*.h5
    run smoke.ini 1 smoke && {
        src=$(newest smoke)
        # (g1) the dataset IS in the file: it must win over [scalar.N] initial
        sed -e 's/^initial = 3.0/initial = 99.0/' -e 's/^field_prefix = smoke/field_prefix = rst_hit/' \
            smoke.ini > .rst_hit.ini
        printf '\n[restart]\nfile = %s\n' "$src" >> .rst_hit.ini
        run .rst_hit.ini 1 rst_hit
        # (g2) the dataset is ABSENT (renamed scalar): warn + reinitialise
        sed -e 's/^initial = 3.0/initial = 99.0/' -e 's/^name = theta/name = theta_new/' \
            -e 's/^field_prefix = smoke/field_prefix = rst_miss/' smoke.ini > .rst_miss.ini
        printf '\n[restart]\nfile = %s\n' "$src" >> .rst_miss.ini
        run .rst_miss.ini 1 rst_miss
        python3 - "$src" "$(newest rst_hit)" "$(newest rst_miss)" <<'PY' | sed 's/^/   /'
import sys, h5py, numpy as np
src, hit, miss = sys.argv[1:4]
s = h5py.File(src)["theta"][...]
h = h5py.File(hit)["theta"][...]
m = h5py.File(miss)["theta_new"][...]
print(f"restart WITH the dataset:    max|s - file| = {np.max(np.abs(h - s)):.3e}"
      f"  (initial = 99 in the ini, value {h.flat[0]})")
print(f"restart WITHOUT the dataset: value = {m.flat[0]} (expected the ini's 99.0)")
ok = np.max(np.abs(h - s)) == 0.0 and np.all(m == 99.0)
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY
        [ ${PIPESTATUS[0]} -eq 0 ] || status=1
        grep -i "warning: restart file has no dataset" rst_miss.log | sed 's/^/   /'
        rm -f .rst_hit.ini .rst_miss.ini
    }
fi

echo
[ $status -eq 0 ] && echo "scalar S1 gates: no failures" || echo "scalar S1 gates: FAILURES"
exit $status
