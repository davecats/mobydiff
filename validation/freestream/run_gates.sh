#!/usr/bin/env bash
# A0 freestream physics gates (docs/next_session_airfoil.md, phase A0).
#
#   ./run_gates.sh              # all gates, sequentially (one job at a time)
#   ./run_gates.sh oblique      # one group: oblique | pois | vortex | ranks | config
#                               #          | low | refined | restart
#
# The last three are the OUTLET-FACE gates (docs/next_session_outlet.md):
#   low      the outlet on a LOW face: a mirrored pair with a non-parallel
#            outflow, and the Poiseuille and vortex cases run in -x (needs the
#            pois and vortex groups' outputs); each must reproduce its x_max
#            image to round-off
#   refined  uniform oblique flow recovered FROM REST through a refined patch
#            that touches nothing / the x_max outlet / the x_min outlet (the
#            exactness gates start from the uniform field and cannot see a
#            halo that is never written), plus 1 == 4 ranks on the second
#   restart  a run restarted at mid-length equals the continuous one
#            (run_restart_outlet.sh)
#
# Environment:  BIN=<solver>   (default ../../build_cpu/main)
#               BIN_GPU=<gpu>  (default ../../build_gpu/main; ranks group)
# Metrics: python3 check_freestream.py (invoked inline below).
set -o pipefail
cd "$(dirname "$0")"
set -u

BIN=${BIN:-../../build_cpu/main}
BIN_GPU=${BIN_GPU:-../../build_gpu/main}
CMP="python3 ../../tools/compare_fields.py"
sel=${1:-all}
status=0
# step 7: prepare the case file the solver will read, when it is missing
PIM="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/tools/prepare_if_missing.sh"

run() { # bin ini ranks log
    local bin=$1 ini=$2 ranks=$3 log=$4
    echo "== $log (ranks $ranks) =="
    "$PIM" "$ranks" "$bin" "$ini" "$log.prep.log" || { echo "   PREPARE FAILED -- see $log.prep.log"; status=1; return 1; }
    if ! mpirun -n "$ranks" "$bin" "$ini" > "$log.log" 2>&1; then
        echo "   FAILED — see $log.log"; status=1; return 1
    fi
    tail -n 2 "$log.log" | sed 's/^/   /'
}

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }

# --- gate (b): uniform oblique freestream preserved exactly ---
if want oblique; then
    rm -f oblique_*.h5
    run "$BIN" oblique.ini 1 oblique && \
    python3 check_freestream.py oblique "$(ls -t oblique_*.h5 | head -1)" \
        --u0 0.9396926207859084 --v0 0.3420201433256687 || status=1
fi

# --- gate (c): inflow/outflow Poiseuille vs the periodic reference ---
if want pois; then
    rm -f pois_ref_*.h5 pois_template_*.h5 pois_io_*.h5 pois_mid_*.h5 IC_pois.h5
    run "$BIN" pois_ref.ini 1 pois_ref && {
        cp "$(ls -t pois_ref_*.h5 | head -1)" pois_ref_final.h5
        sed -e 's|^file = IC_pois.h5|file =|' -e 's/^nsteps.*/nsteps = 1/' \
            -e 's/^field_prefix.*/field_prefix = pois_template/' pois_io.ini > .tmpl.ini
        run "$BIN" .tmpl.ini 1 pois_template && {
            python3 make_freestream_ics.py pois --src pois_ref_final.h5 \
                --template "$(ls -t pois_template_*.h5 | head -1)"
            sed -e 's/^nsteps.*/nsteps = 5000/' -e 's/^field_prefix.*/field_prefix = pois_mid/' \
                pois_io.ini > .mid.ini
            run "$BIN" .mid.ini 1 pois_mid
            run "$BIN" pois_io.ini 1 pois_io
            python3 check_freestream.py pois "$(ls -t pois_io_*.h5 | head -1)" \
                pois_ref_final.h5 --drift "$(ls -t pois_mid_*.h5 | head -1)" || status=1
        }
        rm -f .tmpl.ini .mid.ini
    }
fi

# --- gate (d): Lamb-Oseen vortex exits; report the reflected fraction ---
if want vortex; then
    rm -f vortex_template_*.h5 lamboseen_*.h5 IC_vortex.h5
    sed -e 's|^file = IC_vortex.h5|file =|' -e 's/^nsteps.*/nsteps = 1/' -e 's/^t_final.*/t_final = 0.0/' \
        -e 's/^field_prefix.*/field_prefix = vortex_template/' lamboseen.ini > .tmpl.ini
    run "$BIN" .tmpl.ini 1 vortex_template && {
        python3 make_freestream_ics.py vortex --template "$(ls -t vortex_template_*.h5 | head -1)"
        run "$BIN" lamboseen.ini 1 lamboseen && \
        python3 check_freestream.py vortex $(ls lamboseen_*.h5 | sort -t_ -k2 -n) || status=1
    }
    rm -f .tmpl.ini
fi

# --- gate (e): 1 == 4 ranks EXACT; CPU vs GPU ---
if want ranks; then
    for r in 1 4; do
        sed -e 's/^nsteps.*/nsteps = 200/' -e "s/^field_prefix.*/field_prefix = obl_r${r}/" \
            oblique.ini > .r.ini
        rm -f obl_r${r}_*.h5; run "$BIN" .r.ini $r obl_r$r; rm -f .r.ini
    done
    $CMP "$(ls -t obl_r1_*.h5 | head -1)" "$(ls -t obl_r4_*.h5 | head -1)" --tolerance 0 \
        && echo "oblique 1==4 ranks: EXACT" || { echo "oblique ranks MISMATCH"; status=1; }
    if [ -x "$BIN_GPU" ]; then
        sed -e 's/^nsteps.*/nsteps = 200/' -e 's/^field_prefix.*/field_prefix = obl_gpu/' \
            oblique.ini > .g.ini
        rm -f obl_gpu_*.h5; run "$BIN_GPU" .g.ini 1 obl_gpu; rm -f .g.ini
        $CMP "$(ls -t obl_r1_*.h5 | head -1)" "$(ls -t obl_gpu_*.h5 | head -1)" --tolerance 1e-12 \
            && echo "oblique CPU vs GPU: OK (<=1e-12)" || { echo "oblique CPU/GPU MISMATCH"; status=1; }
    fi
    # inflow/outflow variant (needs the pois group to have built IC_pois.h5)
    if [ -f IC_pois.h5 ]; then
        for r in 1 4; do
            sed -e 's/^nsteps.*/nsteps = 200/' -e "s/^field_prefix.*/field_prefix = pio_r${r}/" \
                pois_io.ini > .r.ini
            rm -f pio_r${r}_*.h5; run "$BIN" .r.ini $r pio_r$r; rm -f .r.ini
        done
        $CMP "$(ls -t pio_r1_*.h5 | head -1)" "$(ls -t pio_r4_*.h5 | head -1)" --tolerance 0 \
            && echo "pois_io 1==4 ranks: EXACT" || { echo "pois_io ranks MISMATCH"; status=1; }
        if [ -x "$BIN_GPU" ]; then
            sed -e 's/^nsteps.*/nsteps = 200/' -e 's/^field_prefix.*/field_prefix = pio_gpu/' \
                pois_io.ini > .g.ini
            rm -f pio_gpu_*.h5; run "$BIN_GPU" .g.ini 1 pio_gpu; rm -f .g.ini
            $CMP "$(ls -t pio_r1_*.h5 | head -1)" "$(ls -t pio_gpu_*.h5 | head -1)" --tolerance 1e-12 \
                && echo "pois_io CPU vs GPU: OK (<=1e-12)" || { echo "pois_io CPU/GPU MISMATCH"; status=1; }
        fi
    else
        echo "pois_io ranks leg SKIPPED (run the pois group first)"; status=1
    fi
fi

# --- gate (f): declared wall == inferred wall bit-exact; contradiction stops ---
if want config; then
    if [ ! -f IC_pois.h5 ]; then
        echo "config group SKIPPED (run the pois group first)"; status=1
    else
        for tag in decl infr; do
            sed -e 's/^nsteps.*/nsteps = 100/' -e "s/^field_prefix.*/field_prefix = cfg_${tag}/" \
                pois_io.ini > .c.ini
            [ $tag = infr ] && sed -i -e '/^y_min_patch/d' -e '/^y_max_patch/d' .c.ini
            rm -f cfg_${tag}_*.h5; run "$BIN" .c.ini 1 cfg_$tag; rm -f .c.ini
        done
        $CMP "$(ls -t cfg_decl_*.h5 | head -1)" "$(ls -t cfg_infr_*.h5 | head -1)" --tolerance 0 \
            && echo "declared wall == inferred: EXACT" || { echo "wall twin MISMATCH"; status=1; }
        for bad in "x_max_p_type = neumann" "x_max_u_type = neumann"; do
            { sed 's/^nsteps.*/nsteps = 1/' pois_io.ini; printf '\n[boundary]\n%s\n' "$bad"; } > .bad.ini
            if mpirun -n 1 "$BIN" .bad.ini > .bad.log 2>&1; then
                echo "CONTRADICTION NOT CAUGHT: $bad"; status=1
            else
                grep -q "contradicts the declared patch type" .bad.log \
                    && echo "contradiction '$bad': error-stops as required" \
                    || { echo "wrong error for: $bad"; tail -3 .bad.log; status=1; }
            fi
            rm -f .bad.ini .bad.log
        done
    fi
fi

# The x-mirror image of a freestream ini: inlet at x_max, outlet at x_min,
# flow in -x (u values and the initial u change sign; v is unchanged).
mirror_ini() {  # in out prefix
    sed -e 's/^x_min_patch = inlet/@IN@/' -e 's/^x_max_patch = outlet/x_min_patch = outlet/' \
        -e 's/^@IN@/x_max_patch = inlet/' \
        -e 's/^x_min_u_value *= *\(.*\)/x_max_u_value = -\1/' \
        -e 's/^x_min_u_profile/x_max_u_profile/' -e 's/^x_min_v_/x_max_v_/' \
        -e 's/^\(y_m[a-z]*_u_value\) *= *\(.*\)/\1 = -\2/' \
        -e 's/^initial_u *= *\(.*\)/initial_u = -\1/' \
        -e "s/^field_prefix.*/field_prefix = $3/" "$1" > "$2"
}

# --- outlet gate: the outlet on a LOW face ---
if want low; then
    rm -f outlet_high_*.h5 outlet_low_*.h5
    run "$BIN" outlet_high.ini 1 outlet_high && run "$BIN" outlet_low.ini 1 outlet_low && \
    python3 check_freestream.py mirror "$(ls -t outlet_high_*.h5 | head -1)" \
        "$(ls -t outlet_low_*.h5 | head -1)" || status=1

    if [ -f pois_ref_final.h5 ] && ls pois_io_*.h5 > /dev/null 2>&1; then
        rm -f pois_low_*.h5 pois_lowmid_*.h5 pois_lowtmpl_*.h5 IC_pois_low.h5
        mirror_ini pois_io.ini .pl.ini pois_low
        sed -i 's|^file = IC_pois.h5|file = IC_pois_low.h5|' .pl.ini
        sed -e 's|^file = IC_pois_low.h5|file =|' -e 's/^nsteps.*/nsteps = 1/' \
            -e 's/^field_prefix.*/field_prefix = pois_lowtmpl/' .pl.ini > .tmpl.ini
        run "$BIN" .tmpl.ini 1 pois_lowtmpl && {
            python3 make_freestream_ics.py pois --mirror --src pois_ref_final.h5 \
                --template "$(ls -t pois_lowtmpl_*.h5 | head -1)" --out IC_pois_low.h5
            sed -e 's/^nsteps.*/nsteps = 5000/' -e 's/^field_prefix.*/field_prefix = pois_lowmid/' \
                .pl.ini > .mid.ini
            run "$BIN" .mid.ini 1 pois_lowmid
            run "$BIN" .pl.ini 1 pois_low
            python3 check_freestream.py pois --mirror "$(ls -t pois_low_*.h5 | head -1)" \
                pois_ref_final.h5 --drift "$(ls -t pois_lowmid_*.h5 | head -1)" || status=1
            python3 check_freestream.py mirror "$(ls -t pois_io_*.h5 | head -1)" \
                "$(ls -t pois_low_*.h5 | head -1)" --tol 1e-10 || status=1
        }
        rm -f .pl.ini .tmpl.ini .mid.ini
    else
        echo "low: Poiseuille leg SKIPPED (run the pois group first)"; status=1
    fi

    if ls lamboseen_*.h5 > /dev/null 2>&1; then
        rm -f vortex_low_*.h5 vortex_lowtmpl_*.h5 IC_vortex_low.h5
        mirror_ini lamboseen.ini .vl.ini vortex_low
        sed -i 's|^file = IC_vortex.h5|file = IC_vortex_low.h5|' .vl.ini
        sed -e 's|^file = IC_vortex_low.h5|file =|' -e 's/^nsteps.*/nsteps = 1/' \
            -e 's/^t_final.*/t_final = 0.0/' -e 's/^field_prefix.*/field_prefix = vortex_lowtmpl/' \
            .vl.ini > .tmpl.ini
        run "$BIN" .tmpl.ini 1 vortex_lowtmpl && {
            python3 make_freestream_ics.py vortex --mirror \
                --template "$(ls -t vortex_lowtmpl_*.h5 | head -1)" --out IC_vortex_low.h5
            run "$BIN" .vl.ini 1 vortex_low && {
                python3 check_freestream.py vortex --mirror \
                    $(ls vortex_low_*.h5 | sort -t_ -k3 -n) || status=1
                python3 check_freestream.py mirror \
                    "$(ls lamboseen_*.h5 | sort -t_ -k2 -n | tail -1)" \
                    "$(ls vortex_low_*.h5 | sort -t_ -k3 -n | tail -1)" --tol 1e-9 || status=1
            }
        }
        rm -f .vl.ini .tmpl.ini
    else
        echo "low: vortex leg SKIPPED (run the vortex group first)"; status=1
    fi
fi

# --- outlet gate: uniform flow recovered from rest through a refined patch ---
if want refined; then
    U0=0.9396926207859084; V0=0.3420201433256687
    zero_ini() {  # out prefix "box"
        sed -e 's/^initial_u.*/initial_u = 0.0/' -e 's/^initial_v.*/initial_v = 0.0/' \
            -e 's/^nsteps.*/nsteps = 2000/' -e 's/^niter.*/niter = 12\naccel = chebyshev/' \
            -e 's/^dtmax.*/dtmax = 2.5e-3/' \
            -e "s/^nb = 8/nb = 8\nrefine = $3\nrefine_levels = 1/" \
            -e "s/^field_prefix.*/field_prefix = $2/" oblique.ini > "$1"
    }
    rm -f zr_inner_*.h5 zr_high_*.h5 zr_low_*.h5 zr_high1_*.h5 zr_high4_*.h5
    zero_ini .zi.ini zr_inner "0.25 0.75 0.25 0.75 0.0 1.0"
    zero_ini .zh.ini zr_high  "0.5 1.0 0.25 0.75 0.0 1.0"
    zero_ini .zt.ini zr_low   "0.0 0.5 0.25 0.75 0.0 1.0"
    mirror_ini .zt.ini .zl.ini zr_low
    run "$BIN" .zi.ini 1 zr_inner && \
        python3 check_freestream.py uniform "$(ls -t zr_inner_*.h5 | head -1)" --u0 $U0 --v0 $V0 || status=1
    run "$BIN" .zh.ini 1 zr_high && \
        python3 check_freestream.py uniform "$(ls -t zr_high_*.h5 | head -1)" --u0 $U0 --v0 $V0 || status=1
    run "$BIN" .zl.ini 1 zr_low && \
        python3 check_freestream.py uniform "$(ls -t zr_low_*.h5 | head -1)" --u0 -$U0 --v0 $V0 || status=1
    for r in 1 4; do
        sed -e 's/^nsteps.*/nsteps = 200/' -e "s/^field_prefix.*/field_prefix = zr_high${r}/" .zh.ini > .r.ini
        run "$BIN" .r.ini $r zr_high$r; rm -f .r.ini
    done
    $CMP "$(ls -t zr_high1_*.h5 | head -1)" "$(ls -t zr_high4_*.h5 | head -1)" --tolerance 0 \
        && echo "refined-on-outlet 1==4 ranks: EXACT" || { echo "refined-on-outlet ranks MISMATCH"; status=1; }
    rm -f .zi.ini .zh.ini .zt.ini .zl.ini
fi

# --- outlet gate: restart == continuous ---
if want restart; then
    SOLVER="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")" ./run_restart_outlet.sh || status=1
fi

echo
echo "freestream gates done (status $status)"
exit $status
