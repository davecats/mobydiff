#!/usr/bin/env bash
# A0 freestream physics gates (docs/next_session_airfoil.md, phase A0).
#
#   ./run_gates.sh              # all gates, sequentially (one job at a time)
#   ./run_gates.sh oblique      # one group: oblique | pois | vortex | ranks | config
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

run() { # bin ini ranks log
    local bin=$1 ini=$2 ranks=$3 log=$4
    echo "== $log (ranks $ranks) =="
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

echo
echo "freestream gates done (status $status)"
exit $status
