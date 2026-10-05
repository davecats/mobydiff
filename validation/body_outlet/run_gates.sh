#!/usr/bin/env bash
# Body-at-outlet gates (docs/next_session_body_at_outlet.md): an immersed
# body crossing an inlet and an outlet plane. Poiseuille flow between two
# immersed slabs (make_slabs.py) in an inflow/outflow channel.
#
#   SOLVER=/abs/build_cpu/moby_solve [RANKS=1] ./run_gates.sh [group]
#   groups: profile | mirror | permute | limits | ranks | restart | gpu | twin | all
#
#   profile  pois_body (inlet u = 1 through the slabs, outlet at x_max) and its
#            periodic twin, 20000 steps (t = 40): the outlet-face profile equals
#            the twin's, p linear and pinned, solid faces at round-off, mass
#            balance, the inlet rows inside the slabs zeroed (B2)
#   mirror   the same with the outlet on the LOW face: x -> lx - x image
#   permute  the stream along y and along z: the x case under permutation
#   limits   the body ending one cell before the plane / starting at it
#   ranks    pois_body on RANKS ranks == 1 rank, max_abs 0 (needs RANKS > 1)
#   restart  a run stopped at mid-length and restarted == the continuous run
#   gpu      GPU_SOLVER=/abs/build_gpu_nofma/moby_solve on 1 rank == the CPU run
#            (give SOLVER the cpu_nofma build for this group; same case file)
#   twin     a case file whose coefficients are ZEROED (the body crossing
#            the outlet with coef = 0) == the body-free run, max_abs 0
#
# Every solve is preceded by tools/prepare_if_missing.sh with the moby_prepare
# beside SOLVER. SHORT=1 cuts every run to 2000 steps (the exactness groups
# do not need t = 40). Numbers are printed; the thresholds are in README.md.
set -u
cd "$(dirname "$0")" || exit 1
ROOT=$(cd ../.. && pwd)
SOLVER=${SOLVER:?absolute path to moby_solve}
RANKS=${RANKS:-1}
SHORT=${SHORT:-0}
sel=${1:-all}
status=0
export OMPI_MCA_hwloc_base_binding_policy=none

steps_of() { [ "$SHORT" = 1 ] && echo 2000 || echo 20000; }

# run <ini> <prefix> <ranks> [steps]
run() {
    local ini=$1 pfx=$2 ranks=$3 steps=${4:-$(steps_of)}
    sed -e "s/^nsteps = .*/nsteps = $steps/" -e "s/^field_prefix = .*/field_prefix = $pfx/" "$ini" > ".$pfx.ini"
    rm -f "${pfx}_"*.h5 "$pfx.case.h5"
    "$ROOT/tools/prepare_if_missing.sh" "$ranks" "$SOLVER" ".$pfx.ini" ".$pfx.prep.log" \
        || { echo "   PREPARE FAILED ($pfx)"; return 2; }
    if ! mpirun -n "$ranks" "$SOLVER" ".$pfx.ini" > ".$pfx.log" 2>&1; then
        echo "   RUN FAILED ($pfx) -- see .$pfx.log"; return 2
    fi
    grep -h 'boundary: non-zero Dirichlet\|removed .* of\|rdenom recomputed' ".$pfx.log" | sed 's/^/   /'
    tail -2 ".$pfx.log" | grep -q 'field filename' || { echo "   did not finish ($pfx)"; return 2; }
}
last() { ls -t "${1}_"*.h5 2>/dev/null | grep -v case | head -1; }
want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }

if want profile; then
    echo "== profile (t = $(steps_of) steps, 1 rank)"
    run pois_body.ini pb_io 1 && run pois_body_periodic.ini pb_per 1 \
        && python3 check_body_outlet.py profile "$(last pb_io)" "$(last pb_per)" pb_io.case.h5 | sed 's/^/   /' \
        && python3 check_body_outlet.py solid "$(last pb_io)" pb_io.case.h5 | sed 's/^/   /' || status=1
fi
if want mirror; then
    echo "== mirror (outlet on the low face)"
    [ -f "$(last pb_io)" ] || run pois_body.ini pb_io 1 || status=1
    run pois_body_low.ini pb_low 1 \
        && python3 check_body_outlet.py mirror "$(last pb_io)" "$(last pb_low)" | sed 's/^/   /' \
        && python3 check_body_outlet.py solid "$(last pb_low)" pb_low.case.h5 | sed 's/^/   /' || status=1
fi
if want permute; then
    echo "== permute (stream along y, along z)"
    [ -f "$(last pb_io)" ] || run pois_body.ini pb_io 1 || status=1
    run pois_body_y.ini pb_y 1 && python3 check_body_outlet.py permute "$(last pb_io)" "$(last pb_y)" y | sed 's/^/   /' || status=1
    run pois_body_z.ini pb_z 1 && python3 check_body_outlet.py permute "$(last pb_io)" "$(last pb_z)" z | sed 's/^/   /' || status=1
fi
if want limits; then
    echo "== limits (the body ends one cell before the plane / starts at it)"
    for c in ends starts; do
        run pois_body_$c.ini pb_$c 1 && python3 check_body_outlet.py solid "$(last pb_$c)" pb_$c.case.h5 | sed 's/^/   /' || status=1
    done
fi
if want ranks; then
    echo "== ranks ($RANKS vs 1, 2000 steps)"
    run pois_body.ini pb_r1 1 2000 && run pois_body.ini pb_rN "$RANKS" 2000 \
        && python3 check_body_outlet.py same "$(last pb_r1)" "$(last pb_rN)" || status=1
fi
if want restart; then
    echo "== restart (1000 + 1000 == 2000, $RANKS rank(s))"
    run pois_body.ini pb_a "$RANKS" 1000 && run pois_body.ini pb_c "$RANKS" 2000 || status=1
    sed -e 's/^nsteps = .*/nsteps = 1000/' -e 's/^field_prefix = .*/field_prefix = pb_b/' \
        -e "s|^\[output\]|[restart]\nfile = $(last pb_a)\n\n[output]|" pois_body.ini > .pb_b.ini
    rm -f pb_b_*.h5 pb_b.case.h5
    "$ROOT/tools/prepare_if_missing.sh" "$RANKS" "$SOLVER" .pb_b.ini .pb_b.prep.log \
        && mpirun -n "$RANKS" "$SOLVER" .pb_b.ini > .pb_b.log 2>&1 \
        && python3 check_body_outlet.py same "$(last pb_c)" "$(last pb_b)" || { echo "   restart leg failed"; status=1; }
fi
if want gpu; then
    echo "== gpu (GPU_SOLVER on 1 rank == CPU 1 rank, 2000 steps, tolerance 0)"
    [ -n "${GPU_SOLVER:-}" ] || { echo "   GPU_SOLVER not set"; status=1; }
    # its OWN CPU leg (pb_gc), so the pair is always the two binaries given
    # here (a pb_r1 left by the ranks group may be a production build); the
    # SAME case file on both sides (the CPU build prepares; a GPU-built
    # moby_prepare is host code too but is not the canonical one), and nofma
    # binaries on both sides: CPU and GPU contract FMAs differently.
    run pois_body.ini pb_gc 1 2000 || status=1
    sed -e 's/^nsteps = .*/nsteps = 2000/' -e 's/^field_prefix = .*/field_prefix = pb_g/' pois_body.ini > .pb_g.ini
    rm -f pb_g_*.h5; cp pb_gc.case.h5 pb_g.case.h5
    mpirun -n 1 "$GPU_SOLVER" .pb_g.ini > .pb_g.log 2>&1 \
        && python3 check_body_outlet.py same "$(last pb_gc)" "$(last pb_g)" || { echo "   gpu leg failed"; status=1; }
fi
if want twin; then
    echo "== twin (zeroed coefficients == body-free, 2000 steps)"
    # body-free: the same ini without the body (no removal either)
    sed -e 's/^enabled = true/enabled = false/' -e '/^stl_file/d' -e 's/^remove_solid = true/remove_solid = false/' \
        pois_body.ini > .pb_free.src.ini
    run .pb_free.src.ini pb_free 1 2000 || status=1
    # the twin: prepared WITH the body, every coefficient then zeroed in the case file
    sed -e 's/^nsteps = .*/nsteps = 2000/' -e 's/^field_prefix = .*/field_prefix = pb_twin/' \
        -e 's/^remove_solid = true/remove_solid = false/' pois_body.ini > .pb_twin.ini
    rm -f pb_twin_*.h5 pb_twin.case.h5
    "$ROOT/tools/prepare_if_missing.sh" 1 "$SOLVER" .pb_twin.ini .pb_twin.prep.log || status=1
    python3 - <<'PY'
import h5py
with h5py.File("pb_twin.case.h5", "r+") as f:
    f["coef_blocks"][...] = 0.0
PY
    mpirun -n 1 "$SOLVER" .pb_twin.ini > .pb_twin.log 2>&1 \
        && python3 check_body_outlet.py same "$(last pb_free)" "$(last pb_twin)" || { echo "   twin leg failed"; status=1; }
fi

echo
[ $status -eq 0 ] && echo "body_outlet gates ($sel): ALL PASS" || echo "body_outlet gates ($sel): FAILURES"
exit $status
