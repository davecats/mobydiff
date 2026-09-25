#!/usr/bin/env bash
# S3 gates: the passive scalars' IMMERSED-BODY coefficients
# (docs/next_session_scalar.md Section 4).
#
#   ./run_gates_s3.sh [solid|conserve|balance|prep|refine|det|missing|cyl|all]
#
# Environment: BIN  (default ../../build_cpu/moby_solve)
#              PREP (default the moby_prepare next to BIN)
#              NBIN/NPREP (nofma CPU pair, det group; default build_cpu_nofma)
#              GBIN (nofma GPU binary, det group; default build_gpu_nofma)
#              RANKS (default 1)
#
# The `cyl` group is the expensive one (2.1 M cells, tens of t.u. from the
# committed steady Re = 40 restart); everything else runs in minutes.
# Group order matters: `conserve` needs `solid`'s snapshot, `balance` needs
# `conserve`'s seeded IC, `prep`/`refine` need their case files.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

BIN=${BIN:-$ROOT/build_cpu/moby_solve}
PREP=${PREP:-$(dirname "$BIN")/moby_prepare}
NBIN=${NBIN:-$ROOT/build_cpu_nofma/moby_solve}
NPREP=${NPREP:-$ROOT/build_cpu_nofma/moby_prepare}
GBIN=${GBIN:-$ROOT/build_gpu_nofma/moby_solve}
RANKS=${RANKS:-1}
PY=${PY:-python3}
CMP="$PY $ROOT/tools/compare_fields.py --tolerance 0"
sel=${1:-all}
status=0

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
run() { echo "   \$ $*"; "$@"; }
report() { if [ "$1" -eq 0 ]; then echo "   PASS"; else echo "   FAIL"; status=1; fi }

# ini variant helpers -------------------------------------------------------
with_file() {  # <base.ini> <case.h5> <prefix> <out.ini>
    sed -e "s|^enabled = true|enabled = true\ncoeff_file = $2|" \
        -e "s|^field_prefix.*|field_prefix = $3|" "$1" > "$4"
}
strip_scalar() {  # <in.ini> <out.ini>: drop every [scalar]/[scalar.N] section
    $PY - "$1" "$2" <<'EOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
out, skip = [], False
for line in open(src):
    s = line.strip()
    if s.startswith("["):
        skip = s.startswith("[scalar")
    if not skip:
        out.append(line)
open(dst, "w").write("".join(out))
EOF
}

# --- (n) dirichlet body mode: solid cells hold ibm_value EXACTLY -----------
# Analytic wavy wall, cold start. theta is penalised toward 1 with the
# coefficient SOLID/Re inside the body, so mu_s underflows and every solid
# cell must read exactly 1.0 -- the sharpest statement the penalization can
# make, and it is an EQUALITY, not a tolerance.
if want solid || want conserve; then
    echo "== (n) dirichlet body mode: solid cell == ibm_value"
    rm -f ibmwavy_200.h5
    run mpirun -n "$RANKS" "$BIN" ibmwavy.ini > ibmwavy.log 2>&1
    run $PY check_scalar_ibm.py solid ibmwavy_200.h5 --scalar theta --value 1.0
    report $?
fi

# --- (o) adiabatic body mode: int phi dV conserved + sealed cells frozen ---
# The scalar is seeded through the restart path with a manufactured
# 2 + sin(k.x) (a NON-ZERO mean, so the integral is a real quantity), x and z
# are periodic and the y walls are Neumann 0: with every solid face masked
# the flux form has no sink at all.
if want conserve; then
    echo "== (o) adiabatic body mode: conservation + sealed cells"
    run $PY make_scalar_ic.py ibmwavy_200.h5 ibmw_ic.h5 --scalar phi \
        --wave 1 0 1 --amp 1.0 --offset 2.0
    sed -e "s|^field_prefix.*|field_prefix = ibmwc|" ibmwavy.ini > .ibmwc.ini
    printf '\n[restart]\nfile = ibmw_ic.h5\n' >> .ibmwc.ini
    rm -f ibmwc_400.h5
    run mpirun -n "$RANKS" "$BIN" .ibmwc.ini > ibmwc.log 2>&1
    run $PY check_scalar_ibm.py conserve ibmw_ic.h5 ibmwc_400.h5 --scalar phi
    report $?
    echo "   (and the dirichlet scalar again, after the restart)"
    run $PY check_scalar_ibm.py solid ibmwc_400.h5 --scalar theta --value 1.0
    report $?
fi

# --- (o2) the dirichlet body's heat release closes the energy budget ------
# The wavy case has NO boundary flux at all (periodic x,z, Neumann-0 y
# walls), so the heat the body releases must be exactly the rate at which
# the domain stores it. Measured over a SHORT window so the trapezoid error
# is negligible. This is also what establishes the Nusselt methodology for
# gate (s): the A2 penalization integral is structurally incomplete for a
# Dirichlet scalar (it sees ~63 % here) -- see check_scalar_ibm.py's header.
if want balance; then
    echo "== (o2) body heat release == d/dt int theta dV (no boundary flux)"
    [ -f ibmw_ic.h5 ] || { echo "   run the conserve group first (needs ibmw_ic.h5)"; status=1; }
    if [ -f ibmw_ic.h5 ]; then
        sed -e "s|^field_prefix.*|field_prefix = ibmbal|" -e "s|^nsteps.*|nsteps = 20|" \
            -e "s|^field_interval.*|field_interval = 10|" ibmwavy.ini > .ibmbal.ini
        printf '\n[restart]\nfile = ibmw_ic.h5\n' >> .ibmbal.ini
        rm -f ibmbal_*.h5
        run mpirun -n "$RANKS" "$BIN" .ibmbal.ini > ibmbal.log 2>&1
        run $PY check_scalar_ibm.py balance ibmbal_210.h5 ibmbal_220.h5 \
            --case ibmwavy_case.h5 --scalar theta --re 100 --pr 0.71 --tol 0.01
        report $?
    fi
fi

# --- (p) moby_prepare: coef_p_blocks is the ONLY change to the case file ---
if want prep; then
    echo "== (p) prepared case file: + coef_p_blocks and nothing else"
    strip_scalar ibmwavy.ini .ibmwavy_noscalar.ini
    rm -f ibmwavy_case.h5 ibmwavy_case_ns.h5 ibmwavy_case_np4.h5
    run mpirun -n 1 "$PREP" ibmwavy.ini ibmwavy_case.h5 > prep_ibmwavy.log 2>&1
    run mpirun -n 1 "$PREP" .ibmwavy_noscalar.ini ibmwavy_case_ns.h5 >> prep_ibmwavy.log 2>&1
    run $PY ../prepare/h5same.py ibmwavy_case_ns.h5 ibmwavy_case.h5 --ignore coef_p_blocks
    report $?
    echo "   coef_p_blocks vs an independent transcription of the graded formula"
    run $PY check_scalar_ibm.py coefp ibmwavy_case.h5 --re 100.0
    report $?
    echo "   prepare 1 rank == 4 ranks (Z-order row split)"
    run mpirun -n 4 --oversubscribe "$PREP" ibmwavy.ini ibmwavy_case_np4.h5 >> prep_ibmwavy.log 2>&1
    run $PY ../prepare/h5same.py ibmwavy_case.h5 ibmwavy_case_np4.h5
    report $?
    echo "   solve from the case file == the inline analytic solve (tolerance 0)"
    with_file ibmwavy.ini ibmwavy_case.h5 ibmwavyf .ibmwavyf.ini
    rm -f ibmwavyf_200.h5
    run mpirun -n "$RANKS" "$BIN" .ibmwavyf.ini > ibmwavyf.log 2>&1
    run $CMP ibmwavy_200.h5 ibmwavyf_200.h5 un vn wn pn theta phi
    report $?
    echo "   solid cells via the FILE path (coef_p tiles, not the analytic indicator)"
    run $PY check_scalar_ibm.py solid ibmwavyf_200.h5 --scalar theta --value 1.0 \
        --case ibmwavy_case.h5 --re 100.0
    report $?
fi

# --- (p2) a case file without coef_p_blocks is a hard error with scalars ---
if want missing; then
    echo "== (p2) [scalar] + a case file with no coef_p_blocks: hard error"
    strip_scalar ibmwavy.ini .ibmwavy_noscalar.ini
    rm -f ibmwavy_case_ns.h5
    run mpirun -n 1 "$PREP" .ibmwavy_noscalar.ini ibmwavy_case_ns.h5 > prep_ns.log 2>&1
    with_file ibmwavy.ini ibmwavy_case_ns.h5 ibmwmiss .ibmwmiss.ini
    if mpirun -n 1 "$BIN" .ibmwmiss.ini > ibmwmiss.log 2>&1; then
        echo "   solver ACCEPTED a case file without coef_p_blocks"; status=1
    else
        grep -i "coef_p_blocks" ibmwmiss.log | head -2 | sed 's/^/   /'
        grep -qi "re-run moby_prepare" ibmwmiss.log
        report $?
    fi
fi

# --- (q) refine_body: per-level cell-centred coefficients ------------------
if want refine; then
    echo "== (q) refine_body: coef_p at every leaf level"
    rm -f ibmwavyr_case.h5 ibmwavyr_50.h5 ibmwavyrf_50.h5
    run mpirun -n 1 "$PREP" ibmwavyr.ini ibmwavyr_case.h5 > prep_ibmwavyr.log 2>&1
    run $PY check_scalar_ibm.py coefp ibmwavyr_case.h5 --re 100.0
    report $?
    echo "   solve from the multi-level case file == the inline analytic solve"
    with_file ibmwavyr.ini ibmwavyr_case.h5 ibmwavyrf .ibmwavyrf.ini
    run mpirun -n "$RANKS" "$BIN" ibmwavyr.ini > ibmwavyr.log 2>&1
    run mpirun -n "$RANKS" "$BIN" .ibmwavyrf.ini > ibmwavyrf.log 2>&1
    run $CMP ibmwavyr_50.h5 ibmwavyrf_50.h5 un vn wn pn theta phi
    report $?
    run $PY check_scalar_ibm.py solid ibmwavyr_50.h5 --scalar theta --value 1.0
    report $?
fi

# --- (r) determinism on a scalar + IBM case (nofma builds) ----------------
if want det; then
    echo "== (r) determinism: 1 == 4 ranks, CPU == GPU (nofma)"
    for tag in r1 r4 gpu; do
        sed -e "s|^field_prefix.*|field_prefix = ibmd_$tag|" \
            -e "s|^nsteps.*|nsteps = 50|" -e "s|^field_interval.*|field_interval = 50|" \
            ibmwavy.ini > ".ibmd_$tag.ini"
        rm -f "ibmd_${tag}_50.h5"
    done
    run mpirun -n 1 "$NBIN" .ibmd_r1.ini > ibmd_r1.log 2>&1
    run mpirun -n 4 --oversubscribe "$NBIN" .ibmd_r4.ini > ibmd_r4.log 2>&1
    run $CMP ibmd_r1_50.h5 ibmd_r4_50.h5 un vn wn pn theta phi
    report $?
    if [ -x "$GBIN" ]; then
        run mpirun -n 1 "$GBIN" .ibmd_gpu.ini > ibmd_gpu.log 2>&1
        run $CMP ibmd_r1_50.h5 ibmd_gpu_50.h5 un vn wn pn theta phi
        report $?
    else
        echo "   GPU binary $GBIN missing -- SKIPPED"; status=1
    fi
    # The three code paths at once: the eddy diffusivity (nut), the body
    # coefficients (coef_p) and the scalar, on both devices.
    echo "   scalar + IBM + LES (WALE, ibm_aware) together"
    for tag in lr1 lgpu; do
        sed -e "s|^field_prefix.*|field_prefix = ibmd_$tag|" -e "s|^nsteps.*|nsteps = 50|" \
            -e "s|^field_interval.*|field_interval = 50|" ibmwavy.ini > ".ibmd_$tag.ini"
        printf '\n[turbulence]\nmodel = les\n\n[les]\nmodel = wale\nibm_aware = true\n' \
            >> ".ibmd_$tag.ini"
        rm -f "ibmd_${tag}_50.h5"
    done
    run mpirun -n 1 "$NBIN" .ibmd_lr1.ini > ibmd_lr1.log 2>&1
    run $PY check_scalar_ibm.py solid ibmd_lr1_50.h5 --scalar theta --value 1.0
    report $?
    if [ -x "$GBIN" ]; then
        run mpirun -n 1 "$GBIN" .ibmd_lgpu.ini > ibmd_lgpu.log 2>&1
        run $CMP ibmd_lr1_50.h5 ibmd_lgpu_50.h5 un vn wn pn nut theta phi
        report $?
    fi

    echo "   the same through the FILE path (coef_p_blocks read on both sides)"
    rm -f ibmwavy_case_nofma.h5
    run mpirun -n 1 "$NPREP" ibmwavy.ini ibmwavy_case_nofma.h5 > prep_nofma.log 2>&1
    for tag in fr1 fgpu; do
        with_file ibmwavy.ini ibmwavy_case_nofma.h5 "ibmd_$tag" ".ibmd_$tag.ini"
        sed -i -e "s|^nsteps.*|nsteps = 50|" -e "s|^field_interval.*|field_interval = 50|" \
            ".ibmd_$tag.ini"
        rm -f "ibmd_${tag}_50.h5"
    done
    run mpirun -n 1 "$NBIN" .ibmd_fr1.ini > ibmd_fr1.log 2>&1
    if [ -x "$GBIN" ]; then
        run mpirun -n 1 "$GBIN" .ibmd_fgpu.ini > ibmd_fgpu.log 2>&1
        run $CMP ibmd_fr1_50.h5 ibmd_fgpu_50.h5 un vn wn pn theta phi
        report $?
    fi
fi

# --- (s) heated cylinder Re 40: Nusselt vs literature and vs the CV flux ---
# ZEROING pn IN THE RESTART IS LOAD-BEARING (../cylinder/README.md, the A2
# pn-drift caveat): the committed steady restart carries a large
# VELOCITY-NEUTRAL mode in its stored pressure (|pn| ~ 1.2e3 here) whose
# spurious grad p is only self-consistent with the coefficient field that
# grew it. Restarting it as-is on the freshly PREPARED coefficients (they
# differ from the retired mobygeom ones in the near-grazing envelope) kicked
# the forces to |C_L| ~ 1.7e2, oscillating with no sign of decay. With pn
# zeroed the SAME restart is clean at the production niter = 6 -- C_D is back
# at 1.69 within 40 steps, C_L ~ 5e-4 -- so no niter = 60 rebuild phase is
# needed here (the A2 momentum cross-check needed one because it reads the
# pressure; this thermal balance never does).
CYLBIN=${CYLBIN:-$ROOT/build_gpu/moby_solve}
if want cyl; then
    echo "== (s) heated cylinder Re 40, Pr 0.71"
    if [ ! -f ../cylinder/cyl_re40_20001.h5 ]; then
        echo "   ../cylinder/cyl_re40_20001.h5 (the steady A2 restart) missing -- SKIPPED"
        status=1
    else
        if [ ! -f cylheat_case.h5 ]; then
            run mpirun -n 4 --oversubscribe "$PREP" cylheat.ini cylheat_case.h5 \
                > prep_cylheat.log 2>&1
        fi
        $PY - <<'EOF'
import shutil, h5py
shutil.copyfile("../cylinder/cyl_re40_20001.h5", "cylheat_ic.h5")
with h5py.File("cylheat_ic.h5", "r+") as f:
    f["pn"][...] = 0.0
print("cylheat_ic.h5: the steady A2 restart with pn zeroed")
EOF
        with_file cylheat.ini cylheat_case.h5 cylheat .cylheat_run.ini
        sed -i -e "/^stl_file/d" \
            -e "s|^force_sample_interval.*|force_sample_interval = 20|" .cylheat_run.ini
        printf '\n[restart]\nfile = cylheat_ic.h5\n' >> .cylheat_run.ini
        rm -f cylheat_2*.h5 forces_cylheat.txt
        run mpirun -n 1 "$CYLBIN" .cylheat_run.ini > cylheat.log 2>&1
        echo "   final force sample: $(tail -1 forces_cylheat.txt)"

        # Nu at every snapshot: its own drift is the convergence measure.
        for h5 in $(ls -tr cylheat_2*.h5); do
            echo "   --- $h5"
            $PY check_scalar_ibm.py surface "$h5" --case cylheat_case.h5 \
                --scalar theta --re 40 --pr 0.71 --value 1.0 | sed 's/^/   /'
        done
        last=$(ls -t cylheat_2*.h5 | head -1)
        prev=$(ls -t cylheat_2*.h5 | head -2 | tail -1)
        run $PY check_scalar_ibm.py surface "$last" --case cylheat_case.h5 \
            --scalar theta --re 40 --pr 0.71 --value 1.0 --band 2.8 4.0
        report $?
        nubody=$($PY check_scalar_ibm.py surface "$last" --case cylheat_case.h5 \
            --scalar theta --re 40 --pr 0.71 --value 1.0 | awk '/^Nu \(surface\)/{print $4}')
        echo "   the INDEPENDENT Gauss/CV border flux (+ storage term)"
        for box in "4 8 6 10" "3 9 5 11" "2 10 4 12"; do
            run $PY check_scalar_ibm.py cv "$last" --re 40 --pr 0.71 --value 1.0 \
                --box $box --prev "$prev" --nu-pen "$nubody" --tol 0.10
            report $?
        done
    fi
fi

echo
[ $status -eq 0 ] && echo "S3 gates ($sel): ALL PASS" || echo "S3 gates ($sel): FAILURES"
exit $status
