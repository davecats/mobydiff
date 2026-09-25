#!/usr/bin/env bash
# S4 gates: the passive scalars' STATISTICS AND TOOLING
# (docs/next_session_scalar.md Section 9, scalar_stats.f90).
#
#   ./run_gates_s4.sh [stats|accum|plane|levels|restart|det|noeffect|heat|adia|cyl|tools|all]
#
# Environment: BIN  (default ../../build_cpu/moby_solve)
#              GBIN (GPU binary, det + cyl groups; default ../../build_gpu/moby_solve)
#              RANKS (default 1)
#
# Every group runs in seconds to a minute: the statistics groups are 10-40
# steps of the S2 LES channel restarted from one of its own snapshots, the
# body groups 200 steps of the S3 wavy case, and `cyl` 20 steps of the S3
# heated cylinder from its converged field.
#
# THE GATE IDEA, once: the solver samples the END-OF-STEP field and the
# snapshot IS the end-of-step field, so check_scalar_stats.py -- recomputing
# the same seven columns from the snapshots the same run wrote -- must
# reproduce the solver's rows to ROUND-OFF. There is no statistical tolerance
# to hide in. Likewise the body heat release is compared with
# check_scalar_ibm.py's `surface`, the Python form S3 validated against the
# full energy budget.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

BIN=${BIN:-$ROOT/build_cpu/moby_solve}
GBIN=${GBIN:-$ROOT/build_gpu/moby_solve}
RANKS=${RANKS:-1}
PY=${PY:-python3}
sel=${1:-all}
status=0

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
run() { echo "   \$ $*"; "$@"; }
report() { if [ "$1" -eq 0 ]; then echo "   PASS"; else echo "   FAIL"; status=1; fi }

# s4stats.ini variants: <sed args...> <out.ini>
variant() {
    local out=${!#}
    sed "${@:1:$#-1}" s4stats.ini > "$out"
}

check_profile() {   # <stats.h5> <fields...>
    local stats=$1; shift
    run $PY check_scalar_stats.py profile "$stats" "$@" \
        --scalar theta --index 1 --pr 0.71 --prt 0.85
    report $?
    run $PY check_scalar_stats.py profile "$stats" "$@" \
        --scalar theta_kc --index 2 --pr 0.71 --prt 0.85 --prt-model kays
    report $?
}

# --- (t) one sample: the solver's rows == the snapshot's, to round-off -----
# Two scalars differing ONLY in the Pr_t model, so this also gates the
# multi-scalar column layout and the Kays-Crawford branch of the statistics
# kernel's face diffusivity.
if want stats; then
    echo "== (t) profile layout, ONE sample: solver == recomputed from the snapshot"
    variant -e "s|^nsteps = 40|nsteps = 10|" \
            -e "s|^stats_file = s4stats.h5|stats_file = s4one.h5|" \
            -e "s|^field_prefix = s4stats|field_prefix = s4one|" .s4one.ini
    rm -f s4one*.h5
    run mpirun -n "$RANKS" "$BIN" .s4one.ini > s4one.log 2>&1
    check_profile s4one.h5 s4one_67610.h5
fi

# --- (t2) four samples: the ACCUMULATION is the same sum -------------------
if want accum || want restart || want tools; then
    echo "== (t2) profile layout, FOUR samples: the accumulated file"
    rm -f s4stats*.h5
    run mpirun -n "$RANKS" "$BIN" s4stats.ini > s4stats.log 2>&1
    if want accum; then check_profile s4stats.h5 s4stats_676[1-4]0.h5; fi
fi

# --- (t3) the (x,y) plane layout (the boundary-layer form) ----------------
# Exercised on the same channel: the plane statistics of a channel are
# perfectly well defined (x-inhomogeneous rows, z averaged), and running a
# boundary-layer campaign here would gate the same kernel far more slowly.
if want plane; then
    echo "== (t3) plane layout: solver == recomputed from the snapshot"
    variant -e "s|^stats_layout = profile|stats_layout = plane|" \
            -e "s|^stats_file = s4stats.h5|stats_file = s4plane.h5|" \
            -e "s|^field_prefix = s4stats|field_prefix = s4plane|" \
            -e "s|^nsteps = 40|nsteps = 10|" .s4plane.ini
    rm -f s4plane*.h5
    run mpirun -n "$RANKS" "$BIN" .s4plane.ini > s4plane.log 2>&1
    run $PY check_scalar_stats.py plane s4plane.h5 s4plane_67610.h5 --scalar theta --index 1
    report $?
    run $PY check_scalar_stats.py plane s4plane.h5 s4plane_67610.h5 \
        --scalar theta_kc --index 2 --prt-model kays
    report $?
fi

# --- (t3b) the 2:1-refined case: one statistics file PER LEVEL ------------
# The row tables are per refinement level (the channel_stats lvlOff layout),
# so the wall-band case writes s4slab.h5 (level 0) and s4slab_l1.h5. The
# `rows` checker accumulates straight from the leaves -- no global box -- so
# it can gate the level split that `profile` cannot reach; the face-flux
# columns stay gated at single level, where the halo values exist.
if want levels; then
    echo "== (t3b) 2:1 wall-band refined: per-level statistics files"
    sed -e "s|^count = 1|count = 1\nstats_sample_interval = 10\nstats_write_interval = 10\nstats_file = s4slab.h5|" \
        -e "s|^field_interval.*|field_interval = 10|" \
        -e "s|^nsteps.*|nsteps = 10|" \
        -e "s|^field_prefix.*|field_prefix = s4slab|" turbslab.ini > .s4slab.ini
    rm -f s4slab*.h5
    run mpirun -n "$RANKS" "$BIN" .s4slab.ini > s4slab.log 2>&1
    for lev in 0 1; do
        f=s4slab.h5; [ "$lev" = 1 ] && f=s4slab_l1.h5
        run $PY check_scalar_stats.py rows "$f" s4slab_[0-9]*.h5 --level "$lev" \
            --scalar theta --index 1
        report $?
    done
fi

# --- (t4) restart continuation: 2 + 2 samples == 4 samples ---------------
if want restart; then
    echo "== (t4) restart: the accumulators continue from the file"
    variant -e "s|^nsteps = 40|nsteps = 20|" \
            -e "s|^stats_file = s4stats.h5|stats_file = s4rst.h5|" \
            -e "s|^field_prefix = s4stats|field_prefix = s4rst|" .s4rst1.ini
    sed -e "s|^file = turbles_67600.h5|file = s4rst_67620.h5|" .s4rst1.ini > .s4rst2.ini
    rm -f s4rst*.h5
    run mpirun -n "$RANKS" "$BIN" .s4rst1.ini > s4rst1.log 2>&1
    run mpirun -n "$RANKS" "$BIN" .s4rst2.ini > s4rst2.log 2>&1
    grep -h "continuing scalar statistics" s4rst2.log
    run $PY check_scalar_stats.py diff s4stats.h5 s4rst.h5 --tol 0
    report $?
fi

# --- (t5) determinism: 1 == 4 ranks, CPU == GPU ---------------------------
# A TOLERANCE, not an equality: the sampling kernel reduces with atomics and
# the write reduces across ranks, so the summation order differs by
# construction. (The FIELDS are the bit-exact ones -- run_bitexact*.sh.)
if want det; then
    echo "== (t5) determinism of the statistics"
    variant -e "s|^nsteps = 40|nsteps = 10|" \
            -e "s|^stats_file = s4stats.h5|stats_file = s4d_r1.h5|" \
            -e "s|^field_prefix = s4stats|field_prefix = s4d_r1|" .s4d_r1.ini
    sed -e "s|^dims = 1 1 1||" -e "s|s4d_r1|s4d_r4|g" .s4d_r1.ini > .s4d_r4.ini
    sed -e "s|s4d_r1|s4d_gpu|g" .s4d_r1.ini > .s4d_gpu.ini
    rm -f s4d_r1*.h5 s4d_r4*.h5 s4d_gpu*.h5
    run mpirun -n 1 "$BIN" .s4d_r1.ini > s4d_r1.log 2>&1
    run mpirun -n 4 --oversubscribe "$BIN" .s4d_r4.ini > s4d_r4.log 2>&1
    run $PY check_scalar_stats.py diff s4d_r1.h5 s4d_r4.h5 --tol 1e-12
    report $?
    if [ -x "$GBIN" ]; then
        run mpirun -n 1 "$GBIN" .s4d_gpu.ini > s4d_gpu.log 2>&1
        run $PY check_scalar_stats.py diff s4d_r1.h5 s4d_gpu.h5 --tol 1e-11
        report $?
    else
        echo "   no GPU binary at $GBIN -- CPU == GPU SKIPPED"
    fi
fi

# --- (t6) the statistics do not touch the solution ------------------------
if want noeffect; then
    echo "== (t6) statistics ON == statistics OFF (tolerance 0, same binary)"
    variant -e "s|^nsteps = 40|nsteps = 10|" \
            -e "s|^stats_file = s4stats.h5|stats_file = s4d_r1.h5|" \
            -e "s|^field_prefix = s4stats|field_prefix = s4d_r1|" .s4d_r1.ini
    sed -e "s|^stats_sample_interval = 10|stats_sample_interval = 0|" \
        -e "s|^stats_write_interval = 10|stats_write_interval = 0|" \
        -e "s|s4d_r1|s4off|g" .s4d_r1.ini > .s4off.ini
    rm -f s4d_r1*.h5 s4off_*.h5
    run mpirun -n "$RANKS" "$BIN" .s4d_r1.ini > s4d_r1.log 2>&1
    run mpirun -n "$RANKS" "$BIN" .s4off.ini > s4off.log 2>&1
    run $PY $ROOT/tools/compare_fields.py s4d_r1_67610.h5 s4off_67610.h5 --tolerance 0
    report $?
fi

# --- (u) the immersed body's heat release --------------------------------
# The runtime samples against check_scalar_ibm.py's `surface` (round-off),
# and then the physics: with no boundary flux anywhere the body's heat
# release IS the domain's storage rate, closed here with the SOLVER's own Q.
if want heat; then
    echo "== (u) body heat release: solver == the validated Python form"
    [ -f ibmwavy_case.h5 ] || echo "   (ibmwavy_case.h5 missing -- run ./run_gates_s3.sh prep)"
    sed -e "s|^count = 2|count = 2\nheat_interval = 10\nheat_file = s4heat.txt|" \
        -e "s|^field_interval = 200|field_interval = 10|" \
        -e "s|^field_prefix = ibmwavy|field_prefix = s4heat|" ibmwavy.ini > .s4heat.ini
    rm -f s4heat_*.h5 s4heat.txt
    run mpirun -n "$RANKS" "$BIN" .s4heat.ini > s4heat.log 2>&1
    run $PY check_scalar_stats.py heat s4heat.txt "s4heat_STEP.h5" --case ibmwavy_case.h5 \
        --scalar theta --re 100 --pr 0.71 --value 1.0 --balance
    report $?
fi

# --- (u2) an adiabatic scalar exchanges nothing: the positive control -----
# The seeded phi field (2 + sin(k.x), manifestly non-zero) with ibm_wall =
# adiabatic must report EXACTLY zero heat -- no penalization is applied and
# every body face is masked, so any non-zero number would be a flux the
# solver never applied.
if want adia; then
    echo "== (u2) adiabatic body mode: the heat columns are exactly zero"
    [ -f ibmw_ic.h5 ] || echo "   (ibmw_ic.h5 missing -- run ./run_gates_s3.sh conserve)"
    sed -e "s|^count = 2|count = 2\nheat_interval = 10\nheat_file = s4adia.txt|" \
        -e "s|^nsteps = 200|nsteps = 10|" \
        -e "s|^field_interval = 200|field_interval = 10|" \
        -e "s|^field_prefix = ibmwavy|field_prefix = s4adia|" ibmwavy.ini > .s4adia.ini
    printf '\n[restart]\nfile = ibmw_ic.h5\n' >> .s4adia.ini
    rm -f s4adia_*.h5 s4adia.txt
    run mpirun -n "$RANKS" "$BIN" .s4adia.ini > s4adia.log 2>&1
    $PY - <<'EOF'
import sys, h5py, numpy as np
row = np.loadtxt("s4adia.txt", ndmin=2)[-1]
with h5py.File("s4adia_210.h5", "r") as f:
    phi = f["phi"][...]
print(f"   phi field range {phi.min():.6f} .. {phi.max():.6f}"
      f"   heat columns {row[5]:.3e} {row[6]:.3e} {row[7]:.3e}")
ok = phi.ptp() > 1.0 and row[5] == 0.0 and row[6] == 0.0 and row[7] == 0.0
sys.exit(0 if ok else 1)
EOF
    report $?
fi

# --- (u3) heated cylinder: the RUNTIME Nusselt number --------------------
# 20 steps from the S3 campaign's converged t = 120 field, so the number is
# directly comparable with the S3 post-processed 3.3655 (and with
# Churchill-Bernstein's 3.35).
if want cyl; then
    echo "== (u3) heated cylinder Re 40: runtime Nusselt"
    if [ ! -f cylheat_24000.h5 ] || [ ! -f cylheat_case.h5 ]; then
        echo "   cylheat_24000.h5 / cylheat_case.h5 missing -- run ./run_gates_s3.sh cyl -- SKIPPED"
    else
        sed -e "s|^enabled = true|enabled = true\ncoeff_file = cylheat_case.h5|" \
            -e "s|^field_prefix.*|field_prefix = s4cyl|" \
            -e "s|^field_interval.*|field_interval = 20|" \
            -e "s|^t_final.*|t_final = 130.0\nnsteps = 20|" \
            -e "s|^count = 1|count = 1\nheat_interval = 20\nheat_file = s4cyl.txt|" \
            -e "s|^runtime_file.*|runtime_file = forces_s4cyl.txt|" cylheat.ini > .s4cyl.ini
        printf '\n[restart]\nfile = cylheat_24000.h5\n' >> .s4cyl.ini
        rm -f s4cyl_*.h5 s4cyl.txt
        cbin=$BIN; [ -x "$GBIN" ] && cbin=$GBIN
        run mpirun -n 1 "$cbin" .s4cyl.ini > s4cyl.log 2>&1
        run $PY check_scalar_stats.py heat s4cyl.txt "s4cyl_STEP.h5" --case cylheat_case.h5 \
            --scalar theta --re 40 --pr 0.71 --value 1.0
        report $?
        run $PY $ROOT/tools/scalar_stats.py heat s4cyl.txt --re 40 --pr 0.71 --lz 0.25 \
            --value 1.0 --band 3.2 3.5
        report $?
    fi
fi

# --- (v) tooling: compare_fields.py dataset discovery ---------------------
if want tools; then
    echo "== (v) compare_fields.py discovers the datasets present in both files"
    $PY - <<'EOF'
import subprocess, sys
out = subprocess.run([sys.executable, "../../tools/compare_fields.py",
                      "s4stats_67610.h5", "s4stats_67620.h5"],
                     capture_output=True, text=True).stdout
found = out.splitlines()[0].removeprefix("datasets: ").split()
want = ["un", "vn", "wn", "pn", "nut", "theta", "theta_kc"]
print("   discovered:", " ".join(found))
sys.exit(0 if found == ["un", "vn", "wn", "pn"] + sorted(set(want) - {"un","vn","wn","pn"})
         else 1)
EOF
    report $?
    echo "   and the profile/heat readers run on the files this session wrote"
    run $PY $ROOT/tools/scalar_stats.py profile s4stats.h5 --pr 0.71 --walls 1 -1 --rows 2
    report $?
fi

exit $status
