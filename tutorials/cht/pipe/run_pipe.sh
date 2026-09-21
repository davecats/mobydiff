#!/usr/bin/env bash
# Conjugate pipe flow vs Neuhauser (NekRS body-fitted DNS) -- the campaign
# driver. docs/next_session_pipe_cht.md Sections 0.1-0.4.
#
#   ./run_pipe.sh [geom|prepare|check|dev|vel|therm|settle|settle2|stats]
#
# Environment: GRID (prod|coarse|fine|zfine -- prod is the result), BIN, PREP,
# RANKS, PY, and the per-leg step counts below.
#
# THE LEGS, and why each one exists:
#   geom     the annulus STL (pure padding outside r = 0.5; the fluid is the
#            hole). Grid-independent, so every grid shares it.
#   prepare  moby_prepare -> the case file. CPU build: the GPU build computes
#            coef on the device and differs by libm ulps.
#   check    the F5/F2 geometry gate on the prepared file -- the level set
#            -0.1 < phi < 0 must BE the shell 0.5 < r < 0.6.
#   dev      leg A: hydrodynamics, NO scalars. If the velocity does not match
#            Neuhauser nothing thermal means anything, and this isolates the
#            immersed pipe from the conjugate scheme. Start it from
#            make_pipe_ic.py (IC=<file>): a plug plus white noise
#            RELAMINARISES, because a plug carries no shear to feed it.
#   vel      leg A's statistics window: snapshots for pipe_stats.py.
#   therm    leg B: the 7 scalars on leg A's velocity, with c4's capacity
#            temporarily at 1 (the accelerator: the solid MEAN is
#            capacity-independent, C1-gated, so it may be reached fast and
#            the real capacity restored afterwards).
#   settle   leg C: c4 back to solid_rhocp = 16. Only theta' in the solid has
#            to re-equilibrate. CHECK the solid mean did not move.
#   settle2  a second settle when the interface heat says the slowest scalar
#            is not converged yet (FROM=<file> overrides what `stats` starts
#            from).
#   stats    leg D: the measurement -- scalar statistics + snapshots.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../../.. && pwd)

GRID=${GRID:-prod}
BIN=${BIN:-$ROOT/build_gpu/moby_solve}
PREP=${PREP:-$ROOT/build_cpu/moby_prepare}
RANKS=${RANKS:-1}
PREP_RANKS=${PREP_RANKS:-4}
PY=${PY:-$HOME/ibmc/bin/python}
sel=${1:-}

CASE=pipe_${GRID}.h5
TAG=${TAG:-p_${GRID}}

# THE PRODUCTION GRID is `prod`: Delta+ 1.84 radially (which resolves the
# SOLID side -- the shell is only 36 wall units thick) with Dz+ 5.06 axially
# (which resolves the fluid, and is the spacing the core variance actually
# converges in). 58.7 M cells, 22.7 GB on the device. The three coarser grids
# are kept because the resolution study in the README quotes them; they are
# not needed to reproduce the result.
case "$GRID" in
    prod)   NX=256 NZ=896 ;;   # Delta+ 1.84, Dz+ 5.06  -- THE RESULT
    coarse) NX=160 NZ=320 ;;   # Delta+ 2.94, Dz+ 14.2
    fine)   NX=256 NZ=448 ;;   # Delta+ 1.84, Dz+ 10.1
    zfine)  NX=160 NZ=896 ;;   # Delta+ 2.94, Dz+ 5.06  -- the decisive test
    *) echo "GRID must be prod, coarse, fine or zfine"; exit 2 ;;
esac

# Step counts. dt is set by the solver's own limiter -- convection in leg A,
# the isc 3 (C_s = 1/16) solid cut cells in the thermal legs -- so these are
# STEPS, not times: read the achieved dt out of the log and convert.
DEV_STEPS=${DEV_STEPS:-30000}
VEL_STEPS=${VEL_STEPS:-20000}
VEL_FIELD=${VEL_FIELD:-500}
THERM_STEPS=${THERM_STEPS:-40000}
SETTLE_STEPS=${SETTLE_STEPS:-20000}
SETTLE2_STEPS=${SETTLE2_STEPS:-50000}
STAT_STEPS=${STAT_STEPS:-200000}
STAT_FIELD=${STAT_FIELD:-4000}
SAMPLE=${SAMPLE:-25}
FLUSH=${FLUSH:-2000}
HEAT=${HEAT:-1000}
# dtmax: a ceiling, not the step. Leg A is convection-limited (~5e-3 coarse);
# the thermal legs are limited far below this by the solid capacity.
DTMAX=${DTMAX:-5.0e-3}
NOISE=${NOISE:-0.3}

# The newest snapshot of a prefix. Step numbers CONTINUE across legs (the
# restart carries its step in), so a leg that runs N steps does not end at
# step N -- leg A's 20000-step statistics window ended at 35000, not 20000.
# Asking the files rather than doing the arithmetic keeps the legs composable.
final_snapshot() {  # <prefix>
    ls -1 "$1"_[0-9]*.h5 2>/dev/null \
        | sed -E 's|.*_([0-9]+)\.h5$|\1 &|' | sort -n | tail -1 | cut -d' ' -f2
}

run() {  # <name> <ini>
    echo "== $1"
    local t0=$SECONDS
    mpirun -n "$RANKS" "$BIN" "$2" > "$1.log" 2>&1 \
        || { echo "   FAILED:"; tail -20 "$1.log"; exit 1; }
    grep -E "seconds_per_step" "$1.log" | tail -1
    echo "   wall $((SECONDS - t0)) s"
}

hydro() {  # <prefix> <nsteps> <write> <restart|->  <noise>
    sed -e "s|@CASE@|$CASE|" -e "s|@PREFIX@|$1|" -e "s|@NSTEPS@|$2|" \
        -e "s|@WRITE@|$3|" -e "s|@RESTART@|$4|" -e "s|@NOISE@|$5|" \
        -e "s|@DTMAX@|$DTMAX|" \
        -e "s|^nx = .*|nx = $NX|" -e "s|^ny = .*|ny = $NX|" \
        -e "s|^nz = .*|nz = $NZ|" \
        pipe_hydro.ini > ".$1.ini"
    [ "$4" = "-" ] && sed -i '/^\[restart\]/,$d' ".$1.ini"
    grep -v '^;' ".$1.ini" | grep -q '@' && { echo "unsubstituted placeholder"; exit 1; }
    return 0
}

# SRCSCALE rescales every shell sink by a measured factor. The sink must
# balance the DISCRETE generation `int w dV`, which is a property of the GRID
# (the coarse one gives 0.9755 of the continuum u_b = 1), and campaign 1
# showed what an uncorrected 2.8 % imbalance does to time-averaged second
# moments. So each grid measures its own factor before its statistics leg.
thermal() {  # <prefix> <nsteps> <write> <restart> <rhocp4> <src4> <sample> <flush> <statsfile>
    sed -e "s|@CASE@|$CASE|" -e "s|@PREFIX@|$1|" -e "s|@NSTEPS@|$2|" \
        -e "s|@WRITE@|$3|" -e "s|@RESTART@|$4|" -e "s|@RHOCP4@|$5|" \
        -e "s|@SRC4@|$6|" -e "s|@SAMPLE@|$7|" -e "s|@STATS@|$8|" \
        -e "s|@STATSFILE@|$9|" -e "s|@HEAT@|$HEAT|" -e "s|@DTMAX@|$DTMAX|" \
        -e "s|^nx = .*|nx = $NX|" -e "s|^ny = .*|ny = $NX|" \
        -e "s|^nz = .*|nz = $NZ|" \
        pipe_thermal.ini > ".$1.ini"
    if [ -n "${SRCSCALE:-}" ]; then
        awk -v f="$SRCSCALE" \
            '/^solid_source = /{printf "solid_source = %.10g\n", $3*f; next} {print}' \
            ".$1.ini" > ".$1.ini.tmp" && mv ".$1.ini.tmp" ".$1.ini"
        echo "   sinks scaled by $SRCSCALE"
    fi
    grep -v '^;' ".$1.ini" | grep -q '@' && { echo "unsubstituted placeholder"; exit 1; }
    return 0
}

case "$sel" in
geom)
    $PY "$ROOT/tools/make_geometry_stl.py" annulus pipe.stl --axis z --centre 0.65 0.65 \
        --r-inner 0.5 --box-half 1.25 --facets 16384 --a0 -1.0 --a1 13.5 \
        --domain-half 0.65
    ;;
prepare)
    sed -e "s|^nx = .*|nx = $NX|" -e "s|^ny = .*|ny = $NX|" \
        -e "s|^nz = .*|nz = $NZ|" pipe_prep.ini > ".prep_$GRID.ini"
    time mpirun -n "$PREP_RANKS" "$PREP" ".prep_$GRID.ini" "$CASE" \
        > "$CASE.prep.log" 2>&1 || { tail -10 "$CASE.prep.log"; exit 1; }
    ls -la "$CASE"
    ;;
check)
    $PY "$ROOT/tools/check_annulus.py" "$CASE" --centre 0.65 0.65 --r-inner 0.5 \
        --facets 16384 --box-half 1.25 --domain-half 0.65 --band-depth 0.1
    ;;
dev)
    # IC=<file> starts from make_pipe_ic.py's field; "-" is the solver's own
    # cold start, which for a pipe RELAMINARISES (see make_pipe_ic.py).
    hydro "${TAG}_dev" "$DEV_STEPS" "${DEV_FIELD:-$((DEV_STEPS / 6))}" "${IC:--}" "$NOISE"
    run "${TAG}_dev" ".${TAG}_dev.ini"
    ;;
vel)
    hydro "${TAG}_vel" "$VEL_STEPS" "$VEL_FIELD" \
        "${FROM:-$(final_snapshot "${TAG}_dev")}" 0.0
    run "${TAG}_vel" ".${TAG}_vel.ini"
    ;;
therm)
    # c4 accelerated: solid_rhocp 1 with the sink rescaled to the SAME power
    # (C_s * solid_source is what the balance fixes).
    thermal "${TAG}_therm" "$THERM_STEPS" "$((THERM_STEPS / 4))" \
        "$(final_snapshot "${TAG}_vel")" 1.0 -2.21110542 0 0 x
    run "${TAG}_therm" ".${TAG}_therm.ini"
    ;;
settle)
    thermal "${TAG}_settle" "$SETTLE_STEPS" "$((SETTLE_STEPS / 2))" \
        "$(final_snapshot "${TAG}_therm")" 16.0 -0.13819409 0 0 x
    run "${TAG}_settle" ".${TAG}_settle.ini"
    ;;
settle2)
    # A SECOND settle, on the same configuration as leg C, for when the
    # interface heat says the slowest scalars are not there yet. The
    # statistics accumulator is cumulative from a run's first step, so a
    # residual transient at the start of leg D would be averaged INTO the
    # product -- settling longer is cheaper than discarding a statistics leg.
    thermal "${TAG}_settle2" "$SETTLE2_STEPS" "$((SETTLE2_STEPS / 2))" \
        "${FROM:-$(final_snapshot "${TAG}_settle")}" 16.0 -0.13819409 0 0 x
    run "${TAG}_settle2" ".${TAG}_settle2.ini"
    ;;
stats)
    rm -f "${TAG}_stats.h5"
    thermal "${TAG}_stat" "$STAT_STEPS" "$STAT_FIELD" \
        "${FROM:-$(final_snapshot "${TAG}_settle2")}" 16.0 -0.13819409 \
        "$SAMPLE" "$FLUSH" "${TAG}_stats.h5"
    run "${TAG}_stat" ".${TAG}_stat.ini"
    ;;
*)
    sed -n '2,34p' "$0"
    exit 2
    ;;
esac
