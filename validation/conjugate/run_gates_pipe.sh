#!/usr/bin/env bash
# Pipe-campaign prerequisites F1 and F5 (docs/next_session_pipe_cht.md).
#
#   F1  [scalar.N] source_dir -- the Kasagi source on a chosen direction.
#   F5  make_geometry_stl.py annulus -- the pipe's immersed body.
#
#   ./run_gates_pipe.sh [source_dir|guard|annulus|ranks|all]
#
# Environment: BIN   (default ../../build_cpu/moby_solve)
#              PREP  (default the moby_prepare next to BIN)
#              RANKS (default 4, prepare only)
#              PY    (default python3; needs h5py, and ~/ibmc/bin/python here)
#
# The source_dir group runs in seconds; the annulus group prepares two 64^2x8
# cases and takes a couple of minutes (the 16384-facet BVH is the cost).
#
# NOT covered here, deliberately: that source_dir OFF is bit-exact. That is
# the standard suite's job -- validation/scalar/run_bitexact{,_s3}.sh against
# a reference binary -- and it is where the "inert by construction" claim is
# actually settled, on 16 cases, CPU and GPU.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

BIN=${BIN:-$ROOT/build_cpu/moby_solve}
PREP=${PREP:-$(dirname "$BIN")/moby_prepare}
RANKS=${RANKS:-4}
PY=${PY:-python3}
sel=${1:-all}
status=0

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }
run() { echo "   \$ $*"; "$@"; }
report() { if [ "$1" -eq 0 ]; then echo "   PASS"; else echo "   FAIL"; status=1; fi }

solve() {  # <template> <prefix> <steps> <sed-expr...>
    local tpl=$1 pre=$2 steps=$3; shift 3
    local ini=".$pre.ini"
    sed -e "s|@PREFIX@|$pre|" -e "s|@STEPS@|$steps|" "$@" "$tpl" > "$ini"
    mpirun -n 1 "$BIN" "$ini" > "$pre.log" 2>&1
}

# --- F1 (1): the closed form, all three branches ---------------------------
# Uniform (u,v,w) = (1,2,3) is an exact steady solution and a uniform scalar
# has no gradient, so theta = source * u_dir * t EXACTLY -- three different
# numbers for the three directions, which is what makes it a direction test.
if want source_dir; then
    echo "== F1 (1) source_dir vs the closed form theta = source * u_dir * t"
    for d in x y z; do
        solve kasagi_uniform.ini "ku_$d" 20 -e "s|@SDIR@|$d|" || { report 1; continue; }
    done
    run $PY ./check_source_dir.py closed-form ku_x_20.h5 ku_y_20.h5 ku_z_20.h5 \
        --components 1.0 2.0 3.0 --time 0.04
    report $?

    # --- F1 (2): the same arithmetic in x and in z, with structure ---------
    echo "== F1 (2) diagonal channel: source_dir x and z are BIT-IDENTICAL"
    for d in x z; do
        solve kasagi_channel.ini "kc_$d" 50 -e "s|@SDIR@|$d|" || { report 1; continue; }
    done
    run $PY ./check_source_dir.py identical kc_x_50.h5 kc_z_50.h5
    report $?
fi

# --- F1 (3): the config guards --------------------------------------------
if want guard; then
    echo "== F1 (3) config guards"
    # NOTE the substitution order: the direction placeholder is filled with
    # the PROBE's value, not filled with a good one and patched afterwards
    # (which silently leaves a valid ini and "passes" by accepting nothing).
    probe_ini() {  # <name> <sdir> [post-sed]
        sed -e "s|@SDIR@|$2|" -e "s|@PREFIX@|kg|" -e "s|@STEPS@|1|" kasagi_uniform.ini \
            | sed "${3:-}" > ".kg_$1.ini"
    }
    probe_ini uniform-source z 's|^source_type = velocity|source_type = uniform|'
    probe_ini bad-direction q
    for name in uniform-source bad-direction; do
        if mpirun -n 1 "$BIN" ".kg_$name.ini" > ".kg_$name.log" 2>&1; then
            echo "   $name: ACCEPTED -- it must be a hard config error"; report 1
        else
            echo "   $name: rejected -- $(grep -m1 'ERROR STOP' ".kg_$name.log")"; report 0
        fi
    done
fi

# --- F5: the annulus, through moby_prepare --------------------------------
# pipe_case <tag> <axis> <outer-flag> <outer-value>
# The grid follows the axis: 8 cells along it over 0.2 (periodic), 64 across
# the 1.3-wide cross-section. Built by substituting the FINAL values per
# direction rather than by permuting placeholders -- sed applies expressions
# in order, so a later expression matching "@N@" never fires once an earlier
# one has already written the number. That is not hypothetical: it is how the
# first version of this function ran the x case on the z case's grid.
pipe_case() {
    local tag=$1 axis=$2 oflag=$3 oval=$4
    local n=(64 64 64) l=(1.3 1.3 1.3) per=(false false false)
    local d; d=$(case $axis in x) echo 0;; y) echo 1;; z) echo 2;; esac)
    n[$d]=8; l[$d]=0.2; per[$d]=true
    $PY ./make_geometry_stl.py annulus ".$tag.stl" --axis "$axis" --centre 0.65 0.65 \
        --r-inner 0.5 "$oflag" "$oval" --facets 16384 --a0 -0.8 --a1 1.0 \
        --domain-half 0.65 || return 1
    sed -e "s|@STL@|.$tag.stl|" -e "s|@CASE@|$tag.h5|" -e "s|@NB@|8|" \
        -e "s|^nx = @N@|nx = ${n[0]}|"  -e "s|^ny = @N@|ny = ${n[1]}|" \
        -e "s|^nz = @NZ@|nz = ${n[2]}|" -e "s|^lx = @L@|lx = ${l[0]}|" \
        -e "s|^ly = @L@|ly = ${l[1]}|"  -e "s|^lz = @LZ@|lz = ${l[2]}|" \
        -e "s|^periodic_x = .*|periodic_x = ${per[0]}|" \
        -e "s|^periodic_y = .*|periodic_y = ${per[1]}|" \
        -e "s|^periodic_z = .*|periodic_z = ${per[2]}|" \
        pipe_geom.ini | sed '/^coeff_file/d' > ".$tag.prep.ini"
    if grep -v '^;' ".$tag.prep.ini" | grep -q '@'; then
        echo "   unsubstituted placeholder in .$tag.prep.ini"; return 1
    fi
    mpirun -n "$RANKS" "$PREP" ".$tag.prep.ini" "$tag.h5" > "$tag.prep.log" 2>&1
}

if want annulus; then
    echo "== F5 (1) annulus about z, square outer surface"
    pipe_case pipe64 z --box-half 1.25 || report 1
    run $PY ./check_annulus.py pipe64.h5 --centre 0.65 0.65 --r-inner 0.5 \
        --facets 16384 --box-half 1.25 --domain-half 0.65
    report $?

    echo "== F5 (2) annulus about x, cylindrical outer surface"
    pipe_case pipex x --r-outer 1.6 || report 1
    run $PY ./check_annulus.py pipex.h5 --axis x --centre 0.65 0.65 --r-inner 0.5 \
        --facets 16384 --r-outer 1.6 --domain-half 0.65
    report $?
fi

if want ranks; then
    echo "== F5 (3) prepared case file is independent of the rank count"
    [ -f pipe64.h5 ] || pipe_case pipe64 z --box-half 1.25 || report 1
    mpirun -n 1 "$PREP" .pipe64.prep.ini pipe64_np1.h5 > pipe64_np1.prep.log 2>&1 || report 1
    run $PY ../prepare/h5same.py pipe64.h5 pipe64_np1.h5
    report $?
fi

echo
if [ $status -eq 0 ]; then echo "pipe prerequisite gates (F1, F5): ALL PASS"
else echo "pipe prerequisite gates (F1, F5): FAILURES"; fi
exit $status
