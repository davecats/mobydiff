#!/usr/bin/env bash
# IDDES T5 (DDES-shielding) gate runner — docs/next_session_iddes.md, phase T5.
# Local-friendly (the channel legs run in minutes-hours on a GPU; CPU works).
#
#   ./run_gates.sh              # all gates, sequentially (one job at a time)
#   ./run_gates.sh iddes180     # one gate group
#   groups: iddes180 (transient+stats legs), fd0 (fd_force=0 == pure WALE,
#           bit-exact), fd1 (fd_force=1 vs the T2 RANS turb180 fixed point),
#           ibm (les_ibm channel stability), ranks (1==4 rank determinism)
#
# Environment:  BIN=<solver>  (default ../../build_gpu/main)
#               RANKS=<n>     (default 1; the ranks gate manages its own)
# Then: python3 check_gates.py
set -o pipefail
cd "$(dirname "$0")"
set -u

BIN=${BIN:-../../build_gpu/main}
RANKS=${RANKS:-1}
sel=${1:-all}
status=0

run() {  # ini ranks logname
    local ini=$1 ranks=$2 log=$3
    echo "== $log (ranks $ranks) =="
    if ! mpirun -n "$ranks" "$BIN" "$ini" > "$log.log" 2>&1; then
        echo "   FAILED — see $log.log"; status=1; return 1
    fi
    tail -n 2 "$log.log" | sed 's/^/   /'
}

want() { [ "$sel" = all ] || [ "$sel" = "$1" ]; }

# --- gates (a)/(b): the WMLES channel, transient then statistics leg ---
if want iddes180; then
    sed -e 's/^stats_sample_interval.*/stats_sample_interval = 0/' \
        -e 's/^t_final.*/t_final = 5.0/' \
        -e 's/^field_prefix.*/field_prefix = iddes180_dev/' \
        -e 's/^field_interval.*/field_interval = 0/' iddes180.ini > .dev.ini
    run .dev.ini "$RANKS" iddes180_dev && {
        devout=$(ls -t iddes180_dev_*.h5 | head -1)
        sed -e "s|^file = .*|file = ${devout}|" iddes180.ini > .stats.ini
        run .stats.ini "$RANKS" iddes180
    }
    rm -f .dev.ini .stats.ini
fi

# --- gate (c), fd = 0 limit: bit-exact vs pure WALE on the same IC ---
if want fd0; then
    sed -e 's/^nsteps.*/nsteps = 20/' -e 's/^t_final.*/t_final = 0.0/' \
        -e 's/^field_interval.*/field_interval = 0/' \
        -e 's/^stats_sample_interval.*/stats_sample_interval = 0/' \
        -e 's/^field_prefix.*/field_prefix = fd0_iddes/' iddes180.ini > .fd0.ini
    printf '\n[turbulence]\nfd_force = 0.0\n' >> .fd0.ini
    run .fd0.ini "$RANKS" fd0_iddes
    sed -e 's/^nsteps.*/nsteps = 20/' -e 's/^t_final.*/t_final = 0.0/' \
        -e 's/^field_interval.*/field_interval = 0/' \
        -e 's/^stats_sample_interval.*/stats_sample_interval = 0/' \
        -e 's/^field_prefix.*/field_prefix = fd0_wale/' \
        -e 's/^model = iddes/model = les/' iddes180.ini > .wale.ini
    run .wale.ini "$RANKS" fd0_wale
    rm -f .fd0.ini .wale.ini
fi

# --- gate (c), fd = 1 limit: hold the T2 RANS turb180 fixed point ---
if want fd1; then
    if [ ! -f ../rans_sst/turb180_132565.h5 ]; then
        echo "== fd1 SKIPPED: ../rans_sst/turb180_132565.h5 missing"; status=1
    else
        sed -e 's/^model = rans/model = iddes/' \
            -e 's/^nsteps.*/nsteps = 2000/' -e 's/^t_final.*/t_final = 0.0/' \
            -e 's/^field_prefix.*/field_prefix = fd1_rans/' ../rans_sst/turb180.ini > .fd1.ini
        printf '\n[turbulence]\nfd_force = 1.0\n[les]\nmodel = wale\n' >> .fd1.ini
        printf '[restart]\nfile = ../rans_sst/turb180_132565.h5\n' >> .fd1.ini
        run .fd1.ini "$RANKS" fd1_rans
        rm -f .fd1.ini
    fi
fi

# --- gate (d): the les_ibm IBM channel runs IDDES stably ---
if want ibm; then
    if [ ! -f ../rans_geometry/ibm_coeff_blocks_l1.h5 ]; then
        echo "== ibm SKIPPED: ../rans_geometry/ibm_coeff_blocks_l1.h5 missing"; status=1
    else
        sed -e 's/^nsteps.*/nsteps = 2000/' -e 's/^t_final.*/t_final = 0.0/' \
            iddes_ibm.ini > .ibm.ini
        run .ibm.ini "$RANKS" iddes_ibm
        rm -f .ibm.ini
    fi
fi

# --- gate (e) iddes-on determinism: 1 rank == 4 ranks, 20 steps ---
if want ranks; then
    for r in 1 4; do
        sed -e 's/^nsteps.*/nsteps = 20/' -e 's/^t_final.*/t_final = 0.0/' \
            -e 's/^field_interval.*/field_interval = 0/' \
            -e 's/^stats_sample_interval.*/stats_sample_interval = 0/' \
            -e "s/^field_prefix.*/field_prefix = ranks${r}/" iddes180.ini > .r$r.ini
        run .r$r.ini "$r" ranks$r
        rm -f .r$r.ini
    done
fi

echo
echo "runs done (status $status); now: python3 check_gates.py"
exit $status
