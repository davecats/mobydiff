#!/usr/bin/env bash
# A3 INCREMENT 2: NACA 0012 aoa sweep (0/4/8 deg), full-turbulent SST,
# sequential GPU runs (one solver job at a time).
#
#   ./run_sweep.sh            # all three angles
#   ./run_sweep.sh 4          # one angle
#
# Environment:  BIN=<solver>  (default ../../build_gpu/main, fall back cpu)
set -o pipefail
cd "$(dirname "$0")"
set -u

if [ -z "${BIN:-}" ]; then
    BIN=../../build_gpu/main
    [ -x "$BIN" ] || BIN=../../build_gpu_nofma/main
fi
angles=${@:-0 4 8}
status=0

for a in $angles; do
    sed -e "s/^aoa = .*/aoa = ${a}.0/" \
        -e "s/^runtime_file = .*/runtime_file = forces_aoa${a}.txt/" \
        -e "s/^field_prefix = .*/field_prefix = n0012_aoa${a}/" \
        aoa4.ini > .aoa${a}.ini
    rm -f forces_aoa${a}.txt
    echo "== aoa = ${a} =="
    if ! mpirun -n 1 "$BIN" .aoa${a}.ini > aoa${a}.log 2>&1; then
        echo "   FAILED — see aoa${a}.log"; status=1
    else
        tail -n 2 aoa${a}.log | sed 's/^/   /'
        tail -n 1 forces_aoa${a}.txt | sed 's/^/   forces: /'
    fi
    rm -f .aoa${a}.ini
done

python3 check_naca.py $(for a in $angles; do echo forces_aoa${a}.txt; done) || status=1
exit $status
