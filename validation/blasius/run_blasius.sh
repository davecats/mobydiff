#!/usr/bin/env bash
# Laminar zero-pressure-gradient Blasius boundary layer -- the laminar
# STRETCHED-LINE gate (one-sided natural y line, Blasius inlet, two
# Dirichlet-p outlets). Mints the restart template, fills it with the
# analytic Blasius field, holds it for 4000 steps (t = 2000 at the dt cap 0.5)
# and compares against the independent Blasius ODE.
#
#     [SOLVER=../../build_gpu/moby_solve] [RANKS=1] ./run_blasius.sh
#
# The template ini is DERIVED from blasius2d.ini here, never edited by hand:
# a restart file carries its boundary rows (bc_type/bc_value) and the ini
# overrides only the rows it sets, so a template minted from a sibling
# variant silently transplants that variant's faces (the outlet-top case
# restarted with a Dirichlet top from another ini diverged by t ~ 14;
# docs/next_session_after_step8.md). Prepare, mint and run share one case
# file and one rank count.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)
SOLVER=${SOLVER:-$ROOT/build_gpu/moby_solve}
RANKS=${RANKS:-1}
[ -x "$SOLVER" ] || { echo "no solver at $SOLVER"; exit 1; }

sed -e 's/^nsteps *=.*/nsteps = 1/' \
    -e 's/^field_interval *=.*/field_interval = 0/' \
    -e 's/^field_prefix *=.*/field_prefix = template/' \
    -e '/^\[restart\]/,/^\[output\]/{/^\[output\]/!d}' blasius2d.ini > template.ini
grep -q 'IC_blasius' template.ini && { echo "template.ini still restarts"; exit 1; }

"$ROOT/tools/prepare_if_missing.sh" "$RANKS" "$SOLVER" blasius2d.ini prepare.log || { echo "PREPARE FAILED (prepare.log)"; exit 1; }
rm -f template_*.h5 blasius2d_*.h5
mpirun -n "$RANKS" "$SOLVER" template.ini > template.log 2>&1 || { echo "TEMPLATE MINT FAILED (template.log)"; exit 1; }
python3 make_blasius_ic.py > make_ic.log 2>&1 || { echo "IC FAILED (make_ic.log)"; exit 1; }
mpirun -n "$RANKS" "$SOLVER" blasius2d.ini > blasius2d.log 2>&1 || { echo "RUN FAILED (blasius2d.log)"; exit 1; }
last=$(ls -t blasius2d_*.h5 | head -1)
python3 compare_blasius.py "$last" --plot blasius.png | tee compare.log
