#!/usr/bin/env bash
# Run the block-overhead timing matrix. Each config gets its own directory under
# runs/ (gitignored) holding runtime.txt and run.log; the unconditional final
# snapshot is deleted afterwards (it is written outside the timed region, but it
# is 2-4 GB per run).
#
#   ./run_overhead.sh                      # all five configs
#   ./run_overhead.sh singleLevel/base_jacobi.ini ...
#
# Environment:
#   BIN                  solver binary       (default ../../../build_gpu/moby_solve)
#   NRANKS               mpirun ranks        (default 1)
#   CUDA_VISIBLE_DEVICES which GPU           (default 1 -- istmcetus GPU 1; GPU 0
#                                             and corax are the production campaign)
#
# On istmcetus: check nvidia-smi first, the machine is shared. Absolute times are
# ~3.1x corax's; only RATIOS transfer between machines.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
BIN="${BIN:-$ROOT/build_gpu/moby_solve}"
NRANKS="${NRANKS:-1}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-1}"

if [ ! -x "$BIN" ]; then
    echo "solver binary not found or not executable: $BIN" >&2
    exit 1
fi

if [ "$#" -gt 0 ]; then
    configs=("$@")
else
    configs=(
        singleLevel/base_redblack.ini
        singleLevel/nb16_redblack.ini
        singleLevel/base_jacobi.ini
        singleLevel/nb16_jacobi.ini
        multiLevel_xz/refined_yp100_jacobi.ini
    )
fi

for cfg in "${configs[@]}"; do
    abs="$(cd "$(dirname "$HERE/$cfg")" && pwd)/$(basename "$cfg")"
    [ -f "$abs" ] || { echo "no such config: $cfg" >&2; exit 1; }
    name="$(basename "$cfg" .ini)"
    run="$HERE/runs/$name"

    echo "=== $name  ($(date '+%F %T'))"
    rm -rf "$run"
    mkdir -p "$run"
    ( cd "$run" && mpirun -n "$NRANKS" "$BIN" "$abs" 2>&1 | tee run.log )
    rm -f "$run"/overhead_*.h5
done

echo
"$HERE/summarise.py"
