#!/usr/bin/env bash
# The gate drivers' prepare step (numerics review step 7: preparing and
# solving are separate executables, the solver refuses a missing case file).
#
#     tools/prepare_if_missing.sh <ranks> <solver-binary> <input.ini> [log]
#
# Resolves the case file the solver will read from the ini -- [case] file,
# else <field_prefix>.case.h5 (the same rule as config.f90 case_file_name) --
# and, when that file does NOT exist, runs
# the moby_prepare that sits next to the solver binary on <ranks> ranks (the
# same rank count as the solve: an unset [blocks] nb is one block per prepare
# rank). An existing file is left alone, whatever produced it: the committed
# and generated reference files the inis point at must never be recomputed
# here, and a stale one is refused by the solver itself with the key named.
# A solver directory without a moby_prepare (a pre-step-7 reference binary,
# which still builds inline) is a no-op too.
#
# Exit 0 when nothing was needed or the prepare succeeded, else the prepare's
# exit code (its output is in [log], default <ini>.prep.log).
set -u
ranks=$1 solver=$2 ini=$3 log=${4:-$3.prep.log}

prep="$(dirname "$solver")/moby_prepare"
[ -x "$prep" ] || exit 0

case_of() {  # the case file the solver will read (fallback transcription of case_file_name)
    local f
    f=$(awk -F'=' '/^[ \t]*\[/{s=$0; gsub(/[ \t]/,"",s)}
                   s=="[case]" && $1 ~ /^[ \t]*file[ \t]*$/ {v=$2; sub(/[;#].*/,"",v); gsub(/^[ \t"]+|[ \t"]+$/,"",v); print v; exit}' "$1")
    if [ -z "$f" ]; then
        local pfx
        pfx=$(awk -F'=' '$1 ~ /^[ \t]*field_prefix[ \t]*$/ {v=$2; sub(/[;#].*/,"",v); gsub(/^[ \t"]+|[ \t"]+$/,"",v); print v; exit}' "$1")
        f="${pfx}.case.h5"   # an unset prefix is "" in config.f90 too
    fi
    printf '%s\n' "$f"
}

# The binary resolves the name (one rule, config.f90 case_file_name); a
# moby_prepare without the flag (before 2026-10-21) prints nothing usable, and
# then the shell transcription above answers.
case_file=$(mpirun -n 1 "$prep" --case-file-name "$ini" 2>/dev/null | tail -n 1)
case "$case_file" in *.h5) ;; *) case_file=$(case_of "$ini") ;; esac
[ -f "$case_file" ] && exit 0
mpirun -n "$ranks" "$prep" "$ini" > "$log" 2>&1
