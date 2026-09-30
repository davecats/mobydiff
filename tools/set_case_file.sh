#!/usr/bin/env bash
# Point an ini at a case file, IN PLACE:
#
#     tools/set_case_file.sh <input.ini> <case.h5>
#
# writes `file = <case.h5>` into the [case] section, adding the section at the
# top when the ini has none (the [case] file key is the one the solver and
# moby_prepare resolve the case file from -- config.f90 case_file_name). The
# gate drivers use it to derive a solve ini from a case ini without a second
# hand-written copy; an existing `file =` line in [case] is replaced.
set -u
ini=$1 case=$2
if grep -q '^\[case\]' "$ini"; then
    # drop a previous [case] file line (only inside that section), then insert
    awk -v f="$case" '
        /^[ \t]*\[/ { insec = ($0 ~ /^[ \t]*\[case\][ \t]*$/) }
        insec && $0 ~ /^[ \t]*file[ \t]*=/ { next }
        { print }
        insec && $0 ~ /^[ \t]*\[case\][ \t]*$/ { print "file = " f }
    ' "$ini" > "$ini.tmp$$" && mv "$ini.tmp$$" "$ini"
else
    printf '[case]\nfile = %s\n\n' "$case" | cat - "$ini" > "$ini.tmp$$" && mv "$ini.tmp$$" "$ini"
fi
