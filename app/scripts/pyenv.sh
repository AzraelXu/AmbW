#!/bin/bash
# Pick a Python that can import parmed + numpy. Source this file, then use
#   "$AMBER_MD_PY" $AMBER_MD_PY_FLAGS script.py ...
# or call amber_md_run_py script.py ...
#
# AMBER_MD_PYTHON overrides the search. The function prefers an interpreter
# that works under -I (isolated), then one that works with the current
# environment (Amber's amber.sh puts parmed on PYTHONPATH).
#
# This matters on machines where Amber ships parmed for a different Python
# than the system python3: `python3 -I` cannot see PYTHONPATH, so the isolated
# check fails even though `python3 script.py` works.

amber_md_pick_python() {
    local cand
    for cand in "${AMBER_MD_PYTHON:-}" python3 python; do
        [ -n "$cand" ] || continue
        command -v "$cand" >/dev/null 2>&1 || continue
        if "$cand" -I -c 'import parmed, numpy' >/dev/null 2>&1; then
            AMBER_MD_PY=$cand
            AMBER_MD_PY_FLAGS="-I"
            export AMBER_MD_PY AMBER_MD_PY_FLAGS
            return 0
        fi
        if "$cand" -c 'import parmed, numpy' >/dev/null 2>&1; then
            AMBER_MD_PY=$cand
            AMBER_MD_PY_FLAGS=""
            export AMBER_MD_PY AMBER_MD_PY_FLAGS
            return 0
        fi
    done
    AMBER_MD_PY=python3
    AMBER_MD_PY_FLAGS=""
    export AMBER_MD_PY AMBER_MD_PY_FLAGS
    return 1
}

amber_md_run_py() {
    # shellcheck disable=SC2086
    "$AMBER_MD_PY" $AMBER_MD_PY_FLAGS "$@"
}

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
    if amber_md_pick_python; then
        echo "$AMBER_MD_PY $AMBER_MD_PY_FLAGS"
    else
        echo "no Python with parmed+numpy found; falling back to python3" >&2
        echo "$AMBER_MD_PY $AMBER_MD_PY_FLAGS"
        exit 1
    fi
fi
