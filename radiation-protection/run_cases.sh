#!/bin/bash
# Run one or more viewer exports. Each positional JSON supplies its complete setup.
set -eu
here="$(cd "$(dirname "$0")" && pwd)"
if [ -n "${GATE_PYTHON:-}" ]; then
    python_bin="$GATE_PYTHON"
elif [ -x "$here/.venv-gate/bin/python" ]; then
    python_bin="$here/.venv-gate/bin/python"
elif [ -x "$here/../.venv-gate/bin/python" ]; then
    python_bin="$here/../.venv-gate/bin/python"
elif [ -n "${VIRTUAL_ENV:-}" ] && [ -x "$VIRTUAL_ENV/bin/python" ]; then
    python_bin="$VIRTUAL_ENV/bin/python"
else
    python_bin=python3
fi
scanner="${SCANNER:-$here/Revelation_Scanner_Inputs/gate_mammography}"
if [ "$#" -gt 0 ] && { [ -f "$1" ] || [[ "$1" == *.[jJ][sS][oO][nN] ]]; }; then
    # Explicit paths always win over preset-related environment variables.
    for config in "$@"; do
        if [ ! -f "$config" ]; then
            echo "Cannot open the selected settings file: $config" >&2
            exit 2
        fi
    done
    for config in "$@"; do
        "$python_bin" "$here/study.py" run --config "$config" --saved-export \
            --scanner "$scanner" --out "${OUTDIR:-$here/results}" \
            --histories "${HISTORIES:-10000}" --batches "${BATCHES:-1}" --jobs "${JOBS:-1}"
    done
    exit 0
fi
if [ -z "${CONFIG:-}" ]; then
    echo 'Usage: bash run_cases.sh "/path/to/your saved viewer.json" [another.json ...]' >&2
    echo 'On macOS, drag the exported JSON from Finder into Terminal after the command.' >&2
    exit 2
fi
config="$CONFIG"
if [ "$#" -eq 0 ] && [ -z "${LAYOUT:-}" ] && [ -z "${PROGRAMME:-}" ]; then
    exec "$python_bin" "$here/study.py" run --config "$config" --saved-export \
        --scanner "$scanner" --out "${OUTDIR:-$here/results}" \
        --histories "${HISTORIES:-10000}" --batches "${BATCHES:-1}" --jobs "${JOBS:-1}"
fi
programme="${PROGRAMME:-viewer}"
layout_args=()
if [ -n "${LAYOUT:-}" ]; then layout_args=(--layouts "$LAYOUT"); fi
study_name="$(basename "$config" .json)"
destination="${OUTDIR:-$here/results/$study_name/$programme}"
if [ "$#" -eq 0 ]; then
    saved_phantom="$("$python_bin" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("revelation_viewer", {}).get("phantom", "L2"))' "$config")"
    set -- "$saved_phantom"
fi
for phantom in "$@"
do
    echo "Running $phantom with $config and $programme"
    "$python_bin" "$here/study.py" run \
        --config "$config" --scanner "$scanner" --phantom "$phantom" ${layout_args[@]+"${layout_args[@]}"} \
        --programme "$programme" --histories "${HISTORIES:-10000}" \
        --batches "${BATCHES:-1}" --jobs "${JOBS:-1}" --out "$destination/$phantom"
done
