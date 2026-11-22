#!/bin/bash
# Thin wrapper: analyze one replica through bin/amber-kit.
# Usage: bash run_analysis_one.sh LABEL REP [--project DIR] [--rna] [--force]
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
LABEL=${1:-}; REP=${2:-}
[ -n "$LABEL" ] && [ -n "$REP" ] || { echo "usage: $0 LABEL REP [--project DIR] [--rna] [--force]" >&2; exit 2; }
shift 2
exec "$KIT_DIR/bin/amber-kit" analyze --label "$LABEL" --rep "$REP" "$@"
