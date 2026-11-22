#!/bin/bash
# Thin wrapper: run MM/GBSA or MM/PBSA for one replica through bin/amber-kit.
# Usage: bash run_mmpbsa_one.sh LABEL REP [--project DIR] [--method gb|pb] [--force]
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
LABEL=${1:-}; REP=${2:-}
[ -n "$LABEL" ] && [ -n "$REP" ] || { echo "usage: $0 LABEL REP [--project DIR] [--method gb|pb] [--force]" >&2; exit 2; }
shift 2
exec "$KIT_DIR/bin/amber-kit" mmpbsa --label "$LABEL" --rep "$REP" "$@"
