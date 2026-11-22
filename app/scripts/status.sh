#!/bin/bash
# Generic Amber MD progress table for a project tree.
#
# Scans directories matching ROOT/--md-glob, finds *.in stages, and reports the
# last stage with an .out file, NSTEP/TIME(PS), ns/day, and completion flags.
# It does not assume a project name, job prefix, or fixed stage list.
#
# Usage:
#   bash status.sh ROOT [--md-glob 'md/*/r*'] [--json] [--queue]
set -u

ROOT=""
MD_GLOB="md/*/r*"
JSON=0
QUEUE=0

while [ $# -gt 0 ]; do
  case "$1" in
    --md-glob) MD_GLOB=$2; shift 2;;
    --json) JSON=1; shift;;
    --queue) QUEUE=1; shift;;
    -h|--help) sed -n '2,11p' "$0"; exit 0;;
    -*) echo "error: unknown option $1" >&2; exit 2;;
    *) ROOT=$1; shift;;
  esac
done

[ -n "$ROOT" ] || { echo "error: ROOT is required" >&2; exit 2; }
[ -d "$ROOT" ] || { echo "error: not a directory: $ROOT" >&2; exit 2; }

json_escape() {
  printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'
}

stage_list() {
  local d=$1
  find "$d" -maxdepth 1 -name '*.in' -printf '%f\n' 2>/dev/null \
    | sed 's/\.in$//' | sort -V
}

nstlim_of() {
  local f=$1
  grep -oE 'nstlim[[:space:]]*=[[:space:]]*[0-9]+' "$f" 2>/dev/null \
    | tail -1 | grep -oE '[0-9]+' || true
}

rows=()
queue_json=""
if [ "$QUEUE" -eq 1 ]; then
  if command -v qstat >/dev/null 2>&1; then
    read -r q_total q_run <<< "$(timeout 10 qstat -a 2>/dev/null | awk '
      NR>5 && $1 ~ /^[0-9]+/ {total++; if ($10 == "R") run++}
      END {print total+0, run+0}')"
  else
    q_total="NA"; q_run="NA"
  fi
  queue_json="{\"total\": \"$q_total\", \"running\": \"$q_run\"}"
fi

found=0
for d in $ROOT/$MD_GLOB; do
  [ -d "$d" ] || continue
  found=$((found + 1))
  label=$(basename "$(dirname "$d")")
  rep=$(basename "$d")
  replica="$label/$rep"

  last_stage="-"; progress="-"; nsday="-"; note=""
  stages=$(stage_list "$d")
  for stem in $stages; do
    if [ -f "$d/$stem.out" ]; then last_stage=$stem; fi
  done

  add_note() {
    if [ -z "$note" ]; then note=$1; else note="$note,$1"; fi
  }

  if [ "$last_stage" != "-" ]; then
    out="$d/$last_stage.out"
    # The last NSTEP line is near the end; tail avoids scanning large .out files.
    line=$(tail -n 2000 "$out" 2>/dev/null | grep -E 'NSTEP[[:space:]]*=' | tail -1 || true)
    val=""
    if [ -n "$line" ]; then
      val=$(awk '{for (i=1;i<=NF;i++) if ($i=="=" && $(i-1)=="NSTEP") {print $(i+1); exit}}' <<< "$line")
      if [ -n "$val" ]; then progress=$val; fi
    fi
    if [ -f "$d/$last_stage.mdinfo" ]; then
      ns=$(grep 'ns/day' "$d/$last_stage.mdinfo" 2>/dev/null | tail -1 | awk '{print $4}' || true)
      if [ -n "$ns" ]; then nsday=$ns; fi
    fi
    # Real fatal markers start a line; "TESTING RELATIVE ERROR" does not.
    if tail -n 300 "$out" 2>/dev/null | grep -qE '^[[:space:]]*(ERROR|FATAL)|(^|[^A-Za-z])NaN([^A-Za-z]|$)'; then
      add_note "ERROR_IN_OUTPUT"
    fi
    if [ ! -f "$d/$last_stage.rst7" ] && [ ! -f "$d/$last_stage.nc" ]; then
      add_note "NO_RESTART_OR_TRAJ"
    fi
    limit=$(nstlim_of "$d/$last_stage.in")
    if [ -n "$limit" ] && [ -n "$val" ] && [ "$val" -ge "$limit" ] 2>/dev/null; then
      add_note "STAGE_COMPLETE"
    fi
  fi

  missing=""
  for stem in $stages; do
    if [ -f "$d/$stem.in" ] && [ ! -f "$d/$stem.nc" ] && [[ "$stem" == 05_prod* ]]; then
      missing="$missing $stem.nc"
    fi
  done

  if [ -f "$d/prod_final.rst7" ]; then
    add_note "PROD_DONE"
  elif [ -f "$d/npt.rst7" ]; then
    add_note "NPT_DONE"
  fi
  if [ -n "$missing" ]; then
    add_note "MISSING:${missing# }"
  fi
  [ -n "$note" ] || note="-"

  rows+=("$replica"$'\t'"$last_stage"$'\t'"$progress"$'\t'"$nsday"$'\t'"$note")
done

if [ "$found" -eq 0 ]; then
  echo "no replica directories matched: $ROOT/$MD_GLOB" >&2
  exit 1
fi

if [ "$JSON" -eq 1 ]; then
  printf '{"replicas": [\n'
  first=1
  for row in "${rows[@]}"; do
    IFS=$'\t' read -r replica stage progress nsday note <<< "$row"
    [ "$first" -eq 1 ] || printf ',\n'
    first=0
    printf '  {"replica": "%s", "stage": "%s", "nstep": "%s", "ns_day": "%s", "note": "%s"}' \
      "$(json_escape "$replica")" "$(json_escape "$stage")" \
      "$(json_escape "$progress")" "$(json_escape "$nsday")" "$(json_escape "$note")"
  done
  printf '\n]'
  [ "$QUEUE" -eq 1 ] && printf ', "queue": %s' "$queue_json"
  printf '}\n'
  exit 0
fi

printf "%-22s %-18s %-14s %-8s %s\n" REPLICA STAGE NSTEP NS_DAY NOTE
for row in "${rows[@]}"; do
  IFS=$'\t' read -r replica stage progress nsday note <<< "$row"
  printf "%-22s %-18s %-14s %-8s %s\n" "$replica" "$stage" "$progress" "$nsday" "$note"
done
if [ "$QUEUE" -eq 1 ]; then
  echo "--- queue: total=$q_total running=$q_run ---"
fi
