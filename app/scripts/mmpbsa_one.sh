#!/bin/bash
# Single-trajectory MM/GBSA or MM/PBSA for one replica.
#
# Splits a SOLVATED topology with ante-MMPBSA.py, validates the split, selects
# a frame window, renders an MMPBSA input, and parses DELTA TOTAL.
#
# Usage:
#   bash mmpbsa_one.sh --top solvated.prmtop --traj dry.nc \
#     --ligand-mask ':27' [--receptor-mask ':1-26'] \
#     --strip-mask ':WAT,K+,Cl-,Na+' \
#     --radii mbondi2 --method gb --igb 5 --saltcon 0.12 \
#     [--start-frame auto|INT] [--end-frame auto|INT] [--interval 100] \
#     [--decomp-residues '1-26,27'] [--istrng 0.12] [--dry-top] \
#     [--outdir DIR] [--check-only] [--keep-files] [--force]
set -euo pipefail

KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TOP=""; TRAJ=""; LIGAND_MASK=""; STRIP_MASK=":WAT,Na+,Cl-,K+,MG,CA,ZN"
RADII="mbondi2"; METHOD="gb"; IGB="5"; SALTCON="0.12"; ISTRNG="0.12"
START="auto"; END="auto"; INTERVAL="100"; DECOMP=""; OUTDIR="mmpbsa"
RECEPTOR_MASK=""
DRY_TOP=0; CHECK_ONLY=0; KEEP_FILES=0; FORCE=0

die() { echo "error: $*" >&2; exit 2; }
note() { echo "[mmpbsa] $*" >&2; }

while [ $# -gt 0 ]; do
  case "$1" in
    --top) TOP=$2; shift 2;;
    --traj) TRAJ=$2; shift 2;;
    --ligand-mask) LIGAND_MASK=$2; shift 2;;
    --receptor-mask) RECEPTOR_MASK=$2; shift 2;;
    --strip-mask) STRIP_MASK=$2; shift 2;;
    --radii) RADII=$2; shift 2;;
    --method) METHOD=$2; shift 2;;
    --igb) IGB=$2; shift 2;;
    --saltcon) SALTCON=$2; shift 2;;
    --istrng) ISTRNG=$2; shift 2;;
    --start-frame) START=$2; shift 2;;
    --end-frame) END=$2; shift 2;;
    --interval) INTERVAL=$2; shift 2;;
    --decomp-residues) DECOMP=$2; shift 2;;
    --outdir) OUTDIR=$2; shift 2;;
    --dry-top) DRY_TOP=1; shift;;
    --check-only) CHECK_ONLY=1; shift;;
    --keep-files) KEEP_FILES=1; shift;;
    --force) FORCE=1; shift;;
    -h|--help) sed -n '2,18p' "$0"; exit 0;;
    *) die "unknown argument: $1";;
  esac
done

[ -n "$TOP" ] || die "--top is required (use the solvated topology)"
[ -n "$TRAJ" ] || die "--traj is required (a single dry trajectory)"
[ -n "$LIGAND_MASK" ] || die "--ligand-mask is required, e.g. ':27' or ':LIG'"
[ -f "$TOP" ] || die "topology not found: $TOP"
[ -f "$TRAJ" ] || die "trajectory not found: $TRAJ"
[ "$METHOD" = gb ] || [ "$METHOD" = pb ] || die "--method must be gb or pb"
[ "$RADII" = mbondi2 ] || [ "$RADII" = mbondi3 ] || [ "$RADII" = mbondi ] || \
  die "--radii must be mbondi, mbondi2, or mbondi3"
[[ "$INTERVAL" =~ ^[1-9][0-9]*$ ]] || die "--interval must be a positive integer"
if [ "$METHOD" = gb ]; then
  case "$IGB" in 1|2|5|7|8) ;; *) die "--igb must be 1, 2, 5, 7, 8";; esac
  if [ "$IGB" = 8 ] && [ "$RADII" = mbondi2 ]; then
    note "WARNING: igb=8 normally requires mbondi3, not mbondi2"
  fi
  if [ "$IGB" = 5 ] && [ "$RADII" = mbondi3 ]; then
    note "WARNING: igb=5 normally requires mbondi2, not mbondi3"
  fi
fi
if ! command -v ante-MMPBSA.py >/dev/null 2>&1 || ! command -v MMPBSA.py >/dev/null 2>&1; then
  if [ -n "${AMBER_MD_AMBER_SH:-}" ] && [ -f "${AMBER_MD_AMBER_SH}" ]; then
    # amber.sh 读取可能未导出的 LD_LIBRARY_PATH/PERL5LIB/PYTHONPATH；
    # 在 set -u 下 source 会直接退出，因此 source 期间临时关闭 nounset。
    # shellcheck disable=SC1090
    set +u
    source "$AMBER_MD_AMBER_SH"
    set -u
  fi
fi
command -v ante-MMPBSA.py >/dev/null 2>&1 || \
  die "ante-MMPBSA.py not in PATH; source \$AMBERHOME/amber.sh"
command -v MMPBSA.py >/dev/null 2>&1 || \
  die "MMPBSA.py not in PATH; source \$AMBERHOME/amber.sh"
command -v cpptraj >/dev/null 2>&1 || die "cpptraj not in PATH"

mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)
TOP=$(cd "$(dirname "$TOP")" && pwd)/$(basename "$TOP")
TRAJ=$(cd "$(dirname "$TRAJ")" && pwd)/$(basename "$TRAJ")

MANAGED=("$OUTDIR/complex.prmtop" "$OUTDIR/receptor.prmtop" "$OUTDIR/ligand.prmtop"
         "$OUTDIR/mmpbsa.in" "$OUTDIR/mmpbsa.dat" "$OUTDIR/mmpbsa.csv"
         "$OUTDIR/summary.tsv")
if [ "$FORCE" -eq 0 ]; then
  for f in "${MANAGED[@]}"; do
    if [ -e "$f" ]; then die "refusing to overwrite $f (use --force)"; fi
  done
else
  # ante-MMPBSA.py and MMPBSA.py refuse to overwrite; clear only our outputs.
  for f in "${MANAGED[@]}"; do rm -f "$f"; done
fi

residue_labels() {
  awk '/^%FLAG RESIDUE_LABEL/{f=1;next} /^%FLAG/{f=0} f{print}' "$1" | \
    fold -w4 | sed 's/ *$//' | grep -v '^$' || true
}

# Is the topology solvated? Check residue labels, not atom names (a protein
# alpha carbon is called CA, which would false-positive an ion grep).
HAS_WATER=0
if residue_labels "$TOP" | grep -qxE 'WAT|HOH|H2O|SOL|OPC|TP3|SPC|T4D|FB3'; then
  HAS_WATER=1
fi
if [ "$DRY_TOP" -eq 1 ]; then
  if [ "$HAS_WATER" -eq 1 ]; then
    die "--dry-top was given but the topology contains solvent; use the solvated topology instead"
  fi
else
  if [ "$HAS_WATER" -eq 0 ]; then
    die "topology looks dry (no water residue names). Use --dry-top only if you intentionally split a dry topology."
  fi
fi

top_atoms() {
  local top=$1 tmp=$OUTDIR/.parminfo.in
  printf 'parminfo\nquit\n' > "$tmp"
  cpptraj -p "$top" -i "$tmp" 2>/dev/null | \
    awk '/contains/ {print $4; exit}'
}
radii_of() {
  grep -A2 '%FLAG RADIUS_SET' "$1" | tail -1
}

# ---- split topology ------------------------------------------------------
ANTE=(ante-MMPBSA.py -p "$TOP" -c complex.prmtop -r receptor.prmtop \
      -l ligand.prmtop -n "$LIGAND_MASK" --radii "$RADII")
if [ "$DRY_TOP" -eq 0 ]; then
  ANTE+=(-s "$STRIP_MASK")
else
  # A no-match mask still makes ante-MMPBSA.py write complex.prmtop with the
  # requested radii; omitting -s leaves the complex topology unwritten.
  ANTE+=(-s '!*')
fi
( cd "$OUTDIR" && "${ANTE[@]}" ) > "$OUTDIR/ante.log" 2>&1 || {
  tail -40 "$OUTDIR/ante.log" >&2; die "ante-MMPBSA.py failed"; }
for f in complex.prmtop receptor.prmtop ligand.prmtop; do
  [ -s "$OUTDIR/$f" ] || die "ante-MMPBSA.py did not produce $f"
done

N_C=$(top_atoms "$OUTDIR/complex.prmtop")
N_R=$(top_atoms "$OUTDIR/receptor.prmtop")
N_L=$(top_atoms "$OUTDIR/ligand.prmtop")
[ -n "$N_C" ] && [ -n "$N_R" ] && [ -n "$N_L" ] || die "cannot read atom counts from split topologies"
[ "$N_C" -eq $((N_R + N_L)) ] || \
  die "split topology atom counts do not add up: complex $N_C, receptor $N_R, ligand $N_L"
if residue_labels "$OUTDIR/complex.prmtop" | grep -qxE 'WAT|HOH|H2O|SOL|OPC|TP3|SPC|T4D|FB3'; then
  die "complex.prmtop still contains solvent; check --strip-mask"
fi
if residue_labels "$OUTDIR/complex.prmtop" | grep -qxE 'Na\+|Cl-|K\+|MG|CA|ZN|LI|RB|CS'; then
  die "complex.prmtop still contains ions; add them to --strip-mask"
fi

for f in complex.prmtop receptor.prmtop ligand.prmtop; do
  R=$(radii_of "$OUTDIR/$f")
  printf '%s' "$R" | grep -qF "($RADII)" || \
    die "RADIUS_SET in $f is '$R', expected $RADII"
done
note "split: complex=$N_C receptor=$N_R ligand=$N_L atoms; radii=$RADII"

# ---- frame window --------------------------------------------------------
N=$(cpptraj -p "$OUTDIR/complex.prmtop" -y "$TRAJ" -tl 2>/dev/null | \
    awk '/^Frames:/{print $2; exit}')
[ -n "$N" ] || die "cpptraj cannot pair the dry trajectory with complex.prmtop; atom counts differ"
[ "$N" -gt 0 ] || die "trajectory has no frames"

if [ "$START" = auto ]; then START=$((N / 2 + 1)); fi
if [ "$END" = auto ]; then END=$N; fi
[[ "$START" =~ ^[0-9]+$ ]] && [[ "$END" =~ ^[0-9]+$ ]] || \
  die "--start-frame/--end-frame must be 'auto' or integers"
[ "$START" -ge 1 ] || die "--start-frame must be >= 1"
[ "$END" -le "$N" ] || die "--end-frame $END exceeds trajectory length $N"
[ "$START" -le "$END" ] || die "--start-frame must be <= --end-frame"
NFRAMES=$(( (END - START) / INTERVAL + 1 ))
if [ "$NFRAMES" -lt 2 ]; then
  note "WARNING: only $NFRAMES frame(s) selected; reduce --interval or widen the window"
fi

# ---- render input --------------------------------------------------------
DECOMP_BLOCK=""
if [ -n "$DECOMP" ]; then
  DECOMP_BLOCK="&decomp
  idecomp=1, print_res=\"$DECOMP\", dec_verbose=1, csv_format=1,
/"
fi
MASKS=""
if [ -n "$RECEPTOR_MASK" ]; then
  # Explicit masks remove MMPBSA.py's ambiguous sequence-based ligand
  # detection, which is essential when the partner is a protein/peptide.
  MASKS="  ligand_mask='$LIGAND_MASK', receptor_mask='$RECEPTOR_MASK',"
fi
TEMPLATE="$KIT_DIR/../templates/mmpbsa/mmpbsa_${METHOD}.in"
[ -f "$TEMPLATE" ] || die "MMPBSA template missing: $TEMPLATE"
RENDER=(python3 "$KIT_DIR/render_template.py" "$TEMPLATE" --out "$OUTDIR/mmpbsa.in" --force
  --set "START=$START" --set "END=$END" --set "INTERVAL=$INTERVAL"
  --set "KEEP_FILES=$KEEP_FILES" --set "DECOMP=$DECOMP_BLOCK"
  --set "MASKS=$MASKS")
if [ "$METHOD" = gb ]; then
  RENDER+=(--set "IGB=$IGB" --set "SALTCON=$SALTCON")
else
  RENDER+=(--set "ISTRNG=$ISTRNG" --set "FILLRATIO=4.0" --set "SPACE=0.5")
fi
"${RENDER[@]}"

MMPBSA_CMD=(MMPBSA.py -O -i mmpbsa.in -o mmpbsa.dat -eo mmpbsa.csv
            -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop -y "$TRAJ")
note "window: frames $START-$END every $INTERVAL ($NFRAMES frames) of $N"
printf 'MMPBSA command: (cd %q && %q' "$OUTDIR" "${MMPBSA_CMD[0]}" >&2
printf ' %q' "${MMPBSA_CMD[@]:1}" >&2
printf ')\n' >&2

if [ "$CHECK_ONLY" -eq 1 ]; then
  note "check-only: topologies and mmpbsa.in are ready; not running MMPBSA.py"
  exit 0
fi

# ---- run -----------------------------------------------------------------
( cd "$OUTDIR" && "${MMPBSA_CMD[@]}" ) > "$OUTDIR/mmpbsa.log" 2>&1 || {
  tail -50 "$OUTDIR/mmpbsa.log" >&2; die "MMPBSA.py failed"; }

DG=$(awk '/DELTA TOTAL/{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9.]+$/){print $i; exit}}' \
     "$OUTDIR/mmpbsa.dat")
SD=$(awk '/DELTA TOTAL/{c=0; for(i=1;i<=NF;i++) if($i ~ /^-?[0-9.]+$/){c++; if(c==2){print $i; exit}}}' \
     "$OUTDIR/mmpbsa.dat")
[ -n "$DG" ] || { tail -40 "$OUTDIR/mmpbsa.dat" >&2; die "no DELTA TOTAL in mmpbsa.dat"; }
[ -n "$SD" ] || SD=NA

printf 'method\tigb\tsaltcon\tstart\tend\tinterval\tframes\tdg_kcal_mol\tsd_kcal_mol\ttop\ttraj\n' \
  > "$OUTDIR/summary.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$METHOD" "$IGB" "$SALTCON" "$START" "$END" "$INTERVAL" "$NFRAMES" \
  "$DG" "$SD" "$TOP" "$TRAJ" >> "$OUTDIR/summary.tsv"

cat "$OUTDIR/summary.tsv"
note "DELTA TOTAL = $DG kcal/mol (SD $SD); wrote $OUTDIR/summary.tsv"
