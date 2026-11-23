#!/bin/bash
# =============================================================================
# One-replica post-processing for amber-auto:
#   cpptraj autoimage + strip -> analyze_replica -> H-bonds -> MM/GBSA
#
# Runs from the replica directory (PBS_O_WORKDIR). Writes a post_status.json
# with per-step return codes so amber-auto status/report can inspect it.
# =============================================================================
set -u
set -o pipefail

die() { echo "error: $*" >&2; exit 2; }
note() { echo "[post] $*" >&2; }

KIT_DIR=""
TOP=""
PREFIX=""
REP="1"
ANALYSIS_ROOT=""
STRIP_MASK=""
LIGAND_RESNAME=""
PARTNER_RESIDUES=""
PARTNER_KIND="protein"
RECEPTOR_MASK=""
DO_MMPBSA=0
MMPBSA_MASK=""
MMPBSA_METHOD="gb"
MMPBSA_IGB="5"
MMPBSA_SALTCON="0.12"
MMPBSA_INTERVAL="100"
DECOMP_RESIDUES=""
RADII="mbondi2"
DRY_TOP=0
PY_ANALYSIS="/opt/software/amber_general_workflow/envs/analysis/bin/python"
HBOND_STRIDE="2"
NO_HBOND=0
FRAME_PS="10"
AMBER_SH="/home/ylqx/amber22/amber.sh"
FORCE=0

while [ $# -gt 0 ]; do
  case "$1" in
    --kit-dir) KIT_DIR=$2; shift 2;;
    --top) TOP=$2; shift 2;;
    --prefix) PREFIX=$2; shift 2;;
    --rep) REP=$2; shift 2;;
    --analysis-root) ANALYSIS_ROOT=$2; shift 2;;
    --strip-mask) STRIP_MASK=$2; shift 2;;
    --ligand-resname) LIGAND_RESNAME=$2; shift 2;;
    --partner-residues) PARTNER_RESIDUES=$2; shift 2;;
    --partner-kind) PARTNER_KIND=$2; shift 2;;
    --receptor-mask) RECEPTOR_MASK=$2; shift 2;;
    --mmpbsa) DO_MMPBSA=$2; shift 2;;
    --mmpbsa-mask) MMPBSA_MASK=$2; shift 2;;
    --mmpbsa-method) MMPBSA_METHOD=$2; shift 2;;
    --igb) MMPBSA_IGB=$2; shift 2;;
    --saltcon) MMPBSA_SALTCON=$2; shift 2;;
    --interval) MMPBSA_INTERVAL=$2; shift 2;;
    --decomp-residues) DECOMP_RESIDUES=$2; shift 2;;
    --radii) RADII=$2; shift 2;;
    --dry-top) DRY_TOP=$2; shift 2;;
    --py-analysis) PY_ANALYSIS=$2; shift 2;;
    --hbond-stride) HBOND_STRIDE=$2; shift 2;;
    --no-hbond) NO_HBOND=$2; shift 2;;
    --frame-ps) FRAME_PS=$2; shift 2;;
    --amber-sh) AMBER_SH=$2; shift 2;;
    --force) FORCE=1; shift;;
    -h|--help) sed -n '2,12p' "$0"; exit 0;;
    *) die "unknown argument: $1";;
  esac
done

[ -n "$KIT_DIR" ] || die "--kit-dir is required"
[ -n "$TOP" ] || die "--top is required"
[ -n "$PREFIX" ] || die "--prefix is required"
[ -n "$ANALYSIS_ROOT" ] || die "--analysis-root is required"
[ -s "$TOP" ] || die "topology not found: $TOP"
[ -s npt.rst7 ] || echo "warning: npt.rst7 missing; post-processing may be early" >&2

if [ -f "$AMBER_SH" ]; then
  # amber.sh 读取可能未导出的 LD_LIBRARY_PATH/PERL5LIB/PYTHONPATH；
  # 在 set -u 下 source 会直接退出，因此 source 期间临时关闭 nounset。
  # shellcheck disable=SC1090
  set +u
  source "$AMBER_SH" || die "cannot source $AMBER_SH"
  set -u
fi

CPPDIR="$ANALYSIS_ROOT/cpptraj"
PYDIR="$ANALYSIS_ROOT/python"
MMDIR="$ANALYSIS_ROOT/mmpbsa"
TABDIR="$ANALYSIS_ROOT/tables"
mkdir -p "$CPPDIR" "$PYDIR" "$MMDIR" "$TABDIR" || die "cannot create analysis dirs"

if [ "$DRY_TOP" -eq 1 ]; then
  DRY_PRM="$CPPDIR/${PREFIX}_dry.prmtop"
  DRY_NC="$CPPDIR/${PREFIX}_dry.nc"
else
  DRY_PRM="$CPPDIR/${PREFIX}_dry.prmtop"
  DRY_NC="$CPPDIR/${PREFIX}_dry.nc"
fi

status_json="$ANALYSIS_ROOT/post_status_${PREFIX}.json"
step_cpptraj=1
step_analyze=1
step_hbond=2   # 2 = skipped, 0 = success, 1 = failure
step_mmpbsa=2  # 2 = skipped, 0 = success, 1 = failure
start_ts=$(date -Is)

# ---- 1) cpptraj: autoimage, strip, dry NetCDF ------------------------------
note "cpptraj prep ($PREFIX)"
prep_args=(--top "$TOP" --traj "$PWD/05_prod_c*.nc"
           --outdir "$CPPDIR" --prefix "$PREFIX" --run)
if [ "$DRY_TOP" -eq 0 ]; then
  prep_args+=(--strip-mask "$STRIP_MASK")
fi
[ "$FORCE" -eq 1 ] && prep_args+=(--force)
if "$PY_ANALYSIS" -I "$KIT_DIR/scripts/make_cpptraj_inputs.py" "${prep_args[@]}"; then
  step_cpptraj=0
else
  note "cpptraj prep failed; stopping"
fi

# ---- 2) structural analysis ----------------------------------------------
if [ "$step_cpptraj" -eq 0 ]; then
  note "analyze_replica ($PREFIX)"
  ana_args=(--prmtop "$DRY_PRM" --traj "$DRY_NC"
            --label "$PREFIX" --rep "$REP"
            --out "$PYDIR/${PREFIX}.npz" --summary "$PYDIR/${PREFIX}.json"
            --frame-ps "$FRAME_PS")
  if [ -n "$PARTNER_RESIDUES" ]; then
    ana_args+=(--partner-residues "$PARTNER_RESIDUES" --partner-kind "$PARTNER_KIND")
    hbond_stride_use=10
  elif [ -n "$LIGAND_RESNAME" ]; then
    ana_args+=(--ligand-resname "$LIGAND_RESNAME")
    hbond_stride_use=${HBOND_STRIDE:-2}
  else
    ana_args+=(--apo)
    hbond_stride_use=0
  fi
  if [ -n "$RECEPTOR_MASK" ]; then
    ana_args+=(--sasa --receptor-mask "$RECEPTOR_MASK")
  fi
  [ "$FORCE" -eq 1 ] && ana_args+=(--force)
  if "$PY_ANALYSIS" -I "$KIT_DIR/scripts/analyze_replica.py" "${ana_args[@]}"; then
    step_analyze=0
  else
    note "analyze_replica failed"
  fi

  # ---- 3) H-bonds ---------------------------------------------------------
  if [ "$NO_HBOND" -eq 0 ] && [ "$hbond_stride_use" -gt 0 ]; then
    note "H-bond analysis ($PREFIX, stride $hbond_stride_use)"
    hb_args=(--prmtop "$DRY_PRM" --traj "$DRY_NC"
             --out "$PYDIR/${PREFIX}_hb.npz" --summary "$PYDIR/${PREFIX}_hb.json"
             --stride "$hbond_stride_use")
    if [ -n "$PARTNER_RESIDUES" ]; then
      hb_args+=(--partner-residues "$PARTNER_RESIDUES" --partner-kind "$PARTNER_KIND")
    else
      hb_args+=(--ligand-resname "$LIGAND_RESNAME")
    fi
    [ "$FORCE" -eq 1 ] && hb_args+=(--force)
    if "$PY_ANALYSIS" -I "$KIT_DIR/scripts/analyze_ligand_hbonds.py" "${hb_args[@]}"; then
      step_hbond=0
    else
      step_hbond=1
      note "H-bond analysis failed (non-fatal)"
    fi
  fi
fi

# ---- 4) MM/GBSA -----------------------------------------------------------
if [ "$step_cpptraj" -eq 0 ] && [ "$DO_MMPBSA" -eq 1 ]; then
  note "MM/GBSA ($PREFIX)"
  mmp_args=(--top "$TOP" --traj "$DRY_NC" --ligand-mask "$MMPBSA_MASK"
            --receptor-mask "$RECEPTOR_MASK"
            --strip-mask "$STRIP_MASK" --radii "$RADII"
            --method "$MMPBSA_METHOD" --igb "$MMPBSA_IGB"
            --saltcon "$MMPBSA_SALTCON" --interval "$MMPBSA_INTERVAL"
            --outdir "$MMDIR/$PREFIX")
  [ -n "$DECOMP_RESIDUES" ] && mmp_args+=(--decomp-residues "$DECOMP_RESIDUES")
  [ "$DRY_TOP" -eq 1 ] && mmp_args+=(--dry-top)
  [ "$FORCE" -eq 1 ] && mmp_args+=(--force)
  if bash "$KIT_DIR/scripts/mmpbsa_one.sh" "${mmp_args[@]}"; then
    step_mmpbsa=0
  else
    step_mmpbsa=1
    note "MM/GBSA failed"
  fi
fi

end_ts=$(date -Is)
cat > "$status_json" <<EOF
{
  "prefix": "$PREFIX",
  "rep": "$REP",
  "top": "$(cd "$(dirname "$TOP")" && pwd)/$(basename "$TOP")",
  "dry_prmtop": "$DRY_PRM",
  "dry_traj": "$DRY_NC",
  "started": "$start_ts",
  "finished": "$end_ts",
  "steps": {
    "cpptraj": $step_cpptraj,
    "analyze": $step_analyze,
    "hbond": $step_hbond,
    "mmpbsa": $step_mmpbsa
  }
}
EOF

if [ "$step_cpptraj" -ne 0 ] || [ "$step_analyze" -ne 0 ]; then
  exit 1
fi
if [ "$DO_MMPBSA" -eq 1 ] && [ "$step_mmpbsa" -ne 0 ]; then
  exit 1
fi
exit 0
