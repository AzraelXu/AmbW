#!/bin/bash
# Build an Amber system with pdb4amber + tleap and QC it.
#
# Covers RNA-apo, protein-apo, RNA-ligand, protein-ligand; explicit solvent or
# implicit (dry) systems. It does not guess ligand protonation or net charge:
# run param_ligand.sh first and pass the resulting mol2/frcmod.
#
# Usage:
#   bash build_system.sh --type {rna-apo,protein-apo,rna-ligand,protein-ligand} \
#     --receptor RECEPTOR.pdb [--ligand lig.mol2 --frcmod lig.frcmod] \
#     [--lig-resname LIG] [--ligand-pdb POSE.pdb] [--gaff gaff2|gaff] \
#     --ff {RNA.YIL,RNA.OL3,RNA.ROC,RNA.LJbb,ff14SB,ff19SB} \
#     --solvent {explicit,implicit} [--water {tip3p,opc,spce,tip4pd,fb3}] \
#     [--box {oct,rect}] [--pad 10.0] [--cation K+] [--anion Cl-] \
#     [--counterion {cation,anion,none}] [--salt-pairs 0] \
#     [--radii {mbondi,mbondi2,mbondi3}] [--igb 8] \
#     --outdir DIR [--top-name NAME] [--extra-leap FILE] \
#     [--restraint-mask MASK] [--no-pdb4amber] [--allow-hetatm] \
#     [--allow-attn] [--no-strict-qc] [--render-only] [--force]
set -euo pipefail

KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TYPE=""; RECEPTOR=""; LIGAND=""; FRCMOD=""; LIG_RESNAME=""; LIG_PDB=""
GAFF="gaff2"; FF=""; SOLVENT="explicit"; WATER=""; BOX="oct"; PAD="10.0"
CATION=""; ANION="Cl-"; COUNTERION=""; SALT_PAIRS=0; RADII=""; IGB="8"
OUTDIR=""; TOP_NAME=""; EXTRA_LEAP=""; RESTRAINT_MASK=""
NO_PDB4AMBER=0; ALLOW_HETATM=0; ALLOW_ATTN=0; STRICT_QC=1
RENDER_ONLY=0; FORCE=0
ORIG_ARGS=("$@")
PROVIDED=$(printf '%s\n' "${ORIG_ARGS[@]}" | sed -n 's/^--\([^=]*\).*/\1/p' || true)
has_provided() { printf '%s\n' "$PROVIDED" | grep -qx "$1"; }

die() { echo "error: $*" >&2; exit 2; }
note() { echo "[build_system] $*" >&2; }
warn() { echo "[build_system] WARNING: $*" >&2; }

while [ $# -gt 0 ]; do
  case "$1" in
    --type) TYPE=$2; shift 2;;
    --receptor) RECEPTOR=$2; shift 2;;
    --ligand) LIGAND=$2; shift 2;;
    --frcmod) FRCMOD=$2; shift 2;;
    --lig-resname) LIG_RESNAME=$2; shift 2;;
    --ligand-pdb) LIG_PDB=$2; shift 2;;
    --gaff) GAFF=$2; shift 2;;
    --ff) FF=$2; shift 2;;
    --solvent) SOLVENT=$2; shift 2;;
    --water) WATER=$2; shift 2;;
    --box) BOX=$2; shift 2;;
    --pad) PAD=$2; shift 2;;
    --cation) CATION=$2; shift 2;;
    --anion) ANION=$2; shift 2;;
    --counterion) COUNTERION=$2; shift 2;;
    --salt-pairs) SALT_PAIRS=$2; shift 2;;
    --radii) RADII=$2; shift 2;;
    --igb) IGB=$2; shift 2;;
    --outdir) OUTDIR=$2; shift 2;;
    --top-name) TOP_NAME=$2; shift 2;;
    --extra-leap) EXTRA_LEAP=$2; shift 2;;
    --restraint-mask) RESTRAINT_MASK=$2; shift 2;;
    --no-pdb4amber) NO_PDB4AMBER=1; shift;;
    --allow-hetatm) ALLOW_HETATM=1; shift;;
    --allow-attn) ALLOW_ATTN=1; shift;;
    --no-strict-qc) STRICT_QC=0; shift;;
    --render-only) RENDER_ONLY=1; shift;;
    --force) FORCE=1; shift;;
    -h|--help) sed -n '2,22p' "$0"; exit 0;;
    *) die "unknown argument: $1";;
  esac
done

case "$TYPE" in
  rna-apo|protein-apo) HOLO=0;;
  rna-ligand|protein-ligand) HOLO=1;;
  *) die "--type must be rna-apo, protein-apo, rna-ligand, or protein-ligand";;
esac
[ -n "$RECEPTOR" ] || die "--receptor is required"
[ -n "$OUTDIR" ] || die "--outdir is required"
[ -f "$RECEPTOR" ] || die "receptor PDB not found: $RECEPTOR"
[[ "$SALT_PAIRS" =~ ^[0-9]+$ ]] || die "--salt-pairs must be a non-negative integer"
[[ "$IGB" =~ ^(1|2|5|7|8)$ ]] || die "--igb must be 1, 2, 5, 7, or 8"
[ "$GAFF" = gaff ] || [ "$GAFF" = gaff2 ] || die "--gaff must be gaff or gaff2"
[ "$SOLVENT" = explicit ] || [ "$SOLVENT" = implicit ] || \
  die "--solvent must be explicit or implicit"
[ "$BOX" = oct ] || [ "$BOX" = rect ] || die "--box must be oct or rect"
awk -v p="$PAD" 'BEGIN{exit(p<=0)}' || die "--pad must be > 0"

# ---- defaults by system type ---------------------------------------------
case "$TYPE" in
  rna-*)
    [ -n "$FF" ] || FF=RNA.YIL
    [ -n "$CATION" ] || CATION=K+
    ;;
  protein-*)
    [ -n "$FF" ] || FF=ff19SB
    [ -n "$CATION" ] || CATION=Na+
    ;;
esac
case "$FF" in
  RNA.YIL|RNA.OL3|RNA.ROC|RNA.LJbb)
    [ "${TYPE%%-*}" = rna ] || die "--ff $FF is an RNA force field but --type is $TYPE"
    LEAPRC_RECEPTOR="leaprc.$FF"
    ;;
  ff14SB|ff19SB)
    [ "${TYPE%%-*}" = protein ] || die "--ff $FF is a protein force field but --type is $TYPE"
    LEAPRC_RECEPTOR="leaprc.protein.$FF"
    ;;
  *) die "--ff must be RNA.YIL, RNA.OL3, RNA.ROC, RNA.LJbb, ff14SB, or ff19SB";;
esac
if ! has_provided water; then
  case "$FF" in
    ff19SB) WATER=opc;;
    *) WATER=tip3p;;
  esac
fi

case "$FF:$WATER" in
  ff19SB:tip3p) warn "ff19SB was parameterized with OPC; TIP3P is not recommended";;
  ff14SB:opc) warn "ff14SB is usually paired with TIP3P; confirm OPC is intended";;
  RNA.*:opc|RNA.*:tip4pd|RNA.*:fb3|RNA.*:spce)
    warn "RNA force fields are usually paired with TIP3P; $WATER is not standard";;
esac

if [ "$SOLVENT" = implicit ]; then
  [ "$SALT_PAIRS" -eq 0 ] || die "--salt-pairs requires explicit solvent"
  has_provided water && warn "--water is ignored for implicit solvent"
  has_provided box && warn "--box is ignored for implicit solvent"
  has_provided pad && warn "--pad is ignored for implicit solvent"
  if [ -n "$COUNTERION" ] && [ "$COUNTERION" != none ]; then
    warn "adding counterions to an implicit/dry system; usually --counterion none is safer"
  else
    [ -n "$COUNTERION" ] || COUNTERION=none
  fi
  [ "$COUNTERION" = none ] || [ "$COUNTERION" = cation ] || [ "$COUNTERION" = anion ] || \
    die "--counterion must be cation, anion, or none"
  LEAPRC_WATER=""; LEAPRC_WATER_BLOCK=""; WATER_BOX=""
  SOLVATE_BLOCK=""
  [ -n "$RADII" ] || { if [ "$IGB" = 8 ]; then RADII=mbondi3; else RADII=mbondi2; fi; }
else
  [ -n "$COUNTERION" ] || COUNTERION=cation
  [ "$COUNTERION" = cation ] || [ "$COUNTERION" = anion ] || [ "$COUNTERION" = none ] || \
    die "--counterion must be cation, anion, or none"
  [ -n "$RADII" ] || RADII=mbondi2
  case "$WATER" in
    tip3p) LEAPRC_WATER="leaprc.water.tip3p"; WATER_BOX="TIP3PBOX";;
    opc) LEAPRC_WATER="leaprc.water.opc"; WATER_BOX="OPCBOX";;
    spce) LEAPRC_WATER="leaprc.water.spce"; WATER_BOX="SPCBOX";;
    tip4pd) LEAPRC_WATER="leaprc.water.tip4pd"; WATER_BOX="TIP4PDBOX";;
    fb3) LEAPRC_WATER="leaprc.water.fb3"; WATER_BOX="FB3BOX";;
    *) die "--water must be tip3p, opc, spce, tip4pd, or fb3";;
  esac
  LEAPRC_WATER_BLOCK="source $LEAPRC_WATER"
fi
case "$RADII" in
  mbondi|mbondi2|mbondi3) ;;
  *) die "--radii must be mbondi, mbondi2, or mbondi3";;
esac
[ -n "$TOP_NAME" ] || { if [ "$HOLO" -eq 1 ]; then TOP_NAME=complex; else TOP_NAME=apo; fi; }
[[ "$TOP_NAME" =~ ^[A-Za-z0-9._-]+$ ]] || die "unsafe --top-name: $TOP_NAME"

has_provided ff || note "default force field: $FF"
has_provided solvent || note "default solvent: $SOLVENT"
if [ "$SOLVENT" = explicit ] && ! has_provided water; then
  note "default water: $WATER"
fi
has_provided cation || note "default cation: $CATION"
has_provided counterion || note "default counterion: $COUNTERION"
has_provided radii || note "default radii: $RADII"

if [ "$HOLO" -eq 1 ]; then
  [ -n "$LIGAND" ] || die "--ligand is required for $TYPE"
  [ -n "$FRCMOD" ] || die "--frcmod is required for $TYPE"
  [ -f "$LIGAND" ] || die "ligand mol2 not found: $LIGAND"
  [ -s "$FRCMOD" ] || die "ligand frcmod missing or empty: $FRCMOD"
  if [ -n "$LIG_PDB" ] && [ ! -f "$LIG_PDB" ]; then die "--ligand-pdb not found: $LIG_PDB"; fi
fi
if [ -n "$EXTRA_LEAP" ] && [ ! -f "$EXTRA_LEAP" ]; then die "--extra-leap not found: $EXTRA_LEAP"; fi
if [ -n "$RESTRAINT_MASK" ] && [ "$RESTRAINT_MASK" = none ]; then
  die "--restraint-mask none is not a build option; omit it or give a real mask"
fi

# ---- Amber environment ---------------------------------------------------
if ! command -v tleap >/dev/null 2>&1; then
  if [ -n "${AMBER_MD_AMBER_SH:-}" ] && [ -f "$AMBER_MD_AMBER_SH" ]; then
    # amber.sh 读取可能未导出的 LD_LIBRARY_PATH/PERL5LIB/PYTHONPATH；
    # 在 set -u 下 source 会直接退出，因此 source 期间临时关闭 nounset。
    # shellcheck disable=SC1090
    set +u
    source "$AMBER_MD_AMBER_SH"
    set -u
  fi
fi
command -v tleap >/dev/null 2>&1 || die "tleap not in PATH; source \$AMBERHOME/amber.sh"
if [ "$NO_PDB4AMBER" -eq 0 ]; then
  command -v pdb4amber >/dev/null 2>&1 || die "pdb4amber not in PATH"
fi
if [ -n "${AMBERHOME:-}" ]; then
  [ -f "$AMBERHOME/dat/leap/cmd/$LEAPRC_RECEPTOR" ] || \
    die "leaprc not found: $AMBERHOME/dat/leap/cmd/$LEAPRC_RECEPTOR"
  if [ -n "$LEAPRC_WATER" ]; then
    [ -f "$AMBERHOME/dat/leap/cmd/$LEAPRC_WATER" ] || \
      die "leaprc not found: $AMBERHOME/dat/leap/cmd/$LEAPRC_WATER"
  fi
fi

mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)
RECEPTOR_ABS=$(cd "$(dirname "$RECEPTOR")" && pwd)/$(basename "$RECEPTOR")

# Managed outputs
MANAGED=("$OUTDIR/leap.in" "$OUTDIR/tleap.out" "$OUTDIR/build_manifest.json"
         "$OUTDIR/build_qc.json" "$OUTDIR/build_qc.tsv" "$OUTDIR/receptor_clean.pdb"
         "$OUTDIR/${TOP_NAME}.prmtop" "$OUTDIR/${TOP_NAME}.inpcrd")
if [ "$FORCE" -eq 0 ]; then
  for f in "${MANAGED[@]}"; do
    if [ -e "$f" ]; then die "refusing to overwrite $f (use --force)"; fi
  done
else
  for f in "${MANAGED[@]}"; do rm -f "$f"; done
fi

# ---- receptor PDB checks -------------------------------------------------
ALL_RES=$(awk '/^(ATOM|HETATM)/{r=substr($0,18,3); gsub(/ /,"",r); print r}' "$RECEPTOR_ABS" | sort -u)
HET_RES=$(awk '/^HETATM/{r=substr($0,18,3); gsub(/ /,"",r); print r}' "$RECEPTOR_ABS" | sort -u)
if [ "$HOLO" -eq 1 ] && [ -n "$LIG_RESNAME" ] && \
   printf '%s\n' "$ALL_RES" | grep -qx "$LIG_RESNAME"; then
  die "receptor PDB already contains residue $LIG_RESNAME; pass a receptor-only PDB"
fi
ALLOWED_HET="WAT HOH H2O SOL NA NA+ K K+ CL CL- MG MG2 CA CA2 ZN ZN2 LI LI+ RB RB+ CS CS+ F F- BR BR- I I-"
UNKNOWN_HET=""
while IFS= read -r r; do
  [ -n "$r" ] || continue
  if printf '%s\n' "$ALLOWED_HET" | tr ' ' '\n' | grep -qx "$r"; then
    continue
  fi
  UNKNOWN_HET="$UNKNOWN_HET $r"
done <<< "$HET_RES"
if [ -n "$UNKNOWN_HET" ]; then
  if [ "$ALLOW_HETATM" -eq 1 ]; then
    warn "keeping non-standard HETATM residues:$UNKNOWN_HET"
  else
    die "non-standard HETATM residues in receptor:$UNKNOWN_HET. Remove them, or rerun with --allow-hetatm after confirming they belong in the system."
  fi
fi
N_MODELS=$(grep -c '^MODEL' "$RECEPTOR_ABS" || true)
if [ "$N_MODELS" -gt 1 ]; then
  warn "receptor PDB has $N_MODELS models; only the first will be used"
fi

# ---- ligand mol2 checks --------------------------------------------------
LIG_ATOMS=0; LIG_CHARGE=0; LIG_HEAVY=0
if [ "$HOLO" -eq 1 ]; then
  read -r LIG_ATOMS LIG_CHARGE LIG_HEAVY LIG_RESNAMES <<< "$(
    awk '/^@<TRIPOS>ATOM/{f=1;next} /^@<TRIPOS>/{f=0}
         f && NF>=9 {n++; q+=$9; if ($6 !~ /^h/) h++; rn[$8]=1}
         END {s=""; for (r in rn) s=s r ","; sub(/,$/,"",s);
              printf "%d %.4f %d %s", n+0, q+0.0, h+0, s}' "$LIGAND")"
  [ "$LIG_ATOMS" -gt 0 ] || die "no atoms found in ligand mol2: $LIGAND"
  if [ -z "$LIG_RESNAME" ]; then
    LIG_RESNAME=$LIG_RESNAMES
    note "ligand residue name from mol2: $LIG_RESNAME"
  fi
  [ "$LIG_RESNAMES" = "$LIG_RESNAME" ] || \
    die "ligand mol2 residue name is '$LIG_RESNAMES', not '$LIG_RESNAME'"
  if printf '%s\n' "$ALL_RES" | grep -qx "$LIG_RESNAME"; then
    die "receptor PDB already contains residue $LIG_RESNAME; pass a receptor-only PDB"
  fi
  note "ligand: $LIG_RESNAME, $LIG_ATOMS atoms ($LIG_HEAVY heavy), charge $LIG_CHARGE, gaff=$GAFF"
fi

# ---- clean receptor ------------------------------------------------------
if [ "$NO_PDB4AMBER" -eq 1 ]; then
  if printf '%s\n' "$ALL_RES" | grep -qxE 'WAT|HOH|H2O|SOL'; then
    die "receptor PDB contains water but --no-pdb4amber was given; remove water first"
  fi
  RECEPTOR_USED=$RECEPTOR_ABS
  note "using receptor PDB without pdb4amber"
else
  ( cd "$OUTDIR" && pdb4amber -i "$RECEPTOR_ABS" -o receptor_clean.pdb \
      --dry --most-populous --nohyd -l pdb4amber.log ) \
      > "$OUTDIR/pdb4amber.stdout" 2>&1 || {
        tail -30 "$OUTDIR/pdb4amber.stdout" >&2; die "pdb4amber failed"; }
  [ -s "$OUTDIR/receptor_clean.pdb" ] || die "pdb4amber produced no output PDB"
  RECEPTOR_USED="$OUTDIR/receptor_clean.pdb"
fi

# ---- build leap.in blocks ------------------------------------------------
RECEPTOR_VAR=mol
if [ "$HOLO" -eq 1 ]; then UNIT=com; else UNIT=mol; fi

if [ "$SOLVENT" = explicit ]; then
  if [ "$BOX" = oct ]; then
    SOLVATE_BLOCK="solvateOct $UNIT $WATER_BOX $PAD"
  else
    SOLVATE_BLOCK="solvateBox $UNIT $WATER_BOX $PAD"
  fi
else
  SOLVATE_BLOCK="# implicit solvent: no solvate command"
fi

NEUTRALIZE_BLOCK="# no counterions requested"
if [ "$COUNTERION" = cation ]; then
  NEUTRALIZE_BLOCK="addIons $UNIT $CATION 0"
elif [ "$COUNTERION" = anion ]; then
  NEUTRALIZE_BLOCK="addIons $UNIT $ANION 0"
fi

SALT_BLOCK="# no added salt"
if [ "$SALT_PAIRS" -gt 0 ]; then
  SALT_BLOCK="addIonsRand $UNIT $CATION $SALT_PAIRS $ANION $SALT_PAIRS"
  if [ "$COUNTERION" = none ]; then
    warn "salt pairs are added without neutralization; confirm the net charge"
  fi
fi

EXTRA_BLOCK="# (no extra leap commands)"
if [ -n "$EXTRA_LEAP" ]; then
  EXTRA_BLOCK=$(cat "$EXTRA_LEAP")
fi

if [ "$HOLO" -eq 1 ]; then
  TEMPLATE="$KIT_DIR/../templates/build/build_${TYPE%%-*}_ligand.in"
else
  TEMPLATE="$KIT_DIR/../templates/build/build_apo.in"
fi
[ -f "$TEMPLATE" ] || die "build template not found: $TEMPLATE"

RENDER=(python3 "$KIT_DIR/render_template.py" "$TEMPLATE" --out "$OUTDIR/leap.in" --force
  --set "LEAPRC_RECEPTOR=$LEAPRC_RECEPTOR"
  --set "LEAPRC_WATER_BLOCK=$LEAPRC_WATER_BLOCK"
  --set "RECEPTOR_VAR=$RECEPTOR_VAR"
  --set "RECEPTOR_PDB=$RECEPTOR_USED"
  --set "SOLVATE_BLOCK=$SOLVATE_BLOCK"
  --set "NEUTRALIZE_BLOCK=$NEUTRALIZE_BLOCK"
  --set "SALT_BLOCK=$SALT_BLOCK"
  --set "RADII=$RADII"
  --set "EXTRA_LEAP=$EXTRA_BLOCK"
  --set "TOP_NAME=$TOP_NAME")
if [ "$HOLO" -eq 1 ]; then
  RENDER+=(--set "LEAPRC_GAFF=leaprc.$GAFF"
           --set "LIG_FRCMOD=$(cd "$(dirname "$FRCMOD")" && pwd)/$(basename "$FRCMOD")"
           --set "LIG_VAR=LIG"
           --set "LIG_MOL2=$(cd "$(dirname "$LIGAND")" && pwd)/$(basename "$LIGAND")")
fi
"${RENDER[@]}"

note "rendered $OUTDIR/leap.in"
if [ "$RENDER_ONLY" -eq 1 ]; then
  note "render-only: stopping before tleap"
  exit 0
fi

# ---- run tleap -----------------------------------------------------------
( cd "$OUTDIR" && tleap -f leap.in ) > "$OUTDIR/tleap.out" 2>&1 || true
[ -s "$OUTDIR/tleap.out" ] || die "tleap produced no output"
TLEAP_LINE=$(grep 'Exiting LEaP' "$OUTDIR/tleap.out" | tail -1 || true)
[ -n "$TLEAP_LINE" ] || { tail -40 "$OUTDIR/tleap.out" >&2; die "tleap did not finish"; }
TLEAP_ERRORS=$(sed -n 's/.*Errors = \([0-9]*\).*/\1/p' <<< "$TLEAP_LINE")
TLEAP_WARN=$(sed -n 's/.*Warnings = \([0-9]*\).*/\1/p' <<< "$TLEAP_LINE")
TLEAP_ATTN=$(grep -c 'ATTN' "$OUTDIR/tleap.out" || true)
if grep -qE 'Could not find (bond|angle|torsion|improper) parameter|FATAL|Error:' "$OUTDIR/tleap.out"; then
  grep -nE 'Could not find|FATAL|Error:' "$OUTDIR/tleap.out" | head -20 >&2
  die "tleap reported a parameter/unit error; see $OUTDIR/tleap.out"
fi
[ "${TLEAP_ERRORS:-1}" -eq 0 ] || die "tleap finished with Errors = ${TLEAP_ERRORS:-?}; see $OUTDIR/tleap.out"
if [ "$TLEAP_ATTN" -gt 0 ] && [ "$ALLOW_ATTN" -eq 0 ]; then
  grep -n 'ATTN' "$OUTDIR/tleap.out" | head -10 >&2
  die "tleap emitted $TLEAP_ATTN ATTN warning(s); inspect them or rerun with --allow-attn"
fi

PRMTOP="$OUTDIR/$TOP_NAME.prmtop"
INPCRD="$OUTDIR/$TOP_NAME.inpcrd"
[ -s "$PRMTOP" ] || die "tleap did not write $PRMTOP"
[ -s "$INPCRD" ] || die "tleap did not write $INPCRD"

RADIUS_LINE=$(grep -A2 '%FLAG RADIUS_SET' "$PRMTOP" | tail -1 || true)
printf '%s' "$RADIUS_LINE" | grep -qF "($RADII)" || \
  die "RADIUS_SET in topology is '$RADIUS_LINE', expected $RADII"

TLEAP_NETQ=$(grep 'Total unperturbed charge:' "$OUTDIR/tleap.out" | tail -1 | sed 's/.*: *//' || true)
[ -n "$TLEAP_NETQ" ] || TLEAP_NETQ=0
if [ "$COUNTERION" != none ]; then
  awk -v q="$TLEAP_NETQ" 'BEGIN{d=q; if(d<0)d=-d; exit(d>0.01)}' || \
    die "net charge after neutralization is $TLEAP_NETQ; the counterion sign may be wrong (try --counterion anion/cation)"
else
  note "net charge (no counterions): $TLEAP_NETQ"
fi

# ---- QC ------------------------------------------------------------------
# shellcheck source=pyenv.sh
source "$KIT_DIR/pyenv.sh"
if ! amber_md_pick_python; then
  warn "no Python with parmed+numpy found; QC mask validation will be skipped"
fi

QC_NETQ=any
if [ "$COUNTERION" != none ]; then QC_NETQ=0.0; fi
QC_WATER=yes
if [ "$SOLVENT" = implicit ]; then QC_WATER=no; fi
QC_ARGS=(--prmtop "$PRMTOP" --inpcrd "$INPCRD" --label "$(basename "$OUTDIR")"
         --expect-net-charge "$QC_NETQ" --expect-water "$QC_WATER"
         --expect-radii "$RADII" --json "$OUTDIR/build_qc.json" --tsv)
if [ "$HOLO" -eq 1 ]; then
  QC_ARGS+=(--ligand-resname "$LIG_RESNAME" --expect-lig-charge "$LIG_CHARGE")
  if [ -n "$LIG_PDB" ]; then
    QC_ARGS+=(--ligand-pdb "$(cd "$(dirname "$LIG_PDB")" && pwd)/$(basename "$LIG_PDB")")
  fi
fi
if [ -n "$RESTRAINT_MASK" ]; then
  QC_ARGS+=(--restraint-mask "$RESTRAINT_MASK")
fi

QC_RC=0
"$AMBER_MD_PY" $AMBER_MD_PY_FLAGS "$KIT_DIR/qc_system.py" "${QC_ARGS[@]}" \
  > "$OUTDIR/build_qc.tsv" 2> "$OUTDIR/build_qc.err" || QC_RC=$?

# ---- manifest ------------------------------------------------------------
export BUILD_TYPE="$TYPE" BUILD_FF="$FF" BUILD_LEAPRC="$LEAPRC_RECEPTOR"
export BUILD_WATER="$WATER" BUILD_SOLVENT="$SOLVENT" BUILD_BOX="$BOX" BUILD_PAD="$PAD"
export BUILD_CATION="$CATION" BUILD_ANION="$ANION" BUILD_COUNTERION="$COUNTERION"
export BUILD_SALT_PAIRS="$SALT_PAIRS" BUILD_RADII="$RADII" BUILD_IGB="$IGB"
export BUILD_TOP_NAME="$TOP_NAME" BUILD_PRMTOP="$PRMTOP" BUILD_INPCRD="$INPCRD"
export BUILD_TLEAP_ERRORS="$TLEAP_ERRORS" BUILD_TLEAP_WARN="$TLEAP_WARN"
export BUILD_TLEAP_ATTN="$TLEAP_ATTN" BUILD_TLEAP_NETQ="$TLEAP_NETQ"
export BUILD_LIG_RESNAME="$LIG_RESNAME" BUILD_LIG_ATOMS="$LIG_ATOMS"
export BUILD_LIG_CHARGE="$LIG_CHARGE" BUILD_LIG_HEAVY="$LIG_HEAVY"
export BUILD_STAGE="qc_failed"
if [ "$QC_RC" -eq 0 ]; then BUILD_STAGE="ok"; fi

"$AMBER_MD_PY" $AMBER_MD_PY_FLAGS - "$OUTDIR/build_manifest.json" "$OUTDIR/build_qc.json" <<'PY'
import json, os, sys
out, qc_path = sys.argv[1], sys.argv[2]
env = os.environ
def num(x, kind=float):
    try:
        return kind(x)
    except (TypeError, ValueError):
        return x
manifest = {
    "generated_at": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
    "generator": "amber_general_workflow/scripts/build_system.sh",
    "type": env["BUILD_TYPE"],
    "force_field": env["BUILD_FF"],
    "leaprc_receptor": env["BUILD_LEAPRC"],
    "solvent": env["BUILD_SOLVENT"],
    "water": env["BUILD_WATER"] or None,
    "box": env["BUILD_BOX"],
    "padding_A": num(env["BUILD_PAD"]),
    "cation": env["BUILD_CATION"],
    "anion": env["BUILD_ANION"],
    "counterion": env["BUILD_COUNTERION"],
    "salt_pairs": num(env["BUILD_SALT_PAIRS"], int),
    "radii": env["BUILD_RADII"],
    "igb": num(env["BUILD_IGB"], int),
    "top_name": env["BUILD_TOP_NAME"],
    "prmtop": env["BUILD_PRMTOP"],
    "inpcrd": env["BUILD_INPCRD"],
    "ligand": {
        "resname": env["BUILD_LIG_RESNAME"] or None,
        "atoms": num(env["BUILD_LIG_ATOMS"], int),
        "heavy_atoms": num(env["BUILD_LIG_HEAVY"], int),
        "charge": num(env["BUILD_LIG_CHARGE"]),
    } if env["BUILD_LIG_RESNAME"] else None,
    "tleap": {
        "errors": num(env["BUILD_TLEAP_ERRORS"], int),
        "warnings": num(env["BUILD_TLEAP_WARN"], int),
        "attn": num(env["BUILD_TLEAP_ATTN"], int),
        "net_charge": num(env["BUILD_TLEAP_NETQ"]),
    },
    "status": env["BUILD_STAGE"],
}
if os.path.exists(qc_path):
    manifest["qc"] = json.load(open(qc_path))
with open(out, "w") as fh:
    json.dump(manifest, fh, indent=2)
    fh.write("\n")
PY

cat "$OUTDIR/build_qc.tsv" 2>/dev/null || true
note "topology: $PRMTOP"
note "tleap: Errors=$TLEAP_ERRORS Warnings=$TLEAP_WARN ATTN=$TLEAP_ATTN net_charge=$TLEAP_NETQ"
note "manifest: $OUTDIR/build_manifest.json"

if [ "$QC_RC" -ne 0 ]; then
  cat "$OUTDIR/build_qc.err" >&2 || true
  if [ "$STRICT_QC" -eq 1 ]; then
    die "QC failed; inspect $OUTDIR/build_qc.tsv and $OUTDIR/build_qc.err"
  else
    warn "QC failed but --no-strict-qc was given; build artifacts are kept"
  fi
fi
note "done"
