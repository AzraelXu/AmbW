#!/bin/bash
# Parameterize a small molecule with antechamber + parmchk2 (GAFF/GAFF2,
# AM1-BCC), or transplant a validated mol2 charge/type set onto a new pose.
#
# The net charge and protonation state are USER INPUTS. This script checks
# them; it does not infer them.
#
# Usage:
#   bash param_ligand.sh --pdb LIG.pdb --charge INT [--resname LIG]
#        [--gaff gaff2|gaff] --outdir DIR [--prefix lig] [--method bcc]
#        [--transplant TEMPLATE.mol2 --onto POSE.pdb]
#        [--allow-attn] [--dry-run] [--force]
set -euo pipefail

KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PDB=""; CHARGE=""; RESNAME="LIG"; GAFF="gaff2"; OUTDIR=""; PREFIX="lig"
METHOD="bcc"; TEMPLATE=""; ONTO=""; ALLOW_ATTN=0; DRY_RUN=0; FORCE=0

die() { echo "error: $*" >&2; exit 2; }
note() { echo "[param_ligand] $*" >&2; }

while [ $# -gt 0 ]; do
  case "$1" in
    --pdb) PDB=$2; shift 2;;
    --charge) CHARGE=$2; shift 2;;
    --resname) RESNAME=$2; shift 2;;
    --gaff) GAFF=$2; shift 2;;
    --outdir) OUTDIR=$2; shift 2;;
    --prefix) PREFIX=$2; shift 2;;
    --method) METHOD=$2; shift 2;;
    --transplant) TEMPLATE=$2; shift 2;;
    --onto) ONTO=$2; shift 2;;
    --allow-attn) ALLOW_ATTN=1; shift;;
    --dry-run) DRY_RUN=1; shift;;
    --force) FORCE=1; shift;;
    -h|--help) sed -n '2,20p' "$0"; exit 0;;
    *) die "unknown argument: $1";;
  esac
done

[ -n "$PDB" ] || die "--pdb is required"
[ -n "$CHARGE" ] || die "--charge is required (user-confirmed net charge)"
[ -n "$OUTDIR" ] || die "--outdir is required"
[[ "$CHARGE" =~ ^-?[0-9]+$ ]] || die "--charge must be an integer, got $CHARGE"
[[ "$RESNAME" =~ ^[A-Za-z][A-Za-z0-9]{0,3}$ ]] || \
  die "--resname must be 1-4 alphanumeric characters starting with a letter"
[[ "$PREFIX" =~ ^[A-Za-z0-9._-]+$ ]] || die "unsafe --prefix: $PREFIX"
[ "$GAFF" = "gaff" ] || [ "$GAFF" = "gaff2" ] || die "--gaff must be gaff or gaff2"
[ "$METHOD" = "bcc" ] || die "only --method bcc is supported"
[ -f "$PDB" ] || die "ligand PDB not found: $PDB"
if [ -n "$TEMPLATE" ] || [ -n "$ONTO" ]; then
  [ -n "$TEMPLATE" ] && [ -n "$ONTO" ] || \
    die "--transplant and --onto must be given together"
  [ -f "$TEMPLATE" ] || die "template mol2 not found: $TEMPLATE"
  [ -f "$ONTO" ] || die "target pose PDB not found: $ONTO"
fi

if ! command -v antechamber >/dev/null 2>&1 || ! command -v parmchk2 >/dev/null 2>&1; then
  if [ -n "${AMBER_MD_AMBER_SH:-}" ] && [ -f "${AMBER_MD_AMBER_SH}" ]; then
    # amber.sh 读取可能未导出的 LD_LIBRARY_PATH/PERL5LIB/PYTHONPATH；
    # 在 set -u 下 source 会直接退出，因此 source 期间临时关闭 nounset。
    # shellcheck disable=SC1090
    set +u
    source "$AMBER_MD_AMBER_SH"
    set -u
  fi
fi
command -v antechamber >/dev/null 2>&1 || \
  die "antechamber not in PATH; source \$AMBERHOME/amber.sh first"
command -v parmchk2 >/dev/null 2>&1 || \
  die "parmchk2 not in PATH; source \$AMBERHOME/amber.sh first"

PY=${AMBER_MD_PYTHON:-python3}
command -v "$PY" >/dev/null 2>&1 || die "python3 not found"

mkdir -p "$OUTDIR"
OUTDIR=$(cd "$OUTDIR" && pwd)
PDB_ABS=$(cd "$(dirname "$PDB")" && pwd)/$(basename "$PDB")

# Managed outputs: refuse to overwrite without --force.
MANAGED=("$OUTDIR/${PREFIX}.mol2" "$OUTDIR/${PREFIX}.frcmod"
         "$OUTDIR/ligand_params.json" "$OUTDIR/input_clean.pdb")
if [ "$DRY_RUN" -eq 0 ] && [ "$FORCE" -eq 0 ]; then
  for f in "${MANAGED[@]}"; do
    if [ -e "$f" ]; then
      die "refusing to overwrite $f (use --force)"
    fi
  done
fi

# First model only, ATOM/HETATM only, no TER/CONECT. This avoids antechamber
# treating a multi-model or multi-chain file as several molecules.
CLEAN="$OUTDIR/input_clean.pdb"
if [ "$DRY_RUN" -eq 0 ]; then
  awk '/^ENDMDL/{exit} /^(ATOM|HETATM)/{print}' "$PDB_ABS" > "$CLEAN"
else
  CLEAN="$PDB_ABS"
fi

N_ATOMS=$(awk '/^(ATOM|HETATM)/{n++} END{print n+0}' "$CLEAN")
N_H=$(awk '/^(ATOM|HETATM)/{e=substr($0,77,2); gsub(/ /,"",e);
  n=substr($0,13,4); gsub(/ /,"",n); nn=n; sub(/^[0-9]+/,"",nn);
  if (e=="H" || e=="D" || (e=="" && (substr(nn,1,1)=="H" || substr(nn,1,1)=="D"))) h++
  } END{print h+0}' "$CLEAN")
N_RES=$(awk '/^(ATOM|HETATM)/{r=substr($0,18,3); gsub(/ /,"",r); print r}' "$CLEAN" | sort -u | wc -l)
PDB_RES=$(awk '/^(ATOM|HETATM)/{r=substr($0,18,3); gsub(/ /,"",r); print r}' "$CLEAN" | sort -u | tr '\n' ',' | sed 's/,$//')

[ "$N_ATOMS" -gt 0 ] || die "no ATOM/HETATM records in $PDB"
[ "$N_RES" -eq 1 ] || die "ligand PDB must contain one residue; found: $PDB_RES"
[ "$N_H" -gt 0 ] || die "ligand PDB has no hydrogens. Protonate it first (e.g. Reduce/OpenBabel) and confirm the net charge; this script does not guess protonation."

note "input: $PDB_RES, $N_ATOMS atoms ($N_H H), target resname $RESNAME, charge $CHARGE, $GAFF"

mol2_stats() {
  awk '/^@<TRIPOS>ATOM/{f=1;next} /^@<TRIPOS>/{f=0}
       f && NF>=9 {n++; q+=$9; if ($6 !~ /^h/) heavy++; rn[$8]++}
       END {printf "%d %.4f %d", n+0, q+0.0, heavy+0}' "$1"
}

# ---- transplant path ------------------------------------------------------
if [ -n "$TEMPLATE" ]; then
  read -r T_ATOMS T_CHARGE T_HEAVY <<< "$(mol2_stats "$TEMPLATE")"
  awk -v q="$T_CHARGE" -v c="$CHARGE" 'BEGIN{d=q-c; if(d<0)d=-d; exit(d>0.01)}' || \
    die "template mol2 charge $T_CHARGE does not match --charge $CHARGE"
  T_RES=$(awk '/^@<TRIPOS>ATOM/{f=1;next} /^@<TRIPOS>/{f=0} f{print $8}' "$TEMPLATE" \
          | sort -u | tr '\n' ',' | sed 's/,$//')
  note "template: $T_RES, $T_ATOMS atoms, charge $T_CHARGE; pose: $ONTO"

  if [ "$DRY_RUN" -eq 1 ]; then
    note "DRY RUN: would transplant $TEMPLATE -> $ONTO -> $OUTDIR/${PREFIX}.mol2"
    note "DRY RUN: would run parmchk2 -s $GAFF"
    exit 0
  fi

  "$PY" "$KIT_DIR/transplant_lig_coords.py" "$TEMPLATE" "$ONTO" \
    "$OUTDIR/${PREFIX}_raw.mol2" > "$OUTDIR/transplant.log" 2>&1 || {
      cat "$OUTDIR/transplant.log" >&2; die "transplant failed"; }

  # Rename the ligand residue/substructure to --resname.
  awk -v rn="$RESNAME" '
    /^@<TRIPOS>ATOM/{f=1; print; next}
    /^@<TRIPOS>SUBSTRUCTURE/{f=2; print; next}
    /^@<TRIPOS>/{f=0; print; next}
    f==1 && NF>=9 {$8=rn; print; next}
    f==2 && NF>=2 {$2=rn; print; next}
    {print}' "$OUTDIR/${PREFIX}_raw.mol2" > "$OUTDIR/${PREFIX}.mol2"
  rm -f "$OUTDIR/${PREFIX}_raw.mol2"
else
  if [ "$DRY_RUN" -eq 1 ]; then
    note "DRY RUN: would run antechamber -c bcc -nc $CHARGE -rn $RESNAME -at $GAFF"
    note "DRY RUN: would run parmchk2 -s $GAFF"
    exit 0
  fi
  ( cd "$OUTDIR" && antechamber -i "$CLEAN" -fi pdb -o "${PREFIX}.mol2" \
      -fo mol2 -c bcc -nc "$CHARGE" -rn "$RESNAME" -at "$GAFF" -pf n ) \
      > "$OUTDIR/antechamber.log" 2>&1 || {
        tail -40 "$OUTDIR/antechamber.log" >&2; die "antechamber failed"; }
fi

[ -s "$OUTDIR/${PREFIX}.mol2" ] || die "antechamber/transplant produced no mol2"
read -r M_ATOMS M_CHARGE M_HEAVY <<< "$(mol2_stats "$OUTDIR/${PREFIX}.mol2")"
[ "$M_ATOMS" -eq "$N_ATOMS" ] || \
  die "atom count changed: input $N_ATOMS, mol2 $M_ATOMS"
awk -v q="$M_CHARGE" -v c="$CHARGE" 'BEGIN{d=q-c; if(d<0)d=-d; exit(d>0.01)}' || \
  die "mol2 charge sum $M_CHARGE != requested $CHARGE"
M_RES=$(awk '/^@<TRIPOS>ATOM/{f=1;next} /^@<TRIPOS>/{f=0} f{print $8}' \
        "$OUTDIR/${PREFIX}.mol2" | sort -u | tr '\n' ',' | sed 's/,$//')
[ "$M_RES" = "$RESNAME" ] || die "mol2 residue name is '$M_RES', expected '$RESNAME'"

ATTN=0
for log in "$OUTDIR/antechamber.log" "$OUTDIR/transplant.log"; do
  if [ -f "$log" ]; then
    ATTN=$((ATTN + $(grep -c 'ATTN' "$log" || true)))
  fi
done
if [ "$ATTN" -gt 0 ] && [ "$ALLOW_ATTN" -eq 0 ]; then
  die "$ATTN ATTN warning(s) in antechamber/transplant log; inspect them or rerun with --allow-attn"
fi

parmchk2 -i "$OUTDIR/${PREFIX}.mol2" -f mol2 -o "$OUTDIR/${PREFIX}.frcmod" \
  -s "$GAFF" > "$OUTDIR/parmchk2.log" 2>&1 || {
    tail -30 "$OUTDIR/parmchk2.log" >&2; die "parmchk2 failed"; }
[ -s "$OUTDIR/${PREFIX}.frcmod" ] || die "parmchk2 produced an empty frcmod"
P_ATTN=$(grep -c 'ATTN' "$OUTDIR/parmchk2.log" || true)
if [ "$P_ATTN" -gt 0 ] && [ "$ALLOW_ATTN" -eq 0 ]; then
  die "$P_ATTN ATTN warning(s) in parmchk2.log; inspect them or rerun with --allow-attn"
fi

# Record everything needed to rebuild or audit the ligand.
TRANSPLANT_JSON="null"
if [ -n "$TEMPLATE" ]; then
  TRANSPLANT_JSON="{\"template\": \"$(cd "$(dirname "$TEMPLATE")" && pwd)/$(basename "$TEMPLATE")\", \"pose\": \"$(cd "$(dirname "$ONTO")" && pwd)/$(basename "$ONTO")\"}"
fi
cat > "$OUTDIR/ligand_params.json" <<EOF
{
  "generated_at": "$(date -Iseconds)",
  "generator": "amber_general_workflow/scripts/param_ligand.sh",
  "input_pdb": "$PDB_ABS",
  "resname": "$RESNAME",
  "charge": $CHARGE,
  "gaff": "$GAFF",
  "charge_method": "$METHOD",
  "atoms": $M_ATOMS,
  "heavy_atoms": $M_HEAVY,
  "hydrogens": $N_H,
  "mol2_charge": $M_CHARGE,
  "mol2": "$OUTDIR/${PREFIX}.mol2",
  "frcmod": "$OUTDIR/${PREFIX}.frcmod",
  "attn": $((ATTN + P_ATTN)),
  "transplant": $TRANSPLANT_JSON
}
EOF

note "wrote $OUTDIR/${PREFIX}.mol2, ${PREFIX}.frcmod, ligand_params.json"
note "mol2: $M_ATOMS atoms ($M_HEAVY heavy), charge $M_CHARGE, resname $RESNAME, gaff=$GAFF"
