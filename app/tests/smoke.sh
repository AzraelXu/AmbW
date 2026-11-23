#!/bin/bash
# Smoke test for amber_general_workflow.
# - syntax checks
# - Python compile
# - tleap render-only (protein/RNA/holo)
# - tiny real tleap builds (protein-apo, rna-apo, protein-ligand)
# - setup_md_inputs --print-only
# - PBS template rendering + bash -n
# - no CUG/CAG hard-coding
# - amber-kit init/check/build --render-only in a temp project
#
# Usage: bash tests/smoke.sh
# Set SKIP_BUILD=1 to skip the real tleap builds (faster).
set -euo pipefail

KIT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
DATA="$KIT/tests/data"
export PYTHONDONTWRITEBYTECODE=1
PY_QC=${PY_QC:-/opt/software/amber_general_workflow/envs/qc/bin/python}
TMP=$(mktemp -d /tmp/agw_smoke.XXXXXX)
# Keep bytecode out of the (possibly read-only) install tree.
export PYTHONPYCACHEPREFIX="$TMP/pycache"
mkdir -p "$PYTHONPYCACHEPREFIX"
cleanup() { rm -rf "$TMP"; }
trap cleanup EXIT
FAIL=0
ok()   { echo "OK   $*"; }
bad()  { echo "FAIL $*" >&2; FAIL=$((FAIL+1)); }

# ---------------------------------------------------------------------------
# 1. shell/python syntax
# ---------------------------------------------------------------------------
for f in "$KIT"/scripts/*.sh "$KIT"/bin/amber-kit; do
  bash -n "$f" && ok "bash -n ${f#$KIT/}" || bad "bash -n ${f#$KIT/}"
done
mkdir -p "$TMP/pycache"
for f in "$KIT"/scripts/*.py; do
  PYTHONPYCACHEPREFIX="$TMP/pycache" python3 -m py_compile "$f" \
    && ok "py_compile ${f#$KIT/}" || bad "py_compile ${f#$KIT/}"
done

# ---------------------------------------------------------------------------
# 2. render-only tleap inputs
# ---------------------------------------------------------------------------
set +u   # amber.sh reads possibly-unset LD_LIBRARY_PATH/PYTHONPATH
source /home/ylqx/amber22/amber.sh >/dev/null 2>&1 || true
set -u
command -v tleap >/dev/null 2>&1 || { echo "tleap not found; source amber.sh first" >&2; exit 2; }

run_render() {
  local name=$1; shift
  local out="$TMP/render_$name"
  if "$KIT/scripts/build_system.sh" "$@" --outdir "$out" --render-only --force \
      >"$TMP/render_$name.log" 2>&1; then
    if [ -s "$out/leap.in" ] && ! grep -q '{{' "$out/leap.in"; then
      ok "render-only $name"
    else
      bad "render-only $name (leap.in missing or has placeholders)"
    fi
  else
    tail -20 "$TMP/render_$name.log" >&2
    bad "render-only $name"
  fi
}

run_render protein_apo --type protein-apo --receptor "$DATA/tiny_peptide.pdb" \
  --ff ff19SB --solvent explicit --water opc --box oct --pad 10.0 --counterion none
run_render rna_apo --type rna-apo --receptor "$DATA/tiny_rna.pdb" \
  --ff RNA.YIL --solvent explicit --water tip3p --box oct --pad 10.0 \
  --counterion cation --cation K+
run_render protein_holo --type protein-ligand --receptor "$DATA/tiny_peptide.pdb" \
  --ligand "$DATA/lig.mol2" --frcmod "$DATA/lig.frcmod" --lig-resname LIG \
  --ligand-pdb "$DATA/ethanol.pdb" --gaff gaff2 \
  --ff ff19SB --solvent explicit --water opc --box oct --pad 10.0 --counterion none
if "$KIT/scripts/build_system.sh" --type protein-ligand --receptor "$DATA/tiny_peptide.pdb" \
    --ligand "$DATA/lig.mol2" --frcmod "$DATA/lig.frcmod" --lig-resname LIG \
    --ligand-pdb "$DATA/ethanol.pdb" --gaff gaff2 \
    --ff ff19SB --solvent explicit --water opc --box oct --pad 10.0 --counterion none \
    --outdir "$TMP/render_hmr" --render-only --hmr --force \
    >"$TMP/render_hmr.log" 2>&1; then
  bad "render-only rejects --hmr"
elif grep -q -- '--hmr requires a real build' "$TMP/render_hmr.log"; then
  ok "render-only rejects --hmr"
else
  bad "render-only rejects --hmr (unexpected error)"
fi

# ---------------------------------------------------------------------------
# 3. tiny real tleap builds
# ---------------------------------------------------------------------------
if [ "${SKIP_BUILD:-0}" != "1" ]; then
  run_build() {
    local name=$1; shift
    local out="$TMP/build_$name"
    if "$KIT/scripts/build_system.sh" "$@" --outdir "$out" --force \
        >"$TMP/build_$name.log" 2>&1; then
      if [ -s "$out/apo.prmtop" ] || [ -s "$out/complex.prmtop" ]; then
        ok "tleap build $name"
      else
        bad "tleap build $name (no prmtop)"
      fi
    else
      tail -30 "$TMP/build_$name.log" >&2
      bad "tleap build $name"
    fi
  }
  run_build protein_apo --type protein-apo --receptor "$DATA/tiny_peptide.pdb" \
    --ff ff19SB --solvent explicit --water opc --box oct --pad 10.0 --counterion none
  run_build rna_apo --type rna-apo --receptor "$DATA/tiny_rna.pdb" \
    --ff RNA.YIL --solvent explicit --water tip3p --box oct --pad 10.0 \
    --counterion cation --cation K+
  run_build protein_holo --type protein-ligand --receptor "$DATA/tiny_peptide.pdb" \
    --ligand "$DATA/lig.mol2" --frcmod "$DATA/lig.frcmod" --lig-resname LIG \
    --ligand-pdb "$DATA/ethanol.pdb" --gaff gaff2 \
    --ff ff19SB --solvent explicit --water opc --box oct --pad 10.0 --counterion none \
    --hmr
  if [ -s "$TMP/build_protein_holo/complex_hmr.prmtop" ] \
     && [ -s "$TMP/build_protein_holo/complex.prmtop.orig" ] \
     && [ -s "$TMP/build_protein_holo/hmr_manifest.json" ]; then
    ok "HMR topology"
  else
    bad "HMR topology"
  fi

  if [ -s "$TMP/build_protein_holo/complex.prmtop" ]; then
    if "$PY_QC" -I "$KIT/scripts/residue_masks.py" \
        --prmtop "$TMP/build_protein_holo/complex.prmtop" --ligand-resname LIG \
        --json "$TMP/build_protein_holo/masks.json" \
        --env "$TMP/build_protein_holo/masks.env" --force >/dev/null 2>&1; then
      grep -q "RESTRAINT_MASK=':1-3 & !@H='" "$TMP/build_protein_holo/masks.env" \
        && ok "residue_masks holo" || bad "residue_masks holo (unexpected mask)"
    else
      bad "residue_masks holo"
    fi
  fi
fi

# ---------------------------------------------------------------------------
# 4. setup_md_inputs --print-only
# ---------------------------------------------------------------------------
if python3 -I "$KIT/scripts/setup_md_inputs.py" \
    --out "$TMP/md_inputs" --seed 12345 --receptor-type rna --solvent explicit \
    --restraint-mask ':1-3 & !@H=' --heat-ns 0.1 --release-stages 2 \
    --release-ns 0.1 --prod-chunks 2 --prod-chunk-ns 0.1 \
    --hmr --dt 0.004 --traj-interval-ps 200 --print-only \
    >"$TMP/setup_md.out" 2>"$TMP/setup_md.err"; then
  grep -q '01_min_restraint.in' "$TMP/setup_md.out" \
    && grep -q '04a.in' "$TMP/setup_md.out" \
    && ok "setup_md_inputs --print-only" || bad "setup_md_inputs output"
else
  cat "$TMP/setup_md.err" >&2
  bad "setup_md_inputs --print-only"
fi

# ---------------------------------------------------------------------------
# 5. PBS rendering
# ---------------------------------------------------------------------------
python3 -I - "$KIT" "$TMP" <<'PY'
import glob, os, re, subprocess, sys
kit, tmp = sys.argv[1], sys.argv[2]
vals = dict(JOB_NAME='test_job', QUEUE='batch', NODES='1', PPN='8', GPUS='1',
            WALLTIME='24:00:00', JOIN='oe', LOG_DIR=tmp + '/logs/', MAIL='n',
            AMBER_SH='/home/ylqx/amber22/amber.sh', MPIRUN='mpirun',
            TOP='complex.prmtop', CRD='complex.inpcrd', MIN1='01_min_restraint',
            MIN2='02_min_free', HEAT='03_heat', RELEASE_STEMS='04a 04b',
            LAST_RELEASE='04b', PROD_STEMS='05_prod_c01 05_prod_c02',
            LAST_PROD='05_prod_c02', CHUNK='05_prod_c01', PREV='npt.rst7',
            REP1_DIR=tmp + '/r1', REP2_DIR=tmp + '/r2', STAGE='05_prod_c01',
            REF='complex.inpcrd', KIT_DIR=kit, LIGAND_PDB=tmp + '/lig.pdb',
            CHARGE='1', RESNAME='LIG', OUTDIR=tmp + '/lig', GAFF='gaff2',
            PROJECT_ROOT=tmp, KIT_BIN=kit + '/bin/amber-kit', LABEL='lig1', REP='1',
            REP_DIRS=tmp + '/r1 ' + tmp + '/r2', GPUS_PER_NODE='2',
            ANALYSIS_ROOT=tmp + '/analysis', PREFIX='lig1_r1',
            STRIP_MASK=':WAT,Na+,Cl-,K+', LIGAND_RESNAME='LIG',
            PARTNER_RESIDUES='', PARTNER_KIND='protein', RECEPTOR_MASK=':1-2',
            MMPBSA='1', MMPBSA_MASK=':3', MMPBSA_METHOD='gb', IGB='5',
            SALTCON='0.12', MMPBSA_INTERVAL='100', DECOMP_RESIDUES='',
            RADII='mbondi2', DRY_TOP='0',
            PY_ANALYSIS='/opt/software/amber_general_workflow/envs/analysis/bin/python',
            HBOND_STRIDE='2', NO_HBOND='0', FRAME_PS='10')
os.makedirs(tmp + '/logs', exist_ok=True)
outdir = tmp + '/pbs'
os.makedirs(outdir, exist_ok=True)
bad = 0
for t in sorted(glob.glob(kit + '/pbs/*.pbs.in')):
    text = open(t).read()
    keys = set(re.findall(r'\{\{([A-Za-z_][A-Za-z0-9_]*)\}\}', text))
    missing = keys - set(vals)
    if missing:
        print('MISSING', os.path.basename(t), sorted(missing)); bad += 1; continue
    for k in keys:
        text = text.replace('{{' + k + '}}', vals[k])
    out = os.path.join(outdir, os.path.basename(t).replace('.pbs.in', '.pbs'))
    open(out, 'w').write(text)
    rc = subprocess.run(['bash', '-n', out]).returncode
    if rc != 0 or '{{' in text:
        print('FAIL', out); bad += 1; continue
    if '#PBS -N' not in text or '#PBS -q' not in text or '#PBS -l walltime=' not in text:
        print('FAIL header', out); bad += 1; continue
    print('OK   pbs', os.path.basename(out))
sys.exit(1 if bad else 0)
PY
[ $? -eq 0 ] && ok "PBS templates" || bad "PBS templates"

# ---------------------------------------------------------------------------
# 6. no CUG/CAG hard-coding in code/templates
# ---------------------------------------------------------------------------
if grep -RInE 'rna/(cug|cag)|lig1_p2|:1-27|A1AZM|A1AZL|MQC' \
    "$KIT/scripts" "$KIT/pbs" "$KIT/bin" >"$TMP/hardcode.txt" 2>/dev/null; then
  cat "$TMP/hardcode.txt" >&2
  bad "hard-coded CUG/CAG strings"
else
  ok "no CUG/CAG hard-coding"
fi

# ---------------------------------------------------------------------------
# 7. amber-kit init/check/build --render-only
# ---------------------------------------------------------------------------
PROJ="$TMP/proj"
"$KIT/bin/amber-kit" init "$PROJ" >/dev/null 2>&1 || bad "amber-kit init"
cp "$KIT/examples/protein-ligand/config.sh" "$PROJ/config/config.sh"
printf 'label\tligand_pdb\tcharge\tresname\tmodels\treplicas\tseed_base\tcontrol\texperimental\texperimental_type\tnotes\nlig1\tinputs/lig.pdb\t0\tLIG\t\t1\t20261101\t0\tNA\tNA\ttest\n' > "$PROJ/config/systems.tsv"
cp "$DATA/tiny_peptide.pdb" "$PROJ/inputs/receptor.pdb"
cp "$DATA/ethanol.pdb" "$PROJ/inputs/lig.pdb"
mkdir -p "$PROJ/ligands/lig1/template"
cp "$DATA/lig.mol2" "$PROJ/ligands/lig1/template/lig.mol2"
cp "$DATA/lig.frcmod" "$PROJ/ligands/lig1/template/lig.frcmod"
if "$PROJ/bin/amber-kit" --project "$PROJ" check --strict >/dev/null 2>&1; then
  ok "amber-kit check"
else
  bad "amber-kit check"
fi
if "$PROJ/bin/amber-kit" --project "$PROJ" build --label lig1 --render-only --force \
    >"$TMP/kit_build.log" 2>&1; then
  [ -s "$PROJ/systems/lig1/leap.in" ] && ok "amber-kit build --render-only" \
    || bad "amber-kit build --render-only (no leap.in)"
else
  tail -20 "$TMP/kit_build.log" >&2
  bad "amber-kit build --render-only"
fi

echo
if [ "$FAIL" -eq 0 ]; then
  echo "SMOKE TEST PASSED"
  exit 0
fi
echo "SMOKE TEST FAILED: $FAIL item(s)" >&2
exit 1
