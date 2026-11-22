#!/bin/bash
# =============================================================================
# amber-auto smoke tests: detection, ligand charge, preparation, build, PBS
# dependency planning, status and report. No jobs are submitted to Torque.
#
#   bash tests/auto_smoke.sh
#   SKIP_BUILD=1 bash tests/auto_smoke.sh     # fast, no tleap/antechamber
# =============================================================================
set -u
KIT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$KIT"
TMP=$(mktemp -d /tmp/amber_auto_smoke.XXXXXX)
trap 'rm -rf "$TMP"' EXIT
# Keep bytecode out of the (possibly read-only) install tree.
export PYTHONPYCACHEPREFIX="$TMP/pycache"
mkdir -p "$PYTHONPYCACHEPREFIX"
PASS=0
FAIL=0
ok() { PASS=$((PASS+1)); echo "OK   $*"; }
bad() { FAIL=$((FAIL+1)); echo "FAIL $*"; }

# 1. syntax / compile -------------------------------------------------------
if python3 -m py_compile bin/amber-auto scripts/inspect_complex.py \
    scripts/prepare_complex.py scripts/ligand_charge.py \
    scripts/partner_selection.py scripts/verify_linkage.py; then
  ok "py_compile amber-auto helpers"
else
  bad "py_compile amber-auto helpers"
fi
if bash -n scripts/postprocess_one.sh; then ok "bash -n postprocess_one.sh"; else bad "bash -n postprocess_one.sh"; fi

# 2. detection unit tests ---------------------------------------------------
# Run the file directly: a regular `tests` package installed in site-packages
# can shadow this directory's namespace package for `-m unittest tests.<mod>`.
if python3 tests/test_detection.py >"$TMP/detection.out" 2>&1; then
  ok "detection unit tests"
else
  cat "$TMP/detection.out" >&2
  bad "detection unit tests"
fi

# 3. ligand charge inference ------------------------------------------------
PY_RDKIT=/opt/software/amber_general_workflow/envs/analysis/bin/python
OBABEL=/opt/schrodinger2021-2/utilities/obabel
if [ -x "$PY_RDKIT" ] && { [ -x "$OBABEL" ] || command -v obabel >/dev/null 2>&1; }; then
  "$PY_RDKIT" -I scripts/ligand_charge.py --pdb tests/data/auto_benzoate_noH.pdb \
    --out-pdb "$TMP/benzoate.pdb" --json "$TMP/benzoate.json" --force \
    --obabel "$OBABEL" >/dev/null 2>"$TMP/benzoate.err"
  if [ $? -eq 0 ] && python3 -c "import json; d=json.load(open('$TMP/benzoate.json')); raise SystemExit(0 if d['charge']==-1 else 1)"; then
    ok "charge inference: benzoate -1"
  else
    cat "$TMP/benzoate.err" >&2
    bad "charge inference: benzoate -1"
  fi
  "$PY_RDKIT" -I scripts/ligand_charge.py --pdb tests/data/auto_neutral_noH.pdb \
    --out-pdb "$TMP/neutral.pdb" --json "$TMP/neutral.json" --force \
    --obabel "$OBABEL" >/dev/null 2>"$TMP/neutral.err"
  if [ $? -eq 0 ] && python3 -c "import json; d=json.load(open('$TMP/neutral.json')); raise SystemExit(0 if d['charge']==0 else 1)"; then
    ok "charge inference: neutral 0"
  else
    cat "$TMP/neutral.err" >&2
    bad "charge inference: neutral 0"
  fi
  if "$PY_RDKIT" -I scripts/ligand_charge.py --pdb tests/data/auto_benzoate_noH.pdb \
      --out-pdb "$TMP/wrong.pdb" --json "$TMP/wrong.json" --force \
      --charge 0 --obabel "$OBABEL" >/dev/null 2>&1; then
    bad "charge inference rejects a wrong --charge"
  else
    ok "charge inference rejects a wrong --charge"
  fi
else
  echo "SKIP charge inference (RDKit/OpenBabel not found)"
fi

# 4. preparation: MSE, disulfide, glycan/metal escape hatches --------------
python3 scripts/inspect_complex.py --pdb tests/data/auto_mse.pdb \
  --json-out "$TMP/mse.json" --quiet
if python3 scripts/prepare_complex.py --pdb tests/data/auto_mse.pdb \
    --components "$TMP/mse.json" --outdir "$TMP/prep_mse" --no-reduce --force \
    >/dev/null 2>&1 && grep -q 'MET' "$TMP/prep_mse/protein_prepared.pdb" \
    && grep -q ' SD ' "$TMP/prep_mse/protein_prepared.pdb"; then
  ok "prepare: MSE -> MET/SD"
else
  bad "prepare: MSE -> MET/SD"
fi

python3 scripts/inspect_complex.py --pdb tests/data/auto_ss.pdb \
  --json-out "$TMP/ss.json" --quiet
if python3 scripts/prepare_complex.py --pdb tests/data/auto_ss.pdb \
    --components "$TMP/ss.json" --outdir "$TMP/prep_ss" --no-reduce --force \
    >/dev/null 2>&1 && grep -q 'bond mol.1.SG mol.2.SG' "$TMP/prep_ss/extra_leap.in"; then
  ok "prepare: CYS -> CYX + disulfide bond"
else
  bad "prepare: CYS -> CYX + disulfide bond"
fi

if ./bin/amber-auto inspect tests/data/auto_glycan.pdb --campaign "$TMP/inspect_gly" >/dev/null 2>&1; then
  bad "amber-auto inspect rejects glycans"
else
  ok "amber-auto inspect rejects glycans"
fi
if ./bin/amber-auto inspect tests/data/auto_zn.pdb --campaign "$TMP/inspect_zn" --drop-metals >/dev/null 2>&1; then
  ok "amber-auto inspect --drop-metals"
else
  bad "amber-auto inspect --drop-metals"
fi

# 5. build + PBS dependency planning ---------------------------------------
if [ "${SKIP_BUILD:-0}" = "1" ]; then
  echo "SKIP build (SKIP_BUILD=1)"
else
  if ./bin/amber-auto build tests/data/auto_pl.pdb --campaign "$TMP/camp" \
      --force --prod-ns 0.1 --prod-chunks 1 --pad 8 \
      >"$TMP/build.out" 2>"$TMP/build.err"; then
    C="$TMP/camp/complexes/auto_pl"
    if [ -s "$C/systems/auto_pl/complex.prmtop" ] \
       && grep -q "LIGAND_MASK=':3'" "$C/systems/auto_pl/masks.env" \
       && python3 -c "import json; d=json.load(open('$C/md/auto_pl/r1/run_config.json')); raise SystemExit(0 if len(d['production']['files'])==1 else 1)"; then
      ok "amber-auto build (protein-ligand)"
    else
      bad "amber-auto build artifacts"
    fi
    if [ -s "$C/md/auto_pl/r1/run_post.pbs" ] \
       && bash -n "$C/md/auto_pl/r1/run_post.pbs"; then
      ok "run_post.pbs rendered"
    else
      bad "run_post.pbs rendered"
    fi
  else
    cat "$TMP/build.err" >&2
    bad "amber-auto build (protein-ligand)"
  fi
  if ./bin/amber-auto resume --campaign "$TMP/camp" --dry-run \
      >"$TMP/resume.out" 2>"$TMP/resume.err"; then
    if [ -s "$TMP/camp/resume_plan.json" ] \
       && python3 - "$TMP/camp/resume_plan.json" <<'PY'
import json, sys
plan = json.load(open(sys.argv[1]))
assert [p["kind"] for p in plan] == ["equil", "prod", "post"]
assert plan[1]["depends_on"] == plan[0]["jobid"]
assert plan[2]["depends_on"] == plan[1]["jobid"]
PY
    then ok "resume --dry-run dependency chain"; else bad "resume --dry-run dependency chain"; fi
  else
    cat "$TMP/resume.err" >&2
    bad "resume --dry-run"
  fi
  if ./bin/amber-auto status --campaign "$TMP/camp" >"$TMP/status.out" 2>&1 \
      && grep -q "auto_pl" "$TMP/status.out"; then
    ok "amber-auto status"
  else
    bad "amber-auto status"
  fi
  if ./bin/amber-auto report --campaign "$TMP/camp" >/dev/null 2>&1 \
      && [ -s "$TMP/camp/report/summary.tsv" ]; then
    ok "amber-auto report"
  else
    bad "amber-auto report"
  fi
  if ./bin/amber-auto build tests/data/auto_pep.pdb --campaign "$TMP/camp_pep" \
      --force --box rect --pad 8 --prod-ns 0.1 --prod-chunks 1 \
      >"$TMP/pep_build.out" 2>"$TMP/pep_build.err"; then
    P="$TMP/camp_pep/complexes/auto_pep"
    if [ -s "$P/systems/auto_pep/complex.prmtop" ] \
       && grep -q "PARTNER_MASK=':61-68'" "$P/systems/auto_pep/masks.env"; then
      ok "amber-auto build (protein-peptide)"
    else
      bad "amber-auto build protein-peptide artifacts"
    fi
  else
    cat "$TMP/pep_build.err" >&2
    bad "amber-auto build (protein-peptide)"
  fi
fi

echo
if [ "$FAIL" -eq 0 ]; then
  echo "AUTO SMOKE TEST PASSED ($PASS checks)"
  exit 0
fi
echo "AUTO SMOKE TEST FAILED ($FAIL failures, $PASS passed)" >&2
exit 1
