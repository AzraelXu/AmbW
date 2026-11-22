#!/usr/bin/env python3
"""Add hydrogens (if needed), perceive bond orders, and infer the net charge
of a small-molecule ligand from its 3D PDB coordinates.

This script must run under a Python that has RDKit. On this machine that is
/opt/software/amber_general_workflow/envs/analysis/bin/python (RDKit 2024.03). OpenBabel is used first when
the ligand has no hydrogens: RDKit's bond-order perception on a heavy-atom-only
PDB can produce chemically absurd charged solutions instead of failing.

The result is never guessed silently: if the charge is ambiguous, the script
exits non-zero and lists the candidates. Users can force a charge with
--charge or validate a protonation/bond-order template with --smiles.

Usage:
  python ligand_charge.py --pdb ligand_input.pdb --out-pdb ligand_prepared.pdb \
      --json ligand_charge.json [--charge INT] [--smiles SMILES] \
      [--ph 7.4|none] [--obabel PATH] [--resname LIG] [--force]
"""
import argparse
import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdDetermineBonds, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
except Exception as exc:  # pragma: no cover - environment check
    print(f"error: RDKit is required for automatic ligand charge inference "
          f"({exc}). Use --charge with another Python or install RDKit.",
          file=sys.stderr)
    sys.exit(2)


SMARTS = {
    "carboxylate": Chem.MolFromSmarts("[CX3](=O)[O-]"),
    "phosphate": Chem.MolFromSmarts("[PX4](=O)([O-])"),
    "sulfate": Chem.MolFromSmarts("[SX4](=O)(=O)[O-]"),
    "amine": Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]"),
    "guanidinium": Chem.MolFromSmarts("[NX3][CX3](=[NX2])[NX3]"),
}


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def warn(msg):
    print(f"warning: {msg}", file=sys.stderr)


def count_h(mol):
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1)


def count_heavy(mol):
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)


def load_pdb(path):
    if not os.path.isfile(path):
        die(f"ligand PDB not found: {path}")
    mol = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
    if mol is None:
        die(f"RDKit cannot read {path}; check the PDB atom names/columns")
    return mol


def add_h_with_obabel(input_pdb, obabel, ph):
    if not obabel:
        return None, None
    if not os.path.isfile(obabel) and not shutil.which(obabel):
        return None, None
    tmpdir = tempfile.mkdtemp(prefix="ligand_charge_")
    out = os.path.join(tmpdir, "ligand_h.pdb")
    cmd = [obabel, input_pdb, "-O", out]
    if ph is None or str(ph).lower() == "none":
        cmd.append("-h")
        provenance = "obabel -h"
    else:
        cmd += ["-p", str(ph)]
        provenance = f"obabel -p {ph}"
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True)
    if proc.returncode != 0 or not os.path.isfile(out):
        warn(f"OpenBabel failed ({provenance}): {proc.stdout.strip()[:400]}")
        return None, None
    mol = Chem.MolFromPDBFile(out, removeHs=False, sanitize=False)
    if mol is None or count_h(mol) == 0:
        warn(f"OpenBabel produced no readable hydrogens ({provenance})")
        return None, None
    return mol, provenance


def obabel_to_mol2(input_pdb, obabel, ph):
    """Let OpenBabel perceive bond orders and add H; return a sanitized RDKit mol.

    OpenBabel's 3D bond-order perception is far more reliable than RDKit's
    charge enumeration on a heavy-atom-only PDB, so the mol2 path is primary
    when hydrogens are absent.
    """
    if not obabel or (not os.path.isfile(obabel) and not shutil.which(obabel)):
        return None, None
    tmpdir = tempfile.mkdtemp(prefix="ligand_charge_")
    attempts = []
    if ph is not None and str(ph).lower() != "none":
        attempts.append((["-p", str(ph)], f"obabel -p {ph}"))
    attempts.append((["-h"], "obabel -h"))
    for idx, (extra, provenance) in enumerate(attempts):
        out = os.path.join(tmpdir, f"ligand_h_{idx}.mol2")
        cmd = [obabel, input_pdb, "-O", out] + extra
        proc = subprocess.run(cmd, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True)
        if proc.returncode != 0 or not os.path.isfile(out):
            warn(f"OpenBabel failed ({provenance}): "
                 f"{proc.stdout.strip()[:400]}")
            continue
        mol = Chem.MolFromMol2File(out, removeHs=False, sanitize=False)
        if mol is None:
            warn(f"RDKit cannot read the OpenBabel mol2 ({provenance})")
            continue
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:
            warn(f"OpenBabel mol2 is chemically invalid ({provenance}): {exc}")
            continue
        if count_h(mol) == 0:
            warn(f"OpenBabel mol2 contains no hydrogens ({provenance})")
            continue
        return mol, provenance
    return None, None


def molecule_summary(mol):
    heavy = Chem.RemoveHs(mol)
    return {
        "charge": sum(a.GetFormalCharge() for a in mol.GetAtoms()),
        "smiles": Chem.MolToSmiles(heavy, canonical=True),
        "smiles_noiso": Chem.MolToSmiles(heavy, canonical=True,
                                          isomericSmiles=False),
        "inchikey": Chem.MolToInchiKey(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
    }


def has_match(mol, key):
    patt = SMARTS[key]
    return patt is not None and mol.HasSubstructMatch(patt)


def score_candidate(mol, charge):
    atoms = mol.GetAtoms()
    sum_fc = sum(abs(a.GetFormalCharge()) for a in atoms)
    n_charged = sum(1 for a in atoms if a.GetFormalCharge() != 0)
    n_arom = sum(1 for a in atoms if a.GetIsAromatic())
    n_rad = sum(a.GetNumRadicalElectrons() for a in atoms)
    n_nonarom_ring = sum(1 for a in atoms
                         if a.IsInRing() and not a.GetIsAromatic())
    score = (3.0 * sum_fc + 2.0 * n_charged - 2.0 * n_arom
             + 100.0 * n_rad + 5.0 * n_nonarom_ring)
    # Charged carbons are almost always a sign of a bad bond-order solution.
    score += 5.0 * sum(1 for a in atoms
                       if a.GetSymbol() == "C" and a.GetFormalCharge() != 0)
    if charge < 0 and (has_match(mol, "carboxylate")
                       or has_match(mol, "phosphate")
                       or has_match(mol, "sulfate")):
        score -= 3.0
    if charge > 0 and (has_match(mol, "amine")
                       or has_match(mol, "guanidinium")):
        score -= 3.0
    return score


def atom_signature(mol):
    """Heavy-atom connectivity signature: atom count and heavy bond pairs."""
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    pairs = set()
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if mol.GetAtomWithIdx(i).GetAtomicNum() > 1 and \
                mol.GetAtomWithIdx(j).GetAtomicNum() > 1:
            pairs.add(tuple(sorted((i, j))))
    return len(heavy), frozenset(pairs)


def _enumerate_candidates_quiet(base, qmin, qmax):
    sig = atom_signature(base)
    n_atoms = base.GetNumAtoms()
    n_h = count_h(base)
    out = []
    for q in sorted(range(qmin, qmax + 1), key=lambda x: (abs(x), x)):
        cand = Chem.Mol(base)
        try:
            rdDetermineBonds.DetermineBondOrders(cand, charge=q)
            Chem.SanitizeMol(cand)
        except Exception:
            continue
        if cand.GetNumAtoms() != n_atoms or count_h(cand) != n_h:
            continue
        if atom_signature(cand) != sig:
            continue
        if len(Chem.GetMolFrags(cand)) != 1:
            continue
        atoms = cand.GetAtoms()
        if sum(a.GetFormalCharge() for a in atoms) != q:
            continue
        if any(abs(a.GetFormalCharge()) > 2 for a in atoms):
            continue
        if any(a.GetNumRadicalElectrons() for a in atoms):
            continue
        try:
            heavy = Chem.RemoveHs(cand)
            smiles = Chem.MolToSmiles(heavy, canonical=True)
            smiles_noiso = Chem.MolToSmiles(heavy, canonical=True,
                                            isomericSmiles=False)
            inchi = Chem.MolToInchiKey(cand)
            formula = rdMolDescriptors.CalcMolFormula(cand)
        except Exception:
            smiles, smiles_noiso, inchi, formula = "", "", "", ""
        out.append({
            "charge": q,
            "score": round(score_candidate(cand, q), 3),
            "smiles": smiles,
            "smiles_noiso": smiles_noiso,
            "inchikey": inchi,
            "formula": formula,
            "mol": cand,
        })
    return out


def enumerate_candidates(base, qmin=-8, qmax=8):
    """Try every charge in [qmin,qmax]; return accepted candidates.

    RDKit's C++ bond-order code can write molecule diagnostics to fd 1 on
    failed charge attempts. Redirect that fd while enumerating so stdout stays
    valid JSON for callers.
    """
    saved = None
    devnull = None
    try:
        saved = os.dup(1)
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, 1)
        with contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            return _enumerate_candidates_quiet(base, qmin, qmax)
    finally:
        if saved is not None:
            os.dup2(saved, 1)
            os.close(saved)
        if devnull is not None:
            os.close(devnull)


def pick_candidate(cands, requested_charge=None, template_smiles=None,
                   template_charge=None):
    if not cands:
        die("cannot assign chemically valid bond orders and a net charge to "
            "this ligand. Candidates tried q=-8..+8. Inspect ligand_input.pdb "
            "or pass --charge/--ligand-smiles.")
    if requested_charge is not None:
        match = [c for c in cands if c["charge"] == requested_charge]
        if not match:
            listing = "; ".join(f"q={c['charge']} score={c['score']} "
                                f"{c['smiles']}" for c in cands)
            die(f"--charge {requested_charge} is not chemically valid for this "
                f"ligand. Valid candidates: {listing}")
        if len(match) > 1:
            match.sort(key=lambda c: c["score"])
            warn(f"multiple bond-order solutions for q={requested_charge}; "
                 f"using the lowest-score one ({match[0]['smiles']})")
        return match[0], "user_charge"
    if template_smiles is not None:
        match = [c for c in cands if c["charge"] == template_charge
                 and c.get("smiles_noiso", c["smiles"]) == template_smiles]
        if not match:
            listing = "; ".join(f"q={c['charge']} {c['smiles']}" for c in cands)
            die("--ligand-smiles does not match the inferred ligand structure. "
                f"Template={template_smiles} (q={template_charge}); "
                f"candidates: {listing}")
        return match[0], "user_smiles"
    cands = sorted(cands, key=lambda c: (c["score"], abs(c["charge"])))
    best = cands[0]
    if len(cands) > 1:
        second = cands[1]
        margin = second["score"] - best["score"]
        if margin < 2.0:
            listing = "; ".join(f"q={c['charge']} score={c['score']} "
                                f"{c['smiles']}" for c in cands)
            die("ligand charge/bond-order assignment is ambiguous "
                f"(score margin {margin:.2f}): {listing}. "
                "Rerun with --charge or --ligand-smiles.")
        if best["charge"] != 0 and margin < 4.0:
            listing = "; ".join(f"q={c['charge']} score={c['score']} "
                                f"{c['smiles']}" for c in cands)
            die(f"inferred charge {best['charge']} is not confidently "
                f"separated from the next candidate (margin {margin:.2f}): "
                f"{listing}. Confirm with --charge.")
    return best, "rdkit_enumeration"


def name_field(name, element):
    name = "".join(ch for ch in name if ch.isalnum() or ch in "'+-")
    if not name:
        name = element
    if len(name) > 4:
        name = name[:4]
    if len(name) == 1 and len(element) == 1:
        # PDB convention: one-letter element names start in column 14
        return (" " + name).ljust(4)
    if name[0].isdigit():
        return name.rjust(4)
    return name.ljust(4)


def write_pdb(mol, path, resname="LIG", chain="L", resid=1):
    conf = mol.GetConformer()
    used = {}
    lines = []
    for i, atom in enumerate(mol.GetAtoms(), 1):
        info = atom.GetPDBResidueInfo()
        raw_name = info.GetName().strip() if info is not None else ""
        if not raw_name:
            raw_name = atom.GetSymbol() + str(i)
        base = raw_name
        if base in used:
            used[base] += 1
            raw_name = f"{base[:3]}{used[base] % 10}"
        else:
            used[base] = 0
        # ensure uniqueness after truncation
        candidate = raw_name[:4]
        k = 1
        while candidate in used and candidate != base:
            candidate = f"{base[:3]}{k % 10}"
            k += 1
        used.setdefault(candidate, 0)
        element = atom.GetSymbol().upper()
        name = name_field(candidate, element)
        pos = conf.GetAtomPosition(i - 1)
        occ = info.GetOccupancy() if info is not None else 1.0
        bfac = info.GetTempFactor() if info is not None else 0.0
        lines.append(
            f"ATOM  {i:5d} {name}{'':1s}{resname:>3s} {chain:1s}"
            f"{resid:4d}{'':1s}   {pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
            f"{occ:6.2f}{bfac:6.2f}          {element:>2s}")
    lines.append("TER")
    lines.append("END")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True, help="single-residue ligand PDB")
    ap.add_argument("--out-pdb", default=None,
                    help="write the protonated ligand PDB here")
    ap.add_argument("--json", default=None, help="write inference JSON here")
    ap.add_argument("--charge", type=int, default=None)
    ap.add_argument("--smiles", default=None,
                    help="template SMILES used to validate bond orders/charge")
    ap.add_argument("--ph", default="7.4",
                    help="pH for OpenBabel protonation; 'none' uses plain -h")
    ap.add_argument("--obabel", default="obabel")
    ap.add_argument("--resname", default="LIG")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if args.out_pdb and os.path.exists(args.out_pdb) and not args.force:
        die(f"refusing to overwrite {args.out_pdb} (use --force)")
    if args.json and os.path.exists(args.json) and not args.force:
        die(f"refusing to overwrite {args.json} (use --force)")

    mol = load_pdb(args.pdb)
    heavy = count_heavy(mol)
    if heavy < 2:
        die(f"ligand has only {heavy} heavy atom(s); not a small molecule")
    if heavy > 200:
        die(f"ligand has {heavy} heavy atoms; amber-auto supports <=200")

    template_charge = None
    template_smiles = None
    if args.smiles:
        template = Chem.MolFromSmiles(args.smiles)
        if template is None:
            die(f"cannot parse --smiles: {args.smiles}")
        template_charge = sum(a.GetFormalCharge() for a in template.GetAtoms())
        template_smiles = Chem.MolToSmiles(template, canonical=True,
                                           isomericSmiles=False)
        if args.charge is not None and args.charge != template_charge:
            die(f"--charge {args.charge} disagrees with the formal charge "
                f"{template_charge} in --smiles")

    warnings = []
    ranked = []
    margin = float("inf")

    if count_h(mol) == 0:
        # Heavy-atom-only crystal ligand: trust OpenBabel's 3D bond-order
        # perception and pH protonation; use RDKit only for validation.
        omol, provenance = obabel_to_mol2(args.pdb, args.obabel, args.ph)
        if omol is not None:
            summary = molecule_summary(omol)
            q = summary["charge"]
            if args.charge is not None and q != args.charge:
                die(f"OpenBabel inferred charge {q} from the ligand geometry, "
                    f"which disagrees with --charge {args.charge}. Provide an "
                    "H-containing ligand with the desired protonation, or use "
                    "--ligand-smiles to document the intended state.")
            if template_smiles is not None and (
                    q != template_charge
                    or summary["smiles_noiso"] != template_smiles):
                die("OpenBabel's inferred structure does not match --ligand-smiles. "
                    f"OpenBabel: {summary['smiles']} (q={q}); "
                    f"template: {template_smiles} (q={template_charge}).")
            best = dict(summary)
            best["score"] = round(score_candidate(omol, q), 3)
            best["mol"] = omol
            ranked = [best]
            method = "openbabel_mol2"
            confidence = "high"
            warnings.append(f"hydrogens added by {provenance}; protonation is "
                            "heuristic and was recorded")
        else:
            if args.charge is None:
                die("ligand has no hydrogens and OpenBabel could not perceive "
                    "bond orders/protonation. Add hydrogens externally or pass "
                    "--charge together with an H-containing ligand.")
            base = Chem.Mol(mol)
            try:
                rdDetermineBonds.DetermineConnectivity(base)
                base = Chem.AddHs(base, addCoords=True)
            except Exception as exc:
                die(f"cannot add hydrogens with RDKit: {exc}")
            if count_h(base) == 0:
                die("RDKit could not add hydrogens; provide an H-containing ligand")
            cands = enumerate_candidates(base)
            best, method = pick_candidate(
                cands, requested_charge=args.charge,
                template_smiles=template_smiles,
                template_charge=template_charge)
            method = "rdkit_addH_user_charge"
            ranked = sorted(cands, key=lambda c: (c["score"], abs(c["charge"])))
            if len(ranked) > 1:
                margin = ranked[1]["score"] - ranked[0]["score"]
            confidence = "low"
            warnings.append("ligand had no hydrogens; hydrogens were added by "
                            "RDKit and the charge came from --charge")
    else:
        provenance = "input hydrogens"
        base = Chem.Mol(mol)
        try:
            rdDetermineBonds.DetermineConnectivity(base)
        except Exception as exc:
            die(f"RDKit DetermineConnectivity failed: {exc}")
        cands = enumerate_candidates(base)
        best, method = pick_candidate(
            cands, requested_charge=args.charge,
            template_smiles=template_smiles,
            template_charge=template_charge)
        ranked = sorted(cands, key=lambda c: (c["score"], abs(c["charge"])))
        if len(ranked) > 1:
            margin = ranked[1]["score"] - ranked[0]["score"]
        confidence = "high" if margin >= 6 else "medium" if margin >= 2 else "low"

    result = {
        "charge": best["charge"],
        "smiles": best["smiles"],
        "inchikey": best["inchikey"],
        "formula": best["formula"],
        "method": method,
        "provenance": provenance,
        "heavy_atoms": count_heavy(best["mol"]),
        "hydrogens": count_h(best["mol"]),
        "confidence": confidence,
        "score_margin": None if margin == float("inf") else round(margin, 3),
        "candidates": [
            {k: c[k] for k in ("charge", "score", "smiles", "smiles_noiso",
                               "inchikey", "formula")}
            for c in ranked
        ],
        "warnings": warnings,
    }

    text = json.dumps(result, indent=2) + "\n"
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            fh.write(text)
    if args.out_pdb:
        write_pdb(best["mol"], args.out_pdb, resname=args.resname)
    print(text, end="")
    print(f"ligand charge: q={best['charge']} method={method} "
          f"confidence={confidence}", file=sys.stderr)
    if args.out_pdb:
        print(f"wrote {args.out_pdb}", file=sys.stderr)


if __name__ == "__main__":
    main()
