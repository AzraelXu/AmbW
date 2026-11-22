#!/usr/bin/env python3
"""Verify that prmtop residue order matches prepare_complex.py's residue map.

The partner masks, disulfide indices, and MM-GBSA ligand mask all depend on
this mapping. A silent reordering by tleap would produce a wrong topology, so
this check is a hard failure.

Usage:
  python3 verify_linkage.py --prmtop complex.prmtop \
      --residue-map prepare/residue_map.tsv --out build_linkage.json \
      [--ligand] [--ligand-resname LIG]
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import qc_system as qc  # noqa: E402
import inspect_complex as ic  # noqa: E402


def canonical(label):
    s = (label or "").strip().upper()
    # tleap terminal names: NALA, CALA, etc.
    if len(s) > 3 and s[0] in ("N", "C") and s[1:] in (
            ic.AA_STANDARD | ic.AA_AMBER | ic.CAPS):
        s = s[1:]
    # histidine variants are chemically the same residue for ordering purposes
    if s in ("HID", "HIE", "HIP", "HIN"):
        s = "HIS"
    if s in ("CYX", "CYM"):
        s = "CYS"
    if s == "MSE":
        s = "MET"
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--residue-map", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ligand", action="store_true")
    ap.add_argument("--ligand-resname", default="LIG")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if os.path.exists(args.out) and not args.force:
        print(f"error: refusing to overwrite {args.out} (use --force)",
              file=sys.stderr)
        sys.exit(2)

    sections = qc.read_prmtop(args.prmtop)
    pointers = qc.numeric(sections, "POINTERS", int)
    if not pointers:
        print(f"error: {args.prmtop} lacks POINTERS", file=sys.stderr)
        sys.exit(2)
    natom = pointers[0]
    labels, ranges = qc.residue_atom_ranges(sections, natom)
    labels = [x.strip() for x in labels]
    solute = [x for x in labels
              if x.upper() not in qc.WATER_NAMES and x.upper() not in qc.ION_NAMES]

    expected = []
    roles = []
    with open(args.residue_map, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            expected.append(cols[idx["prepared_resname"]])
            roles.append(cols[idx["role"]])

    n_expected = len(expected) + (1 if args.ligand else 0)
    mismatches = []
    if len(solute) != n_expected:
        mismatches.append({
            "kind": "residue_count",
            "expected": n_expected,
            "observed": len(solute),
            "message": "solute residue count does not match prepared protein"
                       + (" + ligand" if args.ligand else ""),
        })
    for i, exp in enumerate(expected):
        if i >= len(solute):
            break
        obs = solute[i]
        if canonical(exp) != canonical(obs):
            mismatches.append({
                "kind": "residue_name",
                "prepared_index": i + 1,
                "expected": exp,
                "observed": obs,
                "role": roles[i],
            })
    ligand_index = None
    if args.ligand:
        if solute and canonical(solute[-1]) == canonical(args.ligand_resname):
            ligand_index = len(solute)
        else:
            mismatches.append({
                "kind": "ligand_position",
                "expected": args.ligand_resname,
                "observed": solute[-1] if solute else None,
                "message": "ligand is not the last solute residue",
            })

    result = {
        "status": "ok" if not mismatches else "error",
        "prmtop": os.path.abspath(args.prmtop),
        "residue_map": os.path.abspath(args.residue_map),
        "n_atoms": natom,
        "n_residues": len(labels),
        "n_solute": len(solute),
        "n_expected_protein": len(expected),
        "ligand": bool(args.ligand),
        "ligand_resname": args.ligand_resname if args.ligand else None,
        "ligand_index": ligand_index,
        "solute_labels": solute,
        "expected_labels": expected,
        "mismatches": mismatches,
    }
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2)
        fh.write("\n")
    if mismatches:
        for m in mismatches:
            print(f"error: linkage mismatch: {m}", file=sys.stderr)
        sys.exit(2)
    print(f"linkage OK: {len(expected)} protein residues"
          + (f" + ligand {args.ligand_resname}" if args.ligand else ""))


if __name__ == "__main__":
    main()
