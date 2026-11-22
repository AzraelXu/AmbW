#!/usr/bin/env python3
"""Derive Amber masks from a built topology instead of hard-coding :1-N.

Reads RESIDUE_LABEL / RESIDUE_POINTER, classifies water and ions, and writes:
  restraint_mask  ":<solute ranges> & !@H="
  receptor_mask   ":<solute ranges excluding ligand>"
  ligand_mask     ":<ligand residue>"
  rna_ranges      residue ranges whose names look like RNA/DNA
plus residue lists (1-based) and a sourceable masks.env.

Usage:
  python residue_masks.py --prmtop complex.prmtop --ligand-resname LIG \
      --json masks.json --env masks.env
  python residue_masks.py --prmtop complex.prmtop \
      --partner-residues 121-180,200-210 --partner-kind protein \
      --json masks.json --env masks.env
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qc_system as qc  # noqa: E402

RNA_BASES = {"A", "C", "G", "U", "RA", "RC", "RG", "RU",
             "DA", "DC", "DG", "DT", "ADE", "CYT", "GUA", "URI"}


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def compress(nums):
    """[1,2,3,7,9,10] -> '1-3,7,9-10'."""
    nums = sorted(set(int(x) for x in nums))
    if not nums:
        return ""
    parts = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
            continue
        parts.append(str(start) if start == prev else f"{start}-{prev}")
        start = prev = n
    parts.append(str(start) if start == prev else f"{start}-{prev}")
    return ",".join(parts)


def parse_ranges(text):
    """'1-3,7,9-10' -> [1,2,3,7,9,10]."""
    nums = []
    for part in str(text).replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            nums.extend(range(int(a), int(b) + 1))
        else:
            nums.append(int(part))
    return sorted(set(nums))


def ranges_to_selection(ranges):
    """'5-20,30' -> 'resid 5-20 or resid 30' (MDAnalysis syntax)."""
    parts = [p.strip() for p in str(ranges).split(",") if p.strip()]
    return " or ".join(f"resid {p}" for p in parts)


def looks_rna(label):
    name = label.strip()
    if name in RNA_BASES:
        return True
    # tleap terminal names: G5, A5, C3, U3
    if len(name) == 2 and name[0] in "ACGU" and name[1] in "53":
        return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--inpcrd", help="optional; accepted for CLI symmetry")
    ap.add_argument("--ligand-resname", default="",
                    help="ligand residue name; omit for apo systems")
    ap.add_argument("--ligand-resid", type=int, default=None,
                    help="1-based ligand residue index (alternative to resname)")
    ap.add_argument("--partner-residues", default="",
                    help="peptide/antibody partner residues, 1-based ranges "
                         "such as '121-180,200-210'")
    ap.add_argument("--partner-kind", default="protein",
                    help="partner kind recorded in masks (protein/peptide/antibody)")
    ap.add_argument("--json", help="write JSON here")
    ap.add_argument("--env", help="write sourceable shell variables here")
    ap.add_argument("--print-only", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not os.path.isfile(args.prmtop):
        die(f"topology not found: {args.prmtop}")
    sections = qc.read_prmtop(args.prmtop)
    pointers = qc.numeric(sections, "POINTERS", int)
    if not pointers:
        die("topology lacks POINTERS")
    natom = pointers[0]
    labels, ranges = qc.residue_atom_ranges(sections, natom)
    labels = [x.strip() for x in labels]

    water = [i + 1 for i, x in enumerate(labels)
             if x.upper() in qc.WATER_NAMES]
    ions = [i + 1 for i, x in enumerate(labels)
            if x.upper() in qc.ION_NAMES]
    solute = [i + 1 for i in range(len(labels))
              if (i + 1) not in set(water) | set(ions)]
    if not solute:
        die("no solute residues found")

    lig_name = args.ligand_resname.strip()
    partner_mode = bool(args.partner_residues.strip())
    ligand = []
    if partner_mode:
        try:
            ligand = parse_ranges(args.partner_residues)
        except ValueError as exc:
            die(f"cannot parse --partner-residues: {exc}")
        bad = [r for r in ligand if r < 1 or r > len(labels)]
        if bad:
            die(f"--partner-residues out of range 1..{len(labels)}: {bad}")
        not_solute = [r for r in ligand if r not in set(solute)]
        if not_solute:
            die(f"--partner-residues includes non-solute residues: {not_solute}")
        if not ligand:
            die("--partner-residues selected no residues")
        lig_name = ""
    elif args.ligand_resid is not None:
        if args.ligand_resid < 1 or args.ligand_resid > len(labels):
            die(f"--ligand-resid out of range 1..{len(labels)}")
        if args.ligand_resid not in set(solute):
            die(f"--ligand-resid {args.ligand_resid} is water or an ion")
        ligand = [args.ligand_resid]
        lig_name = ""
    elif lig_name:
        ligand = [i + 1 for i, x in enumerate(labels) if x == lig_name]
        if len(ligand) != 1:
            die(f"expected exactly one residue named {lig_name!r}, "
                f"found {len(ligand)}: {ligand}")
    receptor = [r for r in solute if r not in set(ligand)]
    if ligand and not receptor:
        die("ligand/partner selection leaves no receptor residues")

    rna = [r for r in solute
           if r not in set(ligand) and looks_rna(labels[r - 1])]

    solute_ranges = compress(solute)
    receptor_ranges = compress(receptor)
    ligand_ranges = compress(ligand)
    rna_ranges = compress(rna)

    out = {
        "prmtop": os.path.abspath(args.prmtop),
        "n_atoms": int(natom),
        "n_residues": len(labels),
        "solute_residues": solute,
        "solute_ranges": solute_ranges,
        "restraint_mask": f":{solute_ranges} & !@H=",
        "receptor_residues": receptor,
        "receptor_ranges": receptor_ranges,
        "receptor_mask": f":{receptor_ranges}" if receptor_ranges else "",
        "ligand_residues": ligand,
        "ligand_resid": ligand[0] if (ligand and not partner_mode) else None,
        "ligand_mask": f":{ligand_ranges}" if ligand_ranges else "",
        "ligand_resname_mask": f":{lig_name}" if lig_name else "",
        "rna_residues": rna,
        "rna_ranges": rna_ranges,
        "water_residues": water,
        "ion_residues": ions,
    }
    if partner_mode:
        out.update({
            "partner_residues": ligand,
            "partner_ranges": ligand_ranges,
            "partner_mask": f":{ligand_ranges}" if ligand_ranges else "",
            "partner_kind": args.partner_kind,
            "ligand_resrange": ligand_ranges,
            "ligand_selection": ranges_to_selection(ligand_ranges),
        })
    elif args.ligand_resid is not None:
        out.update({
            "ligand_resrange": ligand_ranges,
            "ligand_selection": ranges_to_selection(ligand_ranges),
        })

    text = json.dumps(out, indent=2) + "\n"
    if args.print_only:
        print(text, end="")
    if args.json:
        if os.path.exists(args.json) and not args.force:
            die(f"refusing to overwrite {args.json} (use --force)")
        with open(args.json, "w", encoding="utf-8") as fh:
            fh.write(text)
    if args.env:
        if os.path.exists(args.env) and not args.force:
            die(f"refusing to overwrite {args.env} (use --force)")
        values = {
            "RESTRAINT_MASK": out["restraint_mask"],
            "RECEPTOR_MASK": out["receptor_mask"],
            "LIGAND_MASK": out["ligand_mask"],
            "LIGAND_RESNAME_MASK": out["ligand_resname_mask"],
            "LIGAND_RESNAME": lig_name,
            "LIGAND_RESID": str(out["ligand_resid"] or ""),
            "SOLUTE_RESIDUES": solute_ranges,
            "LAST_SOLUTE_RESID": str(max(solute)),
            "RNA_RESRANGE": rna_ranges,
        }
        if partner_mode:
            values.update({
                "PARTNER_MASK": out["partner_mask"],
                "PARTNER_RESIDUES": ligand_ranges,
                "PARTNER_RANGES": ligand_ranges,
                "PARTNER_KIND": args.partner_kind,
                "LIGAND_RESRANGE": ligand_ranges,
                "LIGAND_SELECTION": out["ligand_selection"],
            })
        elif args.ligand_resid is not None:
            values.update({
                "LIGAND_RESRANGE": ligand_ranges,
                "LIGAND_SELECTION": out["ligand_selection"],
            })
        with open(args.env, "w", encoding="utf-8") as fh:
            fh.write("# generated by residue_masks.py; source this file\n")
            for key, value in values.items():
                fh.write(f"{key}='{value}'\n")
    if not args.print_only and not args.json and not args.env:
        print(text, end="")


if __name__ == "__main__":
    main()
