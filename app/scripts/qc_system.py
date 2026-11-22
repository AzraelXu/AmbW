#!/usr/bin/env python3
"""Quality-control an Amber topology + coordinates (prmtop + inpcrd/rst7).

Reports atoms, residues, total and ligand charge, box, water model, ion counts
(including K+, which older QC scripts miss), ligand residue number, ligand pose
RMSD vs a reference PDB, minimum ligand-receptor distance, radii set, termini,
and whether a restraint mask hits water/ions. Native prmtop/inpcrd parsing means
the basic report works without parmed; parmed is used only for mask validation.

Usage:
  python3 qc_system.py --prmtop complex.prmtop --inpcrd complex.inpcrd \
      --ligand-resname LIG --ligand-pdb lig.pdb \
      --restraint-mask ':1-20 & !@H=' --strict

  python3 qc_system.py ... --json qc.json --tsv --label model_07
"""
import argparse
import json
import math
import os
import re
import sys

import numpy as np

WATER_NAMES = {
    "WAT", "HOH", "H2O", "SOL", "TIP", "TIP3", "TP3", "T3P", "SPC", "SPCE",
    "OPC", "OP3", "T4P", "TIP4", "T5P", "FB3", "FB", "DOD",
}
ION_NAMES = {
    "NA", "NA+", "K", "K+", "CL", "CL-", "MG", "MG2", "CA", "CA2", "ZN",
    "ZN2", "LI", "LI+", "RB", "RB+", "CS", "CS+", "F", "F-", "BR", "BR-",
    "I", "I-", "BA", "BA2", "SR", "SR2", "MN", "MN2", "FE", "FE2", "FE3",
    "CU", "CU2", "CO", "CO2", "NI", "NI2", "CD", "CD2", "AL", "AL3",
}
TERMINI = {"G5", "A5", "C5", "U5", "G3", "A3", "C3", "U3",
           "ACE", "NME", "NH2", "NHE"}


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def read_prmtop(path):
    sections = {}
    flag = None
    buf = []
    for raw in open(path, encoding="utf-8", errors="replace"):
        line = raw.rstrip("\n")
        if line.startswith("%FLAG"):
            if flag is not None:
                sections[flag] = buf
            flag = line.split()[1]
            buf = []
        elif line.startswith("%FORMAT") or line.startswith("%COMMENT"):
            continue
        elif flag is not None:
            buf.append(line)
    if flag is not None:
        sections[flag] = buf
    return sections


def numeric(sections, name, kind=float):
    if name not in sections:
        return []
    return [kind(x) for x in " ".join(sections[name]).split()]


def text(sections, name, count=None):
    if name not in sections:
        return []
    toks = " ".join(sections[name]).split()
    if count is None or len(toks) == count:
        return toks
    # Fall back to fixed 4-character fields (Amber a4 format).
    chars = "".join(line.ljust(80)[:80] for line in sections[name])
    return [chars[i:i + 4].strip() for i in range(0, 4 * count, 4)]


def read_inpcrd(path):
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    if len(lines) < 2:
        die(f"coordinate file too short: {path}")
    try:
        natom = int(lines[1].split()[0])
    except (IndexError, ValueError):
        die(f"cannot read atom count from {path} line 2")
    data = []
    for line in lines[2:]:
        data.extend(line.split())
    if len(data) < 3 * natom:
        die(f"{path}: expected {3 * natom} coordinate values, found {len(data)}")
    coords = np.array([float(x) for x in data[:3 * natom]], dtype=float).reshape(natom, 3)
    rest = data[3 * natom:]
    box = [float(x) for x in rest[:6]] if rest else None
    return natom, coords, box


def residue_atom_ranges(sections, natom):
    labels = text(sections, "RESIDUE_LABEL")
    ptr = numeric(sections, "RESIDUE_POINTER", int)
    if not labels or not ptr or len(labels) != len(ptr):
        die("topology lacks a usable RESIDUE_LABEL/RESIDUE_POINTER pair")
    ranges = []
    for i, start in enumerate(ptr):
        first = start - 1
        last = ptr[i + 1] - 1 if i + 1 < len(ptr) else natom
        if first < 0 or last > natom or last <= first:
            die(f"bad residue pointer at residue {i + 1}")
        ranges.append((first, last))
    return labels, ranges


def box_volume(box):
    if not box:
        return None
    if len(box) >= 6:
        a, b, c = box[:3]
        al, be, ga = box[3:6]
        # Amber solvateOct: a=b=c, alpha=beta=gamma=109.47 deg. The actual
        # truncated-octahedron volume is a^3*sqrt(2)/2, not the parallelepiped.
        if (abs(a - b) < 1e-3 and abs(b - c) < 1e-3
                and all(abs(x - 109.4712) < 0.01 for x in (al, be, ga))):
            return a ** 3 * math.sqrt(2.0) / 2.0
        al, be, ga = [math.radians(x) for x in (al, be, ga)]
        ca, cb, cg = math.cos(al), math.cos(be), math.cos(ga)
        term = 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg
        if term <= 0:
            return None
        return a * b * c * math.sqrt(term)
    if len(box) >= 3:
        return box[0] * box[1] * box[2]
    return None


def kabsch_rmsd(P, Q):
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    diff = Pc @ R.T - Qc
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def read_ligand_pdb(path):
    coords = {}
    order = []
    model = 0
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("MODEL"):
            model += 1
            if model > 1:
                break
            continue
        if line.startswith("ENDMDL") and model == 1:
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[16:17] not in (" ", "A", ""):
            continue
        name = line[12:16].strip()
        element = line[76:78].strip().upper()
        base = re.sub(r"^\d+", "", name)
        if element in ("H", "D") or (not element and base[:1].upper() in ("H", "D")):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        if name not in coords:
            coords[name] = xyz
            order.append(name)
    return coords, order


def min_distance(A, B):
    if len(A) == 0 or len(B) == 0:
        return float("nan")
    best = float("inf")
    for i in range(0, len(A), 64):
        chunk = A[i:i + 64]
        d = np.linalg.norm(chunk[:, None, :] - B[None, :, :], axis=2)
        best = min(best, float(d.min()))
    return best


def validate_mask_parmed(prmtop, mask, ligand_resname=None):
    result = {"checked": False, "atoms": None,
              "bad_residues": [], "ligand_selected": None}
    try:
        import parmed as pmd
        from parmed.amber import AmberMask
    except Exception:
        return result
    parm = pmd.load_file(prmtop)
    sel = AmberMask(parm, mask).Selection()
    if not isinstance(sel, list):
        sel = list(sel)
    result["checked"] = True
    result["atoms"] = int(sum(1 for x in sel if x))
    bad = set()
    lig_hit = False
    for atom, ok in zip(parm.atoms, sel):
        if not ok:
            continue
        name = atom.residue.name.strip()
        if name.upper() in WATER_NAMES or name.upper() in ION_NAMES:
            bad.add(name)
        if ligand_resname and name == ligand_resname:
            lig_hit = True
    result["bad_residues"] = sorted(bad)
    result["ligand_selected"] = lig_hit if ligand_resname else None
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--inpcrd", required=True)
    ap.add_argument("--ligand-resname")
    ap.add_argument("--ligand-pdb")
    ap.add_argument("--label", default="")
    ap.add_argument("--restraint-mask")
    ap.add_argument("--expect-lig-charge", type=float)
    ap.add_argument("--expect-net-charge", default="0.0",
                    help="expected total charge, or 'any' to skip")
    ap.add_argument("--expect-water", choices=["yes", "no", "any"], default="any")
    ap.add_argument("--expect-radii", help="substring, e.g. mbondi2 or mbondi3")
    ap.add_argument("--max-pose-rmsd", type=float, default=1.0)
    ap.add_argument("--max-min-dist", type=float, default=1.2,
                    help="fail if the closest ligand-receptor contact is below this (A)")
    ap.add_argument("--ion-names", help="comma-separated residue names to count as ions")
    ap.add_argument("--json", metavar="FILE")
    ap.add_argument("--tsv", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    if not os.path.isfile(args.prmtop):
        die(f"prmtop not found: {args.prmtop}")
    if not os.path.isfile(args.inpcrd):
        die(f"inpcrd not found: {args.inpcrd}")

    sections = read_prmtop(args.prmtop)
    pointers = numeric(sections, "POINTERS", int)
    if not pointers:
        die("topology has no POINTERS section")
    natom_top = int(pointers[0])
    natom_crd, coords, box = read_inpcrd(args.inpcrd)
    if natom_top != natom_crd:
        die(f"atom mismatch: topology {natom_top}, coordinates {natom_crd}")

    labels, ranges = residue_atom_ranges(sections, natom_top)
    charges = np.array(numeric(sections, "CHARGE"), dtype=float) / 18.2223
    if len(charges) != natom_top:
        die("topology CHARGE array does not match NATOM")
    masses = np.array(numeric(sections, "MASS"), dtype=float)
    atomic_numbers = numeric(sections, "ATOMIC_NUMBER", int)
    if len(atomic_numbers) == natom_top:
        is_h = np.array([z == 1 for z in atomic_numbers])
    else:
        is_h = masses < 1.5
    heavy = ~is_h

    ion_names = set(ION_NAMES)
    if args.ion_names:
        ion_names = {x.strip().upper() for x in args.ion_names.split(",") if x.strip()}
    water_res = [i for i, name in enumerate(labels)
                 if name.strip().upper() in WATER_NAMES]
    ion_res = [i for i, name in enumerate(labels)
               if name.strip().upper() in ion_names]

    water_atoms = sum(ranges[i][1] - ranges[i][0] for i in water_res)
    per_water = sorted({ranges[i][1] - ranges[i][0] for i in water_res})
    if not per_water:
        water_model = "none"
    elif len(per_water) == 1:
        water_model = {3: "3-site", 4: "4-site", 5: "5-site"}.get(per_water[0],
                                                                   f"{per_water[0]}-atom")
    else:
        water_model = "mixed:" + ",".join(str(x) for x in per_water)

    ion_counts = {}
    for i in ion_res:
        key = labels[i].strip()
        ion_counts[key] = ion_counts.get(key, 0) + 1

    lig_res = []
    if args.ligand_resname:
        lig_res = [i for i, name in enumerate(labels)
                   if name.strip() == args.ligand_resname]
        if len(lig_res) != 1:
            die(f"expected exactly one {args.ligand_resname!r} residue, "
                f"found {len(lig_res)}")
        i = lig_res[0]
        lig_atoms = np.arange(ranges[i][0], ranges[i][1])
    else:
        lig_atoms = np.array([], dtype=int)

    lig_charge = float(charges[lig_atoms].sum()) if len(lig_atoms) else None
    lig_heavy = lig_atoms[heavy[lig_atoms]] if len(lig_atoms) else np.array([], dtype=int)
    excluded = set(lig_atoms.tolist())
    for i in water_res + ion_res:
        excluded.update(range(ranges[i][0], ranges[i][1]))
    recept_atoms = np.array([i for i in range(natom_top) if i not in excluded], dtype=int)
    recept_heavy = recept_atoms[heavy[recept_atoms]] if len(recept_atoms) else np.array([], dtype=int)

    pose_rmsd = None
    if len(lig_heavy) and args.ligand_pdb:
        pdb_coords, pdb_order = read_ligand_pdb(args.ligand_pdb)
        top_names = text(sections, "ATOM_NAME", natom_top)
        by_name = {}
        for atom_idx in lig_heavy:
            name = top_names[atom_idx].strip()
            if name and name not in by_name:
                by_name[name] = atom_idx
        common = [n for n in pdb_order if n in by_name]
        if len(common) < 3:
            print("warning: fewer than 3 ligand heavy atoms matched by name; "
                  "skipping pose RMSD", file=sys.stderr)
        else:
            P = np.array([coords[by_name[n]] for n in common])
            Q = np.array([pdb_coords[n] for n in common])
            pose_rmsd = kabsch_rmsd(P, Q)

    mind = min_distance(coords[lig_heavy], coords[recept_heavy]) if len(lig_heavy) else float("nan")
    volume = box_volume(box)
    density = (float(masses.sum()) / volume * 1.660539) if volume else None
    radii = " ".join(sections.get("RADIUS_SET", [])).strip()
    termini = sorted({labels[i].strip() for i in range(len(labels))
                      if labels[i].strip() in TERMINI})

    mask_info = {"checked": False}
    if args.restraint_mask and args.restraint_mask.lower() != "none":
        mask_info = validate_mask_parmed(args.prmtop, args.restraint_mask,
                                         args.ligand_resname)

    failures = []
    warnings = []

    if str(args.expect_net_charge).lower() != "any":
        try:
            expected_q = float(args.expect_net_charge)
        except ValueError:
            die("--expect-net-charge must be a number or 'any'")
        if abs(float(charges.sum()) - expected_q) > 0.01:
            failures.append(f"net charge {charges.sum():+.3f} != expected {expected_q:+.3f}")

    if args.expect_water == "yes" and not water_res:
        failures.append("expected explicit water, found none")
    if args.expect_water == "no" and water_res:
        failures.append(f"expected no water, found {len(water_res)} residues")

    if args.ligand_resname:
        if args.expect_lig_charge is not None and abs(lig_charge - args.expect_lig_charge) > 0.01:
            failures.append(f"ligand charge {lig_charge:+.3f} != expected "
                            f"{args.expect_lig_charge:+.3f}")
        if pose_rmsd is not None and pose_rmsd > args.max_pose_rmsd:
            failures.append(f"ligand pose RMSD {pose_rmsd:.3f} A > "
                            f"{args.max_pose_rmsd:.3f} A")
        if not math.isnan(mind) and mind < args.max_min_dist:
            failures.append(f"ligand-receptor minimum distance {mind:.3f} A < "
                            f"{args.max_min_dist:.3f} A")
        if mask_info.get("checked") and mask_info.get("atoms") == 0:
            failures.append("restraint mask selects no atoms")
        if mask_info.get("bad_residues"):
            failures.append("restraint mask selects water/ions: "
                            + ",".join(mask_info["bad_residues"]))

    if args.expect_radii:
        # Match the exact "(mbondiN)" token; "mbondi" is a substring of
        # "mbondi2"/"mbondi3" and must not pass by accident.
        exact = f"({args.expect_radii})" in radii
        fallback = (args.expect_radii not in ("mbondi", "mbondi2", "mbondi3")
                    and args.expect_radii in radii)
        if not (exact or fallback):
            failures.append(f"RADIUS_SET {radii!r} does not match "
                            f"{args.expect_radii!r}")
    if not mask_info.get("checked") and args.restraint_mask:
        warnings.append("restraint mask not validated (parmed unavailable)")

    report = {
        "label": args.label,
        "prmtop": os.path.abspath(args.prmtop),
        "inpcrd": os.path.abspath(args.inpcrd),
        "atoms": natom_top,
        "residues": len(labels),
        "total_charge": round(float(charges.sum()), 4),
        "box": box,
        "box_volume_A3": round(volume, 2) if volume else None,
        "density_g_cm3": round(density, 4) if density else None,
        "water_residues": len(water_res),
        "water_atoms": water_atoms,
        "water_model": water_model,
        "ions": ion_counts,
        "ion_residues": len(ion_res),
        "receptor_atoms": int(len(recept_atoms)),
        "receptor_heavy_atoms": int(len(recept_heavy)),
        "ligand": {
            "resname": args.ligand_resname,
            "resid": int(lig_res[0] + 1) if lig_res else None,
            "atoms": int(len(lig_atoms)),
            "heavy_atoms": int(len(lig_heavy)),
            "charge": round(lig_charge, 4) if lig_charge is not None else None,
            "pose_rmsd_A": round(pose_rmsd, 4) if pose_rmsd is not None else None,
            "min_dist_A": round(mind, 4) if not math.isnan(mind) else None,
        },
        "radii": radii,
        "termini": termini,
        "restraint_mask": args.restraint_mask,
        "restraint": mask_info,
        "checks": {"failures": failures, "warnings": warnings},
    }

    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2)
            fh.write("\n")

    if args.tsv:
        cols = ["label", "atoms", "residues", "total_charge", "box_edge",
                "density_g_cm3", "water_residues", "water_model",
                "ion_residues", "lig_resid", "lig_atoms", "lig_charge",
                "lig_pose_rmsd_A", "lig_min_dist_A", "receptor_atoms",
                "radii", "restraint_atoms", "status"]
        print("\t".join(cols))
        box_edge = f"{box[0]:.2f}" if box else "NA"
        status = "FAIL" if failures else ("WARN" if warnings else "OK")
        print("\t".join(str(x) for x in [
            args.label, natom_top, len(labels), f"{charges.sum():.3f}", box_edge,
            f"{density:.4f}" if density else "NA", len(water_res), water_model,
            len(ion_res), lig_res[0] + 1 if lig_res else "NA",
            len(lig_atoms) if len(lig_atoms) else "NA",
            f"{lig_charge:.3f}" if lig_charge is not None else "NA",
            f"{pose_rmsd:.3f}" if pose_rmsd is not None else "NA",
            f"{mind:.3f}" if not math.isnan(mind) else "NA",
            len(recept_atoms), radii,
            mask_info.get("atoms", "NA"), status]))
        for msg in failures:
            print(f"FAIL: {msg}", file=sys.stderr)
        for msg in warnings:
            print(f"WARN: {msg}", file=sys.stderr)
    elif not args.quiet:
        def show(key, value):
            print(f"{key:<24} {value}")
        show("label", args.label or "-")
        show("atoms", natom_top)
        show("residues", len(labels))
        show("total charge", f"{charges.sum():+.4f}")
        show("box", " ".join(f"{x:.3f}" for x in box) if box else "none")
        show("box volume (A^3)", f"{volume:.1f}" if volume else "-")
        show("density (g/cm^3)", f"{density:.4f}" if density else "-")
        show("water", f"{len(water_res)} residues, {water_model}")
        show("ions", ", ".join(f"{k}:{v}" for k, v in sorted(ion_counts.items())) or "none")
        show("receptor atoms", f"{len(recept_atoms)} ({len(recept_heavy)} heavy)")
        if args.ligand_resname:
            show("ligand", f"{args.ligand_resname} residue {lig_res[0] + 1}: "
                           f"{len(lig_atoms)} atoms, {len(lig_heavy)} heavy, "
                           f"charge {lig_charge:+.3f}")
            if pose_rmsd is not None:
                show("ligand pose RMSD", f"{pose_rmsd:.4f} A")
            if not math.isnan(mind):
                show("min ligand-receptor", f"{mind:.3f} A")
        show("RADIUS_SET", radii or "-")
        show("termini", ", ".join(termini) or "none")
        if args.restraint_mask:
            show("restraint mask", args.restraint_mask)
            if mask_info.get("checked"):
                show("restraint atoms", mask_info["atoms"])
                if args.ligand_resname:
                    show("ligand in mask", mask_info.get("ligand_selected"))
        show("checks", "OK" if not failures and not warnings else
             f"{len(failures)} failure(s), {len(warnings)} warning(s)")
        for msg in failures:
            print(f"  FAIL: {msg}")
        for msg in warnings:
            print(f"  WARN: {msg}")

    if args.strict and failures:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
