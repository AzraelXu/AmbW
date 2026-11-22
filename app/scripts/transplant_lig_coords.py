#!/usr/bin/env python3
"""Move validated mol2 charges/atom types onto the coordinates of a new pose.

Atom names are the join key. The template must contain every atom of the
target pose; extra target atoms are reported and ignored. This is the safe way
to reuse one AM1-BCC/RESP parameter set across multiple docking poses without
re-running antechamber and accidentally changing charges.

Usage:
  python3 transplant_lig_coords.py TEMPLATE.mol2 TARGET.pdb OUT.mol2
"""
import sys


def read_pdb(path):
    """Return (name -> xyz, ordered names). First model, blank/A altloc only."""
    coords = {}
    order = []
    model_seen = 0
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("MODEL"):
            model_seen += 1
            if model_seen > 1:
                break
            continue
        if line.startswith("ENDMDL") and model_seen == 1:
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        alt = line[16:17]
        if alt not in (" ", "A", ""):
            continue
        name = line[12:16].strip()
        if not name:
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            sys.exit(f"bad coordinate line in {path}: {line.rstrip()}")
        if name in coords:
            # Keep the first occurrence; alternate locations are already filtered.
            continue
        coords[name] = xyz
        order.append(name)
    if not coords:
        sys.exit(f"no ATOM/HETATM records found in {path}")
    return coords, order


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    template, target, out = sys.argv[1:4]

    coords, order = read_pdb(target)
    lines = open(template, encoding="utf-8").read().splitlines()

    out_lines = []
    in_atoms = False
    transplanted = 0
    missing = []
    seen = set()
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            out_lines.append(line)
            continue
        if line.startswith("@<TRIPOS>"):
            in_atoms = False
            out_lines.append(line)
            continue
        if in_atoms and line.strip():
            fields = line.split()
            if len(fields) < 9:
                sys.exit(f"unexpected mol2 atom line: {line!r}")
            name = fields[1]
            if name not in coords:
                missing.append(name)
                out_lines.append(line)
                continue
            x, y, z = coords[name]
            fields[2] = f"{x:.4f}"
            fields[3] = f"{y:.4f}"
            fields[4] = f"{z:.4f}"
            out_lines.append(" ".join(fields))
            seen.add(name)
            transplanted += 1
        else:
            out_lines.append(line)

    if missing:
        sys.exit(f"template atoms missing from target PDB: {missing}")
    extra = [x for x in order if x not in seen]
    if extra:
        print(f"note: {len(extra)} target atoms not in mol2 (ignored): "
              f"{extra[:10]}", file=sys.stderr)

    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out_lines) + "\n")
    print(f"transplanted {transplanted} atoms from {target} -> {out}")


if __name__ == "__main__":
    main()
