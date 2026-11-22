#!/usr/bin/env python3
"""Split one MODEL of a complex PDB into receptor, ligand, and a reference
complex. The ligand residue is renamed (default LIG), moved to chain L, and
renumbered to residue 1; atom names and coordinates are preserved. Ligand
protonation is NOT changed: the input must already have the desired hydrogens.

Usage:
  python split_complex.py --pdb complex.pdb --model 7 \
      --ligand-resname LIG --receptor-out rec.pdb --ligand-out lig.pdb
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from extract_model import read_models  # noqa: E402


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def pdb_resname(line):
    return line[17:20].strip()


def set_fields(line, resname=None, chain=None, resseq=None):
    rec = line[0:6]
    name = line[12:16]
    alt = line[16]
    rn = f"{resname:>3}" if resname is not None else line[17:20]
    ch = chain if chain is not None else line[21]
    rs = f"{resseq:>4}" if resseq is not None else line[22:26]
    ic = line[26]
    return f"{rec}{line[6:12]}{name}{alt}{rn} {ch}{rs}{ic}{line[27:]}"


def write_pdb(path, atoms, ter_after_chain=True):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        last_chain = None
        for line in atoms:
            if len(line) < 22:
                continue
            ch = line[21]
            if ter_after_chain and last_chain is not None and ch != last_chain:
                fh.write("TER\n")
            fh.write(line.rstrip("\n") + "\n")
            last_chain = ch
        fh.write("TER\nEND\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--model", type=int, default=1)
    ap.add_argument("--ligand-resname", default="LIG",
                    help="comma-separated residue names, e.g. LIG,COF")
    ap.add_argument("--out-resname", default="LIG")
    ap.add_argument("--receptor-out", default="receptor.pdb")
    ap.add_argument("--ligand-out", default="ligand.pdb")
    ap.add_argument("--complex-ref-out", default="complex_ref.pdb")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not os.path.isfile(args.pdb):
        die(f"PDB not found: {args.pdb}")
    if len(args.out_resname.strip()) < 1 or len(args.out_resname.strip()) > 4:
        die("--out-resname must be 1-4 characters")
    lig_names = [x.strip() for x in args.ligand_resname.split(",") if x.strip()]
    if not lig_names:
        die("--ligand-resname is empty")

    header, models = read_models(args.pdb)
    if args.model < 1 or args.model > len(models):
        die(f"--model {args.model} out of range 1..{len(models)}")
    atoms = [line for line in models[args.model - 1]
             if line.startswith(("ATOM", "HETATM"))]
    if not atoms:
        die(f"model {args.model} has no ATOM/HETATM records")

    receptor = []
    ligand = []
    lig_resids = set()
    for line in atoms:
        rn = pdb_resname(line)
        if rn in lig_names:
            ligand.append(set_fields(line, resname=args.out_resname,
                                     chain="L", resseq=1))
            lig_resids.add((line[21], line[22:27]))
        else:
            receptor.append(line)

    if not ligand:
        die(f"no ligand residue {lig_names} found in model {args.model}; "
            "for an apo system use extract_model.py instead")
    if len(lig_resids) != 1:
        die(f"found {len(lig_resids)} ligand residues ({sorted(lig_resids)}); "
            "param_ligand.sh requires exactly one ligand residue. Split "
            "multi-residue ligands manually first.")
    if not receptor:
        die("no receptor atoms found")

    outs = [args.receptor_out, args.ligand_out, args.complex_ref_out]
    if not args.force:
        for path in outs:
            if os.path.exists(path):
                die(f"refusing to overwrite {path} (use --force)")

    write_pdb(args.receptor_out, receptor)
    write_pdb(args.ligand_out, ligand)
    write_pdb(args.complex_ref_out, receptor + ligand)
    print(f"model {args.model}: receptor {len(receptor)} atoms, "
          f"ligand {len(ligand)} atoms ({args.out_resname})")
    print(f"  {args.receptor_out}\n  {args.ligand_out}\n  {args.complex_ref_out}")


if __name__ == "__main__":
    main()
