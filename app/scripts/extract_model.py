#!/usr/bin/env python3
"""Extract one MODEL from a multi-model PDB (NMR ensembles).

No MODEL records -> the file is copied as-is (single model).
The output contains the global header lines (before the first MODEL) plus the
selected model's ATOM/HETATM/TER records and a final END. This is what
pdb4amber/tleap expect.

Usage:
  python extract_model.py --pdb multi.pdb --model 7 --out model_07.pdb
  python extract_model.py --pdb multi.pdb --list
"""
import argparse
import os
import sys


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def read_models(path):
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    has_model = any(line.startswith("MODEL") for line in lines)
    if not has_model:
        return lines, [lines]

    header = []
    models = []
    cur = None
    for line in lines:
        if line.startswith("MODEL"):
            cur = []
            continue
        if line.startswith("ENDMDL"):
            if cur is not None:
                models.append(cur)
            cur = None
            continue
        if cur is None:
            header.append(line)
        else:
            cur.append(line)
    if cur is not None:
        models.append(cur)
    return header, models


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--model", type=int, default=1)
    ap.add_argument("--out")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not os.path.isfile(args.pdb):
        die(f"PDB not found: {args.pdb}")
    header, models = read_models(args.pdb)
    if args.list:
        print(f"{args.pdb}: {len(models)} model(s)")
        return
    if args.model < 1 or args.model > len(models):
        die(f"--model {args.model} out of range 1..{len(models)}")
    if not args.out:
        die("--out is required unless --list is given")
    if os.path.exists(args.out) and not args.force:
        die(f"refusing to overwrite {args.out} (use --force)")

    keep = [line for line in models[args.model - 1]
            if line.startswith(("ATOM", "HETATM", "TER", "ANISOU"))]
    if not keep:
        die(f"model {args.model} has no ATOM/HETATM records")
    out_lines = []
    for line in header:
        if line.startswith(("ATOM", "HETATM", "MODEL", "ENDMDL")):
            continue
        out_lines.append(line)
    out_lines.extend(keep)
    out_lines.append("END")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out_lines) + "\n")
    print(f"{args.pdb}: wrote model {args.model}/{len(models)} -> {args.out}")


if __name__ == "__main__":
    main()
