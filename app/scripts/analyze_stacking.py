#!/usr/bin/env python3
"""Ring-stacking proxies between ligand aromatic rings and receptor rings.

For each frame and each (ligand ring, target ring) pair, computes:
  - centroid distance (A)
  - interplanar angle (deg; 0 = parallel, 90 = perpendicular)
  - projected convex-hull overlap area (A^2)

The ligand ring atom names and target residues are user inputs; nothing is
hard-coded to a particular compound or RNA loop.

Ligand ring syntax:
  "ring1:C1,C2,C3,C4,C5,C6;ring2:C9,C10,C11,C12,C13,C14"

Target residues:
  comma-separated residue numbers, e.g. "7,20,6,8,19,21".
  Ring atoms are chosen by residue name (standard RNA bases / protein aromatics)
  unless --target-atoms is given.

Usage:
  python -I analyze_stacking.py --prmtop dry.prmtop --traj dry.nc \
    --ligand-resname LIG --ligand-ring-atoms "ring1:C1,C2,C3,C4,C5,C6" \
    --target-residues 7,20 --out stack.npz [--stride 10]
"""
import argparse
import json
import os
import sys

import numpy as np
import MDAnalysis as mda
from shapely.geometry import Polygon

RING_ATOMS = {
    "A": ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"],
    "G": ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"],
    "C": ["N1", "C2", "N3", "C4", "C5", "C6"],
    "U": ["N1", "C2", "N3", "C4", "C5", "C6"],
    "T": ["N1", "C2", "N3", "C4", "C5", "C6"],
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def base_name(resname):
    name = resname.strip().upper()
    if len(name) == 2 and name[0] in "ACGU" and name[1] in "53":
        name = name[0]
    if len(name) == 3 and name[0] == "R" and name[1:] in ("A", "C", "G", "U"):
        name = name[1:]
    if len(name) == 2 and name[0] == "D" and name[1] in "ACGT":
        name = name[1]
    return name


def parse_ligand_rings(text):
    rings = {}
    for group in text.split(";"):
        group = group.strip()
        if not group:
            continue
        if ":" not in group:
            die(f"ligand ring must be name:atom1,atom2,...; got {group!r}")
        name, atoms = group.split(":", 1)
        names = [x.strip() for x in atoms.split(",") if x.strip()]
        if not name.strip() or len(names) < 3:
            die(f"ligand ring {group!r} needs a name and >=3 atoms")
        rings[name.strip()] = names
    if not rings:
        die("--ligand-ring-atoms is empty")
    return rings


def plane(P):
    c = P.mean(0)
    _, _, vt = np.linalg.svd(P - c)
    return c, vt[2], vt[0], vt[1]


def overlap(Plig, Pbase):
    c, n, e1, e2 = plane(Pbase)

    def proj(P):
        Q = P - c
        return np.column_stack([Q @ e1, Q @ e2])

    p1 = Polygon(proj(Plig)).convex_hull
    p2 = Polygon(proj(Pbase)).convex_hull
    if not p1.is_valid or not p2.is_valid:
        return 0.0
    return float(p1.intersection(p2).area)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--ligand-resname", default="LIG")
    ap.add_argument("--ligand-ring-atoms", required=True)
    ap.add_argument("--target-residues", required=True)
    ap.add_argument("--target-atoms", default="",
                    help="override ring atom names for all targets (comma list)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--stride", type=int, default=10)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=-1)
    ap.add_argument("--frame-ps", type=float, default=10.0)
    ap.add_argument("--summary", help="optional JSON summary")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if args.stride < 1:
        die("--stride must be >= 1")
    if os.path.exists(args.out) and not args.force:
        die(f"refusing to overwrite {args.out} (use --force)")
    if not os.path.isfile(args.prmtop) or not os.path.isfile(args.traj):
        die("prmtop or trajectory not found")

    u = mda.Universe(args.prmtop, args.traj)
    names = np.array([n.strip() for n in u.atoms.names])
    resids = np.array(u.atoms.resids)
    resnames = np.array([r.strip() for r in u.atoms.resnames])

    lig_rings = {}
    for ring_name, atom_names in parse_ligand_rings(args.ligand_ring_atoms).items():
        idx = []
        for a in atom_names:
            hit = np.where((resnames == args.ligand_resname) & (names == a))[0]
            if len(hit) != 1:
                print(f"warning: ligand ring {ring_name}: atom {a} not found; skipped",
                      file=sys.stderr)
                idx = []
                break
            idx.append(int(hit[0]))
        if idx:
            lig_rings[ring_name] = np.array(idx, dtype=int)

    try:
        targets = [int(x) for x in args.target_residues.replace(",", " ").split()]
    except ValueError:
        die(f"cannot parse --target-residues: {args.target_residues!r}")

    override = [x.strip() for x in args.target_atoms.split(",") if x.strip()]
    tgt = {}
    for r in targets:
        sel = np.where(resids == r)[0]
        if len(sel) == 0:
            print(f"warning: target residue {r} not found; skipped", file=sys.stderr)
            continue
        rn = base_name(str(resnames[sel[0]]))
        atom_names = override or RING_ATOMS.get(rn)
        if not atom_names:
            print(f"warning: no ring atom definition for residue {r} ({rn}); skipped",
                  file=sys.stderr)
            continue
        idx = []
        for a in atom_names:
            hit = sel[names[sel] == a]
            if len(hit) != 1:
                print(f"warning: residue {r} ({rn}) lacks atom {a}; skipped",
                      file=sys.stderr)
                idx = []
                break
            idx.append(int(hit[0]))
        if idx:
            tgt[r] = np.array(idx, dtype=int)

    if not lig_rings or not tgt:
        die("no usable ligand ring / target ring pairs")

    pairs = [(lr, tr) for lr in lig_rings for tr in sorted(tgt)]
    end = u.trajectory.n_frames if args.end < 0 else min(args.end, u.trajectory.n_frames)
    if args.start < 0 or args.start >= end:
        die(f"bad frame range start={args.start} end={end}")
    frames = list(range(args.start, end, args.stride))
    D = np.full((len(frames), len(pairs)), np.nan)
    A = np.full((len(frames), len(pairs)), np.nan)
    O = np.full((len(frames), len(pairs)), np.nan)

    for k, fi in enumerate(frames):
        u.trajectory[fi]
        pos = u.trajectory.ts.positions
        for j, (lr, tr) in enumerate(pairs):
            P1 = pos[lig_rings[lr]]
            P2 = pos[tgt[tr]]
            c1, n1, _, _ = plane(P1)
            c2, n2, _, _ = plane(P2)
            D[k, j] = np.linalg.norm(c1 - c2)
            A[k, j] = np.degrees(np.arccos(np.clip(abs(np.dot(n1, n2)), -1, 1)))
            O[k, j] = overlap(P1, P2)

    labels = np.array([f"{lr}--{tr}" for lr, tr in pairs])
    np.savez_compressed(
        args.out,
        labels=labels,
        frames=np.array(frames, dtype=int),
        time_ns=np.array(frames, dtype=float) * args.frame_ps / 1000.0,
        dist=D,
        angle=A,
        overlap=O,
        ligand_resname=args.ligand_resname,
        target_residues=np.array(sorted(tgt), dtype=int),
        stride=args.stride,
    )

    summary = {
        "prmtop": os.path.abspath(args.prmtop),
        "traj": os.path.abspath(args.traj),
        "ligand_resname": args.ligand_resname,
        "ligand_rings": {k: v.tolist() for k, v in lig_rings.items()},
        "target_residues": sorted(tgt),
        "n_frames": len(frames),
        "stride": args.stride,
        "pairs": [{"label": lab, "dist_mean": float(np.nanmean(D[:, j])),
                   "angle_mean": float(np.nanmean(A[:, j])),
                   "overlap_mean": float(np.nanmean(O[:, j]))}
                  for j, lab in enumerate(labels)],
    }
    if args.summary:
        with open(args.summary, "w", encoding="utf-8") as fh:
            json.dump(summary, fh, indent=2)
            fh.write("\n")

    print(f"{len(frames)} frames x {len(pairs)} ring pairs -> {args.out}")
    for j, lab in enumerate(labels):
        print(f"  {lab:24s} dist {np.nanmean(D[:, j]):5.2f} A  "
              f"angle {np.nanmean(A[:, j]):5.1f} deg  overlap {np.nanmean(O[:, j]):6.1f} A^2")


if __name__ == "__main__":
    main()
