#!/usr/bin/env python3
"""Select representative + diverse models from an NMR ensemble.

This is a generic version of the CUG/CAG-style model selection:
  - core residues define the folded scaffold (Kabsch RMSD matrix)
  - medoid = model with the smallest mean RMSD to all others
  - diverse = farthest-point sampling from the medoid
  - ligand-loop / ligand-terminal distances flag end-stacked poses

The residue lists are user inputs; nothing is hard-coded to a particular
duplex or loop.

Usage:
  python -I select_models.py --pdb multi.pdb --out selected_models.tsv \
    --core-residues 1-6,8-13,14-19,21-26 --loop-residues 7,20 \
    --terminal-residues 1,13,14,26 --ligand-resname LIG --n-select 3
"""
import argparse
import os
import sys

import numpy as np
import MDAnalysis as mda

RNA_BACKBONE = ["P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
PROTEIN_BACKBONE = ["N", "CA", "C", "O"]


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def parse_residues(text):
    out = set()
    for part in text.replace(" ", ",").split(","):
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.update(range(int(a), int(b) + 1))
        else:
            out.add(int(part))
    return out


def kabsch_rmsd(P, Q):
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return float(np.sqrt(np.mean(np.sum((Pc @ R.T - Qc) ** 2, 1))))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--core-residues", required=True)
    ap.add_argument("--loop-residues", default="")
    ap.add_argument("--terminal-residues", default="")
    ap.add_argument("--ligand-resname", default="")
    ap.add_argument("--rna", action="store_true")
    ap.add_argument("--n-select", type=int, default=3)
    ap.add_argument("--end-stacked-loop-dist", type=float, default=15.0)
    ap.add_argument("--end-stacked-min-dist", type=float, default=4.5)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not os.path.isfile(args.pdb):
        die(f"PDB not found: {args.pdb}")
    if os.path.exists(args.out) and not args.force:
        die(f"refusing to overwrite {args.out} (use --force)")
    core_res = parse_residues(args.core_residues)
    loop_res = parse_residues(args.loop_residues) if args.loop_residues else set()
    term_res = parse_residues(args.terminal_residues) if args.terminal_residues else set()

    u = mda.Universe(args.pdb)
    n_models = u.trajectory.n_frames
    if n_models < 1:
        die("no models in PDB")
    if args.n_select < 1 or args.n_select > n_models:
        die(f"--n-select must be 1..{n_models}")

    backbone_names = RNA_BACKBONE if args.rna else PROTEIN_BACKBONE
    name_sel = " or ".join(f"name {n}" for n in backbone_names)
    lig_excl = f" and not resname {args.ligand_resname}" if args.ligand_resname else ""
    core = u.select_atoms(f"({name_sel}) and resid " +
                          " ".join(str(r) for r in sorted(core_res)) + lig_excl)
    if core.n_atoms < 3:
        core = u.select_atoms("not name H* and resid " +
                              " ".join(str(r) for r in sorted(core_res)) + lig_excl)
    if core.n_atoms < 3:
        die("core selection has <3 atoms")

    def atoms_for(residues):
        if not residues:
            return None
        sel = " or ".join(f"resid {r}" for r in sorted(residues))
        ag = u.select_atoms(f"({sel}) and not name H*{lig_excl}")
        return ag if ag.n_atoms else None

    loop = atoms_for(loop_res)
    term = atoms_for(term_res)
    lig = u.select_atoms(f"resname {args.ligand_resname} and not name H*") \
        if args.ligand_resname else None
    if lig is not None and lig.n_atoms == 0:
        lig = None

    coords = []
    for fi in range(n_models):
        u.trajectory[fi]
        coords.append(u.trajectory.ts.positions.copy())
    core_xyz = np.array([c[core.indices] for c in coords])

    rmsd = np.zeros((n_models, n_models))
    for i in range(n_models):
        for j in range(i + 1, n_models):
            r = kabsch_rmsd(core_xyz[i], core_xyz[j])
            rmsd[i, j] = rmsd[j, i] = r

    medoid = int(np.argmin(rmsd.sum(1)))
    selected = [medoid]
    while len(selected) < args.n_select:
        # farthest from the already selected set
        dist = np.min(rmsd[:, selected], axis=1)
        nxt = int(np.argmax(dist))
        if nxt in selected:
            break
        selected.append(nxt)

    rows = []
    for i in range(n_models):
        u.trajectory[i]
        pos = u.trajectory.ts.positions
        min_loop = np.nan
        min_term = np.nan
        end_stacked = "NA"
        if lig is not None and loop is not None:
            from scipy.spatial.distance import cdist
            min_loop = float(cdist(pos[lig.indices], pos[loop.indices]).min())
            if term is not None:
                min_term = float(cdist(pos[lig.indices], pos[term.indices]).min())
            end_stacked = int(min_loop > args.end_stacked_loop_dist or
                              (min_term < min_loop and min_term < args.end_stacked_min_dist))
        rows.append((i + 1, "medoid" if i == medoid else
                     ("diverse" if i in selected else ""),
                     float(rmsd[i, medoid]), min_loop, min_term,
                     end_stacked, int(i in selected)))

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("model\trole\tcore_rmsd_to_medoid_A\tmin_lig_loop_A\t"
                 "min_lig_terminal_A\tend_stacked\tselected\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print(f"{args.pdb}: {n_models} models; medoid {medoid + 1}; "
          f"selected {[i + 1 for i in selected]}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
