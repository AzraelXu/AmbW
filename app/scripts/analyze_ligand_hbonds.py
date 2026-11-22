#!/usr/bin/env python3
"""Bidirectional ligand-receptor hydrogen bonds for a dry Amber trajectory.

Counts:
  - ligand donor (N/O/S-H) -> receptor acceptor (N/O/S)
  - receptor donor (N/O/S-H) -> ligand acceptor (N/O/S)

Criteria (defaults, all configurable):
  H...A < 2.5 A, D...A < 3.5 A, D-H...A > 120 deg.

Donor pairs are assigned once from the first frame by the nearest N/O/S within
1.3 A of each H. This assumes rigid X-H bonds, which is true for SHAKE runs.

Usage:
  python -I analyze_ligand_hbonds.py --prmtop dry.prmtop --traj dry.nc \
    --ligand-resname LIG --out hb.npz [--stride 2] \
    [--receptor-selection "all and not resname LIG"] \
    [--ha-cutoff 2.5] [--da-cutoff 3.5] [--angle 120]

  python -I analyze_ligand_hbonds.py --prmtop dry.prmtop --traj dry.nc \
    --partner-residues 61-68 --out hb.npz --stride 10
"""
import argparse
import json
import os
import sys

import numpy as np
import MDAnalysis as mda
from scipy.spatial.distance import cdist

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from partner_selection import select_partner  # noqa: E402

DONOR_ELEMENTS = {"N", "O", "S"}
ACCEPTOR_ELEMENTS = {"N", "O", "S"}


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def elements_of(ag):
    elems = np.array([(a.element or "").upper() for a in ag.atoms])
    # Amber prmtop names sometimes lack element guesses; fall back to name.
    for i, (e, name) in enumerate(zip(elems, ag.names)):
        if not e:
            n = name.strip().upper()
            if n.startswith(("H", "D")):
                elems[i] = "H"
            elif n.startswith("C"):
                elems[i] = "C"
            elif n.startswith("N"):
                elems[i] = "N"
            elif n.startswith("O"):
                elems[i] = "O"
            elif n.startswith("S"):
                elems[i] = "S"
            elif n.startswith("P"):
                elems[i] = "P"
    return elems


def donor_pairs(ag, elements, max_bond=1.3):
    """Return (h_indices, heavy_indices) in the AtomGroup's local indexing."""
    heavy = np.where(np.isin(elements, list(DONOR_ELEMENTS)))[0]
    hydro = np.where(elements == "H")[0]
    if len(heavy) == 0 or len(hydro) == 0:
        return np.array([], int), np.array([], int)
    # Reference frame coordinates; X-H bonds do not change under SHAKE.
    pos = ag.positions
    h_keep, d_keep = [], []
    for h in hydro:
        d = np.linalg.norm(pos[heavy] - pos[h], axis=1)
        j = int(np.argmin(d))
        if d[j] <= max_bond:
            h_keep.append(h)
            d_keep.append(heavy[j])
    return np.array(h_keep, int), np.array(d_keep, int)


def angle_ok(p, d, h, a, min_angle):
    v1 = p[d] - p[h]
    v2 = p[a] - p[h]
    denom = np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12
    cosang = np.dot(v1, v2) / denom
    return np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))) >= min_angle


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--ligand-resname", default="LIG")
    ap.add_argument("--partner-residues", default="",
                    help="peptide/antibody partner residue ranges, e.g. '121-180'")
    ap.add_argument("--partner-kind", default="protein")
    ap.add_argument("--receptor-selection", default=None,
                    help="MDAnalysis selection; default 'all and not resname <LIG>'")
    ap.add_argument("--ha-cutoff", type=float, default=2.5)
    ap.add_argument("--da-cutoff", type=float, default=3.5)
    ap.add_argument("--angle", type=float, default=120.0)
    ap.add_argument("--summary", help="optional JSON summary")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if args.stride < 1:
        die("--stride must be >= 1")
    if not os.path.isfile(args.prmtop) or not os.path.isfile(args.traj):
        die("prmtop or trajectory not found")
    if os.path.exists(args.out) and not args.force:
        die(f"refusing to overwrite {args.out} (use --force)")

    u = mda.Universe(args.prmtop, args.traj)
    if args.partner_residues.strip():
        try:
            lig, rec, rec_sel = select_partner(
                u, args.prmtop, args.partner_residues)
        except ValueError as exc:
            die(str(exc))
        if args.receptor_selection:
            rec = u.select_atoms(args.receptor_selection)
    else:
        lig = u.select_atoms(f"resname {args.ligand_resname}")
        if lig.n_atoms == 0:
            die(f"no ligand atoms with resname {args.ligand_resname!r}")
        rec_sel = args.receptor_selection or f"all and not resname {args.ligand_resname}"
        rec = u.select_atoms(rec_sel)
    if rec.n_atoms == 0:
        die(f"receptor selection selects no atoms: {rec_sel!r}")

    lig_el = elements_of(lig)
    rec_el = elements_of(rec)
    lig_h, lig_d = donor_pairs(lig, lig_el)
    rec_h, rec_d = donor_pairs(rec, rec_el)

    lig_acc = np.where(np.isin(lig_el, list(ACCEPTOR_ELEMENTS)))[0]
    rec_acc = np.where(np.isin(rec_el, list(ACCEPTOR_ELEMENTS)))[0]

    rec_resids = rec.resids
    rec_resnames = np.array([r.strip() for r in rec.resnames])
    rec_unique = sorted(set(int(x) for x in rec_resids))
    rec_index = {r: i for i, r in enumerate(rec_unique)}
    rec_atom_res = np.array([rec_index[int(r)] for r in rec_resids])

    print(f"ligand: {lig.n_atoms} atoms, donors {len(lig_h)}, acceptors {len(lig_acc)}")
    print(f"receptor: {rec.n_atoms} atoms, donors {len(rec_h)}, acceptors {len(rec_acc)}")
    print(f"frames: {u.trajectory.n_frames}, stride {args.stride}")

    counts, counts_rev = [], []
    phos = []
    donor_occ = np.zeros(len(lig_h), dtype=float)
    donor_occ_rev = np.zeros(len(rec_h), dtype=float)
    res_occ = np.zeros(len(rec_unique), dtype=float)
    nframes = 0

    for fi in range(0, u.trajectory.n_frames, args.stride):
        u.trajectory[fi]
        pos = u.trajectory.ts.positions
        hb = 0
        if len(lig_h) and len(rec_acc):
            d_ha = cdist(pos[lig.atoms.indices[lig_h]], pos[rec.atoms.indices[rec_acc]])
            ii, jj = np.where(d_ha < args.ha_cutoff)
            for i, j in zip(ii, jj):
                hi = lig_h[i]
                ai = rec_acc[j]
                di = lig_d[i]
                p = pos
                if np.linalg.norm(p[lig.atoms.indices[di]] - p[rec.atoms.indices[ai]]) > args.da_cutoff:
                    continue
                if angle_ok(p, lig.atoms.indices[di], lig.atoms.indices[hi],
                            rec.atoms.indices[ai], args.angle):
                    hb += 1
                    donor_occ[i] += 1
                    res_occ[rec_atom_res[ai]] += 1

        hb_rev = 0
        if len(rec_h) and len(lig_acc):
            d_ha = cdist(pos[rec.atoms.indices[rec_h]], pos[lig.atoms.indices[lig_acc]])
            ii, jj = np.where(d_ha < args.ha_cutoff)
            for i, j in zip(ii, jj):
                hi = rec_h[i]
                ai = lig_acc[j]
                di = rec_d[i]
                p = pos
                if np.linalg.norm(p[rec.atoms.indices[di]] - p[lig.atoms.indices[ai]]) > args.da_cutoff:
                    continue
                if angle_ok(p, rec.atoms.indices[di], rec.atoms.indices[hi],
                            lig.atoms.indices[ai], args.angle):
                    hb_rev += 1
                    donor_occ_rev[i] += 1
                    res_occ[rec_atom_res[di]] += 1

        counts.append(hb)
        counts_rev.append(hb_rev)
        nframes += 1

    counts = np.array(counts, dtype=float)
    counts_rev = np.array(counts_rev, dtype=float)
    np.savez_compressed(
        args.out,
        hbond_count=counts,
        hbond_count_rev=counts_rev,
        donor_occ=donor_occ / max(nframes, 1),
        donor_occ_rev=donor_occ_rev / max(nframes, 1),
        ligand_donor_h_idx=lig_h,
        ligand_donor_heavy_idx=lig_d,
        receptor_donor_h_idx=rec_h,
        receptor_donor_heavy_idx=rec_d,
        receptor_resids=np.array(rec_unique),
        receptor_resnames=np.array([rec_resnames[rec_resids == r][0] for r in rec_unique]),
        receptor_res_occ=res_occ / max(nframes, 1),
        ha_cutoff=args.ha_cutoff,
        da_cutoff=args.da_cutoff,
        angle=args.angle,
    )

    summary = {
        "prmtop": os.path.abspath(args.prmtop),
        "traj": os.path.abspath(args.traj),
        "frames": int(nframes),
        "stride": args.stride,
        "ligand_resname": args.ligand_resname,
        "partner_residues": args.partner_residues or None,
        "partner_kind": args.partner_kind if args.partner_residues else None,
        "receptor_selection": rec_sel,
        "ha_cutoff_A": args.ha_cutoff,
        "da_cutoff_A": args.da_cutoff,
        "angle_deg": args.angle,
        "ligand_donors": int(len(lig_h)),
        "receptor_donors": int(len(rec_h)),
        "mean_hbond_lig_to_rec": float(counts.mean()) if nframes else None,
        "mean_hbond_rec_to_lig": float(counts_rev.mean()) if nframes else None,
    }
    if args.summary:
        with open(args.summary, "w", encoding="utf-8") as fh:
            json.dump(summary, fh, indent=2)
            fh.write("\n")

    print(f"{nframes} frames -> {args.out}")
    print(f"mean H-bonds lig->rec {counts.mean():.2f}; rec->lig {counts_rev.mean():.2f}")


if __name__ == "__main__":
    main()
