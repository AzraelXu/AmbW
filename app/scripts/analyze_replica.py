#!/usr/bin/env python3
"""Generic per-replica structural analysis for protein/RNA-ligand MD.

Reads a dry Amber topology/trajectory (cpptraj autoimage + strip output) and
writes a compressed npz plus an optional JSON summary.

Metrics (all optional except receptor RMSD/Rg):
  - receptor backbone RMSD after Kabsch alignment to the first analyzed frame
  - ligand heavy-atom RMSD after receptor alignment (pose stability)
  - ligand internal RMSD (conformational change)
  - receptor Rg, ligand Rg, ligand COM drift
  - ligand-receptor minimum distance, contact atom/residue counts, per-residue
    contact occupancy, buried fraction
  - bidirectional ligand-receptor H-bonds
  - loop-ligand contacts (--loop-residues)
  - receptor LCPO SASA via cpptraj surf (--sasa --receptor-mask; optional)

The Python selections are MDAnalysis selections, not Amber masks. Use
--receptor-mask only for the optional cpptraj SASA calculation.

Usage:
  python -I analyze_replica.py --prmtop dry.prmtop --traj dry.nc \
    --label lig1 --rep 1 --ligand-resname LIG \
    --out analysis/python/lig1_r1.npz --summary analysis/python/lig1_r1.json
  python -I analyze_replica.py --prmtop dry.prmtop --traj dry.nc \
    --label pep1 --rep 1 --partner-residues 61-68 --partner-kind protein-peptide \
    --out analysis/python/pep1_r1.npz --summary analysis/python/pep1_r1.json
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import MDAnalysis as mda
from scipy.spatial.distance import cdist

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from analyze_ligand_hbonds import angle_ok, donor_pairs, elements_of  # noqa: E402
from partner_selection import select_partner  # noqa: E402

RNA_BACKBONE = ["P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
PROTEIN_BACKBONE = ["N", "CA", "C", "O"]


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def kabsch(P, Q):
    """Superpose P onto Q. Returns R, P_center, Q_center, RMSD."""
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    rms = float(np.sqrt(np.mean(np.sum((Pc @ R.T - Qc) ** 2, 1))))
    return R, P.mean(0), Q.mean(0), rms


def radius_of_gyration(positions, masses):
    com = np.average(positions, axis=0, weights=masses)
    return float(np.sqrt(np.sum(masses * np.sum((positions - com) ** 2, axis=1)) / masses.sum()))


def selection_for_names(names):
    return " or ".join(f"name {n}" for n in names)


def hbond_counts(pos, lig_global, rec_global, lig_h, lig_d, lig_acc,
                 rec_h, rec_d, rec_acc, ha_cut, da_cut, min_angle):
    """Count ligand->receptor and receptor->ligand H-bonds for one frame."""
    n_lr = 0
    if len(lig_h) and len(rec_acc):
        d = cdist(pos[lig_global[lig_h]], pos[rec_global[rec_acc]])
        ii, jj = np.where(d < ha_cut)
        for i, j in zip(ii, jj):
            dg = lig_global[lig_d[i]]
            hg = lig_global[lig_h[i]]
            ag = rec_global[rec_acc[j]]
            if np.linalg.norm(pos[dg] - pos[ag]) > da_cut:
                continue
            if angle_ok(pos, dg, hg, ag, min_angle):
                n_lr += 1
    n_rl = 0
    if len(rec_h) and len(lig_acc):
        d = cdist(pos[rec_global[rec_h]], pos[lig_global[lig_acc]])
        ii, jj = np.where(d < ha_cut)
        for i, j in zip(ii, jj):
            dg = rec_global[rec_d[i]]
            hg = rec_global[rec_h[i]]
            ag = lig_global[lig_acc[j]]
            if np.linalg.norm(pos[dg] - pos[ag]) > da_cut:
                continue
            if angle_ok(pos, dg, hg, ag, min_angle):
                n_rl += 1
    return n_lr, n_rl


def cpptraj_sasa(prmtop, traj, mask, stride, n_frames):
    """Return SASA per original frame via cpptraj 'surf' (LCPO), or None."""
    if not shutil.which("cpptraj"):
        print("warning: cpptraj not in PATH; SASA skipped", file=sys.stderr)
        return None
    tmpdir = tempfile.mkdtemp(prefix="agw_sasa_")
    inp = os.path.join(tmpdir, "sasa.in")
    dat = os.path.join(tmpdir, "sasa.dat")
    log = os.path.join(tmpdir, "sasa.log")
    with open(inp, "w", encoding="utf-8") as fh:
        fh.write(f"parm {prmtop}\n")
        fh.write(f"trajin {traj} 1 last {stride}\n")
        fh.write(f"surf receptor '{mask}' out {dat}\n")
        fh.write("run\nquit\n")
    try:
        with open(log, "w", encoding="utf-8") as fh:
            rc = subprocess.run(["cpptraj", "-i", inp], stdout=fh,
                                stderr=subprocess.STDOUT).returncode
        if rc != 0 or not os.path.isfile(dat):
            print("warning: cpptraj surf failed; SASA skipped", file=sys.stderr)
            return None
        out = np.full(n_frames, np.nan)
        row = 0
        with open(dat, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split()
                if len(f) < 2:
                    continue
                try:
                    value = float(f[1])
                except ValueError:
                    continue
                # cpptraj writes output-frame indices 1..N, not original
                # trajectory frame numbers. With "trajin ... 1 last stride",
                # output row i corresponds to original frame i*stride.
                idx = row * stride
                if 0 <= idx < n_frames:
                    out[idx] = value
                row += 1
        return out
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--label", default="replica")
    ap.add_argument("--rep", default="")
    ap.add_argument("--ligand-resname", default="LIG")
    ap.add_argument("--partner-residues", default="",
                    help="peptide/antibody partner residue ranges, e.g. '121-180'")
    ap.add_argument("--partner-kind", default="protein")
    ap.add_argument("--apo", action="store_true",
                    help="no ligand/partner; receptor-only metrics")
    ap.add_argument("--receptor-selection", default=None)
    ap.add_argument("--rna", action="store_true")
    ap.add_argument("--loop-residues", default="")
    ap.add_argument("--contact-cutoff", type=float, default=4.0)
    ap.add_argument("--contact-stride", type=int, default=10)
    ap.add_argument("--hbond-ha-cutoff", type=float, default=2.5)
    ap.add_argument("--hbond-da-cutoff", type=float, default=3.5)
    ap.add_argument("--hbond-angle", type=float, default=120.0)
    ap.add_argument("--hbond-stride", type=int, default=1)
    ap.add_argument("--sasa", action="store_true")
    ap.add_argument("--sasa-stride", type=int, default=10)
    ap.add_argument("--receptor-mask", default="", help="Amber mask for cpptraj surf")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=-1)
    ap.add_argument("--frame-ps", type=float, default=10.0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if args.stride < 1 or args.contact_stride < 1 or args.hbond_stride < 1:
        die("--stride/--contact-stride/--hbond-stride must be >= 1")
    if not os.path.isfile(args.prmtop) or not os.path.isfile(args.traj):
        die("prmtop or trajectory not found")
    if os.path.exists(args.out) and not args.force:
        die(f"refusing to overwrite {args.out} (use --force)")

    u = mda.Universe(args.prmtop, args.traj)
    n_total = u.trajectory.n_frames
    start = args.start
    end = n_total if args.end < 0 else min(args.end, n_total)
    if start < 0 or start >= end:
        die(f"bad frame range start={start} end={end} n_frames={n_total}")
    frame_indices = list(range(start, end, args.stride))

    if args.apo:
        lig = u.select_atoms("resname __AMBER_AUTO_APO__")
        holo = False
        rec_sel = args.receptor_selection or "all"
        rec = u.select_atoms(rec_sel)
    elif args.partner_residues.strip():
        try:
            lig, rec, rec_sel = select_partner(
                u, args.prmtop, args.partner_residues)
        except ValueError as exc:
            die(str(exc))
        holo = lig.n_atoms > 0
        if args.receptor_selection:
            rec = u.select_atoms(args.receptor_selection)
    else:
        lig = u.select_atoms(f"resname {args.ligand_resname}")
        holo = lig.n_atoms > 0
        rec_sel = args.receptor_selection
        if not rec_sel:
            rec_sel = f"all and not resname {args.ligand_resname}" if holo else "all"
        rec = u.select_atoms(rec_sel)
    if rec.n_atoms == 0:
        die(f"receptor selection selects no atoms: {rec_sel!r}")

    lig_el = elements_of(lig) if holo else np.array([])
    rec_el = elements_of(rec)
    lig_heavy = np.where(lig_el != "H")[0] if holo else np.array([], int)
    rec_heavy = np.where(rec_el != "H")[0]

    backbone_names = RNA_BACKBONE if args.rna else PROTEIN_BACKBONE
    backbone = rec.select_atoms(selection_for_names(backbone_names))
    if backbone.n_atoms < 3:
        backbone = rec.select_atoms("not name H*")
    if backbone.n_atoms < 3:
        die("cannot find enough receptor backbone/heavy atoms for alignment")

    loop_atoms = None
    if args.loop_residues.strip():
        try:
            loop_res = [int(x) for x in args.loop_residues.replace(",", " ").split()]
        except ValueError:
            die(f"cannot parse --loop-residues: {args.loop_residues!r}")
        sel = " or ".join(f"resid {r}" for r in loop_res)
        loop_atoms = rec.select_atoms(sel)
        if loop_atoms.n_atoms == 0:
            die(f"--loop-residues selects no receptor atoms: {loop_res}")

    # Global index arrays (positions are full-system coordinates).
    rec_global = rec.atoms.indices
    lig_global = lig.atoms.indices if holo else np.array([], int)
    rec_heavy_global = rec_global[rec_heavy]
    lig_heavy_global = lig_global[lig_heavy] if holo else np.array([], int)
    # backbone.atoms.indices are global indices.
    rec_backbone_global = backbone.atoms.indices

    # Donor/acceptor pairs (assigned on the first frame; X-H bonds are rigid).
    u.trajectory[start]
    lig_h, lig_d, lig_acc = np.array([], int), np.array([], int), np.array([], int)
    rec_h, rec_d, rec_acc = np.array([], int), np.array([], int), np.array([], int)
    if holo:
        lig_h, lig_d = donor_pairs(lig, lig_el)
        lig_acc = np.where(np.isin(lig_el, ["N", "O", "S"]))[0]
    rec_h, rec_d = donor_pairs(rec, rec_el)
    rec_acc = np.where(np.isin(rec_el, ["N", "O", "S"]))[0]

    # Reference frame.
    ref_pos = u.trajectory.ts.positions.copy()
    ref_backbone = ref_pos[rec_backbone_global]
    ref_lig = ref_pos[lig_heavy_global] if holo else None
    ref_lig_com = ref_lig.mean(0) if holo and len(ref_lig) else None

    rec_resids = np.array([int(r) for r in rec.resids])
    rec_resnames = np.array([r.strip() for r in rec.resnames])
    rec_unique = sorted(set(int(x) for x in rec_resids))
    rec_index = {r: i for i, r in enumerate(rec_unique)}
    rec_atom_res = np.array([rec_index[int(r)] for r in rec_resids])

    # Per-frame arrays (NaN = not computed at this frame).
    n = len(frame_indices)
    receptor_rmsd = np.full(n, np.nan)
    ligand_rmsd = np.full(n, np.nan)
    ligand_internal_rmsd = np.full(n, np.nan)
    receptor_rg = np.full(n, np.nan)
    ligand_rg = np.full(n, np.nan)
    ligand_com_drift = np.full(n, np.nan)
    min_dist = np.full(n, np.nan)
    n_contact_atoms = np.full(n, np.nan)
    n_contact_residues = np.full(n, np.nan)
    buried_frac = np.full(n, np.nan)
    loop_min_dist = np.full(n, np.nan)
    loop_contact = np.full(n, np.nan)
    hbond_lr = np.full(n, np.nan)
    hbond_rl = np.full(n, np.nan)
    receptor_sasa = np.full(n, np.nan)
    contact_occ_sum = np.zeros(len(rec_unique), dtype=float)
    min_dist_sum = np.zeros(len(rec_unique), dtype=float)
    n_contact_frames = 0

    rec_masses = rec.atoms.masses[rec_heavy]
    lig_masses = lig.atoms.masses[lig_heavy] if holo else np.array([])

    sasa_full = None
    if args.sasa:
        if not args.receptor_mask:
            print("warning: --sasa needs --receptor-mask; SASA skipped", file=sys.stderr)
        else:
            sasa_full = cpptraj_sasa(args.prmtop, args.traj, args.receptor_mask,
                                     args.sasa_stride, n_total)

    for k, fi in enumerate(frame_indices):
        u.trajectory[fi]
        pos = u.trajectory.ts.positions
        R, Pm, Qm, rrms = kabsch(pos[rec_backbone_global], ref_backbone)
        receptor_rmsd[k] = rrms
        receptor_rg[k] = radius_of_gyration(pos[rec_heavy_global], rec_masses)

        if holo and len(lig_heavy_global):
            lig_aligned = (pos[lig_heavy_global] - Pm) @ R.T + Qm
            ligand_rmsd[k] = float(np.sqrt(np.mean(np.sum((lig_aligned - ref_lig) ** 2, 1))))
            _, _, _, lrms = kabsch(pos[lig_heavy_global], ref_lig)
            ligand_internal_rmsd[k] = lrms
            ligand_rg[k] = radius_of_gyration(pos[lig_heavy_global], lig_masses)
            lig_com = pos[lig_heavy_global].mean(0)
            ligand_com_drift[k] = float(np.linalg.norm(lig_com - ref_lig_com))

        if k % args.contact_stride == 0:
            if holo and len(lig_heavy_global):
                D = cdist(pos[lig_heavy_global], pos[rec_heavy_global])
                min_dist[k] = float(D.min())
                within = D < args.contact_cutoff
                n_contact_atoms[k] = int(within.any(axis=0).sum())
                res_hit = np.zeros(len(rec_unique), dtype=bool)
                res_min = np.full(len(rec_unique), np.inf)
                for j in range(D.shape[1]):
                    col = D[:, j]
                    if col.min() < args.contact_cutoff:
                        res_hit[rec_atom_res[rec_heavy[j]]] = True
                    ridx = rec_atom_res[rec_heavy[j]]
                    if col.min() < res_min[ridx]:
                        res_min[ridx] = col.min()
                n_contact_residues[k] = int(res_hit.sum())
                buried_frac[k] = float(within.any(axis=1).mean())
                contact_occ_sum += res_hit
                min_dist_sum += np.where(np.isfinite(res_min), res_min, 0.0)
                n_contact_frames += 1
                if loop_atoms is not None and loop_atoms.n_atoms:
                    Dl = cdist(pos[lig_heavy_global], pos[loop_atoms.indices])
                    loop_min_dist[k] = float(Dl.min())
                    loop_contact[k] = int((Dl.min(0) < args.contact_cutoff).sum())

        if k % args.hbond_stride == 0:
            if holo:
                lr, rl = hbond_counts(pos, lig_global, rec_global, lig_h, lig_d, lig_acc,
                                      rec_h, rec_d, rec_acc, args.hbond_ha_cutoff,
                                      args.hbond_da_cutoff, args.hbond_angle)
                hbond_lr[k] = lr
                hbond_rl[k] = rl

        if sasa_full is not None:
            receptor_sasa[k] = sasa_full[fi]

    contact_occ = contact_occ_sum / max(n_contact_frames, 1)
    mean_min = min_dist_sum / max(n_contact_frames, 1)

    out = dict(
        label=args.label,
        rep=args.rep,
        frames=np.array(frame_indices, dtype=int),
        time_ns=np.array(frame_indices, dtype=float) * args.frame_ps / 1000.0,
        receptor_rmsd=receptor_rmsd,
        ligand_rmsd=ligand_rmsd,
        ligand_internal_rmsd=ligand_internal_rmsd,
        receptor_rg=receptor_rg,
        ligand_rg=ligand_rg,
        ligand_com_drift=ligand_com_drift,
        min_dist=min_dist,
        n_contact_atoms=n_contact_atoms,
        n_contact_residues=n_contact_residues,
        buried_frac=buried_frac,
        loop_min_dist=loop_min_dist,
        loop_contact=loop_contact,
        hbond_lig_to_rec=hbond_lr,
        hbond_rec_to_lig=hbond_rl,
        receptor_sasa=receptor_sasa,
        receptor_resids=np.array(rec_unique),
        receptor_resnames=np.array([rec_resnames[rec_resids == r][0] for r in rec_unique]),
        contact_occupancy=contact_occ,
        mean_min_dist_per_res=mean_min,
        contact_cutoff=args.contact_cutoff,
        hbond_ha_cutoff=args.hbond_ha_cutoff,
        hbond_da_cutoff=args.hbond_da_cutoff,
        hbond_angle=args.hbond_angle,
        receptor_selection=rec_sel,
        ligand_resname=args.ligand_resname,
        partner_residues=args.partner_residues,
        partner_kind=args.partner_kind,
        apo=args.apo,
        frame_ps=args.frame_ps,
        n_frames_total=n_total,
    )
    np.savez_compressed(args.out, **out)

    def stat(x):
        x = np.asarray(x, dtype=float)
        x = x[np.isfinite(x)]
        if x.size == 0:
            return None
        return {"mean": float(x.mean()), "sd": float(x.std(ddof=1)) if x.size > 1 else 0.0,
                "median": float(np.median(x)), "min": float(x.min()), "max": float(x.max())}

    summary = {
        "label": args.label,
        "rep": args.rep,
        "prmtop": os.path.abspath(args.prmtop),
        "traj": os.path.abspath(args.traj),
        "receptor_selection": rec_sel,
        "ligand_resname": args.ligand_resname,
        "partner_residues": args.partner_residues or None,
        "partner_kind": args.partner_kind if args.partner_residues else None,
        "apo": bool(args.apo),
        "n_frames_total": n_total,
        "n_frames_analyzed": n,
        "frame_range": [start, end, args.stride],
        "contact_cutoff_A": args.contact_cutoff,
        "hbond": {"ha_cutoff_A": args.hbond_ha_cutoff,
                  "da_cutoff_A": args.hbond_da_cutoff,
                  "angle_deg": args.hbond_angle},
        "metrics": {k: stat(out[k]) for k in (
            "receptor_rmsd", "ligand_rmsd", "ligand_internal_rmsd", "receptor_rg",
            "ligand_rg", "ligand_com_drift", "min_dist", "n_contact_atoms",
            "n_contact_residues", "buried_frac", "loop_min_dist", "loop_contact",
            "hbond_lig_to_rec", "hbond_rec_to_lig", "receptor_sasa")},
    }
    if args.summary:
        with open(args.summary, "w", encoding="utf-8") as fh:
            json.dump(summary, fh, indent=2)
            fh.write("\n")

    print(f"{args.label} r{args.rep}: {n}/{n_total} frames -> {args.out}")
    print(f"  receptor RMSD mean {np.nanmean(receptor_rmsd):.2f} A; "
          f"ligand RMSD mean {np.nanmean(ligand_rmsd) if holo else float('nan'):.2f} A; "
          f"min-dist median {np.nanmedian(min_dist) if holo else float('nan'):.2f} A")


if __name__ == "__main__":
    main()
