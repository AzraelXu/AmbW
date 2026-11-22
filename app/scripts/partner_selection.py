#!/usr/bin/env python3
"""Select a peptide/antibody partner by residue range from an Amber topology.

Used by analyze_replica.py and analyze_ligand_hbonds.py when the binding
partner is a protein chain rather than a small molecule. The prmtop residue
order is authoritative: it is the same order produced by prepare_complex.py
and verified by build_linkage.json.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import qc_system as qc  # noqa: E402


def parse_ranges(text):
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


def partner_atom_indices(prmtop, partner_residues):
    sections = qc.read_prmtop(prmtop)
    pointers = qc.numeric(sections, "POINTERS", int)
    if not pointers:
        raise ValueError(f"{prmtop}: missing POINTERS")
    natom = pointers[0]
    labels, ranges = qc.residue_atom_ranges(sections, natom)
    resids = parse_ranges(partner_residues)
    if not resids:
        raise ValueError("partner residue selection is empty")
    bad = [r for r in resids if r < 1 or r > len(ranges)]
    if bad:
        raise ValueError(f"partner residues out of range 1..{len(ranges)}: {bad}")
    idxs = []
    for r in resids:
        first, last = ranges[r - 1]
        idxs.extend(range(first, last))
    return np.array(sorted(set(idxs)), dtype=int), natom, labels


def select_partner(universe, prmtop, partner_residues):
    """Return (ligand_atomgroup, receptor_atomgroup, description)."""
    idxs, natom, _ = partner_atom_indices(prmtop, partner_residues)
    if universe.atoms.n_atoms != natom:
        raise ValueError(
            f"trajectory/topology atom count mismatch: universe {universe.atoms.n_atoms}, "
            f"prmtop {natom}")
    lig = universe.atoms[idxs]
    rec_idx = np.setdiff1d(np.arange(natom), idxs)
    rec = universe.atoms[rec_idx]
    return lig, rec, f"partner residues {partner_residues}"
