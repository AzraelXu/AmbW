#!/usr/bin/env python3
"""Inspect a protein complex PDB and classify its components.

Pure standard library. Reads one MODEL (default 1), applies an altloc policy,
groups residues into chain segments, and reports:

  - protein / peptide / antibody candidates and chain sizes
  - small-molecule ligand candidates with contact-based scores
  - water, monovalent ions, divalent metals, glycans, crystallization additives
  - disulfides (SSBOND and SG-SG distance)
  - chain breaks (C-N distance and numbering gaps)
  - missing side-chain heavy atoms (best effort)
  - histidines, modified residues, altlocs, models

The JSON output is consumed by bin/amber-auto and prepare_complex.py. When
--json-out is given the script exits 0 even if classification.errors is
non-empty, so callers can inspect and report the errors themselves.

Usage:
  python3 inspect_complex.py --pdb complex.pdb --json-out components.json
  python3 inspect_complex.py --pdb complex.pdb --json
  python3 inspect_complex.py --pdb complex.pdb --type protein-peptide \
      --receptor-chains A --partner-chains B
"""
import argparse
import hashlib
import json
import math
import os
import sys

# ---------------------------------------------------------------------------
# residue dictionaries
# ---------------------------------------------------------------------------
AA_STANDARD = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
AA_AMBER = {"HID", "HIE", "HIP", "HIN", "CYX", "CYM", "ASH", "GLH", "LYN"}
AA_ALIASES = {"HSD": "HID", "HSE": "HIE", "HSP": "HIP"}
CAPS = {"ACE", "NME", "NH2", "NHE"}

NUCLEIC = {
    "A", "C", "G", "U", "DA", "DC", "DG", "DT",
    "RA", "RC", "RG", "RU", "ADE", "CYT", "GUA", "URI",
    "A5", "C5", "G5", "U5", "A3", "C3", "G3", "U3",
}

WATER = {
    "WAT", "HOH", "H2O", "SOL", "TIP", "TIP3", "TP3", "T3P", "SPC", "SPCE",
    "OPC", "OP3", "T4P", "TIP4", "T5P", "FB3", "FB", "DOD",
}

ION_MONO = {
    "NA", "NA+", "K", "K+", "CL", "CL-", "LI", "LI+", "RB", "RB+", "CS",
    "CS+", "F", "F-", "BR", "BR-", "I", "I-", "NH4",
}

METAL_DIVALENT = {
    "MG", "MG2", "CA", "CA2", "ZN", "ZN2", "MN", "MN2", "FE", "FE2", "FE3",
    "CU", "CU2", "CO", "CO2", "NI", "NI2", "CD", "CD2", "HG", "HG2", "SR",
    "SR2", "BA", "BA2",
}

GLYCAN = {
    "NAG", "NDG", "NGA", "NBG", "BMA", "MAN", "BGC", "GLC", "GAL", "GLA",
    "GCS", "FUC", "FUL", "XYS", "XYP", "RAM", "RHA", "SIA", "SLT", "KDO",
    "PSE", "ARA", "XYL", "GUL", "IDR", "SGN",
}

ADDITIVE = {
    "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "P6G", "MPD",
    "DMS", "ACT", "ACY", "CIT", "FMT", "TLA", "MES", "EPE", "HEP", "TRS",
    "SCN", "NO3", "BOG", "LDA", "C8E", "XPE", "IPA",
}

MOD_CONVERTIBLE = {"MSE"}
MOD_HARD = {
    "SEP", "TPO", "PTR", "CSO", "CSD", "KCX", "LLP", "MLY", "M3L", "HYP",
    "PCA", "FME", "CME", "OCS", "CAS", "SEC", "PYL", "ASX", "GLX", "PHD",
    "DAL", "NLE", "ORN", "SAR", "ABA",
}

# expected side-chain heavy atoms for missing-atom reporting
AA_HEAVY = {
    "ALA": {"N", "CA", "C", "O", "CB"},
    "ARG": {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "ASN": {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"},
    "ASP": {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"},
    "CYS": {"N", "CA", "C", "O", "CB", "SG"},
    "GLN": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"},
    "GLU": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"},
    "GLY": {"N", "CA", "C", "O"},
    "HIS": {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "ILE": {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"},
    "LEU": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"},
    "LYS": {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"},
    "MET": {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"},
    "PHE": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PRO": {"N", "CA", "C", "O", "CB", "CG", "CD"},
    "SER": {"N", "CA", "C", "O", "CB", "OG"},
    "THR": {"N", "CA", "C", "O", "CB", "OG1", "CG2"},
    "TRP": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2",
            "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ",
            "OH"},
    "VAL": {"N", "CA", "C", "O", "CB", "CG1", "CG2"},
}
AA_HEAVY["HID"] = AA_HEAVY["HIE"] = AA_HEAVY["HIP"] = AA_HEAVY["HIS"]
AA_HEAVY["HIN"] = AA_HEAVY["HIS"]
AA_HEAVY["CYX"] = AA_HEAVY["CYM"] = AA_HEAVY["CYS"]
AA_HEAVY["ASH"] = AA_HEAVY["ASP"]
AA_HEAVY["GLH"] = AA_HEAVY["GLU"]
AA_HEAVY["LYN"] = AA_HEAVY["LYS"]

TWO_LETTER_ELEMENTS = {
    "HE", "LI", "BE", "NE", "NA", "MG", "AL", "SI", "CL", "AR", "CA", "SC",
    "TI", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE",
    "BR", "KR", "RB", "SR", "ZR", "MO", "RU", "RH", "PD", "AG", "CD", "IN",
    "SN", "SB", "TE", "XE", "CS", "BA", "LA", "CE", "PR", "ND", "SM", "EU",
    "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "RE", "OS",
    "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA",
    "AC", "TH", "PA", "U",
}


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def warn(msg):
    print(f"warning: {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------
def _float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def infer_element(name, resname=""):
    """Best-effort element from a PDB atom name."""
    n = (name or "").strip().upper()
    if not n:
        return ""
    while n and (n[0].isdigit() or n[0] in "+-"):
        n = n[1:]
    if not n:
        return ""
    rn = (resname or "").strip().upper()
    protein_like = rn in AA_STANDARD or rn in AA_AMBER or rn in AA_ALIASES or rn in CAPS
    # Protein atom names are one-letter element + remote atom name (CA, CD, ...)
    if protein_like:
        if n.startswith("SE"):
            return "SE"
        return n[0]
    if n[:2] in TWO_LETTER_ELEMENTS:
        return n[:2]
    return n[0]


def parse_atom(line, line_no=0):
    line = line.ljust(80)
    record = line[0:6].strip()
    name = line[12:16].strip()
    altloc = line[16].strip()
    resname = line[17:20].strip()
    chain = line[21].strip()
    resseq = _int(line[22:26], 0)
    icode = line[26].strip()
    x = _float(line[30:38], float("nan"))
    y = _float(line[38:46], float("nan"))
    z = _float(line[46:54], float("nan"))
    occ = _float(line[54:60], 1.0)
    bfac = _float(line[60:66], 0.0)
    element = line[76:78].strip().upper()
    if not element:
        element = infer_element(name, resname)
    return {
        "record": record,
        "serial": _int(line[6:11], 0),
        "name": name,
        "altloc": altloc,
        "resname": resname,
        "chain": chain,
        "resseq": resseq,
        "icode": icode,
        "x": x,
        "y": y,
        "z": z,
        "occ": occ,
        "bfac": bfac,
        "element": element,
        "line": line.rstrip("\n"),
        "line_no": line_no,
    }


def read_pdb(path):
    """Return (models, header). Each model is a list of atom dicts plus TER markers.

    TER is stored as {"record": "TER"} so segment breaks survive parsing.
    """
    if not os.path.isfile(path):
        die(f"PDB not found: {path}")
    if path.lower().endswith((".cif", ".mmcif")):
        die("mmCIF input is not supported; convert it to PDB first "
            "(e.g. gemmi convert input.cif output.pdb)")
    models = []
    header = []
    current = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh, 1):
            rec = line[0:6].strip()
            if rec == "MODEL":
                if current is not None:
                    models.append(current)
                current = []
                continue
            if rec == "ENDMDL":
                if current is not None:
                    models.append(current)
                    current = None
                continue
            if rec in ("ATOM", "HETATM"):
                atom = parse_atom(line, i)
                if current is None:
                    # no MODEL records: keep a single implicit model; atoms
                    # outside MODEL/ENDMDL after a model block form a new model
                    current = []
                current.append(atom)
                continue
            if rec == "TER":
                if current is not None:
                    current.append({"record": "TER", "line_no": i})
                else:
                    header.append(line.rstrip("\n"))
                continue
            if current is None:
                header.append(line.rstrip("\n"))
    if current is not None:
        models.append(current)
    if not models:
        die(f"no ATOM/HETATM records in {path}")
    return models, header


def choose_altlocs(atoms, policy="most-populous"):
    """Filter altlocs residue by residue. Keeps blank altloc plus the chosen one."""
    # Microheterogeneity (different residue names at the same position) cannot
    # be merged safely into one Amber residue.
    resnames_at_pos = {}
    for atom in atoms:
        if atom.get("record") == "TER":
            continue
        pos = (atom["chain"], atom["resseq"], atom["icode"])
        resnames_at_pos.setdefault(pos, set()).add(atom["resname"])
    conflicts = {pos: sorted(names) for pos, names in resnames_at_pos.items()
                 if len(names) > 1}
    if conflicts:
        pos, names = next(iter(conflicts.items()))
        raise ValueError(
            f"alternate conformers at {pos[0]}:{pos[1]}{pos[2]} have different "
            f"residue names {names}; microheterogeneity is not supported")
    # group atom indices by residue key, preserving first-seen order
    order = []
    groups = {}
    for idx, atom in enumerate(atoms):
        if atom.get("record") == "TER":
            continue
        key = (atom["chain"], atom["resseq"], atom["icode"], atom["resname"])
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(idx)

    keep = []
    altloc_report = []
    for key in order:
        idxs = groups[key]
        alt_counts = {}
        alt_occ = {}
        alt_resnames = {}
        for i in idxs:
            a = atoms[i]
            alt = a["altloc"]
            alt_counts[alt] = alt_counts.get(alt, 0) + 1
            alt_occ[alt] = alt_occ.get(alt, 0.0) + a["occ"]
            alt_resnames.setdefault(alt, set()).add(a["resname"])
        nonblank = [x for x in alt_counts if x]
        chosen = ""
        if nonblank:
            if policy in (None, "", "most-populous"):
                chosen = max(nonblank, key=lambda x: (alt_occ[x], alt_counts[x]))
            elif policy in alt_counts:
                chosen = policy
            else:
                chosen = max(nonblank, key=lambda x: (alt_occ[x], alt_counts[x]))
            # altloc resname conflicts are not safe to merge
            names = set()
            for alt in [""] + nonblank:
                names |= alt_resnames.get(alt, set())
            if len(names) > 1:
                raise ValueError(
                    f"alternate conformers of {key[3]} {key[0]}:{key[1]} have "
                    f"different residue names {sorted(names)}; not supported")
            altloc_report.append({
                "chain": key[0], "resseq": key[1], "icode": key[2],
                "resname": key[3], "altlocs": sorted(alt_counts),
                "chosen": chosen or "blank",
            })
        for i in idxs:
            if atoms[i]["altloc"] in ("", chosen):
                keep.append(i)
    kept = [atoms[i] for i in keep]
    # preserve TER markers in their original relative positions
    # (rebuild by original order)
    out = []
    keep_set = set(keep)
    for i, atom in enumerate(atoms):
        if atom.get("record") == "TER":
            out.append(atom)
        elif i in keep_set:
            out.append(atom)
    return out, altloc_report


def group_residues(atoms):
    """Group consecutive ATOM/HETATM records into residues, preserving TER breaks."""
    residues = []
    order = []
    by_key = {}
    ter_pending = False
    for atom in atoms:
        if atom.get("record") == "TER":
            ter_pending = True
            continue
        key = (atom["chain"], atom["resseq"], atom["icode"], atom["resname"])
        if key not in by_key:
            res = {
                "index": len(residues) + 1,
                "chain": atom["chain"],
                "resseq": atom["resseq"],
                "icode": atom["icode"],
                "resname": atom["resname"],
                "atoms": [],
                "ter_before": ter_pending,
                "record_types": set(),
            }
            by_key[key] = res
            residues.append(res)
            order.append(key)
            ter_pending = False
        else:
            if ter_pending:
                # same key after TER: treat as a new segment by marking the
                # previous residue
                residues[-1]["ter_after"] = True
                ter_pending = False
        by_key[key]["atoms"].append(atom)
        by_key[key]["record_types"].add(atom["record"])
    if ter_pending and residues:
        residues[-1]["ter_after"] = True
    return residues


def classify_resname(resname):
    rn = (resname or "").strip().upper()
    if rn in AA_ALIASES or rn in AA_STANDARD or rn in AA_AMBER or rn in CAPS:
        return "protein"
    if rn in MOD_CONVERTIBLE:
        return "modified_convertible"
    if rn in NUCLEIC:
        return "nucleic"
    if rn in WATER:
        return "water"
    if rn in ION_MONO:
        return "ion_mono"
    if rn in METAL_DIVALENT:
        return "metal_divalent"
    if rn in GLYCAN:
        return "glycan"
    if rn in ADDITIVE:
        return "additive"
    if rn in MOD_HARD:
        return "modified_hard"
    return "unknown"


def atom_coord(atom):
    return (atom["x"], atom["y"], atom["z"])


def distance(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def is_heavy(atom):
    return atom["element"].upper() not in ("H", "D")


def atom_by_name(residue, name):
    for atom in residue["atoms"]:
        if atom["name"].upper() == name.upper():
            return atom
    return None


# ---------------------------------------------------------------------------
# chain segmentation
# ---------------------------------------------------------------------------
def segment_residues(residues, cn_break_cutoff=2.0):
    """Assign segment ids, splitting on TER, chain change, C-N break, and gaps."""
    segments = []
    current = None
    prev = None
    for res in residues:
        new_segment = False
        reason = ""
        if current is None:
            new_segment = True
            reason = "first"
        elif res.get("ter_before") or prev.get("ter_after"):
            new_segment = True
            reason = "TER"
        elif res["chain"] != prev["chain"]:
            new_segment = True
            reason = "chain"
        else:
            c_atom = atom_by_name(prev, "C")
            n_atom = atom_by_name(res, "N")
            if c_atom is not None and n_atom is not None:
                d = distance(atom_coord(c_atom), atom_coord(n_atom))
                if d > cn_break_cutoff:
                    new_segment = True
                    reason = f"C-N {d:.2f} A"
            else:
                gap = res["resseq"] - prev["resseq"]
                if gap > 1:
                    new_segment = True
                    reason = "numbering gap with missing C/N"
        if new_segment:
            current = {
                "id": "",
                "chain": res["chain"],
                "residues": [],
                "start_reason": reason,
            }
            segments.append(current)
        current["residues"].append(res)
        prev = res
    # label segments: use chain when unique, otherwise chain#2, chain#3...
    seen = {}
    for seg in segments:
        chain = seg["chain"]
        base = chain if chain else "S"
        seen[base] = seen.get(base, 0) + 1
        seg["id"] = base if seen[base] == 1 else f"{base}#{seen[base]}"
    return segments


def segment_disulfides(segments, ss_cutoff=2.5):
    """Return (pairs, by_segment_count). Pairs are residue index pairs."""
    sg_res = []
    for seg in segments:
        for res in seg["residues"]:
            if res["resname"].upper() in ("CYS", "CYX"):
                sg = atom_by_name(res, "SG")
                if sg is not None and is_heavy(sg):
                    sg_res.append((res, sg, seg["id"]))
    pairs = []
    by_seg = {seg["id"]: 0 for seg in segments}
    for i in range(len(sg_res)):
        for j in range(i + 1, len(sg_res)):
            r1, a1, s1 = sg_res[i]
            r2, a2, s2 = sg_res[j]
            d = distance(atom_coord(a1), atom_coord(a2))
            if d < ss_cutoff:
                pairs.append({
                    "res1": r1["index"], "chain1": r1["chain"],
                    "resseq1": r1["resseq"], "icode1": r1["icode"],
                    "res2": r2["index"], "chain2": r2["chain"],
                    "resseq2": r2["resseq"], "icode2": r2["icode"],
                    "sg_distance_A": round(d, 3), "source": "distance",
                })
                if s1 == s2:
                    by_seg[s1] = by_seg.get(s1, 0) + 1
    return pairs, by_seg


def parse_ssbond(header, residues):
    """Map SSBOND records to residue index pairs."""
    index = {}
    for res in residues:
        index[(res["chain"], res["resseq"], res["icode"])] = res
    pairs = []
    for line in header:
        if not line.startswith("SSBOND"):
            continue
        line = line.ljust(80)
        try:
            chain1 = line[15].strip()
            seq1 = _int(line[17:21], 0)
            icode1 = line[21].strip()
            chain2 = line[29].strip()
            seq2 = _int(line[31:35], 0)
            icode2 = line[35].strip()
        except Exception:
            continue
        r1 = index.get((chain1, seq1, icode1))
        r2 = index.get((chain2, seq2, icode2))
        if r1 is not None and r2 is not None:
            pairs.append({
                "res1": r1["index"], "chain1": chain1, "resseq1": seq1,
                "icode1": icode1, "res2": r2["index"], "chain2": chain2,
                "resseq2": seq2, "icode2": icode2,
                "sg_distance_A": None, "source": "SSBOND",
            })
    return pairs


def merge_disulfides(a, b):
    """Merge SSBOND and distance pairs without duplicates."""
    out = []
    seen = set()
    for pair in a + b:
        key = tuple(sorted((pair["res1"], pair["res2"])))
        if key in seen:
            continue
        seen.add(key)
        out.append(pair)
    return out


def chain_breaks(segments, cn_break_cutoff=2.0):
    breaks = []

    def add_break(prev, cur, segment_id):
        c_atom = atom_by_name(prev, "C")
        n_atom = atom_by_name(cur, "N")
        d = None
        if c_atom is not None and n_atom is not None:
            d = distance(atom_coord(c_atom), atom_coord(n_atom))
        gap = cur["resseq"] - prev["resseq"]
        if (d is not None and d > cn_break_cutoff) or (d is None and gap > 1):
            breaks.append({
                "chain": prev["chain"],
                "segment": segment_id,
                "after_resseq": prev["resseq"],
                "after_icode": prev["icode"],
                "before_resseq": cur["resseq"],
                "before_icode": cur["icode"],
                "cn_distance_A": round(d, 3) if d is not None else None,
                "numbering_gap": gap if gap > 1 else None,
            })

    for seg in segments:
        residues = seg["residues"]
        for prev, cur in zip(residues, residues[1:]):
            add_break(prev, cur, seg["id"])
    # A break that caused segment_residues() to split the same chain into two
    # segments is not visible inside either segment.
    for prev_seg, cur_seg in zip(segments, segments[1:]):
        if prev_seg["chain"] != cur_seg["chain"]:
            continue
        reason = cur_seg.get("start_reason", "")
        if reason.startswith("C-N") or reason.startswith("numbering"):
            add_break(prev_seg["residues"][-1], cur_seg["residues"][0],
                      cur_seg["id"])
    return breaks


def missing_atoms(residues):
    missing = []
    for res in residues:
        rn = AA_ALIASES.get(res["resname"].upper(), res["resname"].upper())
        if res["resname"].upper() == "MSE":
            rn = "MET"
        expected = AA_HEAVY.get(rn)
        if not expected:
            continue
        actual = {a["name"].upper() for a in res["atoms"] if is_heavy(a)}
        if res["resname"].upper() == "MSE":
            # SE is the selenium surrogate for SD; prepare_complex renames it
            if "SE" in actual:
                actual.discard("SE")
                actual.add("SD")
        # terminal OXT is optional; never report it
        lack = sorted(expected - actual)
        if lack:
            missing.append({
                "chain": res["chain"], "resseq": res["resseq"],
                "icode": res["icode"], "resname": res["resname"],
                "missing": lack,
            })
    return missing


# ---------------------------------------------------------------------------
# ligand detection
# ---------------------------------------------------------------------------
def residue_heavy_atoms(residue):
    return [a for a in residue["atoms"] if is_heavy(a)]


def candidate_score(residue, protein_atoms, contact_cutoff, covalent_cutoff):
    lig_atoms = residue_heavy_atoms(residue)
    heavy = len(lig_atoms)
    if heavy == 0:
        return None
    contact_pairs = 0
    contact_res = set()
    buried = 0
    min_dist = float("inf")
    for la in lig_atoms:
        lp = atom_coord(la)
        hit = False
        for pa, pres in protein_atoms:
            d = distance(lp, atom_coord(pa))
            if d < min_dist:
                min_dist = d
            if d < contact_cutoff:
                contact_pairs += 1
                contact_res.add(pres["index"])
                hit = True
        if hit:
            buried += 1
    buried_frac = buried / heavy if heavy else 0.0
    score = (10.0 * min(heavy, 40) + 3.0 * min(contact_pairs, 60)
             + 50.0 * buried_frac + 5.0 * min(len(contact_res), 20))
    return {
        "heavy_atoms": heavy,
        "hydrogens": sum(1 for a in residue["atoms"] if a["element"] in ("H", "D")),
        "elements": _element_counts(residue["atoms"]),
        "contact_pairs": contact_pairs,
        "contact_residues": len(contact_res),
        "buried_fraction": round(buried_frac, 3),
        "min_distance_A": round(min_dist, 3) if min_dist != float("inf") else None,
        "covalent": min_dist < covalent_cutoff,
        "score": round(score, 2),
        "additive": residue["resname"].upper() in ADDITIVE,
    }


def _element_counts(atoms):
    counts = {}
    for atom in atoms:
        el = atom["element"].upper() or "?"
        counts[el] = counts.get(el, 0) + 1
    return dict(sorted(counts.items()))


def parse_ligand_spec(spec):
    parts = spec.split(":")
    if len(parts) > 3:
        raise ValueError(f"bad --ligand spec: {spec} (expected RESNAME[:CHAIN[:RESID]])")
    resname = parts[0].strip()
    chain = parts[1].strip() if len(parts) > 1 and parts[1].strip() else None
    resid = None
    if len(parts) > 2 and parts[2].strip():
        resid = _int(parts[2].strip(), None)
    return resname, chain, resid


# ---------------------------------------------------------------------------
# classification
# ---------------------------------------------------------------------------
def _contact_between_segments(seg_a, seg_b, cutoff=5.0):
    """True if two segments contact.

    CA-CA distances are used first for speed. The CA threshold is expanded to
    8 A because a short peptide can contact a receptor through side chains
    while its CA atoms are farther than 5 A apart. If that still fails and the
    segment pair is small enough, an all-heavy-atom check is done.
    """
    def cas(seg):
        out = []
        for res in seg["residues"]:
            ca = atom_by_name(res, "CA")
            if ca is not None:
                out.append(ca)
        return out
    a_ca, b_ca = cas(seg_a), cas(seg_b)
    if a_ca and b_ca:
        for x in a_ca:
            xp = atom_coord(x)
            for y in b_ca:
                if distance(xp, atom_coord(y)) < max(cutoff, 8.0):
                    return True
        return False
    a = [x for r in seg_a["residues"] for x in residue_heavy_atoms(r)]
    b = [x for r in seg_b["residues"] for x in residue_heavy_atoms(r)]
    if len(a) * len(b) > 4000000:
        return False
    for x in a:
        xp = atom_coord(x)
        for y in b:
            if distance(xp, atom_coord(y)) < cutoff:
                return True
    return False


def classify(segments, residues, ligands, opts, disulfides, by_seg):
    errors = []
    warnings = []
    evidence = []
    type_override = opts.get("type") or "auto"
    receptor_chains = _split_list(opts.get("receptor_chains"))
    partner_chains = _split_list(opts.get("partner_chains"))

    def err(code, message):
        errors.append({"code": code, "message": message})

    if opts.get("nucleic"):
        err("nucleic_unsupported",
            "protein/nucleic-acid complexes are not supported by amber-auto; "
            "use the amber-kit RNA workflow")
    if opts.get("glycans") and not opts.get("drop_glycans"):
        names = ", ".join(f"{x['resname']} {x['chain']}:{x['resseq']}"
                          for x in opts["glycans"][:20])
        err("glycan_unsupported",
            f"N-glycans are not supported: {names}. Remove them or rerun with "
            "--drop-glycans (explicit, logged).")
    if opts.get("metals") and not opts.get("drop_metals"):
        names = ", ".join(f"{x['resname']} {x['chain']}:{x['resseq']}"
                          for x in opts["metals"][:20])
        err("metal_unsupported",
            f"divalent metals are not supported: {names}. Remove them or rerun "
            "with --drop-metals (explicit, logged).")
    if opts.get("modified_hard") and not opts.get("drop_modified"):
        names = ", ".join(f"{x['resname']} {x['chain']}:{x['resseq']}"
                          for x in opts["modified_hard"][:20])
        err("modified_unsupported",
            f"modified residues are not supported: {names}. Only MSE is "
            "converted automatically; rerun with --drop-modified to remove them.")
    if opts.get("unknown_hetatm") and not opts.get("drop_unknown_hetatm"):
        names = ", ".join(f"{x['resname']} {x['chain']}:{x['resseq']}"
                          for x in opts["unknown_hetatm"][:20])
        err("unknown_hetatm",
            f"unknown HETATM residues (not recognized as ligand/water/ion/"
            f"additive): {names}. Specify --ligand or --drop-unknown-hetatm.")

    if not segments:
        err("no_protein", "no protein residues found")
        return {
            "type": "unknown", "confidence": "none",
            "receptor": {"chains": [], "segments": []},
            "partner": {"chains": [], "segments": []},
            "ligand": None, "evidence": evidence,
            "warnings": warnings, "errors": errors,
        }

    # ligand selection
    selected = None
    if ligands:
        spec = opts.get("ligand")
        if spec:
            try:
                rn, ch, rid = parse_ligand_spec(spec)
            except ValueError as exc:
                err("bad_ligand_spec", str(exc))
                rn, ch, rid = "", None, None
            matches = [c for c in ligands
                       if c["resname"].upper() == rn.upper()
                       and (ch is None or c["chain"] == ch)
                       and (rid is None or c["resseq"] == rid)]
            if len(matches) != 1:
                err("ligand_selection_failed",
                    f"--ligand {spec} matched {len(matches)} candidates; "
                    f"available: {_ligand_list(ligands)}")
            else:
                selected = matches[0]
        else:
            ranked = sorted(ligands, key=lambda c: c["score"], reverse=True)
            selected = ranked[0]
            if len(ranked) > 1:
                second = ranked[1]
                if second["score"] >= 0.75 * ranked[0]["score"]:
                    err("ligand_ambiguous",
                        "multiple plausible ligand candidates: "
                        f"{_ligand_list(ranked)}. Choose one with "
                        "--ligand RESNAME[:CHAIN[:RESID]] and add "
                        "--drop-other-ligands.")
                else:
                    warnings.append(
                        "selected ligand %s %s:%s (score %.1f) over %d lower-"
                        "scoring candidate(s)" % (
                            selected["resname"], selected["chain"],
                            selected["resseq"], selected["score"],
                            len(ranked) - 1))
        if selected is not None:
            if selected.get("covalent"):
                err("covalent_ligand",
                    f"ligand {selected['resname']} is covalently bonded to the "
                    f"protein (min distance {selected['min_distance_A']} A); "
                    "covalent ligands are not supported")
            others = [c for c in ligands if c is not selected and not c["additive"]]
            if others and not opts.get("drop_other_ligands"):
                err("other_ligands",
                    f"other non-additive ligand candidates are present: "
                    f"{_ligand_list(others)}. Rerun with --drop-other-ligands "
                    "to keep only the selected ligand.")
            same_res = [c for c in ligands
                        if c["resname"].upper() == selected["resname"].upper()
                        and c is not selected]
            if same_res and not opts.get("drop_other_ligands"):
                err("multiple_ligand_copies",
                    f"ligand {selected['resname']} has multiple copies: "
                    f"{_ligand_list([selected] + same_res)}. Select one with "
                    "--ligand RESNAME:CHAIN:RESID and add --drop-other-ligands.")
            evidence.append(
                f"selected ligand {selected['resname']} {selected['chain']}:"
                f"{selected['resseq']} ({selected['heavy_atoms']} heavy atoms, "
                f"{selected['contact_residues']} contact residues, score "
                f"{selected['score']})")

    def assignment(rec_chains, part_chains, kind, confidence):
        rec_segs = [s for s in segments if _seg_matches(s, rec_chains)]
        part_segs = [s for s in segments if _seg_matches(s, part_chains)]
        if not rec_segs:
            err("receptor_missing", f"receptor chains {rec_chains} matched no segment")
        if part_chains and not part_segs:
            err("partner_missing", f"partner chains {part_chains} matched no segment")
        if rec_segs and part_segs and not _segments_contact(rec_segs, part_segs):
            err("no_interface",
                "receptor and partner do not contact each other (<5 A); "
                "check --receptor-chains/--partner-chains")
        return {
            "type": kind, "confidence": confidence,
            "receptor": {
                "chains": sorted({s["chain"] for s in rec_segs}),
                "segments": [s["id"] for s in rec_segs],
            },
            "partner": {
                "chains": sorted({s["chain"] for s in part_segs}),
                "segments": [s["id"] for s in part_segs],
            },
            "ligand": _ligand_payload(selected) if selected else None,
            "evidence": evidence, "warnings": warnings, "errors": errors,
        }

    # explicit chain roles win
    if receptor_chains or partner_chains:
        if not receptor_chains:
            receptor_chains = [s["id"] for s in segments
                               if not _seg_matches(s, partner_chains)]
        if not partner_chains and not selected:
            partner_chains = [s["id"] for s in segments
                              if not _seg_matches(s, receptor_chains)]
        kind = type_override if type_override != "auto" else (
            "protein-ligand" if selected else "protein-protein")
        if selected:
            rec = [s["id"] for s in segments if _seg_matches(s, receptor_chains)]
            if not rec:
                err("receptor_missing", "no receptor segment after chain assignment")
            return assignment(receptor_chains, [], kind, "user")
        return assignment(receptor_chains, partner_chains, kind, "user")

    if type_override == "protein-ligand" and not selected:
        err("type_mismatch",
            "--type protein-ligand was requested but no small-molecule ligand "
            "candidate was found; specify --ligand or remove --type")
    if type_override == "apo-protein" and selected:
        err("type_mismatch",
            "--type apo-protein was requested but ligand candidates were found; "
            "use --drop-other-ligands or remove --type")

    if selected:
        return assignment([s["id"] for s in segments], [], "protein-ligand", "high")

    # Multiple fragments of the same original chain are chain breaks, not a
    # receptor/partner pair. Treat them as one apo protein unless the user
    # explicitly assigned chain roles above.
    chain_ids = {s["chain"] for s in segments}
    if len(chain_ids) == 1 and len(segments) > 1:
        evidence.append(
            f"{len(segments)} fragments of the same chain "
            f"({', '.join(s['id'] for s in segments)}) come from chain breaks")
        warnings.append("chain breaks split one protein chain into multiple "
                        "fragments; they will be capped with TER")
        return assignment([s["id"] for s in segments], [],
                          type_override if type_override != "auto" else "apo-protein",
                          "medium")

    # protein-only classification
    n = len(segments)
    ig = []
    for seg in segments:
        # Fab chains are ~220 aa; full IgG heavy chains can be ~450 aa.
        if 90 <= len(seg["residues"]) <= 500 and by_seg.get(seg["id"], 0) >= 1:
            ig.append(seg)
    if ig:
        non_ig = [s for s in segments if s not in ig]
        if non_ig and any(_contact_between_segments(a, b) for a in ig for b in non_ig):
            if type_override not in ("auto", "protein-antibody"):
                return assignment([s["id"] for s in non_ig], [s["id"] for s in ig],
                                  type_override, "user")
            evidence.append(
                "Ig-like chains with intrachain disulfides: "
                + ", ".join(f"{s['id']}({len(s['residues'])} aa)" for s in ig))
            return assignment([s["id"] for s in non_ig], [s["id"] for s in ig],
                              "protein-antibody",
                              "high" if len(ig) >= 2 else "medium")

    short = [s for s in segments if len(s["residues"]) <= opts.get("peptide_max", 50)]
    larger = [s for s in segments if s not in short]
    if short and larger and any(_contact_between_segments(a, b)
                                for a in short for b in larger):
        if type_override not in ("auto", "protein-peptide", "protein-protein"):
            return assignment([s["id"] for s in larger], [s["id"] for s in short],
                              type_override, "user")
        evidence.append("short protein chain(s): "
                        + ", ".join(f"{s['id']}({len(s['residues'])} aa)" for s in short))
        return assignment([s["id"] for s in larger], [s["id"] for s in short],
                          "protein-peptide", "high")

    if n == 1:
        if type_override not in ("auto", "apo-protein"):
            return assignment([segments[0]["id"]], [], type_override, "user")
        return assignment([segments[0]["id"]], [], "apo-protein", "high")

    if n == 2:
        lengths = [len(s["residues"]) for s in segments]
        if (type_override == "auto"
                and max(lengths) <= 1.25 * min(lengths)):
            warnings.append(
                "two protein chains of similar length with no clear partner; "
                "treating them as an apo oligomer (use --partner-chains to "
                "run interface analysis/MM-GBSA)")
            return assignment([s["id"] for s in segments], [],
                              "apo-protein", "medium")
        rec = max(segments, key=lambda s: len(s["residues"]))
        part = min(segments, key=lambda s: len(s["residues"]))
        kind = type_override if type_override != "auto" else "protein-protein"
        return assignment([rec["id"]], [part["id"]], kind, "medium")

    lengths = [len(s["residues"]) for s in segments]
    if max(lengths) <= 1.25 * min(lengths):
        warnings.append("multi-chain system with similar chain lengths and no "
                        "clear partner; treating as an apo oligomer")
        return assignment([s["id"] for s in segments], [],
                          type_override if type_override != "auto" else "apo-protein",
                          "medium")
    err("ambiguous_chains",
        "ambiguous multi-chain system: "
        + ", ".join(f"{s['id']}({len(s['residues'])} aa)" for s in segments)
        + ". Specify --receptor-chains and --partner-chains.")
    return {
        "type": "unknown", "confidence": "none",
        "receptor": {"chains": [], "segments": []},
        "partner": {"chains": [], "segments": []},
        "ligand": _ligand_payload(selected) if selected else None,
        "evidence": evidence, "warnings": warnings, "errors": errors,
    }


def _split_list(value):
    if not value:
        return []
    if isinstance(value, list):
        return [str(x).strip() for x in value if str(x).strip()]
    return [x.strip() for x in str(value).replace(";", ",").split(",") if x.strip()]


def _seg_matches(seg, chains):
    if not chains:
        return False
    return seg["id"] in chains or seg["chain"] in chains


def _segments_contact(a_segs, b_segs, cutoff=5.0):
    for a in a_segs:
        for b in b_segs:
            if _contact_between_segments(a, b, cutoff):
                return True
    return False


def _ligand_payload(cand):
    if cand is None:
        return None
    return {k: cand[k] for k in (
        "resname", "chain", "resseq", "icode", "heavy_atoms", "hydrogens",
        "elements", "contact_pairs", "contact_residues", "buried_fraction",
        "min_distance_A", "covalent", "score", "additive")}


def _ligand_list(cands):
    return "; ".join(
        f"{c['resname']} {c['chain']}:{c['resseq']} (heavy {c['heavy_atoms']}, "
        f"contacts {c['contact_residues']}, score {c['score']})"
        for c in cands)


# ---------------------------------------------------------------------------
# main inspection
# ---------------------------------------------------------------------------
def inspect(pdb, opts):
    models, header = read_pdb(pdb)
    model = int(opts.get("model") or 1)
    if model < 1 or model > len(models):
        raise ValueError(f"--model {model} out of range 1..{len(models)}")
    atoms, altloc_report = choose_altlocs(models[model - 1],
                                          opts.get("altloc") or "most-populous")
    residues = group_residues(atoms)

    # classify residues
    protein_residues = []
    nucleic = []
    water = []
    ions = []
    metals = []
    glycans = []
    additives = []
    modified_convertible = []
    modified_hard = []
    unknown_hetatm = []
    for res in residues:
        kind = classify_resname(res["resname"])
        if kind == "protein":
            protein_residues.append(res)
        elif kind == "modified_convertible":
            modified_convertible.append(res)
            protein_residues.append(res)
        elif kind == "nucleic":
            nucleic.append(res)
        elif kind == "water":
            water.append(res)
        elif kind == "ion_mono":
            ions.append(res)
        elif kind == "metal_divalent":
            metals.append(res)
        elif kind == "glycan":
            glycans.append(res)
        elif kind == "additive":
            additives.append(res)
        elif kind == "modified_hard":
            modified_hard.append(res)
        else:
            unknown_hetatm.append(res)

    # segments only over protein residues; non-protein residues do not break
    # chain continuity but do occupy sequence numbers (tleap builds them
    # separately, so they do not affect protein residue indices after cleaning)
    segments = segment_residues(protein_residues,
                                opts.get("cn_break_cutoff", 2.0))
    dis_pairs, by_seg = segment_disulfides(segments, opts.get("ss_cutoff", 2.5))
    dis_pairs = merge_disulfides(parse_ssbond(header, protein_residues), dis_pairs)
    breaks = chain_breaks(segments, opts.get("cn_break_cutoff", 2.0))
    missing = missing_atoms(protein_residues)

    # ligand candidates: unknown HETATM and explicitly requested additives.
    # A residue with a standard AA/nucleic name is never a ligand.
    protein_atoms = []
    for res in protein_residues:
        for atom in res["atoms"]:
            if is_heavy(atom):
                protein_atoms.append((atom, res))
    ligand_candidates = []
    for res in unknown_hetatm:
        heavy = residue_heavy_atoms(res)
        if len(heavy) < opts.get("min_ligand_heavy", 3):
            continue
        if len(heavy) > opts.get("max_ligand_heavy", 200):
            continue
        info = candidate_score(res, protein_atoms,
                               opts.get("contact_cutoff", 4.0),
                               opts.get("covalent_cutoff", 1.9))
        if info is None or info["contact_pairs"] < 1:
            continue
        cand = {
            "resname": res["resname"], "chain": res["chain"],
            "resseq": res["resseq"], "icode": res["icode"],
            "residue_index": res["index"],
        }
        cand.update(info)
        ligand_candidates.append(cand)

    # if the user explicitly names an additive, add it as a candidate
    spec = opts.get("ligand")
    if spec:
        try:
            rn, ch, rid = parse_ligand_spec(spec)
        except ValueError:
            rn = None
        if rn:
            for res in additives:
                if res["resname"].upper() != rn.upper():
                    continue
                if ch is not None and res["chain"] != ch:
                    continue
                if rid is not None and res["resseq"] != rid:
                    continue
                info = candidate_score(res, protein_atoms,
                                       opts.get("contact_cutoff", 4.0),
                                       opts.get("covalent_cutoff", 1.9))
                if info is None:
                    continue
                cand = {
                    "resname": res["resname"], "chain": res["chain"],
                    "resseq": res["resseq"], "icode": res["icode"],
                    "residue_index": res["index"],
                }
                cand.update(info)
                ligand_candidates.append(cand)

    opts2 = dict(opts)
    opts2["nucleic"] = nucleic
    opts2["glycans"] = glycans
    opts2["metals"] = metals
    opts2["modified_hard"] = modified_hard
    selected_indices = {c["residue_index"] for c in ligand_candidates}
    opts2["unknown_hetatm"] = [r for r in unknown_hetatm
                               if r["index"] not in selected_indices]
    classification = classify(segments, residues, ligand_candidates, opts2,
                              dis_pairs, by_seg)

    # warnings for non-fatal conditions
    if len(models) > 1:
        classification["warnings"].append(
            f"input has {len(models)} models; using model {model}")
    if altloc_report:
        classification["warnings"].append(
            f"{len(altloc_report)} residue(s) have alternate locations; "
            "the most populated conformer was kept")
    if modified_convertible:
        classification["warnings"].append(
            "MSE will be converted to MET: "
            + ", ".join(f"{r['chain']}:{r['resseq']}" for r in modified_convertible))
    if breaks:
        classification["warnings"].append(
            f"{len(breaks)} chain break(s) detected; TER records will be inserted "
            "and each fragment will get charged termini")
        if opts.get("strict_gaps"):
            classification["errors"].append({
                "code": "chain_break",
                "message": f"{len(breaks)} chain break(s) detected and "
                           "--strict-gaps is set; fix the model or remove the flag",
            })
    if missing:
        classification["warnings"].append(
            f"{len(missing)} residue(s) have missing side-chain heavy atoms; "
            "tleap will add them")
    if additives:
        classification["warnings"].append(
            "crystallization additives will be dropped: "
            + ", ".join(sorted({r["resname"] for r in additives})))
    if water:
        classification["warnings"].append(f"{len(water)} water residue(s) will be dropped")
    if ions:
        classification["warnings"].append(
            f"{len(ions)} monovalent ion residue(s) will be dropped")

    out = {
        "input": os.path.abspath(pdb),
        "sha256": _sha256(pdb),
        "format": "pdb",
        "models": len(models),
        "model_used": model,
        "chain_count": len({s["chain"] for s in segments}),
        "chains": [
            {
                "id": s["id"], "chain": s["chain"],
                "type": "protein",
                "residues": len(s["residues"]),
                "first": _res_identity(s["residues"][0]),
                "last": _res_identity(s["residues"][-1]),
                "intrachain_disulfides": by_seg.get(s["id"], 0),
            } for s in segments
        ],
        "residues": {
            "standard_protein": {"count": len([r for r in protein_residues
                                               if r not in modified_convertible])},
            "modified_convertible": [
                {"resname": r["resname"], "chain": r["chain"],
                 "resseq": r["resseq"], "icode": r["icode"]}
                for r in modified_convertible],
            "modified_hard": [
                {"resname": r["resname"], "chain": r["chain"],
                 "resseq": r["resseq"], "icode": r["icode"]}
                for r in modified_hard],
            "histidines": [
                {"chain": r["chain"], "resseq": r["resseq"],
                 "resname": r["resname"]}
                for r in protein_residues if r["resname"].upper() in
                ("HIS", "HID", "HIE", "HIP")],
            "cysteines": [
                {"chain": r["chain"], "resseq": r["resseq"],
                 "resname": r["resname"]}
                for r in protein_residues if r["resname"].upper() in ("CYS", "CYX")],
        },
        "water": {"count": len(water),
                  "resnames": sorted({r["resname"] for r in water})},
        "ions": {
            "monovalent": _count_resnames(ions),
            "divalent": _count_resnames(metals),
        },
        "glycans": {
            "detected": bool(glycans),
            "residues": [_res_identity(r) for r in glycans],
        },
        "metals": {"divalent": [_res_identity(r) for r in metals]},
        "ligands": ligand_candidates,
        "disulfides": dis_pairs,
        "chain_breaks": breaks,
        "missing_atoms": missing,
        "altlocs": altloc_report,
        "classification": classification,
    }
    return out


def _res_identity(res):
    return {"resname": res["resname"], "chain": res["chain"],
            "resseq": res["resseq"], "icode": res["icode"]}


def _count_resnames(residues):
    counts = {}
    for res in residues:
        rn = res["resname"]
        counts[rn] = counts.get(rn, 0) + 1
    return counts


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def print_summary(data):
    cls = data["classification"]
    print(f"input: {data['input']}")
    print(f"models: {data['models']} (using {data['model_used']})")
    print("chains:")
    for chain in data["chains"]:
        print(f"  {chain['id']:>4}  {chain['residues']:5d} aa  "
              f"SS={chain['intrachain_disulfides']}  "
              f"{chain['first']['resname']}{chain['first']['resseq']}.."
              f"{chain['last']['resname']}{chain['last']['resseq']}")
    if data["ligands"]:
        print("ligand candidates:")
        for lig in data["ligands"]:
            print(f"  {lig['resname']} {lig['chain']}:{lig['resseq']}  "
                  f"heavy={lig['heavy_atoms']} H={lig['hydrogens']} "
                  f"contacts={lig['contact_residues']} "
                  f"min={lig['min_distance_A']} score={lig['score']}"
                  + (" ADDITIVE" if lig["additive"] else ""))
    print(f"type: {cls['type']} (confidence {cls['confidence']})")
    print(f"receptor: chains={cls['receptor']['chains']} segments={cls['receptor']['segments']}")
    print(f"partner:  chains={cls['partner']['chains']} segments={cls['partner']['segments']}")
    if cls["ligand"]:
        lig = cls["ligand"]
        print(f"ligand:   {lig['resname']} {lig['chain']}:{lig['resseq']}")
    if data["disulfides"]:
        print("disulfides: " + ", ".join(
            f"{d['chain1']}:{d['resseq1']}-{d['chain2']}:{d['resseq2']}"
            + (f" ({d['sg_distance_A']} A)" if d.get("sg_distance_A") else "")
            for d in data["disulfides"]))
    if data["chain_breaks"]:
        print(f"chain breaks: {len(data['chain_breaks'])}")
    for w in cls["warnings"]:
        print(f"warning: {w}")
    for e in cls["errors"]:
        print(f"error: [{e['code']}] {e['message']}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--model", type=int, default=1)
    ap.add_argument("--altloc", default="most-populous")
    ap.add_argument("--type", default="auto",
                    choices=["auto", "protein-ligand", "protein-peptide",
                             "protein-antibody", "protein-protein", "apo-protein"])
    ap.add_argument("--ligand", default=None,
                    help="RESNAME[:CHAIN[:RESID]]")
    ap.add_argument("--receptor-chains", default=None)
    ap.add_argument("--partner-chains", default=None)
    ap.add_argument("--drop-glycans", action="store_true")
    ap.add_argument("--drop-metals", action="store_true")
    ap.add_argument("--drop-modified", action="store_true")
    ap.add_argument("--drop-unknown-hetatm", action="store_true")
    ap.add_argument("--drop-other-ligands", action="store_true")
    ap.add_argument("--strict-gaps", action="store_true")
    ap.add_argument("--peptide-max", type=int, default=50)
    ap.add_argument("--min-ligand-heavy", type=int, default=3)
    ap.add_argument("--max-ligand-heavy", type=int, default=200)
    ap.add_argument("--contact-cutoff", type=float, default=4.0)
    ap.add_argument("--covalent-cutoff", type=float, default=1.9)
    ap.add_argument("--ss-cutoff", type=float, default=2.5)
    ap.add_argument("--cn-break-cutoff", type=float, default=2.0)
    ap.add_argument("--json-out", default=None)
    ap.add_argument("--json", action="store_true", help="print JSON to stdout")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    opts = {
        "model": args.model,
        "altloc": args.altloc,
        "type": args.type,
        "ligand": args.ligand,
        "receptor_chains": args.receptor_chains,
        "partner_chains": args.partner_chains,
        "drop_glycans": args.drop_glycans,
        "drop_metals": args.drop_metals,
        "drop_modified": args.drop_modified,
        "drop_unknown_hetatm": args.drop_unknown_hetatm,
        "drop_other_ligands": args.drop_other_ligands,
        "strict_gaps": args.strict_gaps,
        "peptide_max": args.peptide_max,
        "min_ligand_heavy": args.min_ligand_heavy,
        "max_ligand_heavy": args.max_ligand_heavy,
        "contact_cutoff": args.contact_cutoff,
        "covalent_cutoff": args.covalent_cutoff,
        "ss_cutoff": args.ss_cutoff,
        "cn_break_cutoff": args.cn_break_cutoff,
        "drop_other_ligands": False,
    }
    try:
        data = inspect(args.pdb, opts)
    except ValueError as exc:
        die(str(exc))
    text = json.dumps(data, indent=2) + "\n"
    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as fh:
            fh.write(text)
    if args.json:
        print(text, end="")
    elif not args.quiet:
        print_summary(data)
        if args.json_out:
            print(f"wrote {args.json_out}")
    if not args.quiet:
        for e in data["classification"]["errors"]:
            print(f"error: [{e['code']}] {e['message']}", file=sys.stderr)
    # Machine-readable mode is for callers that will decide what to do.
    if data["classification"]["errors"] and not args.json_out and not args.json:
        sys.exit(2)


if __name__ == "__main__":
    main()
