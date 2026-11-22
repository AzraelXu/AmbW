#!/usr/bin/env python3
"""Prepare a complex PDB for Amber system building.

Reads the components.json produced by inspect_complex.py and writes:

  protein_raw_for_reduce.pdb   protein-only heavy atoms with TER breaks
  protein_reduced.pdb          pdb4amber --reduce output (with H)
  protein_prepared.pdb         heavy atoms, HID/HIE/HIP + flips preserved
  receptor_prepared.pdb        protein-only receptor (protein-ligand)
  complex_prepared.pdb         all protein chains (peptide/antibody/apo)
  ligand_input.pdb             selected ligand, original protonation
  ligand_prepared.pdb          ligand with H, resname LIG
  ligand_charge.json           OpenBabel/RDKit charge inference record
  extra_leap.in                disulfide bond commands
  residue_map.tsv              prepared index -> original identity/role
  prepare_manifest.json        complete preparation record

Important implementation notes:
  * pdb4amber --reduce must NOT be combined with --nohyd: in Amber22 the
    Reduce step runs before assign_histidine(), and --nohyd strips hydrogens
    before histidine assignment, so HID/HIE/HIP would be lost.
  * pdb4amber does not reliably convert CYS to CYX from SSBOND records; this
    script detects SG-SG pairs itself and emits explicit tleap bond commands.
  * All downstream builds use build_system.sh --no-pdb4amber so the prepared
    residue/atom order is authoritative.

Usage:
  python3 prepare_complex.py --pdb complex.pdb --components components.json \
      --outdir prepare [--charge -1|--ligand-smiles ...] [--force]
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import inspect_complex as ic  # noqa: E402


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def warn(msg):
    print(f"warning: {msg}", file=sys.stderr)


def note(msg):
    print(f"[prepare] {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
# PDB output
# ---------------------------------------------------------------------------
def atom_name_field(name, element):
    name = (name or "").strip()
    if not name:
        name = element or "X"
    if len(name) > 4:
        name = name[:4]
    if len(name) == 1 and len(element or "") == 1:
        return (" " + name).ljust(4)
    if name[0].isdigit():
        return name.rjust(4)
    return name.ljust(4)


def format_atom(atom, serial):
    element = (atom.get("element") or "").upper()
    name = atom_name_field(atom.get("name", ""), element)
    resname = (atom.get("resname") or "UNK")[:3]
    chain = (atom.get("chain") or " ")[:1]
    resseq = int(atom.get("resseq") or 1)
    if resseq < 0:
        resseq = 0
    icode = (atom.get("icode") or " ")[:1] or " "
    return (
        f"{atom.get('record', 'ATOM'):<6s}{serial:5d} {name}"
        f"{'':1s}{resname:>3s} {chain}{resseq:4d}{icode}   "
        f"{atom.get('x', 0.0):8.3f}{atom.get('y', 0.0):8.3f}"
        f"{atom.get('z', 0.0):8.3f}{atom.get('occ', 1.0):6.2f}"
        f"{atom.get('bfac', 0.0):6.2f}          {element:>2s}"
    )


def write_residues_pdb(path, residues, include_h=False):
    """Write residues in order, inserting TER at ter_before boundaries."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    lines = []
    serial = 0
    for idx, res in enumerate(residues):
        if idx > 0 and res.get("ter_before"):
            lines.append("TER")
        for atom in res["atoms"]:
            if not include_h and atom.get("element", "").upper() in ("H", "D"):
                continue
            serial += 1
            atom = dict(atom)
            atom["resname"] = res.get("resname", atom.get("resname"))
            atom["chain"] = res.get("chain", atom.get("chain"))
            atom["resseq"] = res.get("resseq", atom.get("resseq"))
            atom["icode"] = res.get("icode", atom.get("icode"))
            lines.append(format_atom(atom, serial))
    lines.append("TER")
    lines.append("END")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def write_atoms_pdb(path, atoms, resname="LIG", chain="L", resid=1):
    """Write one ligand residue, renaming every atom's residue information."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    lines = []
    used = {}
    for i, atom in enumerate(atoms, 1):
        a = dict(atom)
        base = (a.get("name") or a.get("element") or "X").strip()[:4]
        name = base
        k = 1
        while name in used:
            name = f"{base[:3]}{k % 10}"
            k += 1
        used[name] = True
        a["name"] = name
        a["resname"] = resname
        a["chain"] = chain
        a["resseq"] = resid
        a["icode"] = " "
        lines.append(format_atom(a, i))
    lines.append("TER")
    lines.append("END")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# command helpers
# ---------------------------------------------------------------------------
def run_command(cmd, cwd=None, amber_sh=None, log_path=None, quiet=False):
    if amber_sh:
        wrapped = ["bash", "-c", 'source "$1" || exit 1; shift; exec "$@"',
                   "_", amber_sh] + cmd
    else:
        wrapped = cmd
    proc = subprocess.run(wrapped, cwd=cwd, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True)
    if log_path:
        with open(log_path, "w", encoding="utf-8") as fh:
            fh.write(proc.stdout or "")
    if not quiet and proc.returncode != 0:
        tail = "\n".join((proc.stdout or "").splitlines()[-30:])
        warn(f"command failed (rc={proc.returncode}): {' '.join(cmd[:4])}\n{tail}")
    return proc


# ---------------------------------------------------------------------------
# residue helpers
# ---------------------------------------------------------------------------
def residue_key(res):
    return (res["chain"], res["resseq"], res["icode"])


def heavy_atoms(res):
    return [a for a in res["atoms"] if a.get("element", "").upper() not in ("H", "D")]


def normalize_resname(resname):
    rn = (resname or "").strip().upper()
    return ic.AA_ALIASES.get(rn, rn)


def convert_mse(res):
    """Return a copy of the residue with MSE -> MET / SE -> SD."""
    if res["resname"].upper() != "MSE":
        return res, False
    out = dict(res)
    out["resname"] = "MET"
    atoms = []
    for atom in res["atoms"]:
        a = dict(atom)
        if a["name"].upper() == "SE":
            a["name"] = "SD"
            a["element"] = "S"
        if a.get("resname", "").upper() == "MSE":
            a["resname"] = "MET"
        atoms.append(a)
    out["atoms"] = atoms
    return out, True


def allowed_reduced_resname(raw_name, reduced_name):
    """Whether pdb4amber/reduce may rename a residue this way.

    Reduce may reassign histidine tautomers, Asn/Gln flips, and some
    protonation states based on the local hydrogen-bond network. The heavy-atom
    coordinates/names from its output are used, so these renames are safe; a
    change between unrelated residue types is not.
    """
    raw = normalize_resname(raw_name)
    red = normalize_resname(reduced_name)
    if raw == red:
        return True
    his = {"HIS", "HID", "HIE", "HIP", "HIN"}
    if raw in his and red in his:
        return True
    if raw in ("CYS", "CYX", "CYM") and red in ("CYS", "CYX", "CYM"):
        return True
    if raw in ("ASP", "ASH") and red in ("ASP", "ASH"):
        return True
    if raw in ("GLU", "GLH") and red in ("GLU", "GLH"):
        return True
    if raw in ("LYS", "LYN") and red in ("LYS", "LYN"):
        return True
    if raw == "MSE" and red == "MET":
        return True
    return False


def merge_disulfides(pairs):
    out = []
    seen = set()
    for pair in pairs:
        key = tuple(sorted((pair["index1"], pair["index2"])))
        if key in seen:
            continue
        seen.add(key)
        out.append(pair)
    return out


def find_disulfides(residues, header, ss_cutoff):
    """Detect SG-SG pairs and merge SSBOND records; return (pairs, warnings)."""
    warnings = []
    by_key = {residue_key(r): r for r in residues}
    idx_by_key = {residue_key(r): i for i, r in enumerate(residues, 1)}
    pairs = []
    sg = []
    for idx, res in enumerate(residues, 1):
        rn = res["resname"].upper()
        if rn not in ("CYS", "CYX"):
            continue
        s_atom = None
        for atom in heavy_atoms(res):
            if atom["name"].upper() == "SG":
                s_atom = atom
                break
        if s_atom is not None:
            sg.append((idx, res, s_atom))
    partners = {}
    for i in range(len(sg)):
        for j in range(i + 1, len(sg)):
            i1, r1, a1 = sg[i]
            i2, r2, a2 = sg[j]
            d = ic.distance(ic.atom_coord(a1), ic.atom_coord(a2))
            if d < ss_cutoff:
                partners.setdefault(i1, []).append(i2)
                partners.setdefault(i2, []).append(i1)
                pairs.append({
                    "index1": i1, "index2": i2,
                    "chain1": r1["chain"], "resseq1": r1["resseq"],
                    "chain2": r2["chain"], "resseq2": r2["resseq"],
                    "distance_A": round(d, 3), "source": "distance",
                })
    # SSBOND records are authoritative even if the geometry is slightly off.
    for line in header:
        if not line.startswith("SSBOND"):
            continue
        line = line.ljust(80)
        try:
            c1, s1, ic1 = line[15].strip(), int(line[17:21]), line[21].strip()
            c2, s2, ic2 = line[29].strip(), int(line[31:35]), line[35].strip()
        except Exception:
            continue
        r1 = by_key.get((c1, s1, ic1))
        r2 = by_key.get((c2, s2, ic2))
        if r1 is None or r2 is None:
            continue
        i1 = idx_by_key[(c1, s1, ic1)]
        i2 = idx_by_key[(c2, s2, ic2)]
        pairs.append({
            "index1": i1, "index2": i2,
            "chain1": c1, "resseq1": s1, "chain2": c2, "resseq2": s2,
            "distance_A": None, "source": "SSBOND",
        })
    for idx, plist in partners.items():
        if len(set(plist)) > 1:
            res = residues[idx - 1]
            die(f"cysteine {res['chain']}:{res['resseq']} has multiple SG "
                f"partners within {ss_cutoff} A; fix the structure or pass an "
                "explicit disulfide assignment")
    return merge_disulfides(pairs), warnings


def estimate_charge(residues, ligand_charge):
    """Estimate the net charge from residue names and free termini."""
    charge = 0
    for res in residues:
        rn = normalize_resname(res["resname"])
        if rn in ("ASP", "GLU", "CYM"):
            charge -= 1
        elif rn in ("LYS", "ARG", "HIP"):
            charge += 1
    # free termini: split at ter_before
    start = 0
    for idx in range(1, len(residues) + 1):
        if idx == len(residues) or residues[idx].get("ter_before"):
            seg = residues[start:idx]
            if seg:
                first = normalize_resname(seg[0]["resname"])
                last = normalize_resname(seg[-1]["resname"])
                if first not in ("ACE", "NME", "NH2", "NHE"):
                    charge += 1
                if last not in ("ACE", "NME", "NH2", "NHE"):
                    charge -= 1
            start = idx
    if ligand_charge is not None:
        charge += int(ligand_charge)
    return charge


def compress(nums):
    nums = sorted(set(int(x) for x in nums))
    if not nums:
        return ""
    parts = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
            continue
        parts.append(str(start) if start == prev else f"{start}-{prev}")
        start = prev = n
    parts.append(str(start) if start == prev else f"{start}-{prev}")
    return ",".join(parts)


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


# ---------------------------------------------------------------------------
# main preparation
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--components", default=None,
                    help="components.json from inspect_complex.py")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--model", type=int, default=1)
    ap.add_argument("--altloc", default="most-populous")
    ap.add_argument("--type", default="auto")
    ap.add_argument("--ligand", default=None, help="RESNAME[:CHAIN[:RESID]]")
    ap.add_argument("--ligand-resname", default="LIG")
    ap.add_argument("--charge", type=int, default=None)
    ap.add_argument("--ligand-smiles", default=None)
    ap.add_argument("--ph", default="7.4")
    ap.add_argument("--rdkit-python",
                    default="/opt/software/amber_general_workflow/envs/analysis/bin/python")
    ap.add_argument("--obabel", default="/opt/schrodinger2021-2/utilities/obabel")
    ap.add_argument("--amber-sh", default="/opt/software/amber22/amber.sh")
    ap.add_argument("--drop-glycans", action="store_true")
    ap.add_argument("--drop-metals", action="store_true")
    ap.add_argument("--drop-modified", action="store_true")
    ap.add_argument("--drop-unknown-hetatm", action="store_true")
    ap.add_argument("--drop-other-ligands", action="store_true")
    ap.add_argument("--ss-cutoff", type=float, default=2.5)
    ap.add_argument("--cn-break-cutoff", type=float, default=2.0)
    ap.add_argument("--strict-gaps", action="store_true")
    ap.add_argument("--no-reduce", action="store_true")
    ap.add_argument("--build-unit", choices=["auto", "com", "mol"], default="auto")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not os.path.isfile(args.pdb):
        die(f"PDB not found: {args.pdb}")
    os.makedirs(args.outdir, exist_ok=True)

    # ---- components -------------------------------------------------------
    if args.components:
        with open(args.components, encoding="utf-8") as fh:
            components = json.load(fh)
    else:
        opts = {
            "model": args.model, "altloc": args.altloc, "type": args.type,
            "ligand": args.ligand,
            "drop_glycans": args.drop_glycans,
            "drop_metals": args.drop_metals,
            "drop_modified": args.drop_modified,
            "drop_unknown_hetatm": args.drop_unknown_hetatm,
            "drop_other_ligands": args.drop_other_ligands,
            "strict_gaps": args.strict_gaps,
            "min_ligand_heavy": 3, "max_ligand_heavy": 200,
            "contact_cutoff": 4.0, "covalent_cutoff": 1.9,
            "ss_cutoff": args.ss_cutoff, "cn_break_cutoff": args.cn_break_cutoff,
            "peptide_max": 50,
        }
        components = ic.inspect(args.pdb, opts)

    cls = components["classification"]
    if cls.get("errors"):
        for err in cls["errors"]:
            print(f"error: [{err['code']}] {err['message']}", file=sys.stderr)
        die("inspection reported errors; fix them or use the documented escape "
            "hatch before preparing the system")

    # ---- re-parse with the requested model/altloc -------------------------
    models, header = ic.read_pdb(args.pdb)
    if args.model < 1 or args.model > len(models):
        die(f"--model {args.model} out of range 1..{len(models)}")
    atoms, _ = ic.choose_altlocs(models[args.model - 1], args.altloc)
    residues = ic.group_residues(atoms)

    warnings = []
    dropped = {"water": 0, "ions": 0, "additives": 0, "glycans": 0,
               "metals": 0, "modified": 0, "unknown": 0, "other_ligands": 0}
    mse_converted = []
    protein = []
    for res in residues:
        kind = ic.classify_resname(res["resname"])
        if kind == "protein":
            protein.append(res)
        elif kind == "modified_convertible":
            conv, changed = convert_mse(res)
            protein.append(conv)
            if changed:
                mse_converted.append(f"{res['chain']}:{res['resseq']}")
        elif kind == "water":
            dropped["water"] += 1
        elif kind == "ion_mono":
            dropped["ions"] += 1
        elif kind == "additive":
            dropped["additives"] += 1
        elif kind == "glycan":
            dropped["glycans"] += 1
        elif kind == "metal_divalent":
            dropped["metals"] += 1
        elif kind == "modified_hard":
            dropped["modified"] += 1
        elif kind == "nucleic":
            die("nucleic-acid residues are not supported by amber-auto")
        else:
            dropped["unknown"] += 1
    if not protein:
        die("no protein residues remain after preparation")

    # ---- segments / TER ---------------------------------------------------
    segments = ic.segment_residues(protein, args.cn_break_cutoff)
    prepared_raw = []
    seg_of = {}
    for si, seg in enumerate(segments):
        for ri, res in enumerate(seg["residues"]):
            res = dict(res)
            res["ter_before"] = (ri == 0 and si > 0)
            res["segment"] = seg["id"]
            prepared_raw.append(res)
            seg_of[residue_key(res)] = seg["id"]

    breaks = ic.chain_breaks(segments, args.cn_break_cutoff)
    if breaks and args.strict_gaps:
        die(f"{len(breaks)} chain break(s) detected and --strict-gaps is set")
    if breaks:
        warnings.append(
            f"{len(breaks)} chain break(s) will be capped with TER; each fragment "
            "gets charged termini")

    raw_pdb = os.path.join(args.outdir, "protein_raw_for_reduce.pdb")
    write_residues_pdb(raw_pdb, prepared_raw, include_h=False)

    # ---- Reduce -----------------------------------------------------------
    reduced_pdb = os.path.join(args.outdir, "protein_reduced.pdb")
    reduce_used = False
    if args.no_reduce:
        warnings.append("--no-reduce: histidines keep their input name; tleap "
                        "will choose HIE for HIS and Asn/Gln flips are not applied")
        shutil.copyfile(raw_pdb, reduced_pdb)
    else:
        log = os.path.join(args.outdir, "pdb4amber_reduce.log")
        proc = run_command(
            ["pdb4amber", "-i", raw_pdb, "-o", reduced_pdb, "--dry",
             "--most-populous", "--reduce", "-l", log],
            cwd=args.outdir, amber_sh=args.amber_sh, log_path=log,
            quiet=True)
        if proc.returncode != 0 or not os.path.isfile(reduced_pdb):
            die("pdb4amber --reduce failed; see "
                f"{log}. Rerun with --no-reduce to fall back to tleap defaults.")
        reduce_used = True

    # ---- verify Reduce preserved residue order and build prepared PDB -----
    rmodels, _ = ic.read_pdb(reduced_pdb)
    ratoms, _ = ic.choose_altlocs(rmodels[0], args.altloc)
    rresidues = ic.group_residues(ratoms)
    if len(rresidues) != len(prepared_raw):
        die(f"Reduce changed the residue count: {len(prepared_raw)} -> "
            f"{len(rresidues)}; refusing to build a mismatched topology")
    prepared = []
    histidines = []
    for i, (raw, red) in enumerate(zip(prepared_raw, rresidues)):
        if not allowed_reduced_resname(raw["resname"], red["resname"]):
            die(f"Reduce changed residue {i + 1} from {raw['resname']} to "
                f"{red['resname']}; refusing to continue")
        out = {
            "index": i + 1,
            "chain": raw["chain"],
            "resseq": raw["resseq"],
            "icode": raw["icode"],
            "resname": red["resname"].upper() if reduce_used else raw["resname"],
            "ter_before": raw.get("ter_before", False),
            "segment": raw.get("segment", ""),
            "atoms": [],
        }
        # Preserve original residue numbering/chain but take Reduce's heavy
        # atom coordinates, names, and flips.
        for atom in red["atoms"]:
            if atom.get("element", "").upper() in ("H", "D"):
                continue
            a = dict(atom)
            a["chain"] = raw["chain"]
            a["resseq"] = raw["resseq"]
            a["icode"] = raw["icode"]
            a["resname"] = out["resname"]
            out["atoms"].append(a)
        prepared.append(out)
        if out["resname"] in ("HID", "HIE", "HIP", "HIN"):
            histidines.append(f"{out['chain']}:{out['resseq']}={out['resname']}")

    # ---- disulfides -------------------------------------------------------
    disulfides, _ = find_disulfides(prepared, header, args.ss_cutoff)
    paired = set()
    for pair in disulfides:
        paired.add(pair["index1"])
        paired.add(pair["index2"])
    for idx, res in enumerate(prepared, 1):
        if idx in paired and res["resname"].upper() == "CYS":
            res["resname"] = "CYX"
        elif idx not in paired and res["resname"].upper() == "CYX":
            res["resname"] = "CYS"
            warnings.append(f"CYX {res['chain']}:{res['resseq']} has no SG "
                            "partner; renamed to CYS")
    for res in prepared:
        for atom in res["atoms"]:
            atom["resname"] = res["resname"]

    unit = args.build_unit
    if unit == "auto":
        unit = "com" if cls["type"] == "protein-ligand" else "mol"
    extra_leap = os.path.join(args.outdir, "extra_leap.in")
    with open(extra_leap, "w", encoding="utf-8") as fh:
        fh.write("# generated by prepare_complex.py\n")
        fh.write(f"# unit: {unit}\n")
        if disulfides:
            for pair in disulfides:
                fh.write(f"bond {unit}.{pair['index1']}.SG "
                         f"{unit}.{pair['index2']}.SG\n")
        else:
            fh.write("# no disulfide bonds\n")

    protein_pdb = os.path.join(args.outdir, "protein_prepared.pdb")
    write_residues_pdb(protein_pdb, prepared, include_h=False)
    # For protein-ligand builds the receptor is the protein-only PDB; for
    # peptide/antibody/apo builds the whole multi-chain protein is the complex.
    shutil.copyfile(protein_pdb, os.path.join(args.outdir, "receptor_prepared.pdb"))
    shutil.copyfile(protein_pdb, os.path.join(args.outdir, "complex_prepared.pdb"))

    # ---- ligand extraction ------------------------------------------------
    lig_payload = cls.get("ligand")
    ligand_input = os.path.join(args.outdir, "ligand_input.pdb")
    ligand_prepared = os.path.join(args.outdir, "ligand_prepared.pdb")
    ligand_charge_json = os.path.join(args.outdir, "ligand_charge.json")
    ligand_info = None
    lres = None
    if lig_payload:
        target = {
            "resname": lig_payload["resname"],
            "chain": lig_payload["chain"],
            "resseq": lig_payload["resseq"],
            "icode": lig_payload.get("icode", ""),
        }
        for res in residues:
            if (res["resname"] == target["resname"]
                    and res["chain"] == target["chain"]
                    and res["resseq"] == target["resseq"]
                    and res["icode"] == target["icode"]):
                lres = res
                break
        if lres is None:
            die(f"selected ligand {target} not found in the input PDB")
        lig_atoms = list(lres["atoms"])
        write_atoms_pdb(ligand_input, lig_atoms, resname=args.ligand_resname)
        lcmd = [args.rdkit_python, "-I",
                os.path.join(HERE, "ligand_charge.py"),
                "--pdb", ligand_input,
                "--out-pdb", ligand_prepared,
                "--json", ligand_charge_json,
                "--resname", args.ligand_resname,
                "--obabel", args.obabel,
                "--ph", str(args.ph),
                "--force"]
        if args.charge is not None:
            lcmd += ["--charge", str(args.charge)]
        if args.ligand_smiles:
            lcmd += ["--smiles", args.ligand_smiles]
        proc = run_command(lcmd, cwd=args.outdir, quiet=True)
        if proc.returncode != 0 or not os.path.isfile(ligand_prepared):
            tail = "\n".join((proc.stdout or "").splitlines()[-25:])
            die("ligand charge/protonation failed:\n" + tail)
        with open(ligand_charge_json, encoding="utf-8") as fh:
            ligand_info = json.load(fh)
        ligand_info["resname"] = args.ligand_resname
        ligand_info["original_resname"] = lres["resname"]
        ligand_info["original_chain"] = lres["chain"]
        ligand_info["original_resseq"] = lres["resseq"]
        if args.charge is not None and ligand_info["charge"] != args.charge:
            die(f"ligand charge {ligand_info['charge']} does not match "
                f"--charge {args.charge}")
        warnings.extend(ligand_info.get("warnings", []))

    # ---- residue map and roles -------------------------------------------
    rec_chains = set(cls.get("receptor", {}).get("chains", []))
    part_chains = set(cls.get("partner", {}).get("chains", []))
    partner_residues = []
    map_lines = ["prepared_index\trole\torig_chain\torig_resseq\torig_icode\t"
                 "orig_resname\tprepared_resname"]
    for idx, res in enumerate(prepared, 1):
        role = "receptor"
        if res["chain"] in part_chains and res["chain"] not in rec_chains:
            role = "partner"
        elif res["chain"] in part_chains and part_chains:
            role = "partner"
        if role == "partner":
            partner_residues.append(idx)
        map_lines.append(
            f"{idx}\t{role}\t{res['chain']}\t{res['resseq']}\t"
            f"{res['icode']}\t{res['resname']}\t{res['resname']}")
    residue_map = os.path.join(args.outdir, "residue_map.tsv")
    with open(residue_map, "w", encoding="utf-8") as fh:
        fh.write("\n".join(map_lines) + "\n")

    lig_charge = ligand_info["charge"] if ligand_info else None
    estimated = estimate_charge(prepared, lig_charge)
    counterion = "cation" if estimated <= 0 else "anion"

    manifest = {
        "status": "ok",
        "input": os.path.abspath(args.pdb),
        "components": os.path.abspath(args.components) if args.components else None,
        "model": args.model,
        "altloc": args.altloc,
        "type": cls["type"],
        "receptor_chains": sorted(rec_chains),
        "partner_chains": sorted(part_chains),
        "counts": {
            "protein_residues": len(prepared),
            "protein_atoms": sum(len(r["atoms"]) for r in prepared),
            "ligand_atoms": (len(lres["atoms"]) if lres is not None else 0),
        },
        "dropped": dropped,
        "mse_converted": mse_converted,
        "histidines": histidines,
        "reduce_used": reduce_used,
        "disulfides": disulfides,
        "chain_breaks": breaks,
        "missing_atoms": components.get("missing_atoms", []),
        "ligand": ligand_info,
        "partner_residues": partner_residues,
        "partner_ranges": compress(partner_residues),
        "estimated_net_charge": estimated,
        "suggested_counterion": counterion,
        "prepared": {
            "protein_pdb": os.path.abspath(protein_pdb),
            "receptor_pdb": os.path.abspath(
                os.path.join(args.outdir, "receptor_prepared.pdb")),
            "complex_pdb": os.path.abspath(
                os.path.join(args.outdir, "complex_prepared.pdb")),
            "ligand_input": os.path.abspath(ligand_input) if lres else None,
            "ligand_prepared": os.path.abspath(ligand_prepared) if lres else None,
            "ligand_charge_json": (os.path.abspath(ligand_charge_json)
                                   if lres else None),
            "extra_leap": os.path.abspath(extra_leap),
            "residue_map": os.path.abspath(residue_map),
        },
        "warnings": warnings,
    }
    with open(os.path.join(args.outdir, "prepare_manifest.json"), "w",
              encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    # ---- final safety check: prepared protein must only contain Amber names --
    allowed = (ic.AA_STANDARD | ic.AA_AMBER | ic.CAPS
               | {ic.AA_ALIASES.get(x, x) for x in ic.AA_ALIASES})
    for res in prepared:
        if res["resname"].upper() in allowed:
            continue
        if res["resname"].upper() == "MET":
            continue
        die(f"prepared protein contains unsupported residue "
            f"{res['resname']} {res['chain']}:{res['resseq']}")

    note(f"prepared {len(prepared)} protein residues, "
         f"{manifest['counts']['protein_atoms']} heavy atoms"
         + (f", ligand {ligand_info['charge']:+d}" if ligand_info else ""))
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
