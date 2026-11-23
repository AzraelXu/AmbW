#!/usr/bin/env python3
"""Hydrogen mass repartitioning (HMR) for an Amber prmtop.

The script changes each *non-water* hydrogen mass to ``--mass`` (default
3.024 Da) and subtracts the excess from the heavy atom it is bonded to, so the
total mass is conserved.  Water hydrogen masses are left unchanged by ParmEd's
``HMassRepartition`` (the SETTLE algorithm for water constraints is analytical),
which is what OPC/TIP3P/TIP4P-style models expect.

Outputs
-------
* ``--out``                 HMR topology (default ``<stem>_hmr.prmtop``)
* ``--backup``              original topology backup (default ``<prmtop>.orig``)
* ``--json``                machine-readable manifest (default
                            ``<outdir>/hmr_manifest.json``)
* ``--replace``             after the backup, replace ``--prmtop`` with the HMR
                            topology so all downstream stages use it

Safety checks
-------------
* the total mass before/after must agree to 1e-3 Da;
* every atom mass must remain positive (a too-large ``--mass`` can make a
  heavily hydrogenated heavy atom negative -- the script refuses);
* a topology that already looks like HMR (solute H > 2 Da) is only accepted
  when the original backup exists; otherwise the original cannot be recovered.

The script does not repartition water hydrogens.  It is intended for
SHAKE-constrained MD at dt <= 0.004 ps and does not make an unminimized
structure safe to start at 4 fs.
"""
import argparse
import datetime as _dt
import json
import os
import shutil
import sys

WATER_NAMES = {
    "WAT", "HOH", "H2O", "SOL", "TIP", "TIP3", "TP3", "T3P", "SPC", "SPCE",
    "OPC", "OP3", "T4P", "TIP4", "T5P", "FB3", "FB", "DOD",
}
STANDARD_MASS = 1.008


def die(msg: str, code: int = 2) -> None:
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(code)


def warn(msg: str) -> None:
    print(f"warning: {msg}", file=sys.stderr)


def note(msg: str) -> None:
    print(f"[hmr] {msg}", file=sys.stderr)


def is_water(atom) -> bool:
    return str(atom.residue.name).strip().upper() in WATER_NAMES


def hydrogen_groups(parm):
    solute = [a for a in parm.atoms
              if a.atomic_number == 1 and not is_water(a)]
    water = [a for a in parm.atoms
             if a.atomic_number == 1 and is_water(a)]
    return solute, water


def looks_hmr(parm) -> bool:
    solute, _ = hydrogen_groups(parm)
    if not solute:
        return False
    # Standard HMR is 3.024 Da; anything above 2 Da is not a normal H mass.
    return max(a.mass for a in solute) > 2.0


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prmtop", required=True,
                    help="Amber topology to convert")
    ap.add_argument("--out", default=None,
                    help="HMR topology (default: <stem>_hmr.prmtop)")
    ap.add_argument("--backup", default=None,
                    help="original topology backup (default: <prmtop>.orig)")
    ap.add_argument("--mass", type=float, default=3.024,
                    help="new non-water hydrogen mass in Da (default 3.024)")
    ap.add_argument("--json", default=None,
                    help="manifest JSON (default: <outdir>/hmr_manifest.json)")
    ap.add_argument("--replace", action="store_true",
                    help="replace --prmtop with the HMR topology after backup")
    ap.add_argument("--no-water-check", action="store_true",
                    help="skip the water-H unchanged check (not recommended)")
    ap.add_argument("--force", action="store_true",
                    help="overwrite existing outputs")
    args = ap.parse_args()

    if not os.path.isfile(args.prmtop):
        die(f"topology not found: {args.prmtop}")
    if not (args.mass > STANDARD_MASS):
        die(f"--mass must be > {STANDARD_MASS} Da (got {args.mass})")

    stem = args.prmtop[:-len(".prmtop")] if args.prmtop.endswith(".prmtop") \
        else args.prmtop
    out = args.out or f"{stem}_hmr.prmtop"
    backup = args.backup or f"{args.prmtop}.orig"
    manifest_path = args.json or os.path.join(
        os.path.dirname(os.path.abspath(out)), "hmr_manifest.json")

    for path in (out, manifest_path):
        if os.path.exists(path) and not args.force:
            die(f"refusing to overwrite {path} (use --force)")
    if os.path.exists(backup) and not args.force:
        die(f"refusing to overwrite backup {backup} (use --force)")

    try:
        import parmed as pmd
        from parmed.tools import HMassRepartition
    except Exception as exc:  # pragma: no cover - environment dependent
        die(f"ParmEd is required for HMR: {exc}")

    parm_in = pmd.load_file(args.prmtop)
    if looks_hmr(parm_in):
        if not os.path.isfile(backup):
            die(f"input {args.prmtop} already looks like HMR (solute H > 2 Da) "
                f"and no original backup exists at {backup}; rebuild the system "
                f"or restore the original topology")
        parm = pmd.load_file(backup)
        source = backup
        note(f"input already HMR; using original backup {backup}")
    else:
        parm = parm_in
        source = args.prmtop
        os.makedirs(os.path.dirname(os.path.abspath(backup)), exist_ok=True)
        shutil.copyfile(args.prmtop, backup)
        note(f"backup: {backup}")

    solute_before, water_before = hydrogen_groups(parm)
    if not solute_before:
        warn("no non-water hydrogens found; HMR is a no-op for this topology")
    water_masses_before = [a.mass for a in water_before]
    mass_before = sum(a.mass for a in parm.atoms)

    HMassRepartition(parm, f"{args.mass}").execute()

    solute_after, water_after = hydrogen_groups(parm)
    water_masses_after = [a.mass for a in water_after]
    mass_after = sum(a.mass for a in parm.atoms)

    if abs(mass_after - mass_before) > 1e-3:
        die(f"total mass changed by {mass_after - mass_before:.6f} Da; "
            f"refusing to write an inconsistent topology")
    # Zero-mass virtual sites (e.g. the OPC/TIP4P EPW site) are normal and must
    # not be repartitioned; ParmEd itself only rejects non-positive real atoms.
    bad = [a for a in parm.atoms if a.mass <= 0.0 and a.atomic_number > 0]
    if bad:
        a = bad[0]
        die(f"{len(bad)} atom(s) would have non-positive mass after HMR "
            f"(first: {a.name} in {a.residue.name}); reduce --mass")
    if solute_after:
        lo = min(a.mass for a in solute_after)
        hi = max(a.mass for a in solute_after)
        if abs(hi - args.mass) > 1e-3:
            die(f"solute hydrogen mass {hi:.4f} does not match --mass "
                f"{args.mass:.4f}")
    if not args.no_water_check and water_before:
        if (len(water_masses_before) != len(water_masses_after)
                or any(abs(x - y) > 1e-6
                       for x, y in zip(water_masses_before, water_masses_after))):
            die("water hydrogen masses changed unexpectedly; ParmEd should "
                "leave them untouched (use --no-water-check only if the water "
                "model genuinely requires repartitioning)")

    parm.save(out, overwrite=True)
    if args.replace:
        tmp = args.prmtop + ".hmr.tmp"
        shutil.copyfile(out, tmp)
        os.replace(tmp, args.prmtop)

    manifest = {
        "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
        "generator": "amber_general_workflow/scripts/hmr_topology.py",
        "status": "ok",
        "enabled": True,
        "applied": bool(solute_after),
        "input": os.path.abspath(args.prmtop),
        "source": os.path.abspath(source),
        "out": os.path.abspath(out),
        "backup": os.path.abspath(backup),
        "replace": bool(args.replace),
        "hydrogen_mass": args.mass,
        "water_hydrogens_repartitioned": False,
        "atoms": len(parm.atoms),
        "residues": len(parm.residues),
        "hydrogens_total": len(solute_after) + len(water_after),
        "hydrogens_solute": len(solute_after),
        "hydrogens_water": len(water_after),
        "total_mass_before": mass_before,
        "total_mass_after": mass_after,
        "solute_h_mass_min": min((a.mass for a in solute_after), default=None),
        "solute_h_mass_max": max((a.mass for a in solute_after), default=None),
        "water_h_mass_min": min(water_masses_after, default=None),
        "water_h_mass_max": max(water_masses_after, default=None),
        "parmed_version": getattr(pmd, "__version__", "unknown"),
    }
    os.makedirs(os.path.dirname(os.path.abspath(manifest_path)), exist_ok=True)
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    note(f"HMR topology: {out}")
    note(f"mass {mass_before:.4f} -> {mass_after:.4f} Da; "
         f"solute H {len(solute_after)} x {args.mass:.4f} Da; "
         f"water H {len(water_after)} unchanged")
    if not solute_after:
        note("no solute H: HMR did not change the topology; 4 fs is not justified")
    print(json.dumps({
        "status": "ok",
        "out": os.path.abspath(out),
        "backup": os.path.abspath(backup),
        "applied": bool(solute_after),
        "solute_h": len(solute_after),
        "water_h": len(water_after),
        "mass_before": mass_before,
        "mass_after": mass_after,
    }))


if __name__ == "__main__":
    main()
