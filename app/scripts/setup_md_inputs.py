#!/usr/bin/env python3
"""Generate Amber minimization/heating/equilibration/production inputs.

The script is a careful template renderer, not a black box: it refuses illegal
explicit/implicit combinations, validates ns -> nstlim arithmetic, checks the
restraint mask against the topology when parmed is available, and writes a
machine-readable run_config.json plus a human stage table.

Typical explicit run:
  python3 setup_md_inputs.py --out md/complex/r1 --seed 12345 \
    --receptor-type rna --solvent explicit \
    --restraint-mask ':1-20 & !@H=' \
    --prmtop systems/complex.prmtop --ligand-resname LIG

Typical implicit run:
  python3 setup_md_inputs.py --out md/gb/r1 --seed 12345 \
    --receptor-type protein --solvent implicit \
    --restraint-mask ':1-100 & !@H=' --igb 8 --saltcon 0.0

Run order is fixed; production chunks are stateful. See the generated README.md.
"""
import argparse
import datetime as _dt
import json
import math
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_ROOT = os.path.normpath(os.path.join(HERE, "..", "templates", "md"))
WATER_NAMES = {
    "WAT", "HOH", "H2O", "SOL", "TIP", "TIP3", "TP3", "T3P", "SPC", "SPCE",
    "OPC", "OP3", "T4P", "TIP4", "T5P", "FB3", "FB", "DOD",
}
ION_NAMES = {
    "NA", "NA+", "K", "K+", "CL", "CL-", "MG", "MG2", "CA", "CA2", "ZN",
    "ZN2", "LI", "LI+", "RB", "RB+", "CS", "CS+", "F", "F-", "BR", "BR-",
    "I", "I-", "BA", "BA2", "SR", "SR2", "MN", "MN2", "FE", "FE2", "FE3",
    "CU", "CU2", "CO", "CO2", "NI", "NI2", "CD", "CD2", "AL", "AL3",
}
DEFAULT_RELEASE_WEIGHTS = [10.0, 8.0, 6.0, 4.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.0]


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def warn(msg):
    print(f"warning: {msg}", file=sys.stderr)


def steps(ns, dt, what):
    """Convert a length in ns to an integer number of dt-ps steps."""
    if ns is None or ns <= 0:
        die(f"{what} length must be > 0 ns")
    x = ns * 1000.0 / dt
    n = int(round(x))
    if abs(x - n) > 1e-6:
        die(f"{what}: {ns} ns / {dt} ps is not an integer number of steps ({x})")
    if n <= 0:
        die(f"{what}: nstlim must be positive")
    return n


def parse_weights(text):
    parts = [p for p in text.replace(";", ",").replace(" ", ",").split(",") if p]
    try:
        weights = [float(p) for p in parts]
    except ValueError:
        die(f"cannot parse --release-weights: {text!r}")
    if not weights:
        die("--release-weights is empty")
    if any(w < 0 for w in weights):
        die("release weights must be non-negative")
    return weights


def release_weights(args):
    if args.release_weights:
        weights = parse_weights(args.release_weights)
        n = len(weights)
        if args.release_stages not in (None, n):
            die(f"--release-stages {args.release_stages} conflicts with "
                f"{n} --release-weights")
        return weights
    n = args.release_stages
    if n is None:
        n = 10
    if n < 0:
        die("--release-stages must be >= 0")
    if n == 0:
        return []
    if n == 10:
        return list(DEFAULT_RELEASE_WEIGHTS)
    if n == 1:
        return [0.0]
    return [round(10.0 * (1.0 - i / (n - 1)), 4) for i in range(n)]


def restraint_block(mask, wt):
    if mask.lower() == "none":
        return "  ntr=0,"
    return (f"  ntr=1, restraint_wt={wt:.3f},\n"
            f"  restraintmask='{mask}',")


def release_stem(i, n):
    if n <= 26:
        return f"04{chr(ord('a') + i)}"
    return f"04_{i + 1:02d}"


def render(path, values):
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    for key, value in values.items():
        text = text.replace("{{" + key + "}}", str(value))
    if "{{" in text:
        left = sorted({x.split("}}")[0] for x in text.split("{{")[1:]})
        die(f"unresolved placeholders in {path}: {left}")
    return text


def check_prmtop(path, args, mask, solvent):
    """Best-effort topology checks. Returns a dict for the manifest."""
    info = {"checked": False}
    if not path:
        return info
    if not os.path.isfile(path):
        die(f"--prmtop not found: {path}")
    try:
        import parmed as pmd
        from parmed.amber import AmberMask
    except Exception as exc:  # parmed may be unavailable under isolated Python
        warn(f"parmed unavailable ({exc.__class__.__name__}); "
             "topology/mask checks skipped")
        return info

    try:
        parm = pmd.load_file(path)
    except Exception as exc:
        die(f"cannot read topology with parmed: {exc}")

    info["checked"] = True
    info["atoms"] = len(parm.atoms)
    info["residues"] = len(parm.residues)
    info["net_charge"] = round(sum(a.charge for a in parm.atoms), 4)
    try:
        info["radii"] = " ".join(parm.parm_data.get("RADIUS_SET", []))
    except Exception:
        info["radii"] = ""

    water = [r for r in parm.residues if r.name.strip().upper() in WATER_NAMES]
    ions = [r for r in parm.residues if r.name.strip().upper() in ION_NAMES]
    info["water_residues"] = len(water)
    info["ion_residues"] = len(ions)

    if solvent == "explicit" and not water:
        die("--solvent explicit but the topology contains no water residues; "
            "did you build an implicit/dry system?")
    if solvent == "implicit" and water:
        die(f"--solvent implicit but the topology contains {len(water)} water "
            "residues; rebuild without solvation")

    if abs(info["net_charge"]) > 0.01:
        warn(f"topology net charge is {info['net_charge']:+.3f}; "
             "check addIons/neutralization")

    if mask.lower() != "none":
        try:
            sel = AmberMask(parm, mask).Selection()
            if not isinstance(sel, list):
                sel = list(sel)
        except Exception as exc:
            die(f"cannot parse restraint mask {mask!r}: {exc}")
        nsel = sum(1 for x in sel if x)
        info["restraint_atoms"] = nsel
        if nsel == 0:
            die(f"restraint mask {mask!r} selects no atoms")
        bad = sorted({a.residue.name.strip() for a, ok in zip(parm.atoms, sel)
                      if ok and a.residue.name.strip().upper() in
                      (WATER_NAMES | ION_NAMES)})
        if bad:
            die(f"restraint mask {mask!r} selects water/ions: {bad}; "
                "restrain solute heavy atoms only")
        if args.ligand_resname:
            lig_hit = any(ok and a.residue.name.strip() == args.ligand_resname
                          for a, ok in zip(parm.atoms, sel))
            info["ligand_in_restraint_mask"] = bool(lig_hit)
            if not lig_hit:
                warn(f"ligand {args.ligand_resname!r} is not in the restraint "
                     "mask; this is allowed but confirm it is intentional")

    if solvent == "implicit":
        igb = args.igb
        radii = info.get("radii", "")
        if igb == 8 and "mbondi3" not in radii:
            warn(f"igb=8 normally requires mbondi3 radii; topology says "
                 f"{radii!r}. Rebuild with --radii mbondi3.")
        elif igb == 5 and "mbondi2" not in radii:
            warn(f"igb=5 normally requires mbondi2 radii; topology says "
                 f"{radii!r}. Rebuild with --radii mbondi2.")
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", help="replica directory")
    ap.add_argument("--seed", type=int, help="base random seed")
    ap.add_argument("--receptor-type", choices=["rna", "protein"])
    ap.add_argument("--solvent", choices=["explicit", "implicit"])
    ap.add_argument("--restraint-mask",
                    help="Amber mask; use the literal 'none' for unrestrained")
    ap.add_argument("--prmtop", help="topology to validate mask/water/radii")
    ap.add_argument("--ligand-resname", help="ligand residue name, for reporting")

    ap.add_argument("--heat-ns", type=float, default=None)
    ap.add_argument("--release-stages", type=int, default=None)
    ap.add_argument("--release-ns", type=float, default=None)
    ap.add_argument("--release-weights", default=None)
    ap.add_argument("--prod-chunks", type=int, default=None)
    ap.add_argument("--prod-chunk-ns", type=float, default=None)

    ap.add_argument("--dt", type=float, default=None)
    ap.add_argument("--hmr", action="store_true",
                    help="topology uses hydrogen mass repartitioning; "
                         "recorded in run_config.json (does not alter the .in files)")
    ap.add_argument("--hmr-mass", type=float, default=None,
                    help="non-water hydrogen mass in Da (default 3.024)")
    ap.add_argument("--cut", type=float, default=None)
    ap.add_argument("--temp", type=float, default=None)
    ap.add_argument("--pres", type=float, default=None)
    ap.add_argument("--gamma", type=float, default=None)
    ap.add_argument("--taup", type=float, default=None)
    ap.add_argument("--igb", type=int, default=None)
    ap.add_argument("--saltcon", type=float, default=None)
    ap.add_argument("--rgbmax", type=float, default=None)

    ap.add_argument("--min-maxcyc", type=int, default=None)
    ap.add_argument("--min-ncyc", type=int, default=None)
    ap.add_argument("--min-restraint-wt", type=float, default=None)
    ap.add_argument("--traj-interval-ps", type=float, default=None)
    ap.add_argument("--ntpr", type=int, default=None)
    ap.add_argument("--ntwr", type=int, default=None)
    ap.add_argument("--bench-ns", type=float, default=None)
    ap.add_argument("--seed-mode", choices=["stage", "constant"], default=None)
    ap.add_argument("--nmr-restraints", help="Amber NMR restraint file")

    ap.add_argument("--config", help="JSON file with any of the above keys")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--print-only", action="store_true")

    args = ap.parse_args()

    # Config file fills values not explicitly given on the command line.
    if args.config:
        provided = {a.split("=")[0].lstrip("-").replace("-", "_")
                    for a in sys.argv[1:] if a.startswith("--")}
        try:
            cfg = json.load(open(args.config, encoding="utf-8"))
        except Exception as exc:
            die(f"cannot read --config: {exc}")
        if not isinstance(cfg, dict):
            die("--config must contain a JSON object")
        for key, value in cfg.items():
            dest = key.lstrip("-").replace("-", "_")
            if dest in provided:
                continue
            if hasattr(args, dest):
                setattr(args, dest, value)
            else:
                warn(f"unknown config key ignored: {key}")

    # Config files are untyped JSON; normalize before doing arithmetic.
    def as_int(name):
        value = getattr(args, name)
        if value is None:
            return
        try:
            setattr(args, name, int(value))
        except (TypeError, ValueError):
            die(f"{name} must be an integer, got {value!r}")

    def as_float(name):
        value = getattr(args, name)
        if value is None:
            return
        try:
            setattr(args, name, float(value))
        except (TypeError, ValueError):
            die(f"{name} must be a number, got {value!r}")

    for name in ("seed", "release_stages", "prod_chunks", "igb",
                 "min_maxcyc", "min_ncyc", "ntpr", "ntwr"):
        as_int(name)
    for name in ("heat_ns", "release_ns", "prod_chunk_ns", "dt", "cut",
                 "temp", "pres", "gamma", "taup", "saltcon", "rgbmax",
                 "min_restraint_wt", "traj_interval_ps", "bench_ns"):
        as_float(name)

    for name in ("out", "seed", "receptor_type", "solvent", "restraint_mask"):
        if getattr(args, name) is None:
            die(f"--{name.replace('_', '-')} is required (CLI or --config)")
    if not args.out:
        die("--out is empty")
    if args.receptor_type not in ("rna", "protein"):
        die("--receptor-type must be rna or protein")
    if args.solvent not in ("explicit", "implicit"):
        die("--solvent must be explicit or implicit")
    if args.seed_mode is None:
        args.seed_mode = "stage"
    elif args.seed_mode not in ("stage", "constant"):
        die("--seed-mode must be stage or constant")
    if args.release_weights is not None and not isinstance(args.release_weights, str):
        die("--release-weights must be a comma-separated string")

    # Resolved defaults
    args.heat_ns = 2.0 if args.heat_ns is None else args.heat_ns
    args.release_ns = 1.0 if args.release_ns is None else args.release_ns
    args.prod_chunks = 4 if args.prod_chunks is None else args.prod_chunks
    args.prod_chunk_ns = 25.0 if args.prod_chunk_ns is None else args.prod_chunk_ns
    args.dt = 0.002 if args.dt is None else args.dt
    args.temp = 300.0 if args.temp is None else args.temp
    args.pres = 1.0 if args.pres is None else args.pres
    args.gamma = 1.0 if args.gamma is None else args.gamma
    args.taup = 2.0 if args.taup is None else args.taup
    args.igb = 8 if args.igb is None else args.igb
    args.saltcon = 0.0 if args.saltcon is None else args.saltcon
    args.rgbmax = 25.0 if args.rgbmax is None else args.rgbmax
    args.min_maxcyc = 10000 if args.min_maxcyc is None else args.min_maxcyc
    args.min_ncyc = 5000 if args.min_ncyc is None else args.min_ncyc
    args.min_restraint_wt = 10.0 if args.min_restraint_wt is None else args.min_restraint_wt
    args.traj_interval_ps = 10.0 if args.traj_interval_ps is None else args.traj_interval_ps
    args.ntpr = 5000 if args.ntpr is None else args.ntpr
    args.ntwr = 50000 if args.ntwr is None else args.ntwr
    args.bench_ns = 0.2 if args.bench_ns is None else args.bench_ns
    args.seed_mode = "stage" if args.seed_mode is None else args.seed_mode
    if args.cut is None:
        args.cut = 10.0 if args.solvent == "explicit" else 999.0

    # ---- validation -------------------------------------------------------
    if args.seed <= 0:
        die("--seed must be a positive integer")
    if args.dt <= 0:
        die("--dt must be positive")
    hmr_mass = 3.024 if args.hmr_mass is None else args.hmr_mass
    if args.hmr_mass is not None and args.hmr_mass <= 1.008:
        die("--hmr-mass must be > 1.008 Da")
    if args.dt > 0.002:
        if args.hmr:
            warn(f"dt={args.dt} ps with --hmr; confirm the prmtop really has "
                 "repartitioned non-water hydrogen masses")
        else:
            warn(f"dt={args.dt} ps is larger than the SHAKE-safe 0.002 ps; "
                 "pass --hmr and use an HMR topology (otherwise 4 fs is unsafe)")
    if args.heat_ns <= 0:
        die("--heat-ns must be > 0")
    if args.prod_chunks < 1:
        die("--prod-chunks must be >= 1")
    if args.prod_chunk_ns <= 0:
        die("--prod-chunk-ns must be > 0")
    if args.min_maxcyc < 1 or args.min_ncyc < 0 or args.min_ncyc > args.min_maxcyc:
        die("require 1 <= --min-ncyc <= --min-maxcyc")
    if args.ntpr < 1 or args.ntwr < 1:
        die("--ntpr/--ntwr must be positive")
    if args.traj_interval_ps <= 0:
        die("--traj-interval-ps must be positive")
    if args.solvent == "explicit":
        if args.cut > 20.0:
            warn(f"explicit-solvent cut={args.cut} A is unusually large")
        if args.igb not in (8,) or args.saltcon not in (0.0,):
            # igb/saltcon are simply ignored here; warn if they look intentional
            if args.igb != 8 or args.saltcon != 0.0:
                warn("--igb/--saltcon are ignored for explicit solvent")
    else:
        if args.igb not in (1, 2, 5, 7, 8):
            die("--igb must be one of 1, 2, 5, 7, 8")
        if args.cut < 100.0:
            warn(f"GB runs normally use cut=999.0; using cut={args.cut}")
        if args.pres != 1.0:
            warn("--pres is ignored for implicit solvent (no barostat)")

    weights = release_weights(args)
    n_rel = len(weights)
    if n_rel and args.release_ns <= 0:
        die("--release-ns must be > 0 when release stages are requested")
    if weights and weights[-1] != 0.0:
        warn(f"final release weight is {weights[-1]}, not 0; production is "
             "unrestrained but NPT ends restrained")
    if any(w > 0 for w in weights) and args.restraint_mask.lower() == "none":
        warn("release weights > 0 with --restraint-mask none; all stages will "
             "be rendered with ntr=0")

    nstlim_heat = steps(args.heat_ns, args.dt, "heat")
    nstlim_rel = steps(args.release_ns, args.dt, "release") if n_rel else 0
    nstlim_prod = steps(args.prod_chunk_ns, args.dt, "production chunk")
    nstlim_bench = steps(args.bench_ns, args.dt, "benchmark")
    ntwx_f = args.traj_interval_ps / args.dt
    ntwx = int(round(ntwx_f))
    if abs(ntwx_f - ntwx) > 1e-6 or ntwx < 1:
        die(f"--traj-interval-ps {args.traj_interval_ps} / dt {args.dt} is not "
            "a positive integer step count")
    if ntwx > nstlim_prod:
        warn("trajectory interval is longer than a production chunk")

    # Seeds: by default advance per stage so independent runs do not reuse a stream.
    def stage_seed(k):
        return args.seed if args.seed_mode == "constant" else args.seed + k

    n_total = 2 + n_rel + args.prod_chunks + 1
    for k in range(n_total):
        if stage_seed(k) > 2147483647:
            die("seed overflow: keep --seed + number of stages <= 2147483647")

    nmr = None
    if args.nmr_restraints:
        if not os.path.isfile(args.nmr_restraints):
            die(f"--nmr-restraints not found: {args.nmr_restraints}")
        nmr = os.path.abspath(args.nmr_restraints)

    # Topology checks before writing anything.
    top_info = check_prmtop(args.prmtop, args, args.restraint_mask, args.solvent)

    # ---- render -----------------------------------------------------------
    tdir = os.path.join(TEMPLATE_ROOT, args.solvent)
    if not os.path.isdir(tdir):
        die(f"template directory missing: {tdir}")

    common = {
        "CUT": f"{args.cut:g}", "NTPR": args.ntpr, "NTWR": args.ntwr,
        "DT": f"{args.dt:g}", "TEMP": f"{args.temp:g}",
        "GAMMA": f"{args.gamma:g}", "PRES": f"{args.pres:g}",
        "TAUP": f"{args.taup:g}", "NTWX": ntwx,
    }
    if args.solvent == "implicit":
        common.update({"IGB": args.igb, "SALTCON": f"{args.saltcon:g}",
                       "RGBMAX": f"{args.rgbmax:g}"})

    files = {}
    # min 1
    files["01_min_restraint.in"] = render(
        os.path.join(tdir, "min_restrained.in"),
        dict(common, MIN_MAXCYC=args.min_maxcyc, MIN_NCYC=args.min_ncyc,
             RESTRAINT_BLOCK=restraint_block(args.restraint_mask,
                                             args.min_restraint_wt)))
    # min 2
    files["02_min_free.in"] = render(
        os.path.join(tdir, "min_free.in"),
        dict(common, MIN_MAXCYC=args.min_maxcyc, MIN_NCYC=args.min_maxcyc))
    # heat
    files["03_heat.in"] = render(
        os.path.join(tdir, "heat.in"),
        dict(common, HEAT_NSTLIM=nstlim_heat, HEAT_NS=f"{args.heat_ns:g}",
             SEED=stage_seed(0),
             RESTRAINT_BLOCK=restraint_block(args.restraint_mask,
                                             args.min_restraint_wt),
             NMR_LINE=("DISANG=noe.rst" if nmr else "")))
    # release
    rel_seeds = []
    for i, weight in enumerate(weights):
        seed_i = stage_seed(1 + i)
        rel_seeds.append(seed_i)
        block = ("  ntr=0," if weight == 0.0
                 else restraint_block(args.restraint_mask, weight))
        name = release_stem(i, n_rel) + ".in"
        files[name] = render(
            os.path.join(tdir, "npt_release.in" if args.solvent == "explicit"
                         else "release.in"),
            dict(common, RELEASE_NSTLIM=nstlim_rel,
                 RELEASE_NS=f"{args.release_ns:g}", SEED=seed_i,
                 RESTRAINT_BLOCK=block,
                 NMR_CNTRL=("  nmropt=1," if nmr else ""),
                 NMR_TAIL=("&wt type='END', &end\nDISANG=noe.rst" if nmr else "")))
    # production
    prod_seeds = []
    for c in range(1, args.prod_chunks + 1):
        seed_i = stage_seed(1 + n_rel + c - 1)
        prod_seeds.append(seed_i)
        files[f"05_prod_c{c:02d}.in"] = render(
            os.path.join(tdir, "prod.in"),
            dict(common, PROD_NSTLIM=nstlim_prod,
                 PROD_NS=f"{args.prod_chunk_ns:g}", SEED=seed_i,
                 NMR_CNTRL=("  nmropt=1," if nmr else ""),
                 NMR_TAIL=("&wt type='END', &end\nDISANG=noe.rst" if nmr else "")))
    # benchmark
    files["bench.in"] = render(
        os.path.join(tdir, "bench.in"),
        dict(common, BENCH_NSTLIM=nstlim_bench, BENCH_NS=f"{args.bench_ns:g}",
             SEED=stage_seed(1 + n_rel + args.prod_chunks)))

    # ---- stage table ------------------------------------------------------
    restraint_active = args.restraint_mask.lower() != "none"
    min_ntr = "ntr=1" if restraint_active else "ntr=0"
    min_wt = args.min_restraint_wt if restraint_active else "—"
    rows = [("01_min_restraint.in", "min", "—", 0, min_ntr, min_wt,
             "—", "restrained minimization"),
            ("02_min_free.in", "min", "—", 0, "ntr=0", "—", "—",
             "unrestrained minimization"),
            ("03_heat.in", "NVT", f"{args.heat_ns:g}", nstlim_heat,
             min_ntr, min_wt, stage_seed(0), "0 -> T ramp")]
    for i, weight in enumerate(weights):
        active = restraint_active and weight != 0.0
        rows.append((release_stem(i, n_rel) + ".in",
                     "NPT" if args.solvent == "explicit" else "NVT",
                     f"{args.release_ns:g}", nstlim_rel,
                     "ntr=1" if active else "ntr=0",
                     weight if restraint_active else "—",
                     rel_seeds[i], f"release stage {i + 1}/{n_rel}"))
    for c in range(1, args.prod_chunks + 1):
        rows.append((f"05_prod_c{c:02d}.in", "NPT" if args.solvent == "explicit" else "NVT",
                     f"{args.prod_chunk_ns:g}", nstlim_prod, "ntr=0", "—",
                     prod_seeds[c - 1], f"production chunk {c}"))
    rows.append(("bench.in", "NPT" if args.solvent == "explicit" else "NVT",
                 f"{args.bench_ns:g}", nstlim_bench, "ntr=0", "—",
                 stage_seed(1 + n_rel + args.prod_chunks), "benchmark"))

    table = "stage\tfile\tensemble\tlength_ns\tnstlim\tntr\trestraint_wt\tseed\tnotes\n"
    for r in rows:
        table += "\t".join(str(x) for x in r) + "\n"

    if args.print_only:
        print(table)
        for name in files:
            print(f"\n===== {name} =====")
            print(files[name], end="")
        return

    # ---- write ------------------------------------------------------------
    os.makedirs(args.out, exist_ok=True)
    managed = list(files) + ["run_config.json", "stage_table.tsv", "README.md"]
    if nmr:
        managed.append("noe.rst")
    existing = [n for n in managed if os.path.exists(os.path.join(args.out, n))]
    if existing and not args.force:
        die("refusing to overwrite existing files (use --force): "
            + ", ".join(existing))
    # A shorter release/production schedule must not leave stale stages behind.
    stage_like = [n for n in os.listdir(args.out)
                  if n.endswith(".in") and (
                      n.startswith(("01_", "02_", "03_", "04", "05_prod"))
                      or n == "bench.in")]
    stale = sorted(set(stage_like) - set(files))
    if stale:
        if not args.force:
            die("stale generated inputs would remain: " + ", ".join(stale)
                + " (remove them or use --force)")
        for name in stale:
            os.remove(os.path.join(args.out, name))

    for name, text in files.items():
        with open(os.path.join(args.out, name), "w", encoding="utf-8") as fh:
            fh.write(text)
    with open(os.path.join(args.out, "stage_table.tsv"), "w", encoding="utf-8") as fh:
        fh.write(table)

    if nmr:
        target = os.path.join(args.out, "noe.rst")
        if os.path.abspath(nmr) != os.path.abspath(target):
            shutil.copyfile(nmr, target)

    exec_cmd = "pmemd.cuda"
    engine_note = ("" if args.solvent == "explicit" else
                   "\nGB runs also work with `sander.MPI`; `pmemd.cuda` supports GB.")
    readme = f"""# MD run inputs

Generated by `amber_general_workflow/scripts/setup_md_inputs.py`.

- Solvent: **{args.solvent}**
- Restraint mask: `{args.restraint_mask}` (literal `none` = unrestrained)
- Heat: {args.heat_ns} ns; release: {n_rel} x {args.release_ns} ns;
  production: {args.prod_chunks} x {args.prod_chunk_ns} ns; dt {args.dt} ps.
- Seeds: base {args.seed}, mode `{args.seed_mode}` (see stage_table.tsv).

Stage order is mandatory; each stage reads the previous restart.

```bash
# 换成你所在集群的 Amber 初始化脚本
source "$AMBERHOME/amber.sh"

# 0. minimization (restrained min needs -ref)
{exec_cmd} -O -i 01_min_restraint.in -o 01_min_restraint.out \\
  -p complex.prmtop -c complex.inpcrd -ref complex.inpcrd \\
  -r 01_min_restraint.rst7 -inf 01_min_restraint.mdinfo
{exec_cmd} -O -i 02_min_free.in -o 02_min_free.out \\
  -p complex.prmtop -c 01_min_restraint.rst7 \\
  -r 02_min_free.rst7 -inf 02_min_free.mdinfo

# 1. heating (NVT)
{exec_cmd} -O -i 03_heat.in -o 03_heat.out \\
  -p complex.prmtop -c 02_min_free.rst7 -ref 02_min_free.rst7 \\
  -r 03_heat.rst7 -x 03_heat.nc -inf 03_heat.mdinfo

# 2. release stages, in order (03_heat.rst7 -> 04a.rst7 -> ... -> 04j.rst7)
#    Omit -ref if the stage uses ntr=0.
{exec_cmd} -O -i 04a.in -o 04a.out -p complex.prmtop \\
  -c 03_heat.rst7 -ref 03_heat.rst7 -r 04a.rst7 -x 04a.nc -inf 04a.mdinfo
# ... repeat for every 04*.in in sorted order ...

# 3. production chunks, strictly in order
{exec_cmd} -O -i 05_prod_c01.in -o 05_prod_c01.out -p complex.prmtop \\
  -c <last_release>.rst7 -r 05_prod_c01.rst7 -x 05_prod_c01.nc \\
  -inf 05_prod_c01.mdinfo
# chunk 2 reads 05_prod_c01.rst7, chunk 3 reads 05_prod_c02.rst7, ...
```

Submit with the `cluster3-torque` skill; do not run on the login node.{engine_note}
"""
    with open(os.path.join(args.out, "README.md"), "w", encoding="utf-8") as fh:
        fh.write(readme)

    config = {
        "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
        "generator": "amber_general_workflow/scripts/setup_md_inputs.py",
        "out": os.path.abspath(args.out),
        "seed": args.seed,
        "seed_mode": args.seed_mode,
        "receptor_type": args.receptor_type,
        "solvent": args.solvent,
        "restraint_mask": args.restraint_mask,
        "restraint_active": args.restraint_mask.lower() != "none",
        "ligand_resname": args.ligand_resname,
        "prmtop": os.path.abspath(args.prmtop) if args.prmtop else None,
        "topology_qc": top_info,
        "dt_ps": args.dt,
        "cut_A": args.cut,
        "temperature_K": args.temp,
        "pressure_atm": args.pres,
        "gamma_ps": args.gamma,
        "taup_ps": args.taup,
        "igb": args.igb if args.solvent == "implicit" else None,
        "saltcon_M": args.saltcon if args.solvent == "implicit" else None,
        "rgbmax_A": args.rgbmax if args.solvent == "implicit" else None,
        "nmr_restraints": os.path.basename(nmr) if nmr else None,
        "heat": {"ns": args.heat_ns, "nstlim": nstlim_heat,
                 "seed": stage_seed(0)},
        "release": {"n_stages": n_rel, "ns_each": args.release_ns,
                    "nstlim_each": nstlim_rel, "weights": weights,
                    "seeds": rel_seeds,
                    "files": [release_stem(i, n_rel) + ".in"
                              for i in range(n_rel)]},
        "production": {"chunks": args.prod_chunks,
                       "ns_each": args.prod_chunk_ns,
                       "nstlim_each": nstlim_prod,
                       "seeds": prod_seeds,
                       "files": [f"05_prod_c{c:02d}.in"
                                 for c in range(1, args.prod_chunks + 1)]},
        "bench": {"ns": args.bench_ns, "nstlim": nstlim_bench,
                  "seed": stage_seed(1 + n_rel + args.prod_chunks)},
        "trajectory_interval_ps": args.traj_interval_ps,
        "ntwx": ntwx,
        "ntpr": args.ntpr,
        "ntwr": args.ntwr,
        "hmr": {
            "enabled": bool(args.hmr),
            "hydrogen_mass": hmr_mass if args.hmr else None,
            "water_hydrogens_repartitioned": False,
            "topology": (os.path.basename(args.prmtop)
                         if args.hmr and args.prmtop else None),
        },
        "files": list(files),
    }
    with open(os.path.join(args.out, "run_config.json"), "w", encoding="utf-8") as fh:
        json.dump(config, fh, indent=2)
        fh.write("\n")

    print(table, end="")
    print(f"\nWrote {len(files)} Amber inputs to {os.path.abspath(args.out)}")
    print("Read README.md for the mandatory run order.")


if __name__ == "__main__":
    main()
