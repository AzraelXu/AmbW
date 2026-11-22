#!/usr/bin/env python3
"""Generate cpptraj inputs for an Amber trajectory: autoimage, strip, dry
NetCDF, optional RNA nastruct. It writes files by default and only runs
cpptraj with --run, so the user can inspect masks first.

The topology must match the input trajectory atom-for-atom (usually the
solvated topology and the raw production chunks). If the trajectory is already
dry, pass its dry topology and omit --strip-mask.

Usage:
  python3 make_cpptraj_inputs.py --top complex.prmtop \
    --traj 'md/r1/05_prod_c*.nc' --strip-mask ':WAT,K+,Cl-,Na+' \
    --anchor ':1-13' --fixed ':14-26,27' \
    --outdir analysis/cpptraj --prefix lig1_r1 [--nastruct --resrange 1-26] \
    [--run] [--dry-run] [--force]
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def natural_key(path):
    return [int(x) if x.isdigit() else x
            for x in re.split(r"(\d+)", path)]


def expand_trajs(patterns):
    files = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern), key=natural_key)
        if not matches and os.path.isfile(pattern):
            matches = [pattern]
        if not matches:
            die(f"no trajectory matched: {pattern}")
        files.extend(matches)
    seen = set()
    out = []
    for f in files:
        real = os.path.abspath(f)
        if real not in seen:
            seen.add(real)
            out.append(real)
    return out


def cpptraj_frames(top, traj):
    """Return frame count, or die if cpptraj cannot pair top/traj."""
    cmd = ["cpptraj", "-p", top, "-y", traj, "-tl"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True)
    text = proc.stdout
    m = re.search(r"Frames:\s*(\d+)", text)
    if proc.returncode != 0 or not m:
        tail = "\n".join(text.splitlines()[-12:])
        die(f"cpptraj cannot read {traj} with {top}:\n{tail}")
    return int(m.group(1))


def render(template, values):
    text = open(template, encoding="utf-8").read()
    for key, value in values.items():
        text = text.replace("{{" + key + "}}", str(value))
    if "{{" in text:
        die(f"unresolved placeholder in {template}: "
            f"{sorted({x.split('}}')[0] for x in text.split('{{')[1:]})}")
    return text


def write(path, text, force):
    if os.path.exists(path) and not force:
        die(f"refusing to overwrite {path} (use --force)")
    if not text.endswith("\n"):
        text += "\n"
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", required=True)
    ap.add_argument("--traj", action="append", required=True,
                    help="trajectory file or glob; repeatable")
    ap.add_argument("--strip-mask", help="Amber mask, e.g. ':WAT,K+,Cl-,Na+'")
    ap.add_argument("--anchor", help="autoimage anchor mask")
    ap.add_argument("--fixed", help="autoimage fixed mask")
    ap.add_argument("--no-autoimage", action="store_true")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--nastruct", action="store_true")
    ap.add_argument("--resrange", help="nastruct residue range, e.g. 1-26")
    ap.add_argument("--stride", type=int, default=10,
                    help="nastruct frame stride (default 10)")
    ap.add_argument("--rmsd-mask", help="optional RMSD mask for the prep trajectory")
    ap.add_argument("--run", action="store_true", help="run cpptraj after writing")
    ap.add_argument("--dry-run", action="store_true", help="print inputs, write nothing")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not os.path.isfile(args.top):
        die(f"topology not found: {args.top}")
    if not re.match(r"^[A-Za-z0-9._-]+$", args.prefix):
        die(f"unsafe --prefix: {args.prefix}")
    if args.no_autoimage and (args.anchor or args.fixed):
        die("--no-autoimage cannot be combined with --anchor/--fixed")
    if args.stride < 1:
        die("--stride must be >= 1")
    if args.nastruct and not args.resrange:
        print("warning: --nastruct without --resrange may analyze ligand/ions; "
              "confirm the default", file=sys.stderr)

    trajs = expand_trajs(args.traj)
    for f in trajs:
        if not os.path.isfile(f) or os.path.getsize(f) == 0:
            die(f"trajectory missing or empty: {f}")

    frames = [cpptraj_frames(args.top, f) for f in trajs]
    total_frames = sum(frames)

    outdir = os.path.abspath(args.outdir)
    if not args.dry_run:
        os.makedirs(outdir, exist_ok=True)
    top_abs = os.path.abspath(args.top)
    dry_top = f"{args.prefix}_dry.prmtop"
    dry_nc = f"{args.prefix}_dry.nc"
    prep_path = os.path.join(outdir, f"{args.prefix}_prep.in")

    trajin_lines = "\n".join(f"trajin {f}" for f in trajs)
    autoimage_line = ""
    if not args.no_autoimage:
        if args.anchor and args.fixed:
            autoimage_line = f"autoimage anchor '{args.anchor}' fixed '{args.fixed}'"
        elif args.anchor:
            autoimage_line = f"autoimage anchor '{args.anchor}'"
        elif args.fixed:
            autoimage_line = f"autoimage fixed '{args.fixed}'"
        else:
            autoimage_line = "autoimage"

    strip_line = ""
    if args.strip_mask:
        strip_line = f"strip '{args.strip_mask}' parmout {dry_top}"
        if not os.path.isdir(outdir):
            pass
        # If the prefix contains a directory, cpptraj runs in outdir; the
        # parmout/trajout names are deliberately bare filenames.
    else:
        strip_line = f"parmwrite out {dry_top}"

    rms_line = ""
    if args.rmsd_mask:
        rms_line = f"rms '{args.rmsd_mask}' first mass out {args.prefix}_rmsd.dat"

    prep = "\n".join(x for x in [
        f"parm {top_abs}",
        trajin_lines,
        autoimage_line,
        strip_line,
        rms_line,
        f"trajout {dry_nc} netcdf",
        "run",
        "quit",
        "",
    ] if x != "")

    nstruct = None
    if args.nastruct:
        nstruct_path = os.path.join(outdir, f"{args.prefix}_nastruct.in")
        resrange = f"resrange {args.resrange}" if args.resrange else ""
        nstruct = "\n".join(x for x in [
            f"parm {dry_top}",
            f"trajin {dry_nc} 1 last {args.stride}",
            f"nastruct naout {args.prefix}_na {resrange} "
            "sscalc calcnohb groovecalc 3dna",
            "run",
            "quit",
            "",
        ] if x != "")

    summary = (f"topology: {top_abs}\n"
               f"chunks:   {len(trajs)} ({', '.join(str(x) for x in frames)} frames; "
               f"total {total_frames})\n"
               f"strip:    {args.strip_mask or 'none'}\n"
               f"autoimage:{' off' if args.no_autoimage else ' on'}\n"
               f"output:   {dry_top}, {dry_nc}\n")

    if args.dry_run:
        print(summary)
        print(f"===== {prep_path} =====")
        print(prep)
        if nstruct:
            print(f"===== {args.prefix}_nastruct.in =====")
            print(nstruct)
        return

    write(prep_path, prep, args.force)
    if nstruct:
        write(os.path.join(outdir, f"{args.prefix}_nastruct.in"), nstruct, args.force)

    print(summary)
    print(f"wrote {prep_path}")
    if nstruct:
        print(f"wrote {os.path.join(outdir, f'{args.prefix}_nastruct.in')}")

    if not args.run:
        print("not running cpptraj; inspect the input, then rerun with --run "
              "or run: cpptraj -i " + prep_path)
        return

    if not shutil.which("cpptraj"):
        die("cpptraj not in PATH; source $AMBERHOME/amber.sh")
    log = os.path.join(outdir, f"{args.prefix}_prep.log")
    with open(log, "w", encoding="utf-8") as fh:
        proc = subprocess.run(["cpptraj", "-i", prep_path], cwd=outdir,
                              stdout=fh, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        die(f"cpptraj prep failed; see {log}")
    for needed in (dry_top, dry_nc):
        path = os.path.join(outdir, needed)
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            die(f"cpptraj did not produce {path}; see {log}")
    print(f"cpptraj prep OK: {os.path.join(outdir, dry_nc)}")

    if nstruct:
        nlog = os.path.join(outdir, f"{args.prefix}_nastruct.log")
        with open(nlog, "w", encoding="utf-8") as fh:
            proc = subprocess.run(["cpptraj", "-i",
                                   os.path.join(outdir, f"{args.prefix}_nastruct.in")],
                                  cwd=outdir, stdout=fh, stderr=subprocess.STDOUT)
        if proc.returncode != 0:
            die(f"cpptraj nastruct failed; see {nlog}")
        print(f"cpptraj nastruct OK: {os.path.join(outdir, args.prefix + '_na')}")


if __name__ == "__main__":
    main()
