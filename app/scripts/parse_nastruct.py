#!/usr/bin/env python3
"""Summarize cpptraj nastruct BP/BPstep output.

The cpptraj files are named:
  BP.<prefix>_na
  BPstep.<prefix>_na
(e.g. prefix = lig1_r1 for analysis/cpptraj/lig1_r1_na).

Outputs:
  - a per-prefix summary CSV (mean buckle/propeller/opening/rise/twist/...)
  - an optional per-base-pair CSV with occupancy

Usage:
  python -I parse_nastruct.py --prefix analysis/cpptraj/lig1_r1 \
      --label lig1 --rep 1 --out analysis/tables/nastruct_summary.csv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

BP_COLS = "Frame Base1 Base2 Shear Stretch Stagger Buckle Propeller Opening BP HB".split()
STEP_COLS = "Frame BP1 BP2 Shift Slide Rise Tilt Roll Twist Zp Major Minor".split()


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def read_table(path, names):
    if not os.path.isfile(path):
        die(f"file not found: {path}")
    df = pd.read_csv(path, sep=r"\s+", comment="#", names=names)
    df = df.replace("----", np.nan)
    for col in df.columns:
        if col not in ("Frame", "Base1", "Base2", "BP1", "BP2"):
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prefix", required=True,
                    help="prefix before _na, e.g. analysis/cpptraj/lig1_r1")
    ap.add_argument("--label", default="")
    ap.add_argument("--rep", default="")
    ap.add_argument("--out", default="nastruct_summary.csv")
    ap.add_argument("--pairs-out", default="")
    ap.add_argument("--last-fraction", type=float, default=1.0,
                    help="use the last fraction of frames (default 1.0 = all)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not (0.0 < args.last_fraction <= 1.0):
        die("--last-fraction must be in (0,1]")
    # cpptraj names BP.<prefix>_na, BPstep.<prefix>_na
    bp_path = os.path.join(os.path.dirname(args.prefix),
                           "BP." + os.path.basename(args.prefix) + "_na")
    st_path = os.path.join(os.path.dirname(args.prefix),
                           "BPstep." + os.path.basename(args.prefix) + "_na")
    bp = read_table(bp_path, BP_COLS)
    st = read_table(st_path, STEP_COLS)

    if args.last_fraction < 1.0:
        def tail(df):
            fmax = df.Frame.max()
            fmin = df.Frame.min()
            cut = fmax - args.last_fraction * (fmax - fmin + 1)
            return df[df.Frame >= cut]
        bp, st = tail(bp), tail(st)

    row = {
        "label": args.label,
        "replica": args.rep,
        "prefix": args.prefix,
        "n_frames": int(bp.Frame.nunique()),
        "buckle": bp.Buckle.mean(),
        "propeller": bp.Propeller.mean(),
        "opening": bp.Opening.mean(),
        "shear": bp.Shear.mean(),
        "stretch": bp.Stretch.mean(),
        "stagger": bp.Stagger.mean(),
        "rise": st.Rise.mean(),
        "twist": st.Twist.mean(),
        "roll": st.Roll.mean(),
        "tilt": st.Tilt.mean(),
        "major_groove": pd.to_numeric(st.Major, errors="coerce").mean(),
        "minor_groove": pd.to_numeric(st.Minor, errors="coerce").mean(),
    }

    # Per-pair occupancy: BP flag is 1 for a paired frame, 0 otherwise.
    pairs = []
    for (b1, b2), group in bp.groupby(["Base1", "Base2"], sort=False):
        bp_flag = pd.to_numeric(group.BP, errors="coerce")
        pairs.append({
            "label": args.label,
            "replica": args.rep,
            "base1": b1,
            "base2": b2,
            "n_frames": int(len(group)),
            "bp_occupancy": float(bp_flag.mean()) if len(group) else np.nan,
            "hb_mean": float(pd.to_numeric(group.HB, errors="coerce").mean()),
            "buckle": group.Buckle.mean(),
            "propeller": group.Propeller.mean(),
            "opening": group.Opening.mean(),
        })
    pairs_df = pd.DataFrame(pairs)

    if os.path.exists(args.out) and not args.force:
        die(f"refusing to overwrite {args.out} (use --force)")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    pd.DataFrame([row]).to_csv(args.out, index=False, float_format="%.4f")
    print(f"wrote {args.out}")
    print(pd.DataFrame([row]).to_string(index=False))

    if args.pairs_out:
        if os.path.exists(args.pairs_out) and not args.force:
            die(f"refusing to overwrite {args.pairs_out} (use --force)")
        pairs_df.to_csv(args.pairs_out, index=False, float_format="%.4f")
        print(f"wrote {args.pairs_out}")


if __name__ == "__main__":
    main()
