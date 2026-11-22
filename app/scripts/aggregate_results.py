#!/usr/bin/env python3
"""Aggregate per-replica structural metrics across replicas and systems.

Reads analysis/python/<label>_r<rep>.npz written by analyze_replica.py and
analysis/python/<label>_r<rep>_hb.npz written by analyze_ligand_hbonds.py.

Outputs:
  per_replica_metrics.csv   one row per replica; mean/SD/median over the window
  system_summary.csv        one row per label; mean/SD/SEM across replicas

Usage:
  python -I aggregate_results.py --root . --systems config/systems.tsv \
    --window last50 --out analysis/tables
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

METRICS = [
    "receptor_rmsd", "ligand_rmsd", "ligand_internal_rmsd",
    "receptor_rg", "ligand_rg", "ligand_com_drift",
    "min_dist", "n_contact_atoms", "n_contact_residues", "buried_frac",
    "loop_min_dist", "hbond_lig_to_rec", "hbond_rec_to_lig", "receptor_sasa",
]


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def read_systems(path):
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str,
                     keep_default_na=False)
    required = {"label", "replicas"}
    missing = required - set(df.columns)
    if missing:
        die(f"{path}: missing columns {sorted(missing)}")
    return df


def window_mask(time_ns, window):
    t = np.asarray(time_ns, dtype=float)
    if t.size == 0:
        return np.array([], dtype=bool)
    tmin, tmax = np.nanmin(t), np.nanmax(t)
    if window == "all":
        return np.ones(t.shape, dtype=bool)
    if window == "last50":
        return t >= tmin + 0.5 * (tmax - tmin)
    if window == "last25":
        return t >= tmin + 0.75 * (tmax - tmin)
    die(f"unknown --window {window!r}")


def summarize(arr, mask):
    x = np.asarray(arr, dtype=float)[mask]
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan, np.nan
    sd = float(np.std(x, ddof=1)) if x.size > 1 else 0.0
    return float(np.mean(x)), sd, float(np.median(x))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--systems", required=True)
    ap.add_argument("--window", choices=["last50", "last25", "all"], default="last50")
    ap.add_argument("--out", default="analysis/tables")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    systems = read_systems(args.systems)
    os.makedirs(args.out, exist_ok=True)

    rows = []
    for _, row in systems.iterrows():
        label = str(row["label"]).strip()
        reps = str(row.get("replicas", "1")).strip()
        if not reps or reps == "nan":
            reps = "1"
        for rep in [x.strip() for x in reps.split(",") if x.strip()]:
            prefix = f"{label}_r{rep}"
            npz_path = os.path.join(root, "analysis", "python", f"{prefix}.npz")
            if not os.path.isfile(npz_path):
                print(f"warning: missing {npz_path}", file=sys.stderr)
                continue
            d = np.load(npz_path, allow_pickle=True)
            mask = window_mask(d["time_ns"], args.window)
            out = {"label": label, "replica": rep, "n_frames": int(mask.sum())}
            for m in METRICS:
                if m not in d:
                    continue
                mean, sd, med = summarize(d[m], mask)
                out[f"{m}_mean"] = mean
                out[f"{m}_sd"] = sd
                out[f"{m}_median"] = med
            hb_path = os.path.join(root, "analysis", "python", f"{prefix}_hb.npz")
            if os.path.isfile(hb_path):
                h = np.load(hb_path, allow_pickle=True)
                for key in ("hbond_count", "hbond_count_rev"):
                    if key in h:
                        x = np.asarray(h[key], dtype=float)
                        x = x[np.isfinite(x)]
                        out[f"{key}_mean"] = float(x.mean()) if x.size else np.nan
            rows.append(out)

    if not rows:
        die("no per-replica npz files found")

    per = pd.DataFrame(rows).sort_values(["label", "replica"])
    per_path = os.path.join(args.out, "per_replica_metrics.csv")
    if os.path.exists(per_path) and not args.force:
        die(f"refusing to overwrite {per_path} (use --force)")
    per.to_csv(per_path, index=False, float_format="%.4f")

    # System summary: mean of replica means, SD/SEM across replicas.
    value_cols = [c for c in per.columns
                  if c.endswith("_mean") and c not in ("n_frames_mean",)]
    summary = []
    for label, group in per.groupby("label", sort=False):
        row = {"label": label, "n_replicas": int(len(group))}
        for col in value_cols:
            x = group[col].to_numpy(dtype=float)
            x = x[np.isfinite(x)]
            row[col.replace("_mean", "_system_mean")] = float(x.mean()) if x.size else np.nan
            row[col.replace("_mean", "_system_sd")] = (
                float(np.std(x, ddof=1)) if x.size > 1 else 0.0)
            row[col.replace("_mean", "_system_sem")] = (
                float(np.std(x, ddof=1) / np.sqrt(x.size)) if x.size > 1 else np.nan)
        summary.append(row)
    summ = pd.DataFrame(summary).sort_values("label")
    summ_path = os.path.join(args.out, "system_summary.csv")
    if os.path.exists(summ_path) and not args.force:
        die(f"refusing to overwrite {summ_path} (use --force)")
    summ.to_csv(summ_path, index=False, float_format="%.4f")

    print(per.to_string(index=False))
    print(f"\nwrote {per_path}\nwrote {summ_path}")


if __name__ == "__main__":
    main()
