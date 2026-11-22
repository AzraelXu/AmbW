#!/usr/bin/env python3
"""Collect per-replica MM/GBSA DeltaG and rank ligands.

Reads analysis/mmpbsa/<label>_r<rep>/mmpbsa.dat and config/systems.tsv.
Each replica is one sample; the system mean is the ranking quantity.
Bootstrap resamples replicas; the permutation test uses the calibration
systems with experimental values (control systems are excluded).

Output columns:
  label, control, experimental, experimental_type, n,
  dg_mean, dg_sd, dg_sem, dg_min, dg_max,
  rank_all (includes control), rank_calibration (experimental only),
  exp_rank, p_rank1_bootstrap

Usage:
  python -I rank_ligands.py --root . --systems config/systems.tsv \
    --out analysis/tables/ranking.csv
"""
import argparse
import os
import re
import sys
from itertools import permutations

import numpy as np
import pandas as pd


def die(msg):
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(2)


def parse_delta(path):
    text = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", text)
    if not m:
        raise ValueError(f"DELTA TOTAL not found in {path}")
    return float(m.group(1)), float(m.group(2)), float(m.group(3))


def spearman(x, y):
    return float(pd.Series(x).rank().corr(pd.Series(y).rank()))


def exact_perm_p(observed, x, y):
    vals = [abs(spearman(x, perm)) for perm in permutations(y)]
    return float(np.mean(np.array(vals) >= abs(observed) - 1e-12))


def read_systems(path):
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str,
                     keep_default_na=False)
    for col in ("label", "replicas", "control", "experimental",
                "experimental_type"):
        if col not in df.columns:
            df[col] = ""
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".")
    ap.add_argument("--systems", required=True)
    ap.add_argument("--out", default="analysis/tables/ranking.csv")
    ap.add_argument("--n-boot", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=20261022)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    systems = read_systems(args.systems)
    rows = []
    meta = {}
    for _, srow in systems.iterrows():
        label = str(srow["label"]).strip()
        if not label:
            continue
        reps = str(srow.get("replicas", "")).strip() or "1"
        control = str(srow.get("control", "")).strip()
        experimental = str(srow.get("experimental", "")).strip()
        exptype = str(srow.get("experimental_type", "")).strip() or "NA"
        meta[label] = {
            "control": 1 if control in ("1", "true", "TRUE", "yes", "Y") else 0,
            "experimental": np.nan if experimental in ("", "NA", "nan", "-") else float(experimental),
            "experimental_type": exptype,
            "replicas": [x.strip() for x in reps.split(",") if x.strip()],
        }
        for rep in meta[label]["replicas"]:
            path = os.path.join(root, "analysis", "mmpbsa",
                                f"{label}_r{rep}", "mmpbsa.dat")
            if not os.path.isfile(path):
                print(f"warning: missing {path}", file=sys.stderr)
                continue
            try:
                dg, sd, sem = parse_delta(path)
            except Exception as exc:
                print(f"warning: parse failed {path}: {exc}", file=sys.stderr)
                continue
            rows.append({"label": label, "replica": rep, "dg": dg,
                         "frame_sd": sd, "frame_sem": sem})

    if not rows:
        die("no mmpbsa.dat files found")

    df = pd.DataFrame(rows)
    stats = (df.groupby("label")["dg"]
               .agg(["mean", "std", "count", "min", "max"]).reset_index()
               .rename(columns={"mean": "dg_mean", "std": "dg_sd",
                                "count": "n", "min": "dg_min", "max": "dg_max"}))
    stats["dg_sem"] = stats["dg_sd"] / np.sqrt(stats["n"])
    stats["control"] = stats["label"].map(lambda x: meta.get(x, {}).get("control", 0))
    stats["experimental"] = stats["label"].map(lambda x: meta.get(x, {}).get("experimental", np.nan))
    stats["experimental_type"] = stats["label"].map(
        lambda x: meta.get(x, {}).get("experimental_type", "NA"))
    stats["rank_all"] = stats["dg_mean"].rank(method="min").astype(int)

    calibration = stats[(stats["control"] == 0) & stats["experimental"].notna()].copy()
    calibration["exp_rank"] = calibration["experimental"].rank(method="min").astype(int)
    calibration["rank_calibration"] = calibration["dg_mean"].rank(method="min").astype(int)
    stats = stats.merge(calibration[["label", "exp_rank", "rank_calibration"]],
                        on="label", how="left")
    stats = stats.sort_values("dg_mean").reset_index(drop=True)

    # Bootstrap over replicas (each label independently).
    rng = np.random.default_rng(args.seed)
    labels_cal = calibration["label"].tolist()
    dg_by_label = {lab: df[df.label == lab]["dg"].to_numpy() for lab in labels_cal}
    boot_rank1 = {lab: [] for lab in labels_cal}
    if labels_cal and args.n_boot > 0:
        for _ in range(args.n_boot):
            means = {}
            for lab in labels_cal:
                x = dg_by_label[lab]
                means[lab] = float(np.mean(rng.choice(x, size=len(x), replace=True)))
            order = sorted(labels_cal, key=lambda lab: means[lab])
            for rank, lab in enumerate(order, 1):
                boot_rank1[lab].append(rank)
    stats["p_rank1_bootstrap"] = stats["label"].map(
        {lab: (float(np.mean(np.array(boot_rank1[lab]) == 1)) if boot_rank1[lab] else np.nan)
         for lab in labels_cal})

    rho = pval = np.nan
    if len(calibration) >= 3:
        x = calibration["dg_mean"].to_numpy(dtype=float)
        y = calibration["experimental"].to_numpy(dtype=float)
        rho = spearman(x, y)
        pval = exact_perm_p(rho, x, y)

    cols = ["label", "control", "experimental", "experimental_type", "n",
            "dg_mean", "dg_sd", "dg_sem", "dg_min", "dg_max",
            "rank_all", "rank_calibration", "exp_rank", "p_rank1_bootstrap"]
    for col in cols:
        if col not in stats:
            stats[col] = np.nan
    stats = stats[cols]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    if os.path.exists(args.out) and not args.force:
        die(f"refusing to overwrite {args.out} (use --force)")
    stats.to_csv(args.out, index=False, float_format="%.4f")

    print("Per-replica MM/GBSA DeltaG (kcal/mol):")
    print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("\nSystem summary (sorted by DeltaG):")
    print(stats.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    if len(calibration) >= 3:
        print(f"\nSpearman rho (DeltaG vs experimental, calibration only): {rho:.3f}")
        print(f"Exact permutation p (|rho|, {len(calibration)} points): {pval:.3f}")
    else:
        print("\nNot enough calibration systems for a correlation test.")
    n_cal = len(calibration)
    print(f"\nCalibration systems: {n_cal}")
    if n_cal < 5:
        print("WARNING: fewer than 5 compounds; do not claim generalizability.")
    if n_cal >= 1:
        import math
        print(f"Perfect-rank probability by chance: 1/{math.factorial(n_cal)} "
              f"= {1.0 / math.factorial(n_cal):.3f}")
    print("MM/GBSA is a relative proxy, not a Kd. Do not tune parameters to match "
          "the experimental order.")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
