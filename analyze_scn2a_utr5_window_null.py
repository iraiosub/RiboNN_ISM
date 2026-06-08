#!/usr/bin/env python3
"""
analyze_scn2a_utr5_window_null.py

Build a deletion-window null distribution for the SCN2A 5'UTR.

Default comparison:
    target = the final 15-nt window immediately before the canonical AUG
             offsets -15..-1
    null   = every other 15-nt sliding-window deletion across the 5'UTR
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


WINDOW_RE = re.compile(
    r"^delwin(?P<size>\d+)_pos(?P<start>\d+)_(?P<end>\d+)"
    r"_off(?P<start_offset>[+-]\d+)_(?P<end_offset>[+-]\d+)$"
)


def load_results(path):
    df = pd.read_csv(path, sep="\t")
    if "tx_id" not in df.columns:
        df = pd.read_csv(path, sep="\t", index_col=0)
        df.index.name = "tx_id"
        df = df.reset_index()

    if "mean_predicted_TE" not in df.columns:
        predicted_cols = [col for col in df.columns if col.startswith("predicted_")]
        if not predicted_cols:
            raise ValueError(
                "Prediction output is missing mean_predicted_TE and has no predicted_* columns."
            )
        group_cols = [
            col
            for col in ("tx_id", "utr5_sequence", "cds_sequence", "utr3_sequence")
            if col in df.columns
        ]
        df = df.groupby(group_cols, as_index=False)[predicted_cols].mean()
        df["mean_predicted_TE"] = df[predicted_cols].mean(axis=1)
    return df


def parse_window_id(tx_id):
    match = WINDOW_RE.match(tx_id)
    if not match:
        return None
    parsed = {key: int(value) for key, value in match.groupdict().items()}
    return {
        "window_size": parsed["size"],
        "start_utr5_pos": parsed["start"],
        "end_utr5_pos": parsed["end"],
        "start_offset": parsed["start_offset"],
        "end_offset": parsed["end_offset"],
    }


def build_window_frame(df, te_col, ref_te):
    rows = []
    for _, row in df.iterrows():
        parsed = parse_window_id(row["tx_id"])
        if parsed is None:
            continue
        parsed["tx_id"] = row["tx_id"]
        parsed["te"] = float(row[te_col])
        parsed["delta_TE"] = float(row[te_col]) - ref_te
        rows.append(parsed)
    return pd.DataFrame(rows)


def centered_empirical_p(background_values, target_value):
    bg_mean = float(background_values.mean())
    centered_bg = np.abs(background_values - bg_mean)
    centered_target = abs(target_value - bg_mean)
    return float((centered_bg >= centered_target).mean())


def plot_null(windows, background, target, bg_stats, outdir, window_size):
    bg_vals = background["delta_TE"].values
    target_delta = float(target["delta_TE"])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    bins = min(50, max(10, int(np.sqrt(len(bg_vals)))))
    ax.hist(bg_vals, bins=bins, color="#b7c9d9", edgecolor="white", alpha=0.85)
    ax.axvline(bg_stats["mean"], color="black", lw=1.0, ls="--", label="null mean")
    ax.axvline(target_delta, color="#c43c39", lw=2.0, label="target -15..-1")
    ax.set_xlabel("Delta mean predicted TE")
    ax.set_ylabel("Window count")
    ax.set_title(f"Null distribution from other {window_size}-nt deletions")
    ax.legend(frameon=False, fontsize=8)

    ax2 = axes[1]
    ordered = windows.sort_values("start_utr5_pos")
    ax2.plot(
        ordered["start_utr5_pos"],
        ordered["delta_TE"],
        color="#446c8a",
        lw=1.0,
        marker="o",
        ms=2.5,
        alpha=0.75,
    )
    ax2.scatter(
        [target["start_utr5_pos"]],
        [target_delta],
        color="#c43c39",
        edgecolor="black",
        linewidth=0.5,
        s=70,
        zorder=3,
    )
    ax2.axhline(bg_stats["mean"], color="black", lw=0.8, ls="--", alpha=0.6)
    ax2.set_xlabel("Deleted window start position in 5'UTR (1-based)")
    ax2.set_ylabel("Delta mean predicted TE")
    ax2.set_title("Sliding-window deletion scan")

    fig.tight_layout()
    path = outdir / f"window{window_size}_deletion_null_plot.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved {path}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="results/human/scn2a_utr5_delwin15_prediction_output.txt",
        help="RiboNN prediction output TSV",
    )
    parser.add_argument("--te-col", default="mean_predicted_TE")
    parser.add_argument("--window-size", type=int, default=15)
    parser.add_argument("--target-start-offset", type=int, default=-15)
    parser.add_argument("--target-end-offset", type=int, default=-1)
    parser.add_argument("--outdir", default="plots_ism_scn2a_window15")
    return parser.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading predictions from {args.input} ...")
    df = load_results(args.input)
    if args.te_col not in df.columns:
        raise ValueError(f"Column not found: {args.te_col}")

    ref_rows = df[df["tx_id"] == "reference"]
    if ref_rows.empty:
        raise ValueError("No reference row found in prediction output.")
    ref_te = float(ref_rows.iloc[0][args.te_col])
    print(f"Reference {args.te_col}: {ref_te:.6f}")

    windows = build_window_frame(df, args.te_col, ref_te)
    windows = windows[windows["window_size"] == args.window_size].copy()
    if windows.empty:
        raise ValueError(f"No delwin{args.window_size} rows found in prediction output.")

    is_target = (
        (windows["start_offset"] == args.target_start_offset)
        & (windows["end_offset"] == args.target_end_offset)
    )
    if not is_target.any():
        examples = windows[["tx_id", "start_offset", "end_offset"]].tail(5)
        raise ValueError(
            "Target window not found. Last parsed windows were:\n"
            f"{examples.to_string(index=False)}"
        )
    if is_target.sum() > 1:
        raise ValueError("Multiple target windows matched; check window labels.")

    windows["is_target"] = is_target
    target = windows.loc[is_target].iloc[0]
    background = windows.loc[~is_target].copy()
    if background.empty:
        raise ValueError("No background windows remain after excluding the target.")

    bg_vals = background["delta_TE"].values
    bg_mean = float(bg_vals.mean())
    bg_std = float(bg_vals.std(ddof=1)) if len(bg_vals) > 1 else float("nan")
    target_delta = float(target["delta_TE"])
    target_z = (target_delta - bg_mean) / bg_std if bg_std > 0 else float("nan")
    emp_p_two_tailed = centered_empirical_p(bg_vals, target_delta)
    emp_p_upper = float((bg_vals >= target_delta).mean())
    emp_p_lower = float((bg_vals <= target_delta).mean())
    percentile = float((bg_vals <= target_delta).mean() * 100.0)

    windows["z_vs_background"] = (
        (windows["delta_TE"] - bg_mean) / bg_std if bg_std > 0 else np.nan
    )
    all_windows_path = outdir / f"window{args.window_size}_deletion_all_windows.tsv"
    windows.sort_values("start_utr5_pos").to_csv(
        all_windows_path,
        sep="\t",
        index=False,
        float_format="%.6f",
    )

    summary = pd.DataFrame(
        [
            {
                "target_tx_id": target["tx_id"],
                "window_size": args.window_size,
                "target_start_offset": args.target_start_offset,
                "target_end_offset": args.target_end_offset,
                "target_delta_TE": target_delta,
                "background_n": len(background),
                "background_mean_delta_TE": bg_mean,
                "background_std_delta_TE": bg_std,
                "target_z_score": target_z,
                "emp_p_two_tailed_centered": emp_p_two_tailed,
                "emp_p_upper": emp_p_upper,
                "emp_p_lower": emp_p_lower,
                "target_percentile": percentile,
            }
        ]
    )
    summary_path = outdir / f"window{args.window_size}_deletion_null_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False, float_format="%.6f")

    print("\nTarget window comparison:")
    print(f"  target tx_id       : {target['tx_id']}")
    print(f"  target delta TE    : {target_delta:+.6f}")
    print(f"  background windows : {len(background)}")
    print(f"  background mean    : {bg_mean:+.6f}")
    print(f"  background std     : {bg_std:.6f}")
    print(f"  z score            : {target_z:+.3f}")
    print(f"  empirical p        : {emp_p_two_tailed:.4f}")
    print(f"  percentile         : {percentile:.1f}")
    print(f"\nAll windows TSV: {all_windows_path}")
    print(f"Summary TSV:     {summary_path}")

    bg_stats = {"mean": bg_mean, "std": bg_std}
    print("Generating plot ...")
    plot_null(windows, background, target, bg_stats, outdir, args.window_size)
    print(f"\nDone. Outputs in {outdir}/")


if __name__ == "__main__":
    main()
