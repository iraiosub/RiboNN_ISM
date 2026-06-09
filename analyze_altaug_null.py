#!/usr/bin/env python3
"""
analyze_altaug_null.py – Position-matched substitution null for altAUG positions.

Uses the SNV saturation scan from a RiboNN ISM run to test whether ΔTE at the
specified altAUG positions (-8, -7, -6 by default) are outliers relative to the
distribution of ΔTE across all other scanned positions.

Method (hybrid null, standard in deep-learning ISM work):
    background  =  ΔTE of every SNV at every position NOT in the altAUG set
    signal      =  ΔTE of every SNV at the altAUG positions
    statistics  =  z-score from background mean / std
                   empirical two-tailed p-value  (fraction of |background ΔTE| ≥ |signal ΔTE|)
                   directional empirical p-values (upper / lower tail)

Output:
    <outdir>/altaug_null_summary.tsv  – per-SNV stats table
    <outdir>/altaug_null_plot.png     – background distribution + signal overlay

Usage:
    python analyze_altaug_null.py \\
        [--input results/human/prediction_output.txt] \\
        [--altaug-positions -8 -7 -6] \\
        [--te-col mean_predicted_TE] \\
        [--outdir plots_ism_scn2a] \\
        [--gene-name SCN2A] [--species human]
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sp_stats


# ── helpers (mirrors plot_te_changes.py without importing it) ─────────────────

def load_results(path):
    df = pd.read_csv(path, sep="\t")
    if "tx_id" not in df.columns:
        df = pd.read_csv(path, sep="\t", index_col=0)
        df.index.name = "tx_id"
        df = df.reset_index()

    if "mean_predicted_TE" not in df.columns:
        predicted_cols = [c for c in df.columns if c.startswith("predicted_")]
        if not predicted_cols:
            raise ValueError(
                "Prediction output is missing mean_predicted_TE "
                "and has no predicted_* columns to aggregate."
            )
        group_cols = [
            c for c in ("tx_id", "utr5_sequence", "cds_sequence", "utr3_sequence")
            if c in df.columns
        ]
        df = df.groupby(group_cols, as_index=False)[predicted_cols].mean()
        df["mean_predicted_TE"] = df[predicted_cols].mean(axis=1)
    return df


def snv_label_to_coords(label):
    """Parse '-3_A>C' → (offset=-3, ref='A', alt='C') or None."""
    m = re.match(r"^([+-]\d+)_([ACGT])>([ACGT])$", label)
    if not m:
        return None
    return int(m.group(1)), m.group(2), m.group(3)


# ── core analysis ─────────────────────────────────────────────────────────────

def build_snv_frame(df, te_col, ref_te):
    """Extract all SNV rows, parse labels, compute ΔTE."""
    rows = []
    for _, row in df.iterrows():
        parsed = snv_label_to_coords(row["tx_id"])
        if parsed is None:
            continue
        offset, ref_base, alt = parsed
        rows.append({
            "tx_id":    row["tx_id"],
            "offset":   offset,
            "ref_base": ref_base,
            "alt_base": alt,
            "te":       row[te_col],
            "delta_te": row[te_col] - ref_te,
        })
    return pd.DataFrame(rows)


def annotate_uatg(signal_df):
    """
    For each SNV, check whether the substitution creates (or destroys) an ATG
    at the altAUG triplet formed by the three scanned positions.

    Requires exactly 3 consecutive positions (e.g. -8, -7, -6).  If the
    positions don't form a single codon the column is set to "n/a".
    """
    positions = sorted(signal_df["offset"].unique())
    if len(positions) != 3 or positions[2] - positions[0] != 2:
        signal_df = signal_df.copy()
        signal_df["creates_uATG"]  = "n/a"
        signal_df["destroys_uATG"] = "n/a"
        return signal_df

    o1, o2, o3 = positions

    # reference base at each position (from any row's ref_base field)
    ref_at = {
        o: signal_df.loc[signal_df["offset"] == o, "ref_base"].iloc[0]
        for o in positions
    }
    ref_triplet = ref_at[o1] + ref_at[o2] + ref_at[o3]

    def triplet_after(row):
        bases = dict(ref_at)
        bases[row["offset"]] = row["alt_base"]
        return bases[o1] + bases[o2] + bases[o3]

    signal_df = signal_df.copy()
    signal_df["triplet_after"]  = signal_df.apply(triplet_after, axis=1)
    signal_df["creates_uATG"]   = signal_df["triplet_after"] == "ATG"
    signal_df["destroys_uATG"]  = (ref_triplet == "ATG") & (signal_df["triplet_after"] != "ATG")
    signal_df["ref_triplet"]    = ref_triplet
    return signal_df


def position_null_test(snv_df, altaug_positions):
    """
    Partition SNVs into background (non-altAUG offsets) and signal (altAUG offsets).
    Return (background_df, signal_df, bg_stats_dict).
    """
    altaug_set  = set(altaug_positions)
    is_signal   = snv_df["offset"].isin(altaug_set)
    background  = snv_df[~is_signal].copy()
    signal      = snv_df[is_signal].copy()

    bg_vals = background["delta_te"].values
    bg_mean = float(bg_vals.mean())
    bg_std  = float(bg_vals.std(ddof=1))
    bg_abs  = np.abs(bg_vals)

    # z-score relative to background
    signal = signal.copy()
    signal["z_score"] = (signal["delta_te"] - bg_mean) / bg_std

    # empirical p-values
    signal["emp_p_two_tailed"] = signal["delta_te"].apply(
        lambda x: float((bg_abs >= abs(x)).mean())
    )
    signal["emp_p_upper"] = signal["delta_te"].apply(
        lambda x: float((bg_vals >= x).mean())
    )
    signal["emp_p_lower"] = signal["delta_te"].apply(
        lambda x: float((bg_vals <= x).mean())
    )
    signal["pct_rank"] = signal["delta_te"].apply(
        lambda x: float((bg_vals <= x).mean()) * 100.0
    )

    bg_stats = {
        "n":    int(len(bg_vals)),
        "mean": bg_mean,
        "std":  bg_std,
        "min":  float(bg_vals.min()),
        "max":  float(bg_vals.max()),
    }
    return background, signal, bg_stats


# ── plotting ──────────────────────────────────────────────────────────────────

def title_label(gene_name, species=None):
    if species:
        return f"{species} {gene_name}"
    return gene_name


def plot_null_distribution(background, signal, bg_stats, te_col, altaug_positions, outdir, label):
    bg_vals   = background["delta_te"].values
    positions = sorted(signal["offset"].unique())
    palette   = plt.cm.tab10.colors
    pos_color = {p: palette[i % 10] for i, p in enumerate(positions)}

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5),
                             gridspec_kw={"width_ratios": [2.2, 1]})

    # ── Panel 1: background histogram + KDE + signal rug ─────────────────────
    ax = axes[0]
    ax.hist(bg_vals, bins=40, color="#b0c4de", alpha=0.75,
            edgecolor="white", density=True, label=f"background ({bg_stats['n']} SNVs)")

    xs = np.linspace(bg_vals.min() - 0.06, bg_vals.max() + 0.06, 500)
    try:
        kde = sp_stats.gaussian_kde(bg_vals, bw_method="scott")
        ax.plot(xs, kde(xs), color="#2166ac", lw=1.8, zorder=3)
    except Exception:
        pass

    # ±2 σ shading
    lo = bg_stats["mean"] - 2 * bg_stats["std"]
    hi = bg_stats["mean"] + 2 * bg_stats["std"]
    ax.axvspan(lo, hi, color="grey", alpha=0.07, label="±2 σ background")
    ax.axvline(bg_stats["mean"], color="black", lw=0.9, ls="--", alpha=0.55,
               label=f"background mean ({bg_stats['mean']:+.4f})")

    # altAUG SNVs as vertical lines, grouped by position
    for pos in positions:
        subset = signal[signal["offset"] == pos]
        for _, row in subset.iterrows():
            ax.axvline(
                row["delta_te"], color=pos_color[pos], lw=1.4, alpha=0.9,
                label=f"{row['tx_id']}  z={row['z_score']:+.2f}  p={row['emp_p_two_tailed']:.3f}",
            )

    ax.set_xlabel(f"ΔTE ({te_col})", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title(
        f"Background SNV ΔTE distribution  (n={bg_stats['n']})\n"
        f"μ = {bg_stats['mean']:+.4f},  σ = {bg_stats['std']:.4f},  "
        f"range [{bg_stats['min']:+.4f}, {bg_stats['max']:+.4f}]",
        fontsize=9,
    )
    ax.legend(fontsize=6.5, loc="upper left", framealpha=0.8, ncol=1)

    # ── Panel 2: strip plot background vs altAUG positions ───────────────────
    ax2 = axes[1]
    rng = np.random.default_rng(42)

    # background
    jb = rng.uniform(-0.18, 0.18, len(bg_vals))
    ax2.scatter(jb, bg_vals, color="#cccccc", s=7, alpha=0.45, zorder=1, label="background")

    # signal per position
    x_offset = 0.65
    xticks, xlabels = [0.0], ["background"]
    for i, pos in enumerate(positions):
        xc = x_offset + i * 0.55
        xticks.append(xc)
        xlabels.append(f"offset {pos:+d}")
        subset = signal[signal["offset"] == pos].copy()
        js = rng.uniform(-0.12, 0.12, len(subset))
        ax2.scatter(
            js + xc, subset["delta_te"],
            color=pos_color[pos], s=65, zorder=3,
            edgecolors="black", linewidths=0.5,
            label=f"offset {pos:+d}",
        )
        for j, (_, row) in enumerate(subset.iterrows()):
            creates_label = ""
            if row.get("creates_uATG") is True:
                creates_label = " ★"
            ax2.annotate(
                row["alt_base"] + creates_label,
                xy=(js[j] + xc, row["delta_te"]),
                xytext=(0, 5),
                textcoords="offset points",
                fontsize=7, ha="center", color=pos_color[pos],
            )

    ax2.axhline(bg_stats["mean"], color="black", lw=0.8, ls="--", alpha=0.5)
    ax2.axhline(hi, color="grey", lw=0.5, ls=":", alpha=0.6)
    ax2.axhline(lo, color="grey", lw=0.5, ls=":", alpha=0.6)
    ax2.text(xticks[-1] + 0.3, hi, "+2σ", fontsize=6, va="bottom", color="grey")
    ax2.text(xticks[-1] + 0.3, lo, "−2σ", fontsize=6, va="top",    color="grey")

    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xlabels, rotation=35, ha="right", fontsize=8)
    ax2.set_ylabel(f"ΔTE ({te_col})", fontsize=10)
    ax2.set_title(
        f"altAUG positions\n({', '.join(str(p) for p in altaug_positions)})\n"
        f"★ = creates uATG",
        fontsize=9,
    )

    fig.suptitle(
        f"{label} 5′UTR ISM: position-matched substitution null for altAUG positions",
        fontsize=11, y=1.01,
    )
    fig.tight_layout()
    path = outdir / "altaug_null_plot.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", default="results/human/prediction_output.txt",
        help="RiboNN prediction output TSV (default: results/human/prediction_output.txt)",
    )
    parser.add_argument(
        "--altaug-positions", nargs="+", type=int, default=[-8, -7, -6],
        metavar="N",
        help="Offsets from canonical AUG considered the altAUG signal (default: -8 -7 -6)",
    )
    parser.add_argument(
        "--te-col", default="mean_predicted_TE",
        help="TE column to use for ΔTE (default: mean_predicted_TE)",
    )
    parser.add_argument(
        "--outdir", default="plots_ism_scn2a",
        help="Output directory (default: plots_ism_scn2a)",
    )
    parser.add_argument(
        "--gene-name", default="SCN2A",
        help="Gene label for plot titles (default: SCN2A)",
    )
    parser.add_argument(
        "--species", choices=("human", "mouse"), default=None,
        help="Species label for plot titles",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading predictions from {args.input} ...")
    df = load_results(args.input)
    print(f"  {len(df)} rows loaded.")

    ref_rows = df[df["tx_id"] == "reference"]
    if ref_rows.empty:
        raise ValueError("No 'reference' row found in prediction output.")
    if args.te_col not in df.columns:
        raise ValueError(
            f"Column '{args.te_col}' not in prediction output. "
            f"Available columns: {list(df.columns)}"
        )
    ref_te = float(ref_rows.iloc[0][args.te_col])
    print(f"  Reference {args.te_col} = {ref_te:.4f}")

    print("\nBuilding SNV ΔTE frame ...")
    snv_df = build_snv_frame(df, args.te_col, ref_te)
    print(f"  {len(snv_df)} SNV variants parsed.")

    print(f"\nRunning position-matched null test ...")
    print(f"  altAUG positions : {args.altaug_positions}")
    background, signal, bg_stats = position_null_test(snv_df, args.altaug_positions)

    if len(signal) == 0:
        print(
            "[WARN] No SNVs found at the specified altAUG positions. "
            "Check --altaug-positions against your upstream-bases scan range."
        )
        return

    signal = annotate_uatg(signal)

    print(f"  Background       : {bg_stats['n']} SNVs")
    print(f"    mean  = {bg_stats['mean']:+.4f}")
    print(f"    std   = {bg_stats['std']:.4f}")
    print(f"    range = [{bg_stats['min']:+.4f}, {bg_stats['max']:+.4f}]")
    print(f"  Signal           : {len(signal)} SNVs at altAUG positions")

    # ── console summary table ─────────────────────────────────────────────────
    col_w = {"tx_id": 22, "offset": 8, "delta_te": 11, "z": 9,
             "p2t": 10, "pct": 8, "uATG": 14}
    hdr = (
        f"  {'tx_id':<{col_w['tx_id']}} {'offset':>{col_w['offset']}} "
        f"{'ΔTE':>{col_w['delta_te']}} {'z_score':>{col_w['z']}} "
        f"{'emp_p_2t':>{col_w['p2t']}} {'pct_rank':>{col_w['pct']}} "
        f"{'creates_uATG':>{col_w['uATG']}}"
    )
    print(f"\naltAUG SNV summary (sorted by offset, then alt base):")
    print(hdr)
    print("  " + "-" * (sum(col_w.values()) + 14))
    for _, row in signal.sort_values(["offset", "alt_base"]).iterrows():
        creates = str(row.get("creates_uATG", "?"))
        print(
            f"  {row['tx_id']:<{col_w['tx_id']}} {row['offset']:>+{col_w['offset']}} "
            f"{row['delta_te']:>+{col_w['delta_te']}.4f} "
            f"{row['z_score']:>+{col_w['z']}.2f} "
            f"{row['emp_p_two_tailed']:>{col_w['p2t']}.4f} "
            f"{row['pct_rank']:>{col_w['pct']}.1f} "
            f"{creates:>{col_w['uATG']}}"
        )

    # per-position aggregate (max |z|, min emp_p within each position)
    print("\nPer-position aggregate (all alt bases at that offset):")
    agg = (
        signal
        .groupby("offset")
        .agg(
            n_snvs=("tx_id", "count"),
            max_abs_z=("z_score", lambda s: s.abs().max()),
            min_emp_p=("emp_p_two_tailed", "min"),
            delta_te_range_min=("delta_te", "min"),
            delta_te_range_max=("delta_te", "max"),
        )
        .reset_index()
        .sort_values("offset")
    )
    print(f"  {'offset':>8} {'n_snvs':>8} {'max|z|':>9} {'min_emp_p':>11} "
          f"{'ΔTE_min':>10} {'ΔTE_max':>10}")
    print("  " + "-" * 60)
    for _, row in agg.iterrows():
        print(
            f"  {row['offset']:>+8} {row['n_snvs']:>8} "
            f"{row['max_abs_z']:>9.2f} {row['min_emp_p']:>11.4f} "
            f"{row['delta_te_range_min']:>+10.4f} {row['delta_te_range_max']:>+10.4f}"
        )

    # ── save outputs ──────────────────────────────────────────────────────────
    tsv_path = outdir / "altaug_null_summary.tsv"
    signal.to_csv(tsv_path, sep="\t", index=False, float_format="%.6f")
    print(f"\nSummary TSV → {tsv_path}")

    print("Generating plot ...")
    label = title_label(args.gene_name, args.species)
    plot_null_distribution(background, signal, bg_stats, args.te_col,
                           args.altaug_positions, outdir, label)

    print(f"\nDone. Outputs in {outdir}/")


if __name__ == "__main__":
    main()
