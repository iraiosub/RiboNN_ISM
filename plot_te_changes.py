#!/usr/bin/env python3
"""
plot_te_changes.py – Visualise 5'UTR ISM results from RiboNN predictions.

Reads results/human/prediction_output.txt (or --input) and the reference row
from data/prediction_input.txt (or --variants-input), then produces:

  1. SNV heatmap          – positions × alt bases, coloured by ΔTE
  2. Single-del heatmap   – one row per deleted position
  3. Waterfall chart      – all variants ranked by ΔTE
  4. Growing-deletion trend – ΔTE vs deletion length
  5. Neuron panel         – neuronal cell types vs ΔTE for top variants

Output: plots_ism_scn2a/  (PNG files)

Usage:
    python plot_te_changes.py \\
        [--input results/human/prediction_output.txt] \\
        [--variants-input data/prediction_input.txt] \\
        [--outdir plots_ism_scn2a] \\
        [--upstream-bases 15] \\
        [--gene-name SCN2A] [--species human]
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns

# Neuronal / brain cell types present in RiboNN output. The prediction wrapper
# writes these as predicted_TE_* columns; raw TE_* names are accepted too.
HUMAN_NEURON_TYPES = [
    "TE_neurons",
    "TE_early_neurons",
    "TE_neuronal_precursor_cells",
    "TE_SH.SY5Y",
    "TE_normal_brain_tissue",
    "TE_human_brain_tumor",
]

MOUSE_NEURON_TYPES = [
    "TE_Cerebellum",
    "TE_DRG_neuronal_culture",
    "TE_Dorsal_section_of_lumbar_spinal_cord",
    "TE_ES_cell_derived_neurons",
    "TE_Fetal_cortex",
    "TE_Forebrain",
    "TE_NSC",
    "TE_Neuro2a",
    "TE_Neurons_(DIV_8)_derived_from_CGR8_ES_cells",
    "TE_Primary_cortical_neurons",
    "TE_Striatal_cells",
    "TE_brain",
    "TE_dentate_gyrus",
    "TE_hippocampal",
    "TE_mouse_eye",
    "TE_neural_tube",
]

BASES = ("A", "C", "G", "T")


# ── helpers ──────────────────────────────────────────────────────────────────

def load_results(path):
    df = pd.read_csv(path, sep="\t")
    if "tx_id" not in df.columns:
        # main.py output uses tx_id as first column but may be unnamed
        df = pd.read_csv(path, sep="\t", index_col=0)
        df.index.name = "tx_id"
        df = df.reset_index()

    if "mean_predicted_TE" not in df.columns:
        predicted_cols = [col for col in df.columns if col.startswith("predicted_")]
        if not predicted_cols:
            raise ValueError(
                "Prediction output is missing mean_predicted_TE and has no predicted_* columns to aggregate."
            )
        group_cols = [
            col
            for col in ("tx_id", "utr5_sequence", "cds_sequence", "utr3_sequence")
            if col in df.columns
        ]
        df = df.groupby(group_cols, as_index=False)[predicted_cols].mean()
        df["mean_predicted_TE"] = df[predicted_cols].mean(axis=1)
    return df


def get_ref_te(df):
    ref = df[df["tx_id"] == "reference"]
    if ref.empty:
        raise ValueError("No 'reference' row found in prediction output.")
    return ref.iloc[0]


def delta_te(df, ref_row, col="mean_predicted_TE"):
    ref_val = ref_row[col]
    return df[col] - ref_val


def title_label(gene_name, species=None):
    if species:
        return f"{species} {gene_name}"
    return gene_name


def _predicted_name(te_name):
    return f"predicted_{te_name}" if te_name.startswith("TE_") else te_name


def neuron_columns_for_species(species, columns):
    species = (species or "").lower()
    if species == "human":
        expected = HUMAN_NEURON_TYPES
    elif species == "mouse":
        expected = MOUSE_NEURON_TYPES
    else:
        expected = HUMAN_NEURON_TYPES + MOUSE_NEURON_TYPES

    seen = set()
    available = []
    for raw_name in expected:
        for candidate in (_predicted_name(raw_name), raw_name):
            if candidate in columns and candidate not in seen:
                available.append(candidate)
                seen.add(candidate)
    return available


def clean_te_label(column):
    return column.replace("predicted_TE_", "").replace("TE_", "")


def snv_label_to_coords(label):
    """Parse e.g. '-3_A>C' → (offset=-3, ref='A', alt='C')."""
    m = re.match(r"^([+-]\d+)_([ACGT])>([ACGT])$", label)
    if not m:
        return None
    return int(m.group(1)), m.group(2), m.group(3)


def del1_label_to_coords(label):
    """Parse 'del1_-3_A' → (offset=-3, ref='A')."""
    m = re.match(r"^del1_([+-]\d+)_([ACGT])$", label)
    if not m:
        return None
    return int(m.group(1)), m.group(2)


def deln_label_to_len(label):
    """Parse 'del_5_to_-1' → 5."""
    m = re.match(r"^del_(\d+)_to_-1$", label)
    if not m:
        return None
    return int(m.group(1))


# ── plots ─────────────────────────────────────────────────────────────────────

def plot_snv_heatmap(df, ref_row, upstream_bases, outdir, label):
    offsets_all = list(range(-upstream_bases, 0))
    offsets = sorted(set(offsets_all))
    alts    = list(BASES)

    matrix   = np.full((len(BASES), len(offsets)), np.nan)
    ref_seq  = {}

    for _, row in df.iterrows():
        parsed = snv_label_to_coords(row["tx_id"])
        if parsed is None:
            continue
        offset, ref_base, alt = parsed
        if offset not in offsets:
            continue
        oi = offsets.index(offset)
        ai = alts.index(alt)
        matrix[ai, oi] = row["mean_predicted_TE"] - ref_row["mean_predicted_TE"]
        ref_seq[offset] = ref_base

    xlabels = [f"{o}\n({ref_seq.get(o,'?')})" for o in offsets]

    vmax = np.nanmax(np.abs(matrix))
    if np.isnan(vmax) or vmax == 0:
        vmax = 0.1

    fig, ax = plt.subplots(figsize=(max(8, upstream_bases * 0.6 + 1), 3.5))
    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(offsets)))
    ax.set_xticklabels(xlabels, fontsize=8)
    ax.set_yticks(range(len(alts)))
    ax.set_yticklabels(alts)
    ax.set_xlabel("Position relative to AUG (ref base in parentheses)")
    ax.set_ylabel("Alternate base")
    ax.set_title(f"{label}: SNV effect on mean predicted TE (ΔTE vs reference)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("ΔTE")
    fig.tight_layout()
    path = outdir / "snv_heatmap.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved {path}")


def plot_del1_heatmap(df, ref_row, upstream_bases, outdir, label):
    offsets = sorted(range(-upstream_bases, 0))
    delta_by_offset = {}
    ref_by_offset = {}
    for _, row in df.iterrows():
        parsed = del1_label_to_coords(row["tx_id"])
        if parsed is None:
            continue
        offset, ref_base = parsed
        if offset in offsets:
            delta_by_offset[offset] = row["mean_predicted_TE"] - ref_row["mean_predicted_TE"]
            ref_by_offset[offset] = ref_base

    if not delta_by_offset:
        print("  No single-deletion variants found; skipping del1 heatmap.")
        return

    delta = [delta_by_offset.get(offset, np.nan) for offset in offsets]
    labels = [f"{offset} ({ref_by_offset.get(offset, '?')})" for offset in offsets]
    matrix = np.array(delta).reshape(1, -1)
    vmax   = max(np.nanmax(np.abs(matrix)), 0.01)

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 0.6 + 1), 1.8))
    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8, rotation=45, ha="right")
    ax.set_yticks([0])
    ax.set_yticklabels(["del1"])
    ax.set_title(f"{label}: single-base deletion ΔTE (per position)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("ΔTE")
    fig.tight_layout()
    path = outdir / "del1_heatmap.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved {path}")


def plot_waterfall(df, ref_row, outdir, label, top_n=40):
    d = delta_te(df, ref_row).copy()
    df2 = df.copy()
    df2["delta_TE"] = d
    df2 = df2[df2["tx_id"] != "reference"].copy()
    df2 = df2.sort_values("delta_TE")

    if len(df2) > top_n * 2:
        df2 = pd.concat([df2.head(top_n), df2.tail(top_n)]).drop_duplicates()

    colors = ["#d73027" if v < 0 else "#1a9850" for v in df2["delta_TE"]]

    fig, ax = plt.subplots(figsize=(max(10, len(df2) * 0.3 + 2), 5))
    ax.bar(range(len(df2)), df2["delta_TE"], color=colors, width=0.8)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(range(len(df2)))
    ax.set_xticklabels(df2["tx_id"], rotation=90, fontsize=7)
    ax.set_ylabel("ΔTE (vs reference)")
    ax.set_title(f"Waterfall: effect of all {label} 5'UTR variants on mean predicted TE")
    fig.tight_layout()
    path = outdir / "waterfall.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved {path}")


def plot_growing_deletion_trend(df, ref_row, outdir, label):
    rows = []
    for _, row in df.iterrows():
        n = deln_label_to_len(row["tx_id"])
        if n is not None:
            rows.append({"del_len": n, "delta_TE": row["mean_predicted_TE"] - ref_row["mean_predicted_TE"]})
    if not rows:
        print("  No growing-deletion variants found; skipping trend plot.")
        return

    rows_df = pd.DataFrame(rows).sort_values("del_len")

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(rows_df["del_len"], rows_df["delta_TE"], "o-", color="#2166ac", lw=2, ms=6)
    ax.axhline(0, color="grey", lw=0.8, ls="--")
    ax.set_xlabel("Bases deleted from 5'UTR end (before AUG)")
    ax.set_ylabel("ΔTE (vs reference)")
    ax.set_title(f"Effect of growing 5'UTR deletion on {label} predicted TE")
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    fig.tight_layout()
    path = outdir / "growing_deletion_trend.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved {path}")


def plot_neuron_panel(df, ref_row, outdir, label, species=None, top_n=20):
    # Find which neuron cols are present
    avail = neuron_columns_for_species(species, df.columns)
    if not avail:
        print("  No neuronal cell-type columns found; skipping neuron panel.")
        return

    df2 = df[df["tx_id"] != "reference"].copy()
    ref_vals = {c: ref_row[c] for c in avail}

    delta_mat = pd.DataFrame(
        {c: (df2[c] - ref_vals[c]).to_numpy() for c in avail},
        index=df2["tx_id"].values
    )

    # Select top_n variants by |mean ΔTE| across neuron types
    delta_mat["abs_mean"] = delta_mat.abs().mean(axis=1)
    top_idx = delta_mat["abs_mean"].abs().nlargest(top_n).index
    delta_mat = delta_mat.loc[top_idx, avail].rename(columns=clean_te_label)

    vmax = np.nanmax(np.abs(delta_mat.values))
    if vmax == 0:
        vmax = 0.1

    fig, ax = plt.subplots(figsize=(len(avail) * 1.2 + 2, max(6, len(top_idx) * 0.35 + 2)))
    sns.heatmap(
        delta_mat,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.3,
        annot=len(top_idx) <= 30,
        fmt=".2f",
        annot_kws={"size": 7},
        cbar_kws={"label": "ΔTE"},
    )
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Variant")
    ax.set_title(f"{label}: top {len(top_idx)} variants — neuronal cell-type ΔTE")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)
    fig.tight_layout()
    path = outdir / "neuron_panel.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved {path}")


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",           default="results/human/prediction_output.txt")
    parser.add_argument("--variants-input",  default="data/prediction_input.txt",
                        help="Only used to confirm reference row is present; not strictly required")
    parser.add_argument("--outdir",          default="plots_ism_scn2a")
    parser.add_argument("--upstream-bases",  type=int, default=15)
    parser.add_argument("--gene-name",       default="SCN2A")
    parser.add_argument("--species",         choices=("human", "mouse"), default=None)
    parser.add_argument("--top-waterfall",   type=int, default=40,
                        help="Number of top/bottom variants shown in waterfall (default 40)")
    parser.add_argument("--top-neuron",      type=int, default=20,
                        help="Number of variants shown in neuron panel (default 20)")
    return parser.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading predictions from {args.input} ...")
    df = load_results(args.input)
    print(f"  {len(df)} rows loaded.")

    ref_row = get_ref_te(df)
    print(f"  Reference mean_predicted_TE = {ref_row['mean_predicted_TE']:.4f}")

    label = title_label(args.gene_name, args.species)
    print("Generating plots ...")
    plot_snv_heatmap(df, ref_row, args.upstream_bases, outdir, label)
    plot_del1_heatmap(df, ref_row, args.upstream_bases, outdir, label)
    plot_waterfall(df, ref_row, outdir, label, top_n=args.top_waterfall)
    plot_growing_deletion_trend(df, ref_row, outdir, label)
    plot_neuron_panel(df, ref_row, outdir, label, species=args.species, top_n=args.top_neuron)

    print(f"\nAll plots saved to {outdir}/")


if __name__ == "__main__":
    main()
