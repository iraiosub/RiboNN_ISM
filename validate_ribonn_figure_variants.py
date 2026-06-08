#!/usr/bin/env python3
"""
validate_ribonn_figure_variants.py - Prepare and validate selected RiboNN figure SNVs.

Offsets are transcript coordinates relative to the canonical AUG:
  -1 is the base immediately before the A in ATG.

The script only creates the five requested SNVs, plus one matched reference row
per gene so predicted TE deltas can be calculated.

Usage:
    # Build the prediction input and verify the reference bases.
    python validate_ribonn_figure_variants.py

    # Run RiboNN and validate expected directions in one GPU job.
    python validate_ribonn_figure_variants.py --run-predict

    # Or validate an existing RiboNN output file.
    python validate_ribonn_figure_variants.py \\
        --predictions results/human/ribonn_figure_variants_output.txt
"""

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
import pyfaidx

from prepare_scn2a_ism import (
    BASES,
    DEFAULT_FASTA,
    DEFAULT_GTF,
    MAX_UTR5_LEN,
    _interval_len,
    _merge,
    _open,
    _parse_gtf_attrs,
    check_and_truncate,
    split_transcript,
    transcript_regions,
    transcript_sequence,
)


REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = Path("data/ribonn_figure_variants_input.txt")
DEFAULT_OUTPUT = Path("results/human/ribonn_figure_variants_output.txt")
DEFAULT_SUMMARY = Path("results/human/ribonn_figure_variants_summary.txt")
DEFAULT_PLOT = Path("results/human/ribonn_figure_variants_direction_plot.png")


@dataclass(frozen=True)
class FigureVariant:
    gene: str
    gene_names: tuple[str, ...]
    offset: int
    ref_base: str
    alt_base: str
    expected_direction: str

    @property
    def variant_id(self) -> str:
        return f"{self.gene}_{self.offset:+d}_{self.ref_base}>{self.alt_base}"

    @property
    def reference_id(self) -> str:
        return f"{self.gene}_reference"


@dataclass(frozen=True)
class PreparedVariant:
    spec: FigureVariant
    transcript: dict
    utr5_ref: str
    cds: str
    utr3: str
    utr5_mut: str
    transcript_pos: int
    warning: str


class ReferenceBaseMismatch(ValueError):
    """Raised when the requested ref base is absent in candidate transcripts."""


FIGURE_VARIANTS = (
    FigureVariant("ADAM32", ("ADAM32",), -61, "C", "T", "decreased"),
    FigureVariant("NUMA1", ("NUMA1",), -36, "G", "T", "decreased"),
    FigureVariant("COMT", ("COMT",), -212, "G", "A", "decreased"),
    FigureVariant("QARS", ("QARS", "QARS1"), -14, "C", "T", "increased"),
    FigureVariant("AKT3", ("AKT3",), -75, "G", "A", "increased"),
)


def load_gene_transcripts(gtf_path: str, gene_names: Iterable[str], transcript_id: str | None = None):
    """Return candidate transcripts for any gene_name in gene_names."""
    names = set(gene_names)
    transcripts = {}
    keep = {"exon", "CDS", "start_codon", "stop_codon"}

    with _open(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _src, feature, start, end, _sc, strand, _fr, attrs = fields
            if feature not in keep:
                continue

            parsed = _parse_gtf_attrs(attrs)
            if parsed.get("gene_name") not in names:
                continue

            tid = parsed.get("transcript_id")
            if not tid:
                continue
            if transcript_id and tid != transcript_id and tid.split(".")[0] != transcript_id.split(".")[0]:
                continue

            tx = transcripts.setdefault(
                tid,
                {
                    "chrom": chrom,
                    "strand": strand,
                    "gene_name": parsed.get("gene_name", ""),
                    "transcript_id": tid,
                    "transcript_name": parsed.get("transcript_name", ""),
                    "exons": [],
                    "cds": [],
                    "start_codons": [],
                    "stop_codons": [],
                },
            )

            s0, e0 = int(start) - 1, int(end)
            if feature == "exon":
                tx["exons"].append((s0, e0))
            elif feature == "CDS":
                tx["cds"].append((s0, e0))
            elif feature == "start_codon":
                tx["start_codons"].append((s0, e0))
            elif feature == "stop_codon":
                tx["stop_codons"].append((s0, e0))

    candidates = []
    for tx in transcripts.values():
        for key in ("exons", "cds", "start_codons", "stop_codons"):
            tx[key] = _merge(tx[key])
        if tx["exons"] and tx["cds"]:
            candidates.append(tx)

    candidates.sort(key=lambda t: (_interval_len(t["cds"]), _interval_len(t["exons"])), reverse=True)
    return candidates


def parse_transcript_overrides(values: list[str]) -> dict[str, str]:
    overrides = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"Invalid --transcript-id value '{value}'. Use GENE=ENST...")
        gene, transcript_id = value.split("=", 1)
        gene = gene.strip()
        transcript_id = transcript_id.strip()
        if not gene or not transcript_id:
            raise ValueError(f"Invalid --transcript-id value '{value}'. Use GENE=ENST...")
        overrides[gene] = transcript_id
    return overrides


def transcript_override_for(spec: FigureVariant, overrides: dict[str, str]) -> str | None:
    for gene_name in (spec.gene, *spec.gene_names):
        if gene_name in overrides:
            return overrides[gene_name]
    return None


def observed_base_at_offset(fasta, transcript: dict, offset: int):
    tx_seq, coords = transcript_sequence(fasta, transcript)
    regions = transcript_regions(transcript, coords)
    canonical_start = regions["canonical_start"]
    transcript_pos = canonical_start + offset

    if transcript_pos < 0 or transcript_pos >= canonical_start:
        return tx_seq, regions, transcript_pos, None

    return tx_seq, regions, transcript_pos, tx_seq[transcript_pos].upper()


def prepare_variant(fasta, gtf_path: str, spec: FigureVariant, transcript_id: str | None, truncate_utr3: bool):
    candidates = load_gene_transcripts(gtf_path, spec.gene_names, transcript_id=transcript_id)
    if not candidates:
        suffix = f" transcript {transcript_id}" if transcript_id else ""
        raise ValueError(f"No candidate transcript found for {spec.gene}{suffix} in {gtf_path}")

    observations = []
    for transcript in candidates:
        tx_seq, regions, transcript_pos, observed = observed_base_at_offset(fasta, transcript, spec.offset)
        if observed is None:
            observations.append(
                f"{transcript['transcript_id']}({transcript['transcript_name']}): offset outside 5'UTR"
            )
            continue

        observations.append(
            f"{transcript['transcript_id']}({transcript['transcript_name']}): {observed}"
        )
        if observed != spec.ref_base:
            continue

        utr5_ref, cds, utr3 = split_transcript(
            tx_seq,
            regions["canonical_start"],
            regions["canonical_stop_end"],
        )
        utr5_ref, cds, utr3, warning = check_and_truncate(
            utr5_ref,
            cds,
            utr3,
            truncate_utr3,
        )

        if len(utr5_ref) > MAX_UTR5_LEN:
            raise ValueError(
                f"{spec.gene} {transcript['transcript_id']} 5'UTR is {len(utr5_ref)} nt, "
                f"which exceeds RiboNN limit {MAX_UTR5_LEN}."
            )
        if not cds.startswith("ATG"):
            raise ValueError(
                f"{spec.gene} {transcript['transcript_id']} CDS starts with {cds[:3]}, not ATG."
            )
        if spec.alt_base not in BASES:
            raise ValueError(f"{spec.variant_id} alt base is not one of {BASES}")

        utr5_mut = utr5_ref[:transcript_pos] + spec.alt_base + utr5_ref[transcript_pos + 1:]
        return PreparedVariant(
            spec=spec,
            transcript=transcript,
            utr5_ref=utr5_ref,
            cds=cds,
            utr3=utr3,
            utr5_mut=utr5_mut,
            transcript_pos=transcript_pos,
            warning=warning,
        )

    checked = "; ".join(observations[:12])
    extra = "" if len(observations) <= 12 else f"; ... {len(observations) - 12} more"
    raise ReferenceBaseMismatch(
        f"{spec.variant_id}: expected reference base {spec.ref_base} at offset {spec.offset}, "
        f"but no candidate transcript matched. Checked {checked}{extra}"
    )


def write_prediction_input(prepared: list[PreparedVariant], output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["tx_id", "tx_sequence", "utr5_size", "cds_size"])
        for item in prepared:
            writer.writerow([
                item.spec.reference_id,
                item.utr5_ref + item.cds + item.utr3,
                len(item.utr5_ref),
                len(item.cds),
            ])
            writer.writerow([
                item.spec.variant_id,
                item.utr5_mut + item.cds + item.utr3,
                len(item.utr5_mut),
                len(item.cds),
            ])


def aggregate_predictions(df: pd.DataFrame) -> pd.DataFrame:
    if "tx_id" not in df.columns:
        unnamed_cols = [col for col in df.columns if str(col).startswith("Unnamed:")]
        if unnamed_cols:
            df = df.rename(columns={unnamed_cols[0]: "tx_id"})
        else:
            df = df.reset_index().rename(columns={"index": "tx_id"})

    predicted_cols = [col for col in df.columns if col.startswith("predicted_")]
    value_cols = predicted_cols.copy()
    if "mean_predicted_TE" in df.columns:
        value_cols.append("mean_predicted_TE")

    if "mean_predicted_TE" in df.columns and not df["tx_id"].duplicated().any():
        return df

    if not value_cols:
        raise ValueError("Prediction file has no mean_predicted_TE or predicted_* columns.")

    group_cols = [
        col
        for col in ("tx_id", "utr5_sequence", "cds_sequence", "utr3_sequence")
        if col in df.columns
    ]
    aggregated = df.groupby(group_cols, as_index=False)[value_cols].mean()
    if predicted_cols:
        aggregated["mean_predicted_TE"] = aggregated[predicted_cols].mean(axis=1)
    return aggregated


def run_predictions(input_path: Path, output_path: Path, top_k: int, batch_size: int, num_workers: int):
    from run_ribonn_predict import download_weights, load_run_df
    from src.predict import predict_using_nested_cross_validation_models

    download_weights()
    run_df = load_run_df()

    raw = predict_using_nested_cross_validation_models(
        input_path=str(input_path),
        species="human",
        run_df=run_df,
        top_k_models_to_use=top_k,
        batch_size=batch_size,
        num_workers=num_workers,
    )
    aggregated = aggregate_predictions(raw)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    aggregated.to_csv(output_path, sep="\t", index=False)
    return aggregated


def load_predictions(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    return aggregate_predictions(df)


def status_symbol(passed: bool) -> str:
    return "\N{CHECK MARK}" if passed else "x"


def plot_direction_summary(summary: pd.DataFrame, output_path: Path):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] matplotlib is not installed; skipping direction plot.")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)

    plot_df = summary.copy()
    plot_df["variant_label"] = plot_df.apply(
        lambda row: (
            f"{row['gene']} {row['offset_from_aug']:+d} "
            f"{row['ref_base']}>{row['alt_base']}"
        ),
        axis=1,
    )
    plot_df["detail_label"] = plot_df.apply(
        lambda row: (
            f"expected {row['expected_direction']}, "
            f"observed {row['observed_direction']}"
        ),
        axis=1,
    )

    y = range(len(plot_df))
    deltas = plot_df["delta_mean_predicted_TE"].astype(float)
    max_abs = max(float(deltas.abs().max()), 1e-6)
    colors = ["#218c74" if passed else "#b33939" for passed in plot_df["pass"]]

    fig_height = max(3.8, 0.62 * len(plot_df) + 1.8)
    fig, ax = plt.subplots(figsize=(10, fig_height))
    ax.barh(y, deltas, color=colors, height=0.5)
    ax.axvline(0, color="#333333", lw=1.0)
    ax.grid(axis="x", color="#d9d9d9", lw=0.7, alpha=0.8)
    ax.set_axisbelow(True)

    ax.set_yticks(list(y))
    ax.set_yticklabels(plot_df["variant_label"], fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel("Delta mean predicted TE vs matched reference")
    ax.set_title("RiboNN Figure Variant Direction Check", fontsize=13, pad=12)

    x_left = -max_abs * 1.35
    x_right = max_abs * 2.25
    ax.set_xlim(x_left, x_right)
    detail_x = max_abs * 1.1
    status_x = max_abs * 2.0

    for i, row in plot_df.iterrows():
        symbol = status_symbol(bool(row["pass"]))
        ax.text(
            status_x,
            i,
            symbol,
            va="center",
            ha="center",
            fontsize=18,
            fontweight="bold",
            color="#218c74" if row["pass"] else "#b33939",
        )
        ax.text(
            detail_x,
            i,
            row["detail_label"],
            va="center",
            ha="left",
            fontsize=9,
            color="#333333",
        )

    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    print(f"Direction plot written to {output_path}")


def validate_predictions(predictions: pd.DataFrame, prepared: list[PreparedVariant], summary_path: Path, tolerance: float):
    by_id = predictions.set_index("tx_id", drop=False)
    rows = []

    for item in prepared:
        ref_id = item.spec.reference_id
        var_id = item.spec.variant_id
        missing = [tx_id for tx_id in (ref_id, var_id) if tx_id not in by_id.index]
        if missing:
            raise ValueError(
                f"Prediction output is missing rows for: {', '.join(missing)}. "
                "They may have been filtered for sequence length."
            )

        ref_te = float(by_id.loc[ref_id, "mean_predicted_TE"])
        var_te = float(by_id.loc[var_id, "mean_predicted_TE"])
        delta = var_te - ref_te

        if delta > tolerance:
            observed = "increased"
        elif delta < -tolerance:
            observed = "decreased"
        else:
            observed = "unchanged"

        passed = observed == item.spec.expected_direction
        rows.append(
            {
                "status": status_symbol(passed),
                "gene": item.spec.gene,
                "variant": item.spec.variant_id,
                "transcript_id": item.transcript["transcript_id"],
                "transcript_name": item.transcript["transcript_name"],
                "gtf_gene_name": item.transcript["gene_name"],
                "offset_from_aug": item.spec.offset,
                "transcript_pos_0based": item.transcript_pos,
                "ref_base": item.spec.ref_base,
                "alt_base": item.spec.alt_base,
                "expected_direction": item.spec.expected_direction,
                "observed_direction": observed,
                "delta_mean_predicted_TE": delta,
                "reference_mean_predicted_TE": ref_te,
                "variant_mean_predicted_TE": var_te,
                "pass": passed,
            }
        )

    summary = pd.DataFrame(rows)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, sep="\t", index=False)

    print("\nDirection validation:")
    for row in rows:
        print(
            f"  {row['status']} {row['variant']}: expected {row['expected_direction']}, "
            f"observed {row['observed_direction']} "
            f"(delta mean_predicted_TE={row['delta_mean_predicted_TE']:.6g})"
        )
    print(f"\nSummary written to {summary_path}")

    return summary


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", default=DEFAULT_FASTA, help="GRCh38 primary assembly FASTA")
    parser.add_argument("--gtf", default=DEFAULT_GTF, help="GENCODE annotation GTF, optionally gzipped")
    parser.add_argument("--output-input", default=str(DEFAULT_INPUT), help="Prediction input TSV to write")
    parser.add_argument("--predictions", default=None, help="Existing prediction output TSV to validate")
    parser.add_argument("--output-predictions", default=str(DEFAULT_OUTPUT), help="Prediction output TSV for --run-predict")
    parser.add_argument("--summary", default=str(DEFAULT_SUMMARY), help="Validation summary TSV")
    parser.add_argument("--plot", default=str(DEFAULT_PLOT), help="Direction-check plot PNG")
    parser.add_argument("--run-predict", action="store_true", help="Run RiboNN predictions after writing the input")
    parser.add_argument("--top-k", type=int, default=5, help="Top-k models per fold for --run-predict")
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--num-workers", type=int, default=4)
    parser.add_argument("--tolerance", type=float, default=0.0, help="Absolute delta treated as unchanged")
    parser.add_argument("--no-truncate-utr3", action="store_true", help="Do not truncate 3'UTR to fit RiboNN limits")
    parser.add_argument(
        "--transcript-id",
        action="append",
        default=[],
        metavar="GENE=ENST...",
        help="Force a transcript for a gene; may be repeated",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.chdir(REPO_ROOT)

    transcript_overrides = parse_transcript_overrides(args.transcript_id)
    output_input = Path(args.output_input)
    output_predictions = Path(args.output_predictions)
    summary_path = Path(args.summary)
    plot_path = Path(args.plot)

    print("Preparing requested RiboNN figure variants ...")
    fasta = pyfaidx.Fasta(args.fasta)
    try:
        prepared = []
        skipped = []
        for spec in FIGURE_VARIANTS:
            transcript_id = transcript_override_for(spec, transcript_overrides)
            try:
                item = prepare_variant(
                    fasta=fasta,
                    gtf_path=args.gtf,
                    spec=spec,
                    transcript_id=transcript_id,
                    truncate_utr3=not args.no_truncate_utr3,
                )
            except ReferenceBaseMismatch as exc:
                skipped.append(spec.variant_id)
                print(f"  [WARN] Skipping {spec.variant_id}: {exc}")
                continue
            prepared.append(item)
            print(
                f"  {spec.variant_id}: {item.transcript['transcript_id']} "
                f"({item.transcript['transcript_name']}), "
                f"5'UTR[{item.transcript_pos}]={spec.ref_base}, "
                f"expected {spec.expected_direction}"
            )
            if item.warning:
                print(item.warning.strip())
    finally:
        fasta.close()

    write_prediction_input(prepared, output_input)
    print(f"\nWrote {len(prepared) * 2} rows to {output_input}")
    if skipped:
        print(f"Skipped {len(skipped)} variant(s) due to reference-base mismatch: {', '.join(skipped)}")

    if not prepared:
        print("No variants remain after reference-base checks; nothing to predict or plot.")
        return

    predictions = None
    if args.run_predict:
        print("\nRunning RiboNN predictions ...")
        predictions = run_predictions(
            input_path=output_input,
            output_path=output_predictions,
            top_k=args.top_k,
            batch_size=args.batch_size,
            num_workers=args.num_workers,
        )
        print(f"Predictions written to {output_predictions}")
    elif args.predictions:
        predictions = load_predictions(Path(args.predictions))

    if predictions is None:
        print("\nNext steps:")
        print(
            f"  python run_ribonn_predict.py --input {output_input} "
            f"--output {output_predictions} --top-k {args.top_k} "
            f"--batch-size {args.batch_size} --num-workers {args.num_workers}"
        )
        print(f"  python {Path(__file__).name} --predictions {output_predictions}")
        return

    summary = validate_predictions(predictions, prepared, summary_path, args.tolerance)
    plot_direction_summary(summary, plot_path)


if __name__ == "__main__":
    main()
