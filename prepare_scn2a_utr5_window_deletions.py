#!/usr/bin/env python3
"""
prepare_scn2a_utr5_window_deletions.py

Generate RiboNN prediction input for a 5'UTR sliding-window deletion experiment.

The default experiment deletes every 15-nt window across the SCN2A 5'UTR:
positions 1-15, then 2-16, then 3-17, through the final window immediately
before the canonical AUG.

Output columns match RiboNN:
    tx_id | utr5_sequence | cds_sequence | utr3_sequence
"""

import argparse
import csv
from pathlib import Path

import pyfaidx

from prepare_scn2a_ism import (
    DEFAULT_FASTA,
    DEFAULT_GTF,
    check_and_truncate,
    load_gene_transcript,
    split_transcript,
    transcript_regions,
    transcript_sequence,
)


def window_tx_id(window_size, start0, end0, utr5_len):
    """Return a parseable tx_id for a deleted UTR5 window."""
    start_1 = start0 + 1
    end_1 = end0
    start_offset = start0 - utr5_len
    end_offset = end0 - 1 - utr5_len
    return (
        f"delwin{window_size}_pos{start_1:04d}_{end_1:04d}"
        f"_off{start_offset:+d}_{end_offset:+d}"
    )


def generate_window_deletions(utr5, cds, utr3, window_size, stride):
    """Yield reference plus every sliding-window deletion row."""
    utr5_len = len(utr5)
    yield "reference", utr5, cds, utr3

    for start0 in range(0, utr5_len - window_size + 1, stride):
        end0 = start0 + window_size
        tx_id = window_tx_id(window_size, start0, end0, utr5_len)
        yield tx_id, utr5[:start0] + utr5[end0:], cds, utr3


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", default=DEFAULT_FASTA, help="GRCh38 primary assembly FASTA")
    parser.add_argument("--gtf", default=DEFAULT_GTF, help="GENCODE annotation GTF, optionally gzipped")
    parser.add_argument("--gene-name", default="SCN2A")
    parser.add_argument("--transcript-id", default=None, help="Specific transcript ID")
    parser.add_argument("--window-size", type=int, default=15, help="Deleted window size in nt")
    parser.add_argument("--stride", type=int, default=1, help="Window step size in nt")
    parser.add_argument("--truncate-utr3", action="store_true")
    parser.add_argument(
        "--output",
        default="data/scn2a_utr5_delwin15_input.txt",
        help="Output TSV path",
    )
    parser.add_argument("--audit", action="store_true", help="Print experiment info and exit")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.window_size <= 0:
        raise ValueError("--window-size must be positive")
    if args.stride <= 0:
        raise ValueError("--stride must be positive")

    print(f"Loading transcript for {args.gene_name} ...")
    fasta = pyfaidx.Fasta(args.fasta)
    transcript = load_gene_transcript(args.gtf, args.gene_name, args.transcript_id)
    tx_seq, coords = transcript_sequence(fasta, transcript)
    regions = transcript_regions(transcript, coords)
    fasta.close()

    utr5_ref, cds, utr3 = split_transcript(
        tx_seq,
        regions["canonical_start"],
        regions["canonical_stop_end"],
    )
    utr5_ref, cds, utr3, length_warning = check_and_truncate(
        utr5_ref,
        cds,
        utr3,
        args.truncate_utr3,
    )

    utr5_len = len(utr5_ref)
    if args.window_size > utr5_len:
        raise ValueError(
            f"--window-size {args.window_size} exceeds 5'UTR length {utr5_len}"
        )

    n_windows = ((utr5_len - args.window_size) // args.stride) + 1
    final_start0 = utr5_len - args.window_size
    final_end0 = utr5_len

    print(f"Transcript : {transcript['transcript_id']} ({transcript['transcript_name']})")
    print(f"Strand     : {transcript['strand']}  Chrom: {transcript['chrom']}")
    print(f"5'UTR      : {utr5_len} nt")
    print(f"CDS        : {len(cds)} nt")
    print(f"3'UTR      : {len(utr3)} nt")
    print(f"Window size: {args.window_size} nt")
    print(f"Stride     : {args.stride} nt")
    print(f"Windows    : {n_windows}")
    print(
        "Target final window before AUG: "
        f"{window_tx_id(args.window_size, final_start0, final_end0, utr5_len)}"
    )

    if length_warning:
        print()
        print(length_warning.strip())

    if args.audit:
        print("\n--audit flag set; exiting without writing output.")
        return

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["tx_id", "utr5_sequence", "cds_sequence", "utr3_sequence"])
        for row in generate_window_deletions(
            utr5_ref,
            cds,
            utr3,
            args.window_size,
            args.stride,
        ):
            writer.writerow(row)
            written += 1

    print(f"\nWrote {written} rows to {out_path}")


if __name__ == "__main__":
    main()
