#!/usr/bin/env python3
"""Prepare validated transcripts and balanced shards for all-5'UTR mutagenesis."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import pandas as pd

from workflow_common import (
    BASES,
    DEFAULT_REFS,
    MAX_CDS_UTR3_LEN,
    MAX_UTR5_LEN,
    STOP_CODONS,
    load_gtf_transcripts,
    lookup_transcript_sequence,
    open_text,
    read_transcript_fasta,
    transcript_region_bounds,
    transcript_sequence_from_genome,
    weighted_contiguous_shards,
    write_dict_rows,
    write_json,
)


CATALOG_FIELDS = [
    "tx_index",
    "transcript_id",
    "gene_id",
    "gene_name",
    "transcript_name",
    "chrom",
    "strand",
    "tx_sequence",
    "utr5_size",
    "cds_size",
    "utr3_size",
    "original_utr3_size",
    "utr3_truncated",
    "screen_mode",
    "target_position_count",
    "target_positions_1based",
    "orf_target_count",
    "orf_targets_json",
    "variant_count",
]

AUDIT_FIELDS = [
    "transcript_id",
    "gene_id",
    "gene_name",
    "status",
    "reason",
    "tx_size",
    "utr5_size",
    "cds_size",
    "utr3_size",
]

ORF_AUDIT_FIELDS = [
    "transcript_id",
    "gene_id",
    "gene_name",
    "orf_start_1based",
    "orf_stop_1based",
    "orf_frame",
    "annotated",
    "target_codon",
    "in_utr5",
    "status",
    "reason",
]


def normalize_row(
    tx_sequence: str,
    utr5_size: int,
    cds_size: int,
    metadata: dict,
    truncate_utr3: bool,
) -> tuple[dict | None, dict]:
    transcript_id = str(metadata.get("transcript_id", ""))
    audit = {
        "transcript_id": transcript_id,
        "gene_id": metadata.get("gene_id", ""),
        "gene_name": metadata.get("gene_name", ""),
        "status": "excluded",
        "reason": "",
        "tx_size": len(tx_sequence),
        "utr5_size": utr5_size,
        "cds_size": cds_size,
        "utr3_size": max(0, len(tx_sequence) - utr5_size - cds_size),
    }

    sequence = str(tx_sequence).strip().upper().replace("U", "T")
    if not transcript_id:
        audit["reason"] = "missing_transcript_id"
        return None, audit
    if utr5_size <= 0:
        audit["reason"] = "empty_utr5"
        return None, audit
    if cds_size <= 0 or utr5_size + cds_size > len(sequence):
        audit["reason"] = "invalid_region_lengths"
        return None, audit

    utr5 = sequence[:utr5_size]
    cds = sequence[utr5_size : utr5_size + cds_size]
    utr3 = sequence[utr5_size + cds_size :]
    audit.update(
        tx_size=len(sequence),
        utr5_size=len(utr5),
        cds_size=len(cds),
        utr3_size=len(utr3),
    )

    if any(base not in BASES for base in sequence):
        audit["reason"] = "non_acgt_sequence"
        return None, audit
    if len(utr5) > MAX_UTR5_LEN:
        audit["reason"] = f"utr5_exceeds_{MAX_UTR5_LEN}"
        return None, audit
    if len(cds) % 3 != 0:
        audit["reason"] = "cds_length_not_multiple_of_3"
        return None, audit
    if not cds.startswith("ATG"):
        audit["reason"] = "cds_does_not_start_ATG"
        return None, audit
    if cds[-3:] not in STOP_CODONS:
        audit["reason"] = "cds_missing_terminal_stop"
        return None, audit

    original_utr3_size = len(utr3)
    utr3_truncated = False
    if len(cds) + len(utr3) > MAX_CDS_UTR3_LEN:
        if not truncate_utr3:
            audit["reason"] = f"cds_utr3_exceeds_{MAX_CDS_UTR3_LEN}"
            return None, audit
        allowed = MAX_CDS_UTR3_LEN - len(cds)
        if allowed < 0:
            audit["reason"] = f"cds_exceeds_{MAX_CDS_UTR3_LEN}"
            return None, audit
        utr3 = utr3[:allowed]
        utr3_truncated = True

    normalized_sequence = utr5 + cds + utr3
    record = {
        "transcript_id": transcript_id,
        "gene_id": metadata.get("gene_id", ""),
        "gene_name": metadata.get("gene_name", ""),
        "transcript_name": metadata.get("transcript_name", ""),
        "chrom": metadata.get("chrom", ""),
        "strand": metadata.get("strand", ""),
        "tx_sequence": normalized_sequence,
        "utr5_size": len(utr5),
        "cds_size": len(cds),
        "utr3_size": len(utr3),
        "original_utr3_size": original_utr3_size,
        "utr3_truncated": int(utr3_truncated),
        "screen_mode": "all_utr5",
        "target_position_count": len(utr5),
        "target_positions_1based": "",
        "orf_target_count": 0,
        "orf_targets_json": "[]",
        "variant_count": 1 + 4 * len(utr5),
    }
    audit.update(
        status="included",
        reason="",
        tx_size=len(normalized_sequence),
        utr3_size=len(utr3),
    )
    return record, audit


def records_from_gtf(args) -> tuple[list[dict], list[dict]]:
    transcripts = load_gtf_transcripts(args.gtf)
    transcript_sequences = None
    genome = None
    if args.transcript_fasta:
        print(f"Loading full-transcript FASTA: {args.transcript_fasta}")
        transcript_sequences = read_transcript_fasta(args.transcript_fasta)
    else:
        print(f"Opening genome FASTA: {args.genome_fasta}")
        import pyfaidx

        genome = pyfaidx.Fasta(args.genome_fasta)

    records: list[dict] = []
    audit_rows: list[dict] = []
    try:
        for tx in transcripts:
            if args.limit_transcripts and len(audit_rows) >= args.limit_transcripts:
                break
            try:
                start, stop_end = transcript_region_bounds(tx)
                if transcript_sequences is not None:
                    sequence = lookup_transcript_sequence(
                        transcript_sequences, tx["transcript_id"]
                    )
                    if sequence is None:
                        raise ValueError("transcript_not_found_in_fasta")
                else:
                    sequence = transcript_sequence_from_genome(genome, tx)
                if stop_end > len(sequence):
                    raise ValueError("annotation_bounds_exceed_sequence")
                record, audit = normalize_row(
                    sequence,
                    start,
                    stop_end - start,
                    tx,
                    args.truncate_utr3,
                )
            except Exception as exc:
                record = None
                audit = {
                    "transcript_id": tx["transcript_id"],
                    "gene_id": tx.get("gene_id", ""),
                    "gene_name": tx.get("gene_name", ""),
                    "status": "excluded",
                    "reason": str(exc),
                    "tx_size": "",
                    "utr5_size": "",
                    "cds_size": "",
                    "utr3_size": "",
                }
            audit_rows.append(audit)
            if record is not None:
                records.append(record)
    finally:
        if genome is not None:
            genome.close()
    return records, audit_rows


def finite_int(value, column_name: str, transcript_id: str) -> int:
    if pd.isna(value):
        raise ValueError(f"{transcript_id}: missing {column_name}")
    numeric = float(value)
    if not math.isfinite(numeric) or not numeric.is_integer():
        raise ValueError(f"{transcript_id}: non-integer {column_name}={value}")
    return int(numeric)


def infer_table_sep(path: str | Path) -> str:
    suffixes = "".join(Path(path).suffixes).lower()
    if ".tsv" in suffixes or ".txt" in suffixes:
        return "\t"
    return ","


def first_existing_column(columns: Iterable[str], aliases: tuple[str, ...]) -> str | None:
    column_set = set(columns)
    for alias in aliases:
        if alias in column_set:
            return alias
    return None


def parse_orf_id(value: str) -> dict[str, str]:
    parts = str(value).rsplit("_", 3)
    if len(parts) != 4:
        raise ValueError(
            f"Could not parse orf_id={value!r}; expected "
            "<transcript_id>_<orf_start>_<orf_stop>_<orf_frame>"
        )
    transcript_id, start, stop, frame_value = parts
    return {
        "transcript_id": transcript_id,
        "orf_start": start,
        "orf_stop": stop,
        "orf_frame": frame_value,
    }


def inferred_annotated_value(row, annotated_col: str | None, transcript_id: str) -> int:
    if annotated_col:
        return finite_int(row[annotated_col], annotated_col, transcript_id)
    label = str(row.get("orf_label", row.get("model_class", ""))).lower()
    if label and "theoretical" not in label:
        return 1
    return 0


def load_orf_predictions(path: str | Path) -> dict[str, list[dict]]:
    frame = pd.read_csv(path, sep=infer_table_sep(path), dtype=str)
    aliases = {
        "transcript_id": ("transcript_id", "tx_id"),
        "orf_start": ("orf_start", "orf_start_1based"),
        "orf_stop": ("orf_stop", "orf_stop_1based"),
        "orf_frame": ("orf_frame", "frame"),
        "annotated": ("annotated", "orf_annotated", "is_annotated"),
    }
    selected = {
        key: first_existing_column(frame.columns, value)
        for key, value in aliases.items()
    }
    needs_orf_id = any(
        selected[key] is None
        for key in ("transcript_id", "orf_start", "orf_stop", "orf_frame")
    )
    if needs_orf_id and "orf_id" not in frame.columns:
        accepted = {
            key: list(value)
            for key, value in aliases.items()
            if key != "annotated"
        }
        raise ValueError(
            "ORF prediction table is missing coordinate columns. Expected "
            f"aliases {accepted}, or an orf_id formatted as "
            "<transcript_id>_<orf_start>_<orf_stop>_<orf_frame>."
        )

    predictions: dict[str, list[dict]] = {}
    seen: set[tuple[str, int, int, int]] = set()
    for _, row in frame.iterrows():
        parsed_orf_id = parse_orf_id(row["orf_id"]) if needs_orf_id else {}
        transcript_id = str(
            row[selected["transcript_id"]]
            if selected["transcript_id"]
            else parsed_orf_id["transcript_id"]
        )
        if not transcript_id or transcript_id == "nan":
            continue
        start_source = (
            row[selected["orf_start"]]
            if selected["orf_start"]
            else parsed_orf_id["orf_start"]
        )
        stop_source = (
            row[selected["orf_stop"]]
            if selected["orf_stop"]
            else parsed_orf_id["orf_stop"]
        )
        frame_source = (
            row[selected["orf_frame"]]
            if selected["orf_frame"]
            else parsed_orf_id["orf_frame"]
        )
        start = finite_int(start_source, "orf_start", transcript_id)
        stop = finite_int(stop_source, "orf_stop", transcript_id)
        frame_value = finite_int(frame_source, "orf_frame", transcript_id)
        annotated = inferred_annotated_value(row, selected["annotated"], transcript_id)
        key = (transcript_id, start, stop, frame_value)
        if key in seen:
            continue
        seen.add(key)
        predictions.setdefault(transcript_id, []).append(
            {
                "orf_start_1based": start,
                "orf_stop_1based": stop,
                "orf_frame": frame_value,
                "annotated": annotated,
            }
        )
    return predictions


def apply_orf_start_targets(
    records: list[dict],
    orf_predictions_path: str | Path,
    whole_atg_deletion: bool = False,
) -> tuple[list[dict], list[dict], dict[str, tuple[str, str]], int]:
    """Restrict the screen to ATG start codons fully inside the 5'UTR."""
    predictions_by_tx = load_orf_predictions(orf_predictions_path)
    filtered_records: list[dict] = []
    orf_audit_rows: list[dict] = []
    transcript_status: dict[str, tuple[str, str]] = {}
    non_atg_count = 0

    for record in records:
        transcript_id = record["transcript_id"]
        sequence = record["tx_sequence"]
        utr5_size = int(record["utr5_size"])
        predictions = predictions_by_tx.get(transcript_id, [])
        if not predictions:
            transcript_status[transcript_id] = ("excluded", "no_orf_prediction")
            continue

        target_positions: set[int] = set()
        target_by_start: dict[int, dict] = {}
        duplicate_counts: dict[int, int] = {}

        for prediction in predictions:
            start = int(prediction["orf_start_1based"])
            stop = int(prediction["orf_stop_1based"])
            frame_value = int(prediction["orf_frame"])
            annotated = int(prediction["annotated"])
            in_bounds = 1 <= start and start + 2 <= len(sequence)
            codon = sequence[start - 1 : start + 2] if in_bounds else ""
            in_utr5 = 1 <= start and start + 2 <= utr5_size
            status = "excluded"
            reason = ""

            if not in_bounds:
                reason = "orf_start_out_of_transcript_bounds"
            elif codon != "ATG":
                reason = "orf_start_not_ATG"
                non_atg_count += 1
            elif not in_utr5:
                reason = "not_fully_in_utr5"
            else:
                status = "included"
                reason = ""
                duplicate_counts[start] = duplicate_counts.get(start, 0) + 1
                if start not in target_by_start:
                    positions = [start, start + 1, start + 2]
                    target_by_start[start] = {
                        "orf_start_1based": start,
                        "orf_stop_1based": stop,
                        "orf_frame": frame_value,
                        "annotated": annotated,
                        "target_codon": codon,
                        "positions_1based": positions,
                        "n_orfs_same_start": 1,
                    }
                    target_positions.update(positions)
                else:
                    existing = target_by_start[start]
                    existing["annotated"] = max(int(existing["annotated"]), annotated)
                    existing["n_orfs_same_start"] = duplicate_counts[start]

            orf_audit_rows.append(
                {
                    "transcript_id": transcript_id,
                    "gene_id": record.get("gene_id", ""),
                    "gene_name": record.get("gene_name", ""),
                    "orf_start_1based": start,
                    "orf_stop_1based": stop,
                    "orf_frame": frame_value,
                    "annotated": annotated,
                    "target_codon": codon,
                    "in_utr5": int(in_utr5),
                    "status": status,
                    "reason": reason,
                }
            )

        if not target_by_start:
            transcript_status[transcript_id] = (
                "excluded",
                "no_valid_utr5_atg_orf_start",
            )
            continue

        targets = [target_by_start[start] for start in sorted(target_by_start)]
        positions = sorted(target_positions)
        if whole_atg_deletion:
            screen_mode = "whole_atg_deletion"
            target_positions_encoded = ";".join(
                str(target["orf_start_1based"]) for target in targets
            )
            target_position_count = len(targets)
            variant_count = 1 + len(targets)
            status_reason = "whole_atg_deletion_screen"
        else:
            screen_mode = "orf_start_codon"
            target_positions_encoded = ";".join(str(pos) for pos in positions)
            target_position_count = len(positions)
            variant_count = 1 + 4 * len(positions)
            status_reason = "orf_start_codon_screen"
        record = dict(record)
        record["screen_mode"] = screen_mode
        record["target_position_count"] = target_position_count
        record["target_positions_1based"] = target_positions_encoded
        record["orf_target_count"] = len(targets)
        record["orf_targets_json"] = json.dumps(targets, separators=(",", ":"))
        record["variant_count"] = variant_count
        filtered_records.append(record)
        transcript_status[transcript_id] = ("included", status_reason)

    return filtered_records, orf_audit_rows, transcript_status, non_atg_count


def records_from_table(args) -> tuple[list[dict], list[dict]]:
    frame = pd.read_csv(args.input_table, sep="\t", dtype={"tx_id": str})
    if "tx_id" not in frame.columns:
        raise ValueError("--input-table must contain tx_id")

    split_columns = {"utr5_sequence", "cds_sequence", "utr3_sequence"}
    full_columns = {"tx_sequence", "utr5_size", "cds_size"}
    if split_columns.issubset(frame.columns):
        frame = frame.fillna(
            {"utr5_sequence": "", "cds_sequence": "", "utr3_sequence": ""}
        )
        frame["tx_sequence"] = (
            frame["utr5_sequence"].astype(str)
            + frame["cds_sequence"].astype(str)
            + frame["utr3_sequence"].astype(str)
        )
        frame["utr5_size"] = frame["utr5_sequence"].astype(str).str.len()
        frame["cds_size"] = frame["cds_sequence"].astype(str).str.len()
    elif not full_columns.issubset(frame.columns):
        raise ValueError(
            "--input-table needs either utr5_sequence/cds_sequence/utr3_sequence "
            "or tx_sequence/utr5_size/cds_size"
        )

    records: list[dict] = []
    audit_rows: list[dict] = []
    for row_number, row in frame.iterrows():
        if args.limit_transcripts and row_number >= args.limit_transcripts:
            break
        metadata = {
            "transcript_id": row["tx_id"],
            "gene_id": row.get("gene_id", ""),
            "gene_name": row.get("gene_name", row.get("SYMBOL", "")),
            "transcript_name": row.get("transcript_name", ""),
            "chrom": row.get("chrom", ""),
            "strand": row.get("strand", ""),
        }
        record, audit = normalize_row(
            str(row["tx_sequence"]),
            int(row["utr5_size"]),
            int(row["cds_size"]),
            metadata,
            args.truncate_utr3,
        )
        audit_rows.append(audit)
        if record is not None:
            records.append(record)
    return records, audit_rows


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", choices=sorted(DEFAULT_REFS), default="human")
    parser.add_argument("--gtf", default=None)
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--transcript-fasta",
        help="Full transcript FASTA; GTF supplies UTR/CDS boundaries",
    )
    source.add_argument(
        "--genome-fasta",
        help="Genome FASTA used with the GTF to reconstruct transcripts",
    )
    source.add_argument(
        "--input-table",
        help="Existing RiboNN-format transcript TSV instead of FASTA+GTF",
    )
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--num-shards", type=int, default=128)
    parser.add_argument(
        "--orf-predictions",
        default=None,
        help=(
            "CSV/CSV.GZ with transcript_id,orf_start,orf_stop,orf_frame,annotated. "
            "When provided, only ATG start-codon bases fully inside the 5'UTR "
            "are mutated."
        ),
    )
    parser.add_argument(
        "--whole-atg-deletion",
        action="store_true",
        help=(
            "With --orf-predictions, predict only one full 3-base ATG deletion "
            "per retained ORF start instead of per-base substitutions/deletions."
        ),
    )
    parser.add_argument(
        "--no-truncate-utr3",
        dest="truncate_utr3",
        action="store_false",
        help="Exclude rather than truncate transcripts above the CDS+3'UTR limit",
    )
    parser.set_defaults(truncate_utr3=True)
    parser.add_argument(
        "--limit-transcripts",
        type=int,
        default=None,
        help="Testing aid: inspect only the first N annotated/input transcripts",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.num_shards <= 0:
        raise ValueError("--num-shards must be positive")
    if args.whole_atg_deletion and not args.orf_predictions:
        raise ValueError("--whole-atg-deletion requires --orf-predictions")

    defaults = DEFAULT_REFS[args.species]
    args.gtf = args.gtf or str(defaults["gtf"])
    if not args.input_table:
        args.genome_fasta = args.genome_fasta or str(defaults["genome_fasta"])
        if not Path(args.gtf).is_file():
            raise FileNotFoundError(f"GTF not found: {args.gtf}")
        source_path = args.transcript_fasta or args.genome_fasta
        if not Path(source_path).is_file():
            raise FileNotFoundError(f"FASTA not found: {source_path}")
        records, audit_rows = records_from_gtf(args)
        source = {
            "mode": "transcript_fasta" if args.transcript_fasta else "genome_fasta",
            "fasta": str(source_path),
            "gtf": str(args.gtf),
        }
    else:
        if not Path(args.input_table).is_file():
            raise FileNotFoundError(f"Input table not found: {args.input_table}")
        records, audit_rows = records_from_table(args)
        source = {"mode": "input_table", "input_table": str(args.input_table)}

    records.sort(key=lambda row: row["transcript_id"])
    seen: set[str] = set()
    unique_records = []
    for record in records:
        if record["transcript_id"] in seen:
            audit_rows.append(
                {
                    "transcript_id": record["transcript_id"],
                    "gene_id": record.get("gene_id", ""),
                    "gene_name": record.get("gene_name", ""),
                    "status": "excluded",
                    "reason": "duplicate_transcript_id",
                }
            )
            continue
        seen.add(record["transcript_id"])
        unique_records.append(record)
    records = unique_records
    if not records:
        raise ValueError("No transcripts passed validation; inspect catalog_audit.tsv.gz")

    orf_audit_rows: list[dict] = []
    transcript_target_status: dict[str, tuple[str, str]] = {}
    non_atg_orf_starts = 0
    screen_mode = "all_utr5"
    if args.orf_predictions:
        if not Path(args.orf_predictions).is_file():
            raise FileNotFoundError(f"ORF prediction CSV not found: {args.orf_predictions}")
        records, orf_audit_rows, transcript_target_status, non_atg_orf_starts = (
            apply_orf_start_targets(
                records,
                args.orf_predictions,
                whole_atg_deletion=args.whole_atg_deletion,
            )
        )
        screen_mode = "whole_atg_deletion" if args.whole_atg_deletion else "orf_start_codon"
        for audit in audit_rows:
            if audit.get("status") != "included":
                continue
            status, reason = transcript_target_status.get(
                audit["transcript_id"],
                ("excluded", "no_orf_prediction"),
            )
            audit["status"] = status
            audit["reason"] = "" if status == "included" else reason
        if not records:
            raise ValueError(
                "No transcripts retained after ORF-start filtering; inspect "
                "catalog_audit.tsv.gz and orf_target_audit.tsv.gz"
            )

    for tx_index, record in enumerate(records):
        record["tx_index"] = tx_index

    outdir = Path(args.outdir)
    catalog_dir = outdir / "catalog"
    shard_dir = catalog_dir / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)

    write_dict_rows(catalog_dir / "transcripts.tsv.gz", CATALOG_FIELDS, records)
    write_dict_rows(catalog_dir / "catalog_audit.tsv.gz", AUDIT_FIELDS, audit_rows)
    if args.orf_predictions:
        write_dict_rows(
            catalog_dir / "orf_target_audit.tsv.gz",
            ORF_AUDIT_FIELDS,
            orf_audit_rows,
        )
        if non_atg_orf_starts:
            raise ValueError(
                f"ORF start ATG check failed for {non_atg_orf_starts:,} "
                "prediction rows. Inspect catalog/orf_target_audit.tsv.gz."
            )

    shards = weighted_contiguous_shards(records, args.num_shards)
    shard_manifest_rows = []
    for shard_id, shard_records in enumerate(shards):
        shard_path = shard_dir / f"shard_{shard_id:03d}.transcripts.tsv.gz"
        write_dict_rows(shard_path, CATALOG_FIELDS, shard_records)
        shard_manifest_rows.append(
            {
                "shard_id": shard_id,
                "transcript_count": len(shard_records),
                "position_count": sum(
                    int(row["target_position_count"]) for row in shard_records
                ),
                "orf_target_count": sum(
                    int(row["orf_target_count"]) for row in shard_records
                ),
                "variant_count": sum(int(row["variant_count"]) for row in shard_records),
                "first_transcript_id": shard_records[0]["transcript_id"],
                "last_transcript_id": shard_records[-1]["transcript_id"],
            }
        )
    write_dict_rows(
        catalog_dir / "shard_manifest.tsv",
        [
            "shard_id",
            "transcript_count",
            "position_count",
            "orf_target_count",
            "variant_count",
            "first_transcript_id",
            "last_transcript_id",
        ],
        shard_manifest_rows,
    )

    included_utr5_positions = sum(int(row["utr5_size"]) for row in records)
    included_target_positions = sum(
        int(row["target_position_count"]) for row in records
    )
    included_orf_targets = sum(int(row["orf_target_count"]) for row in records)
    manifest = {
        "species": args.species,
        "source": source,
        "screen_mode": screen_mode,
        "orf_predictions": str(args.orf_predictions) if args.orf_predictions else None,
        "whole_atg_deletion": args.whole_atg_deletion,
        "truncate_utr3": args.truncate_utr3,
        "requested_shards": args.num_shards,
        "effective_shards": len(shards),
        "included_transcripts": len(records),
        "excluded_transcripts": sum(
            row.get("status") == "excluded" for row in audit_rows
        ),
        "included_utr5_positions": included_utr5_positions,
        "included_target_positions": included_target_positions,
        "included_orf_targets": included_orf_targets,
        "total_variants": sum(int(row["variant_count"]) for row in records),
        "score_definition": (
            "mean predicted TE across all model output cell types and test folds"
        ),
    }
    write_json(outdir / "workflow_manifest.json", manifest)

    print(f"Included transcripts : {manifest['included_transcripts']:,}")
    print(f"Excluded transcripts : {manifest['excluded_transcripts']:,}")
    print(f"Screen mode          : {manifest['screen_mode']}")
    print(f"5'UTR positions      : {manifest['included_utr5_positions']:,}")
    print(f"Target positions     : {manifest['included_target_positions']:,}")
    if args.orf_predictions:
        print(f"ORF start codons     : {manifest['included_orf_targets']:,}")
    print(f"Prediction variants  : {manifest['total_variants']:,}")
    print(f"Effective shards     : {manifest['effective_shards']}")
    print(f"Catalog              : {catalog_dir / 'transcripts.tsv.gz'}")
    print(f"Audit                : {catalog_dir / 'catalog_audit.tsv.gz'}")
    print(f"Manifest             : {outdir / 'workflow_manifest.json'}")


if __name__ == "__main__":
    main()
