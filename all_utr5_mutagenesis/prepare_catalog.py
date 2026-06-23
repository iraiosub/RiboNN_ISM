#!/usr/bin/env python3
"""Prepare validated transcripts and balanced shards for all-5'UTR mutagenesis."""

from __future__ import annotations

import argparse
from pathlib import Path

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

    for tx_index, record in enumerate(records):
        record["tx_index"] = tx_index

    outdir = Path(args.outdir)
    catalog_dir = outdir / "catalog"
    shard_dir = catalog_dir / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)

    write_dict_rows(catalog_dir / "transcripts.tsv.gz", CATALOG_FIELDS, records)
    write_dict_rows(catalog_dir / "catalog_audit.tsv.gz", AUDIT_FIELDS, audit_rows)

    shards = weighted_contiguous_shards(records, args.num_shards)
    shard_manifest_rows = []
    for shard_id, shard_records in enumerate(shards):
        shard_path = shard_dir / f"shard_{shard_id:03d}.transcripts.tsv.gz"
        write_dict_rows(shard_path, CATALOG_FIELDS, shard_records)
        shard_manifest_rows.append(
            {
                "shard_id": shard_id,
                "transcript_count": len(shard_records),
                "position_count": sum(int(row["utr5_size"]) for row in shard_records),
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
            "variant_count",
            "first_transcript_id",
            "last_transcript_id",
        ],
        shard_manifest_rows,
    )

    included_positions = sum(int(row["utr5_size"]) for row in records)
    manifest = {
        "species": args.species,
        "source": source,
        "truncate_utr3": args.truncate_utr3,
        "requested_shards": args.num_shards,
        "effective_shards": len(shards),
        "included_transcripts": len(records),
        "excluded_transcripts": sum(
            row.get("status") == "excluded" for row in audit_rows
        ),
        "included_utr5_positions": included_positions,
        "total_variants": sum(int(row["variant_count"]) for row in records),
        "score_definition": (
            "mean predicted TE across all model output cell types and test folds"
        ),
    }
    write_json(outdir / "workflow_manifest.json", manifest)

    print(f"Included transcripts : {manifest['included_transcripts']:,}")
    print(f"Excluded transcripts : {manifest['excluded_transcripts']:,}")
    print(f"5'UTR positions      : {manifest['included_utr5_positions']:,}")
    print(f"Prediction variants  : {manifest['total_variants']:,}")
    print(f"Effective shards     : {manifest['effective_shards']}")
    print(f"Catalog              : {catalog_dir / 'transcripts.tsv.gz'}")
    print(f"Audit                : {catalog_dir / 'catalog_audit.tsv.gz'}")
    print(f"Manifest             : {outdir / 'workflow_manifest.json'}")


if __name__ == "__main__":
    main()
