#!/usr/bin/env python3
"""Shared helpers for the transcriptome-wide 5'UTR mutagenesis workflow."""

from __future__ import annotations

import csv
import gzip
import heapq
import json
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np


BASES = ("A", "C", "G", "T")
STOP_CODONS = {"TAA", "TAG", "TGA"}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

MAX_UTR5_LEN = 1_381
MAX_CDS_UTR3_LEN = 11_937
MAX_TX_LEN = MAX_UTR5_LEN + MAX_CDS_UTR3_LEN

DEFAULT_REF_ROOT = Path("/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref")
DEFAULT_REFS = {
    "human": {
        "genome_fasta": DEFAULT_REF_ROOT
        / "human"
        / "GRCh38.primary_assembly.genome.fa",
        "gtf": DEFAULT_REF_ROOT
        / "human"
        / "gencode.v44.primary_assembly.annotation.longest_cds_transcripts.gtf.gz",
        "orf_predictions": DEFAULT_REF_ROOT
        / "human"
        / "orfs"
        / "gencode.v44.primary_assembly.annotation.longest_cds.transcript_info.orf_predictions.csv.gz",
    },
    "mouse": {
        "genome_fasta": DEFAULT_REF_ROOT
        / "mouse"
        / "GRCm39.primary_assembly.genome.fa",
        "gtf": DEFAULT_REF_ROOT
        / "mouse"
        / "gencode.vM33.primary_assembly.annotation.longest_cds_transcripts.gtf.gz",
        "orf_predictions": DEFAULT_REF_ROOT
        / "mouse"
        / "orfs"
        / "gencode.vM33.primary_assembly.annotation.longest_cds.transcript_info.orf_predictions.csv.gz",
    },
}


def open_text(path: str | Path, mode: str = "rt", compresslevel: int = 1):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode, compresslevel=compresslevel)
    return path.open(mode)


def parse_gtf_attrs(value: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in value.strip().rstrip(";").split(";"):
        item = item.strip()
        if item and " " in item:
            key, val = item.split(" ", 1)
            attrs[key] = val.strip().strip('"')
    return attrs


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    intervals = sorted(intervals)
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def load_gtf_transcripts(gtf_path: str | Path) -> list[dict]:
    """Load transcript structures needed to identify UTR/CDS boundaries."""
    keep = {"exon", "CDS", "start_codon", "stop_codon"}
    transcripts: dict[str, dict] = {}
    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, attrs = fields
            if feature not in keep:
                continue
            parsed = parse_gtf_attrs(attrs)
            transcript_id = parsed.get("transcript_id")
            if not transcript_id:
                continue
            tx = transcripts.setdefault(
                transcript_id,
                {
                    "transcript_id": transcript_id,
                    "transcript_name": parsed.get("transcript_name", ""),
                    "gene_id": parsed.get("gene_id", ""),
                    "gene_name": parsed.get("gene_name", ""),
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [],
                    "cds": [],
                    "start_codons": [],
                    "stop_codons": [],
                },
            )
            start0, end0 = int(start) - 1, int(end)
            if feature == "exon":
                tx["exons"].append((start0, end0))
            elif feature == "CDS":
                tx["cds"].append((start0, end0))
            elif feature == "start_codon":
                tx["start_codons"].append((start0, end0))
            elif feature == "stop_codon":
                tx["stop_codons"].append((start0, end0))

    output = []
    for tx in transcripts.values():
        for key in ("exons", "cds", "start_codons", "stop_codons"):
            tx[key] = merge_intervals(tx[key])
        if tx["exons"] and tx["cds"]:
            output.append(tx)
    return sorted(output, key=lambda row: row["transcript_id"])


def transcript_coords(transcript: dict) -> np.ndarray:
    if transcript["strand"] == "+":
        exons = sorted(transcript["exons"])
        chunks = [np.arange(start, end, dtype=np.int64) for start, end in exons]
    else:
        exons = sorted(transcript["exons"], reverse=True)
        chunks = [
            np.arange(end - 1, start - 1, -1, dtype=np.int64)
            for start, end in exons
        ]
    return np.concatenate(chunks)


def interval_mask(coords: np.ndarray, intervals: Iterable[tuple[int, int]]) -> np.ndarray:
    mask = np.zeros(len(coords), dtype=bool)
    for start, end in intervals:
        mask |= (coords >= start) & (coords < end)
    return mask


def transcript_region_bounds(transcript: dict) -> tuple[int, int]:
    """Return transcript-coordinate CDS start and stop-codon end."""
    coords = transcript_coords(transcript)
    cds_idx = np.flatnonzero(interval_mask(coords, transcript["cds"]))
    if len(cds_idx) == 0:
        raise ValueError("no CDS bases")
    start_idx = np.flatnonzero(interval_mask(coords, transcript["start_codons"]))
    stop_idx = np.flatnonzero(interval_mask(coords, transcript["stop_codons"]))
    canonical_start = int(start_idx.min()) if len(start_idx) else int(cds_idx.min())
    canonical_stop_end = (
        int(stop_idx.max() + 1) if len(stop_idx) else int(cds_idx.max() + 1)
    )
    return canonical_start, canonical_stop_end


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1].upper()


def transcript_sequence_from_genome(genome, transcript: dict) -> str:
    chrom = transcript["chrom"]
    if transcript["strand"] == "+":
        exons = sorted(transcript["exons"])
        chunks = [genome[chrom][start:end].seq.upper() for start, end in exons]
    else:
        exons = sorted(transcript["exons"], reverse=True)
        chunks = [
            reverse_complement(genome[chrom][start:end].seq)
            for start, end in exons
        ]
    return "".join(chunks)


def fasta_aliases(header: str) -> list[str]:
    token = header.strip().split()[0]
    first_pipe = token.split("|", 1)[0]
    aliases = [token, first_pipe]
    aliases.extend(alias.split(".", 1)[0] for alias in list(aliases))
    return list(dict.fromkeys(alias for alias in aliases if alias))


def read_transcript_fasta(path: str | Path) -> dict[str, str]:
    """Read a full-transcript FASTA and index common Ensembl header aliases."""
    sequences: dict[str, str] = {}
    with open_text(path) as handle:
        header = None
        chunks: list[str] = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(chunks).upper().replace("U", "T")
                    for alias in fasta_aliases(header):
                        sequences.setdefault(alias, sequence)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            sequence = "".join(chunks).upper().replace("U", "T")
            for alias in fasta_aliases(header):
                sequences.setdefault(alias, sequence)
    return sequences


def lookup_transcript_sequence(
    sequence_index: dict[str, str], transcript_id: str
) -> str | None:
    for alias in (transcript_id, transcript_id.split(".", 1)[0]):
        if alias in sequence_index:
            return sequence_index[alias]
    return None


def variant_id(
    tx_index: int,
    kind: str,
    position_1based: int | None = None,
    ref: str | None = None,
    alt: str | None = None,
) -> str:
    if kind == "wt":
        return f"t{tx_index}|wt"
    parts = [f"t{tx_index}", kind, str(position_1based), str(ref)]
    if alt is not None:
        parts.append(str(alt))
    return "|".join(parts)


def parse_variant_id(value: str) -> dict:
    parts = value.split("|")
    if len(parts) < 2 or not parts[0].startswith("t"):
        raise ValueError(f"Invalid variant ID: {value}")
    parsed = {"tx_index": int(parts[0][1:]), "kind": parts[1]}
    if parsed["kind"] == "wt":
        if len(parts) != 2:
            raise ValueError(f"Invalid WT variant ID: {value}")
        return parsed
    if parsed["kind"] == "sub" and len(parts) == 5:
        parsed.update(position_1based=int(parts[2]), ref=parts[3], alt=parts[4])
        return parsed
    if parsed["kind"] == "del" and len(parts) == 4:
        parsed.update(position_1based=int(parts[2]), ref=parts[3])
        return parsed
    raise ValueError(f"Invalid variant ID: {value}")


def weighted_contiguous_shards(
    records: list[dict], requested_shards: int
) -> list[list[dict]]:
    """Split ordered transcripts into near-equal, deterministic variant loads."""
    if requested_shards <= 0:
        raise ValueError("requested_shards must be positive")
    if not records:
        return []
    shard_count = min(requested_shards, len(records))
    shards: list[list[dict]] = []
    start = 0
    remaining_weight = sum(int(row["variant_count"]) for row in records)
    for shard_index in range(shard_count):
        remaining_shards = shard_count - shard_index
        if remaining_shards == 1:
            shards.append(records[start:])
            break
        target = remaining_weight / remaining_shards
        shard: list[dict] = []
        shard_weight = 0
        max_end = len(records) - (remaining_shards - 1)
        while start < max_end:
            next_weight = int(records[start]["variant_count"])
            if shard and shard_weight + next_weight > target:
                break
            shard.append(records[start])
            shard_weight += next_weight
            start += 1
        if not shard:
            shard.append(records[start])
            shard_weight += int(records[start]["variant_count"])
            start += 1
        shards.append(shard)
        remaining_weight -= shard_weight
    return shards


def write_json(path: str | Path, value: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def read_json(path: str | Path) -> dict:
    with Path(path).open() as handle:
        return json.load(handle)


def write_dict_rows(
    path: str | Path,
    fieldnames: list[str],
    rows: Iterable[dict],
    compresslevel: int = 1,
) -> int:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open_text(path, "wt", compresslevel=compresslevel) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
            count += 1
    return count
