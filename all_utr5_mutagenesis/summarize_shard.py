#!/usr/bin/env python3
"""Collapse one variant-score shard to one row per transcript 5'UTR position."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from workflow_common import open_text, parse_variant_id


OUTPUT_FIELDS = [
    "tx_index",
    "transcript_id",
    "gene_id",
    "gene_name",
    "transcript_name",
    "utr5_position_1based",
    "offset_from_cds_start",
    "ref_base",
    "wt_mean_predicted_TE",
    "substitution_mean_predicted_TE",
    "substitution_te_change_mean",
    "substitution_te_change_min",
    "substitution_te_change_max",
    "substitution_direction_vs_wt",
    "substitution_effect_pattern",
    "deletion_mean_predicted_TE",
    "deletion_te_change",
    "deletion_direction_vs_wt",
    "n_substitutions",
]


def direction(delta: float, tolerance: float) -> str:
    if delta > tolerance:
        return "increase"
    if delta < -tolerance:
        return "decrease"
    return "unchanged"


def substitution_pattern(deltas: list[float], tolerance: float) -> str:
    labels = {direction(delta, tolerance) for delta in deltas}
    if labels == {"increase"}:
        return "all_increase"
    if labels == {"decrease"}:
        return "all_decrease"
    if labels == {"unchanged"}:
        return "all_unchanged"
    return "mixed"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--scores", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--direction-tolerance", type=float, default=1e-9)
    return parser.parse_args()


def main():
    args = parse_args()
    effects: dict[int, dict] = {}
    score_rows = 0
    with open_text(args.scores) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"variant_id", "mean_predicted_TE"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Score table must contain {sorted(required)}")
        for row in reader:
            score_rows += 1
            parsed = parse_variant_id(row["variant_id"])
            tx_index = parsed["tx_index"]
            score = float(row["mean_predicted_TE"])
            tx_effects = effects.setdefault(
                tx_index, {"wt": None, "sub": {}, "del": {}}
            )
            if parsed["kind"] == "wt":
                tx_effects["wt"] = score
            elif parsed["kind"] == "sub":
                tx_effects["sub"].setdefault(parsed["position_1based"], []).append(
                    score
                )
            elif parsed["kind"] == "del":
                tx_effects["del"][parsed["position_1based"]] = score

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    position_rows = 0
    transcript_rows = 0
    with open_text(args.catalog) as catalog_handle, open_text(
        output, "wt", compresslevel=1
    ) as output_handle:
        reader = csv.DictReader(catalog_handle, delimiter="\t")
        writer = csv.DictWriter(
            output_handle,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for tx in reader:
            transcript_rows += 1
            tx_index = int(tx["tx_index"])
            utr5_size = int(tx["utr5_size"])
            sequence = tx["tx_sequence"]
            if tx_index not in effects or effects[tx_index]["wt"] is None:
                raise ValueError(f"Missing WT prediction for tx_index={tx_index}")
            wt = float(effects[tx_index]["wt"])

            for pos1, ref in enumerate(sequence[:utr5_size], start=1):
                substitution_scores = effects[tx_index]["sub"].get(pos1, [])
                deletion_score = effects[tx_index]["del"].get(pos1)
                if len(substitution_scores) != 3:
                    raise ValueError(
                        f"tx_index={tx_index}, position={pos1}: expected 3 "
                        f"substitutions, found {len(substitution_scores)}"
                    )
                if deletion_score is None:
                    raise ValueError(
                        f"tx_index={tx_index}, position={pos1}: missing deletion"
                    )
                substitution_deltas = [score - wt for score in substitution_scores]
                substitution_mean = sum(substitution_scores) / len(
                    substitution_scores
                )
                substitution_delta_mean = sum(substitution_deltas) / len(
                    substitution_deltas
                )
                deletion_delta = deletion_score - wt
                writer.writerow(
                    {
                        "tx_index": tx_index,
                        "transcript_id": tx["transcript_id"],
                        "gene_id": tx["gene_id"],
                        "gene_name": tx["gene_name"],
                        "transcript_name": tx["transcript_name"],
                        "utr5_position_1based": pos1,
                        "offset_from_cds_start": pos1 - utr5_size - 1,
                        "ref_base": ref,
                        "wt_mean_predicted_TE": f"{wt:.10g}",
                        "substitution_mean_predicted_TE": f"{substitution_mean:.10g}",
                        "substitution_te_change_mean": (
                            f"{substitution_delta_mean:.10g}"
                        ),
                        "substitution_te_change_min": (
                            f"{min(substitution_deltas):.10g}"
                        ),
                        "substitution_te_change_max": (
                            f"{max(substitution_deltas):.10g}"
                        ),
                        "substitution_direction_vs_wt": direction(
                            substitution_delta_mean, args.direction_tolerance
                        ),
                        "substitution_effect_pattern": substitution_pattern(
                            substitution_deltas, args.direction_tolerance
                        ),
                        "deletion_mean_predicted_TE": f"{deletion_score:.10g}",
                        "deletion_te_change": f"{deletion_delta:.10g}",
                        "deletion_direction_vs_wt": direction(
                            deletion_delta, args.direction_tolerance
                        ),
                        "n_substitutions": len(substitution_scores),
                    }
                )
                position_rows += 1

    print(f"Variant score rows : {score_rows:,}")
    print(f"Transcripts        : {transcript_rows:,}")
    print(f"Position rows      : {position_rows:,}")
    print(f"Position summary   : {output}")


if __name__ == "__main__":
    main()

