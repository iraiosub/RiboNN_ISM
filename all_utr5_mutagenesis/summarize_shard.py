#!/usr/bin/env python3
"""Collapse one variant-score shard to one row per transcript 5'UTR position."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from workflow_common import open_text, parse_variant_id


OUTPUT_FIELDS = [
    "tx_index",
    "transcript_id",
    "gene_id",
    "gene_name",
    "transcript_name",
    "orf_start_1based",
    "orf_stop_1based",
    "orf_frame",
    "orf_annotated",
    "orf_position_in_start_codon",
    "target_codon",
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

ORF_OUTPUT_FIELDS = [
    "tx_index",
    "transcript_id",
    "gene_id",
    "gene_name",
    "transcript_name",
    "orf_start_1based",
    "orf_stop_1based",
    "orf_frame",
    "orf_annotated",
    "target_codon",
    "wt_mean_predicted_TE",
    "orf_substitution_mean_predicted_TE_3nt",
    "orf_substitution_te_change_mean_3nt",
    "orf_substitution_direction_vs_wt",
    "orf_deletion_mean_predicted_TE_3nt",
    "orf_deletion_te_change_mean_3nt",
    "orf_deletion_direction_vs_wt",
    "n_positions",
    "n_orfs_same_start",
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
    parser.add_argument(
        "--orf-output",
        default=None,
        help="Optional per-ORF-start 3-position summary table",
    )
    parser.add_argument("--direction-tolerance", type=float, default=1e-9)
    return parser.parse_args()


def target_positions(row: dict, utr5_size: int) -> list[int]:
    encoded = str(row.get("target_positions_1based", "") or "").strip()
    if encoded:
        return [int(value) for value in encoded.split(";") if value]
    return list(range(1, utr5_size + 1))


def load_orf_targets(row: dict) -> list[dict]:
    raw = str(row.get("orf_targets_json", "") or "").strip()
    if not raw:
        return []
    targets = json.loads(raw)
    if not isinstance(targets, list):
        raise ValueError("orf_targets_json must encode a list")
    return targets


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
    orf_output = Path(args.orf_output) if args.orf_output else None
    if orf_output is not None:
        orf_output.parent.mkdir(parents=True, exist_ok=True)
    position_rows = 0
    orf_rows = 0
    transcript_rows = 0
    with open_text(args.catalog) as catalog_handle, open_text(
        output, "wt", compresslevel=1
    ) as output_handle:
        orf_output_handle = (
            open_text(orf_output, "wt", compresslevel=1)
            if orf_output is not None
            else None
        )
        reader = csv.DictReader(catalog_handle, delimiter="\t")
        writer = csv.DictWriter(
            output_handle,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        orf_writer = None
        if orf_output_handle is not None:
            orf_writer = csv.DictWriter(
                orf_output_handle,
                fieldnames=ORF_OUTPUT_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            orf_writer.writeheader()
        try:
            for tx in reader:
                transcript_rows += 1
                tx_index = int(tx["tx_index"])
                utr5_size = int(tx["utr5_size"])
                sequence = tx["tx_sequence"]
                positions = target_positions(tx, utr5_size)
                orf_targets = load_orf_targets(tx)
                position_to_orf: dict[int, dict] = {}
                for target in orf_targets:
                    start = int(target["orf_start_1based"])
                    for pos1 in target.get("positions_1based", []):
                        position_meta = dict(target)
                        position_meta["orf_position_in_start_codon"] = (
                            int(pos1) - start + 1
                        )
                        position_to_orf[int(pos1)] = position_meta

                if tx_index not in effects or effects[tx_index]["wt"] is None:
                    raise ValueError(f"Missing WT prediction for tx_index={tx_index}")
                wt = float(effects[tx_index]["wt"])

                position_summaries: dict[int, dict] = {}
                for pos1 in positions:
                    ref = sequence[pos1 - 1]
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
                    substitution_deltas = [
                        score - wt for score in substitution_scores
                    ]
                    substitution_mean = sum(substitution_scores) / len(
                        substitution_scores
                    )
                    substitution_delta_mean = sum(substitution_deltas) / len(
                        substitution_deltas
                    )
                    deletion_delta = deletion_score - wt
                    orf_meta = position_to_orf.get(pos1, {})
                    row_out = {
                        "tx_index": tx_index,
                        "transcript_id": tx["transcript_id"],
                        "gene_id": tx["gene_id"],
                        "gene_name": tx["gene_name"],
                        "transcript_name": tx["transcript_name"],
                        "orf_start_1based": orf_meta.get("orf_start_1based", ""),
                        "orf_stop_1based": orf_meta.get("orf_stop_1based", ""),
                        "orf_frame": orf_meta.get("orf_frame", ""),
                        "orf_annotated": orf_meta.get("annotated", ""),
                        "orf_position_in_start_codon": orf_meta.get(
                            "orf_position_in_start_codon", ""
                        ),
                        "target_codon": orf_meta.get("target_codon", ""),
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
                    writer.writerow(row_out)
                    position_summaries[pos1] = {
                        "substitution_mean": substitution_mean,
                        "substitution_delta_mean": substitution_delta_mean,
                        "deletion_score": deletion_score,
                        "deletion_delta": deletion_delta,
                    }
                    position_rows += 1

                if orf_writer is not None:
                    for target in orf_targets:
                        positions_1based = [
                            int(pos) for pos in target["positions_1based"]
                        ]
                        summaries = [
                            position_summaries[pos] for pos in positions_1based
                        ]
                        if len(summaries) != 3:
                            raise ValueError(
                                f"tx_index={tx_index}, orf_start="
                                f"{target['orf_start_1based']}: expected 3 "
                                f"position summaries, found {len(summaries)}"
                            )
                        substitution_mean_3nt = sum(
                            item["substitution_mean"] for item in summaries
                        ) / len(summaries)
                        substitution_delta_3nt = sum(
                            item["substitution_delta_mean"] for item in summaries
                        ) / len(summaries)
                        deletion_mean_3nt = sum(
                            item["deletion_score"] for item in summaries
                        ) / len(summaries)
                        deletion_delta_3nt = sum(
                            item["deletion_delta"] for item in summaries
                        ) / len(summaries)
                        orf_writer.writerow(
                            {
                                "tx_index": tx_index,
                                "transcript_id": tx["transcript_id"],
                                "gene_id": tx["gene_id"],
                                "gene_name": tx["gene_name"],
                                "transcript_name": tx["transcript_name"],
                                "orf_start_1based": target["orf_start_1based"],
                                "orf_stop_1based": target["orf_stop_1based"],
                                "orf_frame": target["orf_frame"],
                                "orf_annotated": target["annotated"],
                                "target_codon": target["target_codon"],
                                "wt_mean_predicted_TE": f"{wt:.10g}",
                                "orf_substitution_mean_predicted_TE_3nt": (
                                    f"{substitution_mean_3nt:.10g}"
                                ),
                                "orf_substitution_te_change_mean_3nt": (
                                    f"{substitution_delta_3nt:.10g}"
                                ),
                                "orf_substitution_direction_vs_wt": direction(
                                    substitution_delta_3nt,
                                    args.direction_tolerance,
                                ),
                                "orf_deletion_mean_predicted_TE_3nt": (
                                    f"{deletion_mean_3nt:.10g}"
                                ),
                                "orf_deletion_te_change_mean_3nt": (
                                    f"{deletion_delta_3nt:.10g}"
                                ),
                                "orf_deletion_direction_vs_wt": direction(
                                    deletion_delta_3nt,
                                    args.direction_tolerance,
                                ),
                                "n_positions": len(summaries),
                                "n_orfs_same_start": target.get(
                                    "n_orfs_same_start", 1
                                ),
                            }
                        )
                        orf_rows += 1
        finally:
            if orf_output_handle is not None:
                orf_output_handle.close()

    print(f"Variant score rows : {score_rows:,}")
    print(f"Transcripts        : {transcript_rows:,}")
    print(f"Position rows      : {position_rows:,}")
    print(f"Position summary   : {output}")
    if orf_output is not None:
        print(f"ORF rows           : {orf_rows:,}")
        print(f"ORF summary        : {orf_output}")


if __name__ == "__main__":
    main()
