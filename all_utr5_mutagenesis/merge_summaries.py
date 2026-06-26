#!/usr/bin/env python3
"""Merge ordered shard summaries into the final transcript-position table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from workflow_common import open_text, read_json, write_json


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def merge_ordered_shards(
    summary_dir: Path,
    final_path: Path,
    expected_shards: int,
    suffix: str,
    required: bool,
) -> tuple[int, set[str]]:
    expected_header = None
    row_count = 0
    transcript_ids: set[str] = set()
    found_any = False
    with open_text(final_path, "wt", compresslevel=1) as output_handle:
        writer = None
        for shard_id in range(expected_shards):
            shard_path = summary_dir / f"shard_{shard_id:03d}.{suffix}.tsv.gz"
            if not shard_path.is_file():
                if required:
                    raise FileNotFoundError(f"Missing shard summary: {shard_path}")
                continue
            found_any = True
            with open_text(shard_path) as shard_handle:
                reader = csv.DictReader(shard_handle, delimiter="\t")
                header = reader.fieldnames
                if expected_header is None:
                    expected_header = header
                    writer = csv.DictWriter(
                        output_handle,
                        fieldnames=expected_header,
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                elif header != expected_header:
                    raise ValueError(f"Header mismatch in {shard_path}")
                for row in reader:
                    writer.writerow(row)
                    row_count += 1
                    transcript_ids.add(row["transcript_id"])

    if not found_any and not required:
        final_path.unlink(missing_ok=True)
    return row_count, transcript_ids


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    manifest_path = outdir / "workflow_manifest.json"
    manifest = read_json(manifest_path)
    expected_shards = int(manifest["effective_shards"])
    screen_mode = manifest.get("screen_mode", "all_utr5")
    summary_dir = outdir / "summaries"
    final_dir = outdir / "final"
    final_dir.mkdir(parents=True, exist_ok=True)
    if screen_mode == "whole_atg_deletion":
        output_path = final_dir / "whole_atg_deletion_scores.tsv.gz"
        orf_output_path = final_dir / "whole_atg_deletion_orf_scores.tsv.gz"
    else:
        output_path = final_dir / "all_utr5_position_scores.tsv.gz"
        orf_output_path = final_dir / "orf_start_codon_scores.tsv.gz"

    row_count, transcript_ids = merge_ordered_shards(
        summary_dir,
        output_path,
        expected_shards,
        "positions",
        required=True,
    )

    expected_positions = int(
        manifest.get("included_target_positions", manifest["included_utr5_positions"])
    )
    if row_count != expected_positions:
        raise ValueError(
            f"Merged {row_count:,} positions, expected "
            f"{expected_positions:,}"
        )
    if len(transcript_ids) != int(manifest["included_transcripts"]):
        raise ValueError(
            f"Merged {len(transcript_ids):,} transcripts, expected "
            f"{int(manifest['included_transcripts']):,}"
        )

    expected_orfs = int(manifest.get("included_orf_targets", 0))
    orf_row_count = 0
    orf_transcript_ids: set[str] = set()
    if expected_orfs:
        orf_row_count, orf_transcript_ids = merge_ordered_shards(
            summary_dir,
            orf_output_path,
            expected_shards,
            "orfs",
            required=True,
        )
        if orf_row_count != expected_orfs:
            raise ValueError(
                f"Merged {orf_row_count:,} ORF starts, expected "
                f"{expected_orfs:,}"
            )

    completed = {
        **manifest,
        "merged_position_rows": row_count,
        "merged_transcripts": len(transcript_ids),
        "final_position_table": str(output_path),
        "merged_orf_rows": orf_row_count,
        "final_orf_table": str(orf_output_path) if expected_orfs else None,
    }
    write_json(final_dir / "run_summary.json", completed)
    print(f"Merged transcripts : {len(transcript_ids):,}")
    print(f"Merged positions   : {row_count:,}")
    print(f"Final table        : {output_path}")
    if expected_orfs:
        print(f"Merged ORF starts  : {orf_row_count:,}")
        print(f"Final ORF table    : {orf_output_path}")
    print(f"Run summary        : {final_dir / 'run_summary.json'}")


if __name__ == "__main__":
    main()
