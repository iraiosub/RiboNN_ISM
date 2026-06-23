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


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    manifest_path = outdir / "workflow_manifest.json"
    manifest = read_json(manifest_path)
    expected_shards = int(manifest["effective_shards"])
    summary_dir = outdir / "summaries"
    final_dir = outdir / "final"
    final_dir.mkdir(parents=True, exist_ok=True)
    output_path = final_dir / "all_utr5_position_scores.tsv.gz"

    expected_header = None
    row_count = 0
    transcript_ids: set[str] = set()
    with open_text(output_path, "wt", compresslevel=1) as output_handle:
        writer = None
        for shard_id in range(expected_shards):
            shard_path = summary_dir / f"shard_{shard_id:03d}.positions.tsv.gz"
            if not shard_path.is_file():
                raise FileNotFoundError(f"Missing shard summary: {shard_path}")
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

    if row_count != int(manifest["included_utr5_positions"]):
        raise ValueError(
            f"Merged {row_count:,} positions, expected "
            f"{int(manifest['included_utr5_positions']):,}"
        )
    if len(transcript_ids) != int(manifest["included_transcripts"]):
        raise ValueError(
            f"Merged {len(transcript_ids):,} transcripts, expected "
            f"{int(manifest['included_transcripts']):,}"
        )

    completed = {
        **manifest,
        "merged_position_rows": row_count,
        "merged_transcripts": len(transcript_ids),
        "final_position_table": str(output_path),
    }
    write_json(final_dir / "run_summary.json", completed)
    print(f"Merged transcripts : {len(transcript_ids):,}")
    print(f"Merged positions   : {row_count:,}")
    print(f"Final table        : {output_path}")
    print(f"Run summary        : {final_dir / 'run_summary.json'}")


if __name__ == "__main__":
    main()

