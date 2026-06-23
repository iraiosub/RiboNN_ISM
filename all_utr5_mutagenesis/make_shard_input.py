#!/usr/bin/env python3
"""Materialize one compressed RiboNN input shard just before prediction."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from workflow_common import BASES, open_text, variant_id


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", required=True, help="Per-shard transcript catalog")
    parser.add_argument("--output", required=True, help="Output TSV or TSV.GZ")
    parser.add_argument("--compresslevel", type=int, default=1)
    return parser.parse_args()


def target_positions(row: dict, utr5_size: int) -> list[int]:
    encoded = str(row.get("target_positions_1based", "") or "").strip()
    if encoded:
        positions = [int(value) for value in encoded.split(";") if value]
    else:
        positions = list(range(1, utr5_size + 1))
    for pos1 in positions:
        if pos1 < 1 or pos1 > utr5_size:
            raise ValueError(
                f"{row.get('transcript_id', row.get('tx_index'))}: target "
                f"position {pos1} is outside 5'UTR length {utr5_size}"
            )
    return positions


def main():
    args = parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    transcript_count = 0
    variant_count = 0
    with open_text(args.catalog) as source, open_text(
        output, "wt", compresslevel=args.compresslevel
    ) as destination:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow(["tx_id", "tx_sequence", "utr5_size", "cds_size"])

        for row in reader:
            transcript_count += 1
            tx_index = int(row["tx_index"])
            sequence = row["tx_sequence"]
            utr5_size = int(row["utr5_size"])
            cds_size = int(row["cds_size"])
            positions = target_positions(row, utr5_size)

            writer.writerow(
                [variant_id(tx_index, "wt"), sequence, utr5_size, cds_size]
            )
            variant_count += 1

            for pos1 in positions:
                pos0 = pos1 - 1
                ref = sequence[pos0]
                for alt in BASES:
                    if alt == ref:
                        continue
                    mutated = sequence[:pos0] + alt + sequence[pos0 + 1 :]
                    writer.writerow(
                        [
                            variant_id(tx_index, "sub", pos1, ref, alt),
                            mutated,
                            utr5_size,
                            cds_size,
                        ]
                    )
                    variant_count += 1

                deleted = sequence[:pos0] + sequence[pos0 + 1 :]
                writer.writerow(
                    [
                        variant_id(tx_index, "del", pos1, ref),
                        deleted,
                        utr5_size - 1,
                        cds_size,
                    ]
                )
                variant_count += 1

    print(f"Transcripts materialized: {transcript_count:,}")
    print(f"Variants materialized   : {variant_count:,}")
    print(f"Compressed input        : {output}")


if __name__ == "__main__":
    main()
