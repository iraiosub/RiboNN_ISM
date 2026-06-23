#!/usr/bin/env python3
"""Predict one shard and retain only one mean predicted-TE score per variant."""

from __future__ import annotations

import argparse
import gc
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", choices=("human", "mouse"), default="human")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--repo-root", default=None)
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--num-workers", type=int, default=4)
    return parser.parse_args()


def main():
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    repo_root = (
        Path(args.repo_root).resolve()
        if args.repo_root
        else script_dir.parent.resolve()
    )
    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    os.chdir(repo_root)
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    import torch

    from run_ribonn_predict import download_weights, load_run_df
    from src.data import RiboNNDataModule
    from src.predict import predict_using_models_trained_in_one_fold
    from src.utils.helpers import extract_config

    models_dir = repo_root / "models" / args.species
    download_weights(models_dir)
    run_df = load_run_df(args.species, models_dir)
    if run_df.empty:
        raise ValueError("No model runs were found")

    config = extract_config(run_df, run_df.run_id.iloc[0])
    config["species"] = args.species
    config["max_utr5_len"] = 1_381
    config["max_cds_utr3_len"] = 11_937
    config["tx_info_path"] = str(input_path)
    config["num_workers"] = args.num_workers
    config["test_batch_size"] = args.batch_size
    config["remove_extreme_txs"] = False
    config["target_column_pattern"] = None
    dm = RiboNNDataModule(config)

    if dm.df.empty:
        raise ValueError("RiboNN removed every input variant")
    if dm.df["tx_id"].duplicated().any():
        raise ValueError("Variant IDs are not unique within the shard")

    accumulated = np.zeros(len(dm.df), dtype=np.float64)
    fold_count = 0
    folds = np.sort(run_df["params.test_fold"].unique())
    print(
        f"Predicting {len(dm.df):,} variants across {len(folds)} folds; "
        f"top_k={args.top_k}, batch_size={args.batch_size}"
    )
    for test_fold in folds:
        fold_string = str(test_fold)
        sub_run_df = run_df.query(
            "`params.test_fold` == @fold_string or `params.test_fold` == @test_fold"
        ).reset_index(drop=True)
        prediction_df = predict_using_models_trained_in_one_fold(
            sub_run_df,
            config,
            dm,
            top_k_models_to_use=args.top_k,
        )
        predicted_columns = [
            column
            for column in prediction_df.columns
            if column.startswith("predicted_")
        ]
        if not predicted_columns:
            raise ValueError(f"Fold {test_fold} produced no predicted_* columns")
        accumulated += prediction_df[predicted_columns].mean(axis=1).to_numpy()
        fold_count += 1
        del prediction_df
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        print(f"Completed test fold {test_fold}")

    mean_scores = accumulated / fold_count
    output = pd.DataFrame(
        {
            "variant_id": dm.df["tx_id"].astype(str).to_numpy(),
            "mean_predicted_TE": mean_scores,
        }
    )
    output.to_csv(
        output_path,
        sep="\t",
        index=False,
        float_format="%.10g",
        compression={"method": "gzip", "compresslevel": 1}
        if output_path.suffix == ".gz"
        else None,
    )
    print(f"Lean variant scores written: {output_path} ({len(output):,} rows)")


if __name__ == "__main__":
    main()

