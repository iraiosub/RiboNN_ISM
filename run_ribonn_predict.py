#!/usr/bin/env python3
"""
run_ribonn_predict.py – Run RiboNN predictions on an ISM input TSV.

Bypasses the main.py hardcoded filename by calling predict_using_nested_cross_validation_models
directly. Writes results in the same format as `make predict_human` or
`make predict_mouse`.

Models are downloaded automatically from Zenodo if models/<species>/ is not present.
(Same logic as the Makefile.)

Usage:
    python run_ribonn_predict.py [--species human] [--input data/prediction_input.txt]
    python run_ribonn_predict.py --species mouse --input data/mouse_scn2a_prediction_input.txt
    python run_ribonn_predict.py --download-only   # just download weights and exit
"""

import argparse
import os
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd

ZENODO_URL  = "https://zenodo.org/records/17258709/files/weights.zip"
DEFAULT_IN  = Path("data/prediction_input.txt")


def weights_ready(models_dir: Path) -> bool:
    return (
        (models_dir / "runs.csv").is_file()
        and any(models_dir.glob("*/state_dict.pth"))
    )


def download_weights(models_dir: Path):
    if weights_ready(models_dir):
        print(f"Model weights already present in {models_dir}/")
        return

    print(f"Downloading model weights from {ZENODO_URL} ...")
    import urllib.request
    Path("tmp").mkdir(parents=True, exist_ok=True)
    Path("models").mkdir(parents=True, exist_ok=True)
    zip_path = Path("tmp/weights.zip")
    urllib.request.urlretrieve(ZENODO_URL, zip_path)
    print("Extracting weights ...")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall("models")
    zip_path.unlink()
    if not weights_ready(models_dir):
        raise FileNotFoundError(
            f"Downloaded weights but could not find runs.csv and state_dict.pth files in {models_dir}."
        )
    print(f"Weights extracted to {models_dir}/")


def load_run_df(species: str, models_dir: Path):
    """Load the run_df that maps fold/model indices to checkpoint paths."""
    runs_csv = models_dir / "runs.csv"
    if runs_csv.is_file():
        return pd.read_csv(runs_csv)

    try:
        from src.predict import get_run_df
        return get_run_df(species)
    except ImportError:
        pass
    try:
        from src.utils import get_run_df
        return get_run_df(species)
    except (ImportError, AttributeError):
        pass
    # Fallback: build run_df by scanning models/<species>/
    rows = []
    for model_dir in sorted(models_dir.glob("fold_*")):
        fold = int(model_dir.name.split("_")[1])
        for ckpt in sorted(model_dir.glob("*.pt")):
            rows.append({"fold": fold, "model_path": str(ckpt)})
    if not rows:
        raise FileNotFoundError(
            f"No model checkpoints found in {models_dir}. "
            f"Expected models/{species}/runs.csv and run_id/state_dict.pth files. "
            "Run with --download-only first."
        )
    return pd.DataFrame(rows)


def aggregate_predictions(df: pd.DataFrame) -> pd.DataFrame:
    if "mean_predicted_TE" in df.columns and not df["tx_id"].duplicated().any():
        return df

    predicted_cols = [col for col in df.columns if col.startswith("predicted_")]
    if not predicted_cols:
        raise ValueError("Prediction output has no predicted_* columns to aggregate.")

    group_cols = [
        col
        for col in ("tx_id", "utr5_sequence", "cds_sequence", "utr3_sequence")
        if col in df.columns
    ]
    if "tx_id" not in group_cols:
        raise ValueError("Prediction output is missing tx_id.")

    aggregated = df.groupby(group_cols, as_index=False)[predicted_cols].mean()
    aggregated["mean_predicted_TE"] = aggregated[predicted_cols].mean(axis=1)
    return aggregated


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--species",       choices=("human", "mouse"), default="human", help="RiboNN model species (default: human)")
    parser.add_argument("--input",         default=str(DEFAULT_IN),  help=f"Prediction input TSV (default: {DEFAULT_IN})")
    parser.add_argument("--output",        default=None, help="Output file (default: results/<species>/prediction_output.txt)")
    parser.add_argument("--top-k",         type=int, default=5,       help="Top-k models per fold (default: 5)")
    parser.add_argument("--batch-size",    type=int, default=1024)
    parser.add_argument("--num-workers",   type=int, default=4)
    parser.add_argument("--download-only", action="store_true",       help="Download weights and exit")
    args = parser.parse_args()

    models_dir = Path("models") / args.species
    output_path = Path(args.output) if args.output else Path("results") / args.species / "prediction_output.txt"

    download_weights(models_dir)

    if args.download_only:
        print("Weights ready. Exiting.")
        return

    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"[ERROR] Input file not found: {input_path}")

    # Ensure we're running from the repo root so relative imports work
    repo_root = Path(__file__).parent.resolve()
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from src.predict import predict_using_nested_cross_validation_models

    run_df = load_run_df(args.species, models_dir)

    print(f"Running RiboNN predictions on {input_path} ...")
    print(f"  species={args.species}")
    print(f"  top_k={args.top_k}  batch_size={args.batch_size}  num_workers={args.num_workers}")

    raw_results_df = predict_using_nested_cross_validation_models(
        input_path=str(input_path),
        species=args.species,
        run_df=run_df,
        top_k_models_to_use=args.top_k,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
    )
    results_df = aggregate_predictions(raw_results_df)

    out_path = output_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(out_path, sep="\t", index=False)
    print(f"Results written to {out_path}  ({len(results_df)} rows)")


if __name__ == "__main__":
    main()
