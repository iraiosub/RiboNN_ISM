#!/bin/bash
# submit_ism_scn2a.sh – SLURM job: SCN2A 5'UTR ISM with RiboNN
#
# Runs the full pipeline:
#   1. prepare_scn2a_ism.py   – extract sequences, generate variants
#   2. run_ribonn_predict.py  – RiboNN predictions (GPU)
#   3. plot_te_changes.py     – generate plots
#
# Usage:
#   sbatch submit_ism_scn2a.sh [upstream_bases] [max_deletion]
#
# Defaults:
#   upstream_bases=15   (mutate/delete the 15 bases immediately before AUG)
#   max_deletion=15     (growing deletions up to 15 bp)

#SBATCH --job-name=ribonn_ism_scn2a
#SBATCH --output=logs/ribonn_ism_scn2a_%j.out
#SBATCH --error=logs/ribonn_ism_scn2a_%j.err
#SBATCH --partition=ga100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/environment.yml" ]]; then
    REPO_ROOT="${SCRIPT_DIR}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]] && [[ -f "${SLURM_SUBMIT_DIR}/environment.yml" ]]; then
    REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
elif [[ -f "${PWD}/environment.yml" ]]; then
    REPO_ROOT="${PWD}"
else
    echo "[ERROR] Could not find environment.yml. Submit this job from the RiboNN repo root." >&2
    exit 1
fi

cd "${REPO_ROOT}"

# UPSTREAM=9999 → prepare_scn2a_ism.py auto-caps to the actual UTR length,
# giving a full saturation scan of every base.  Pass a smaller number to restrict.
UPSTREAM="${1:-9999}"
MAX_DEL="${2:-15}"

FASTA="/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref/human/GRCh38.primary_assembly.genome.fa"
GTF="/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref/human/gencode.v44.primary_assembly.annotation.longest_cds_transcripts.gtf.gz"

# ── environment ──────────────────────────────────────────────────────────────
module load CUDA/12.1.1 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ribonn
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

# PyTorch 1.13.1 does not support newer allocator options such as
# expandable_segments; inherited values can make CUDA init fail.
unset PYTORCH_CUDA_ALLOC_CONF
export PYTHONWARNINGS="${PYTHONWARNINGS:+${PYTHONWARNINGS},}ignore:pkg_resources is deprecated as an API:UserWarning"

echo "=== RiboNN ISM: SCN2A 5'UTR ==="
echo "Host:              $(hostname)"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"
echo "Upstream bases:    ${UPSTREAM}"
echo "Max deletion:      ${MAX_DEL}"
nvidia-smi | head -20

# Quick torch check
python -c "
import importlib.util
import torch, sys
print('python:', sys.executable)
print('torch:', torch.__version__)
print('torch cuda build:', torch.version.cuda)
print('cuda available:', torch.cuda.is_available())
if importlib.util.find_spec('pkg_resources') is None:
    sys.exit('[ERROR] Missing pkg_resources. Run: conda install -n ribonn -c conda-forge \"setuptools<81\"')
if not torch.cuda.is_available():
    sys.exit('[ERROR] PyTorch cannot access CUDA. Run bash setup_hpc.sh --recreate to rebuild the CUDA-compatible env.')
if torch.cuda.is_available():
    print('GPU:', torch.cuda.get_device_name(0))
"

mkdir -p logs data results/human plots_ism_scn2a

# ── step 1: generate variants ─────────────────────────────────────────────────
echo ""
echo "=== Step 1: Generating variant sequences ==="
python prepare_scn2a_ism.py \
    --fasta "${FASTA}" \
    --gtf   "${GTF}" \
    --upstream-bases "${UPSTREAM}" \
    --max-deletion   "${MAX_DEL}" \
    --truncate-utr3 \
    --output data/prediction_input.txt

# ── step 2: predict ───────────────────────────────────────────────────────────
echo ""
echo "=== Step 2: Running RiboNN predictions ==="
python run_ribonn_predict.py \
    --input  data/prediction_input.txt \
    --output results/human/prediction_output.txt \
    --top-k  5 \
    --batch-size 1024 \
    --num-workers 4

# ── step 3: plot ──────────────────────────────────────────────────────────────
echo ""
echo "=== Step 3: Generating plots ==="
python plot_te_changes.py \
    --input          results/human/prediction_output.txt \
    --outdir         plots_ism_scn2a \
    --upstream-bases "${UPSTREAM}"

echo ""
echo "=== Done. Plots in plots_ism_scn2a/ ==="
