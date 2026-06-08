#!/bin/bash
# submit_validate_ribonn_figures.sh - SLURM job for RiboNN figure SNV validation.
#
# Runs:
#   python validate_ribonn_figure_variants.py --run-predict
#
# Usage:
#   sbatch submit_validate_ribonn_figures.sh
#   sbatch submit_validate_ribonn_figures.sh --transcript-id QARS=ENST...

#SBATCH --job-name=ribonn_figval
#SBATCH --output=logs/ribonn_figval_%j.out
#SBATCH --error=logs/ribonn_figval_%j.err
#SBATCH --partition=ga100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00

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

echo "=== RiboNN figure-variant validation ==="
echo "Host:              $(hostname)"
echo "Repo:              ${REPO_ROOT}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"
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
print('GPU:', torch.cuda.get_device_name(0))
"

mkdir -p logs data results/human

echo ""
echo "=== Running requested figure-variant validation ==="
python validate_ribonn_figure_variants.py \
    --run-predict \
    --fasta "${FASTA}" \
    --gtf "${GTF}" \
    --batch-size 1024 \
    --num-workers "${SLURM_CPUS_PER_TASK:-4}" \
    "$@"

echo ""
echo "=== Done ==="
echo "Input:   data/ribonn_figure_variants_input.txt"
echo "Output:  results/human/ribonn_figure_variants_output.txt"
echo "Summary: results/human/ribonn_figure_variants_summary.txt"
echo "Plot:    results/human/ribonn_figure_variants_direction_plot.png"
