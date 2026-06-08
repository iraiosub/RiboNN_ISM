#!/bin/bash
# submit_scn2a_utr5_window15.sh
#
# Submit the SCN2A 5'UTR 15-nt sliding-window deletion experiment.
#
# This is a separate side experiment:
#   Job 1 (CPU): prepare reference plus every 15-nt UTR5 window deletion
#   Job 2 (GPU): run RiboNN predictions on those deleted-window transcripts
#   Job 3 (CPU): build the deletion-window null and compare the final -15..-1 window
#
# Usage:
#   bash submit_scn2a_utr5_window15.sh [window_size] [stride]
#
# Defaults:
#   window_size=15
#   stride=1

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/environment.yml" ]]; then
    REPO_ROOT="${SCRIPT_DIR}"
elif [[ -f "${PWD}/environment.yml" ]]; then
    REPO_ROOT="${PWD}"
else
    echo "[ERROR] Run this script from the RiboNN repo root." >&2
    exit 1
fi

WINDOW_SIZE="${1:-15}"
STRIDE="${2:-1}"

FASTA="/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref/human/GRCh38.primary_assembly.genome.fa"
GTF="/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref/human/gencode.v44.primary_assembly.annotation.longest_cds_transcripts.gtf.gz"

INPUT="data/scn2a_utr5_delwin${WINDOW_SIZE}_input.txt"
OUTPUT="results/human/scn2a_utr5_delwin${WINDOW_SIZE}_prediction_output.txt"
OUTDIR="plots_ism_scn2a_window${WINDOW_SIZE}"

mkdir -p "${REPO_ROOT}/logs"

echo "=== SCN2A UTR5 sliding-window deletion submission ==="
echo "Repo root  : ${REPO_ROOT}"
echo "Window size: ${WINDOW_SIZE} nt"
echo "Stride     : ${STRIDE} nt"
echo "Input      : ${INPUT}"
echo "Output     : ${OUTPUT}"
echo "Outdir     : ${OUTDIR}"
echo ""

JOB1=$(sbatch --parsable \
    --job-name=scn2a_win${WINDOW_SIZE}_1prep \
    --output="${REPO_ROOT}/logs/scn2a_win${WINDOW_SIZE}_1_prepare_%j.out" \
    --error="${REPO_ROOT}/logs/scn2a_win${WINDOW_SIZE}_1_prepare_%j.err" \
    --partition=ncpu \
    --cpus-per-task=4 \
    --mem=16G \
    --time=01:00:00 \
    --chdir="${REPO_ROOT}" \
    << EOF
#!/bin/bash
set -eo pipefail
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate ribonn
export PYTHONWARNINGS="\${PYTHONWARNINGS:+\${PYTHONWARNINGS},}ignore:pkg_resources is deprecated as an API:UserWarning"

echo "Host: \$(hostname)"
mkdir -p logs data results/human "${OUTDIR}"

for ref in "${FASTA}" "${GTF}"; do
    if [[ ! -r "\${ref}" ]]; then
        echo "[ERROR] Reference file is not readable on \$(hostname): \${ref}" >&2
        exit 1
    fi
done

echo "=== Step 1: Preparing ${WINDOW_SIZE}-nt sliding-window deletions ==="
python prepare_scn2a_utr5_window_deletions.py --fasta "${FASTA}" --gtf "${GTF}" --window-size "${WINDOW_SIZE}" --stride "${STRIDE}" --truncate-utr3 --output "${INPUT}"
EOF
)
echo "Job 1 submitted (prepare): ${JOB1}"

JOB2=$(sbatch --parsable \
    --job-name=scn2a_win${WINDOW_SIZE}_2pred \
    --output="${REPO_ROOT}/logs/scn2a_win${WINDOW_SIZE}_2_predict_%j.out" \
    --error="${REPO_ROOT}/logs/scn2a_win${WINDOW_SIZE}_2_predict_%j.err" \
    --partition=ga100 \
    --gres=gpu:1 \
    --cpus-per-task=4 \
    --mem=32G \
    --time=04:00:00 \
    --dependency=afterok:${JOB1} \
    --chdir="${REPO_ROOT}" \
    << EOF
#!/bin/bash
set -eo pipefail
module load CUDA/12.1.1 2>/dev/null || true
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate ribonn
export LD_LIBRARY_PATH="\${CONDA_PREFIX}/lib:\${LD_LIBRARY_PATH}"
unset PYTORCH_CUDA_ALLOC_CONF
export PYTHONWARNINGS="\${PYTHONWARNINGS:+\${PYTHONWARNINGS},}ignore:pkg_resources is deprecated as an API:UserWarning"

echo "Host: \$(hostname)"
nvidia-smi | head -20

python -c "import torch, sys; print('torch:', torch.__version__, ' cuda:', torch.cuda.is_available()); sys.exit('[ERROR] No CUDA. Rebuild env with: bash setup_hpc.sh --recreate') if not torch.cuda.is_available() else print('GPU:', torch.cuda.get_device_name(0))"

echo "=== Step 2: Running RiboNN predictions ==="
python run_ribonn_predict.py --input "${INPUT}" --output "${OUTPUT}" --top-k 5 --batch-size 1024 --num-workers "\${SLURM_CPUS_PER_TASK:-4}"
EOF
)
echo "Job 2 submitted (predict): ${JOB2} [afterok:${JOB1}]"

JOB3=$(sbatch --parsable \
    --job-name=scn2a_win${WINDOW_SIZE}_3null \
    --output="${REPO_ROOT}/logs/scn2a_win${WINDOW_SIZE}_3_null_%j.out" \
    --error="${REPO_ROOT}/logs/scn2a_win${WINDOW_SIZE}_3_null_%j.err" \
    --partition=ncpu \
    --cpus-per-task=2 \
    --mem=8G \
    --time=00:30:00 \
    --dependency=afterok:${JOB2} \
    --chdir="${REPO_ROOT}" \
    << EOF
#!/bin/bash
set -eo pipefail
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate ribonn
export PYTHONWARNINGS="\${PYTHONWARNINGS:+\${PYTHONWARNINGS},}ignore:pkg_resources is deprecated as an API:UserWarning"

echo "Host: \$(hostname)"
echo "=== Step 3: Window-deletion null test ==="
python analyze_scn2a_utr5_window_null.py --input "${OUTPUT}" --window-size "${WINDOW_SIZE}" --target-start-offset "-${WINDOW_SIZE}" --target-end-offset -1 --outdir "${OUTDIR}"
EOF
)
echo "Job 3 submitted (null): ${JOB3} [afterok:${JOB2}]"

echo ""
echo "Pipeline submitted. Monitor with:"
echo "  squeue -j ${JOB1},${JOB2},${JOB3}"
echo "  sacct  -j ${JOB1},${JOB2},${JOB3} --format=JobID,JobName,State,ExitCode,Elapsed,Reason"
echo ""
echo "Outputs:"
echo "  ${INPUT}"
echo "  ${OUTPUT}"
echo "  ${OUTDIR}/window${WINDOW_SIZE}_deletion_null_summary.tsv"
echo "  ${OUTDIR}/window${WINDOW_SIZE}_deletion_null_plot.png"
