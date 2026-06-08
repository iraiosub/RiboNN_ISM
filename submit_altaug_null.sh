#!/bin/bash
# submit_altaug_null.sh – SLURM job: position-matched substitution null for altAUG positions
#
# Runs analyze_altaug_null.py on already-computed ISM predictions to test whether
# ΔTE at the altAUG positions (-8/-7/-6 by default) are outliers relative to the
# background distribution of all other scanned positions.
#
# Prerequisite: run submit_ism_scn2a.sh first (or pass --dependency below).
#
# Usage:
#   sbatch submit_altaug_null.sh [altaug_pos1 altaug_pos2 altaug_pos3]
#   sbatch --dependency=afterok:<ISM_JOB_ID> submit_altaug_null.sh
#
# Defaults:
#   altaug positions = -8 -7 -6
#   input            = results/human/prediction_output.txt
#   outdir           = plots_ism_scn2a

#SBATCH --job-name=ribonn_altaug_null
#SBATCH --output=logs/ribonn_altaug_null_%j.out
#SBATCH --error=logs/ribonn_altaug_null_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00

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

# altAUG positions: taken from positional arguments or defaults
ALTAUG_POS="${1:--8} ${2:--7} ${3:--6}"

INPUT="results/human/prediction_output.txt"
OUTDIR="plots_ism_scn2a"

# ── environment ──────────────────────────────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ribonn
export PYTHONWARNINGS="${PYTHONWARNINGS:+${PYTHONWARNINGS},}ignore:pkg_resources is deprecated as an API:UserWarning"

echo "=== RiboNN ISM: altAUG position-matched null test ==="
echo "Host:            $(hostname)"
echo "Repo:            ${REPO_ROOT}"
echo "Input:           ${INPUT}"
echo "altAUG positions: ${ALTAUG_POS}"
echo "Output dir:      ${OUTDIR}"

if [[ ! -f "${INPUT}" ]]; then
    echo "[ERROR] Prediction output not found: ${INPUT}" >&2
    echo "        Run submit_ism_scn2a.sh first." >&2
    exit 1
fi

mkdir -p logs "${OUTDIR}"

python analyze_altaug_null.py \
    --input              "${INPUT}" \
    --altaug-positions   ${ALTAUG_POS} \
    --outdir             "${OUTDIR}"

echo ""
echo "=== Done ==="
echo "Summary: ${OUTDIR}/altaug_null_summary.tsv"
echo "Plot:    ${OUTDIR}/altaug_null_plot.png"
