#!/bin/bash
# submit_all_ism_scn2a.sh – Submit the full SCN2A ISM pipeline as chained SLURM jobs.
#
# Submits 4 jobs with automatic dependencies:
#   Job 1 (CPU):  prepare_scn2a_ism.py   – generate variant sequences
#   Job 2 (GPU):  run_ribonn_predict.py   – RiboNN predictions      [afterok:job1]
#   Job 3 (CPU):  plot_te_changes.py      – ISM heatmaps / waterfall [afterok:job2]
#   Job 4 (CPU):  analyze_altaug_null.py  – altAUG null test         [afterok:job2]
#
# Usage:
#   bash submit_all_ism_scn2a.sh [upstream_bases] [altaug_pos1 altaug_pos2 altaug_pos3]
#
# Examples:
#   bash submit_all_ism_scn2a.sh                    # full UTR, altAUG at -8 -7 -6
#   bash submit_all_ism_scn2a.sh 9999               # same, explicit
#   bash submit_all_ism_scn2a.sh 9999 -10 -9 -8    # custom altAUG positions

set -eo pipefail

# ── repo root ─────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/environment.yml" ]]; then
    REPO_ROOT="${SCRIPT_DIR}"
elif [[ -f "${PWD}/environment.yml" ]]; then
    REPO_ROOT="${PWD}"
else
    echo "[ERROR] Run this script from the RiboNN repo root." >&2
    exit 1
fi

# ── parameters ────────────────────────────────────────────────────────────────
UPSTREAM="${1:-9999}"     # 9999 → auto-capped to full UTR length by prepare script
ALTAUG_P1="${2:--8}"
ALTAUG_P2="${3:--7}"
ALTAUG_P3="${4:--6}"

FASTA="/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref/human/GRCh38.primary_assembly.genome.fa"
GTF="/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref/human/gencode.v44.primary_assembly.annotation.longest_cds_transcripts.gtf.gz"

mkdir -p "${REPO_ROOT}/logs"

echo "=== SCN2A ISM pipeline submission ==="
echo "Repo root       : ${REPO_ROOT}"
echo "Upstream bases  : ${UPSTREAM} (capped to actual UTR length)"
echo "altAUG positions: ${ALTAUG_P1} ${ALTAUG_P2} ${ALTAUG_P3}"
echo ""

# ── Job 1: prepare variant sequences (CPU) ────────────────────────────────────
JOB1=$(sbatch --parsable \
    --job-name=ism_scn2a_1_prepare \
    --output="${REPO_ROOT}/logs/ism_scn2a_1_prepare_%j.out" \
    --error="${REPO_ROOT}/logs/ism_scn2a_1_prepare_%j.err" \
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
mkdir -p logs data results/human plots_ism_scn2a

for ref in "${FASTA}" "${GTF}"; do
    if [[ ! -r "\${ref}" ]]; then
        echo "[ERROR] Reference file is not readable on \$(hostname): \${ref}" >&2
        exit 1
    fi
done

echo "=== Step 1: Generating variant sequences ==="
python prepare_scn2a_ism.py --fasta "${FASTA}" --gtf "${GTF}" --upstream-bases "${UPSTREAM}" --max-deletion 15 --truncate-utr3 --output data/prediction_input.txt
EOF
)
echo "Job 1 submitted  (prepare)  : ${JOB1}"

# ── Job 2: RiboNN predictions (GPU) ───────────────────────────────────────────
JOB2=$(sbatch --parsable \
    --job-name=ism_scn2a_2_predict \
    --output="${REPO_ROOT}/logs/ism_scn2a_2_predict_%j.out" \
    --error="${REPO_ROOT}/logs/ism_scn2a_2_predict_%j.err" \
    --partition=ga100 \
    --gres=gpu:1 \
    --cpus-per-task=4 \
    --mem=32G \
    --time=06:00:00 \
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

python -c "
import importlib.util, torch, sys
print('torch:', torch.__version__, '  cuda:', torch.cuda.is_available())
if not torch.cuda.is_available():
    sys.exit('[ERROR] No CUDA. Rebuild env with: bash setup_hpc.sh --recreate')
print('GPU:', torch.cuda.get_device_name(0))
"

echo "=== Step 2: Running RiboNN predictions ==="
python run_ribonn_predict.py --input data/prediction_input.txt --output results/human/prediction_output.txt --top-k 5 --batch-size 1024 --num-workers "\${SLURM_CPUS_PER_TASK:-4}"
EOF
)
echo "Job 2 submitted  (predict)   : ${JOB2}  [afterok:${JOB1}]"

# ── Job 3: ISM plots (CPU) ────────────────────────────────────────────────────
JOB3=$(sbatch --parsable \
    --job-name=ism_scn2a_3_plot \
    --output="${REPO_ROOT}/logs/ism_scn2a_3_plot_%j.out" \
    --error="${REPO_ROOT}/logs/ism_scn2a_3_plot_%j.err" \
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

# Detect actual UTR length scanned from the prediction input
ACTUAL_UPSTREAM=\$(python3 -c "
import re, sys
max_off = 0
with open('data/prediction_input.txt') as fh:
    next(fh)
    for line in fh:
        tid = line.split('\t')[0]
        m = re.match(r'^([+-]\d+)_[ACGT]>[ACGT]$', tid)
        if m:
            max_off = max(max_off, abs(int(m.group(1))))
print(max_off if max_off else 15)
")
echo "Detected upstream bases: \${ACTUAL_UPSTREAM}"

echo "=== Step 3: Generating ISM plots ==="
python plot_te_changes.py --input results/human/prediction_output.txt --outdir plots_ism_scn2a --upstream-bases "\${ACTUAL_UPSTREAM}"

echo "Plots saved to plots_ism_scn2a/"
EOF
)
echo "Job 3 submitted  (plot)      : ${JOB3}  [afterok:${JOB2}]"

# ── Job 4: altAUG null test (CPU) ─────────────────────────────────────────────
JOB4=$(sbatch --parsable \
    --job-name=ism_scn2a_4_null \
    --output="${REPO_ROOT}/logs/ism_scn2a_4_null_%j.out" \
    --error="${REPO_ROOT}/logs/ism_scn2a_4_null_%j.err" \
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

echo "=== Step 4: altAUG position-matched null test ==="
python analyze_altaug_null.py --input results/human/prediction_output.txt --altaug-positions ${ALTAUG_P1} ${ALTAUG_P2} ${ALTAUG_P3} --outdir plots_ism_scn2a

echo "Summary: plots_ism_scn2a/altaug_null_summary.tsv"
echo "Plot:    plots_ism_scn2a/altaug_null_plot.png"
EOF
)
echo "Job 4 submitted  (null test) : ${JOB4}  [afterok:${JOB2}]"

echo ""
echo "Pipeline submitted. Monitor with:"
echo "  squeue -j ${JOB1},${JOB2},${JOB3},${JOB4}"
echo "  sacct  -j ${JOB1},${JOB2},${JOB3},${JOB4} --format=JobID,JobName,State,ExitCode,Elapsed,Reason"
echo ""
echo "If a dependency is never satisfied, inspect the upstream job first:"
echo "  tail -n 80 ${REPO_ROOT}/logs/ism_scn2a_1_prepare_${JOB1}.err"
echo "  tail -n 80 ${REPO_ROOT}/logs/ism_scn2a_1_prepare_${JOB1}.out"
