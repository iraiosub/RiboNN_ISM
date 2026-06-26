#!/bin/bash
# Submit transcriptome-wide 5'UTR substitutions + single-base deletions.
#
# The GPU stage is a weighted, restartable SLURM array. Each task creates only
# its own compressed variant input, predicts a lean mean-TE score, summarizes
# one row per position, and removes the large reproducible intermediates.

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

is_repo_root() {
    local candidate="$1"
    [[ -f "${candidate}/run_ribonn_predict.py" ]] &&
        [[ -f "${candidate}/all_utr5_mutagenesis/prepare_catalog.py" ]] &&
        [[ -f "${candidate}/all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh" ]]
}

find_repo_root() {
    local candidate abs git_root dir
    local candidates=()

    # If this launcher is accidentally submitted with sbatch, SLURM may copy it
    # into /tmp/slurmd before execution. In that case BASH_SOURCE points at the
    # spool copy, so prefer explicit/user submission locations when available.
    [[ -n "${RIBONN_ISM_REPO_ROOT:-}" ]] && candidates+=("${RIBONN_ISM_REPO_ROOT}")
    [[ -n "${SLURM_SUBMIT_DIR:-}" ]] && candidates+=("${SLURM_SUBMIT_DIR}")
    candidates+=("${SCRIPT_DIR}/.." "${PWD}")

    for candidate in "${candidates[@]}"; do
        [[ -d "${candidate}" ]] || continue
        abs="$(cd "${candidate}" 2>/dev/null && pwd -P)" || continue

        if is_repo_root "${abs}"; then
            echo "${abs}"
            return 0
        fi

        if command -v git >/dev/null 2>&1; then
            git_root="$(git -C "${abs}" rev-parse --show-toplevel 2>/dev/null || true)"
            if [[ -n "${git_root}" ]] && is_repo_root "${git_root}"; then
                echo "${git_root}"
                return 0
            fi
        fi

        dir="${abs}"
        while [[ "${dir}" != "/" ]]; do
            if is_repo_root "${dir}"; then
                echo "${dir}"
                return 0
            fi
            dir="$(dirname "${dir}")"
        done
    done

    return 1
}

REPO_ROOT="$(find_repo_root || true)"

usage() {
    cat << 'EOF'
Usage:
  bash all_utr5_mutagenesis/submit_all_utr5_mutagenesis.sh [options]

Input options:
  --repo-root PATH          RiboNN_ISM checkout root; normally auto-detected
  --species human|mouse
  --orf-predictions PATH    ORF prediction/master table; defaults are
                            species-specific
  --all-utr5                Mutate every retained 5'UTR base instead of ORF starts
  --transcript-fasta PATH   Full transcript FASTA; requires --gtf
  --genome-fasta PATH       Reconstruct transcripts from genome FASTA + GTF
  --gtf PATH
  --input-table PATH        Existing RiboNN-format transcript TSV

Scaling options:
  --num-shards N            Weighted array shards (default: 128)
  --max-concurrent N        Concurrent GPUs (default: 8)
  --batch-size N            RiboNN inference batch (default: 256)
  --top-k N                 Models per fold (default: 5)
  --te-column NAME|mean_predicted_TE
                            Single predicted_TE_* column to use instead of the
                            repo-standard mean_predicted_TE aggregation
                            (default: human uses predicted_TE_normal_brain_tissue;
                            mouse uses mean_predicted_TE)
  --partition NAME          GPU partition (default: ga100)
  --gpu-time HH:MM:SS       Per-array-task limit (default: 12:00:00)
  --gpu-mem SIZE            Per-array-task RAM (default: 32G)
  --outdir PATH
  --keep-intermediates      Keep shard variant inputs and lean variant scores
  --dry-run                 Print resolved settings without submitting

The default screen mutates only ORF start codons fully inside the 5'UTR:
orf_start, orf_start+1, and orf_start+2. The default sequence source is the
same species-specific genome FASTA + longest-CDS GTF used by the existing
SCN2A RiboNN analyses.
EOF
}

SPECIES="${RIBONN_SPECIES:-human}"
REF_ROOT="${RIBONN_REF_ROOT:-/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref}"
ANALYSIS_RESULTS_ROOT="${RIBONN_ANALYSIS_RESULTS_ROOT:-/camp/lab/ulej/home/shared/oscar_ira_riboloco/analysis_results}"
REPO_ROOT_OVERRIDE=""
TRANSCRIPT_FASTA=""
GENOME_FASTA=""
GTF=""
ORF_PREDICTIONS=""
ALL_UTR5=0
INPUT_TABLE=""
OUTDIR=""
NUM_SHARDS=128
MAX_CONCURRENT=8
BATCH_SIZE=256
TOP_K=5
if [[ -n "${RIBONN_TE_COLUMN+x}" ]]; then
    TE_COLUMN="${RIBONN_TE_COLUMN}"
    TE_COLUMN_SET=1
else
    TE_COLUMN=""
    TE_COLUMN_SET=0
fi
GPU_PARTITION="ga100"
GPU_TIME="12:00:00"
GPU_MEM="32G"
KEEP_INTERMEDIATES=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo-root) REPO_ROOT_OVERRIDE="$2"; shift 2 ;;
        --species) SPECIES="$2"; shift 2 ;;
        --orf-predictions) ORF_PREDICTIONS="$2"; shift 2 ;;
        --all-utr5) ALL_UTR5=1; shift ;;
        --transcript-fasta) TRANSCRIPT_FASTA="$2"; shift 2 ;;
        --genome-fasta|--fasta) GENOME_FASTA="$2"; shift 2 ;;
        --gtf) GTF="$2"; shift 2 ;;
        --input-table) INPUT_TABLE="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --num-shards) NUM_SHARDS="$2"; shift 2 ;;
        --max-concurrent) MAX_CONCURRENT="$2"; shift 2 ;;
        --batch-size) BATCH_SIZE="$2"; shift 2 ;;
        --top-k) TOP_K="$2"; shift 2 ;;
        --te-column) TE_COLUMN="$2"; TE_COLUMN_SET=1; shift 2 ;;
        --partition) GPU_PARTITION="$2"; shift 2 ;;
        --gpu-time) GPU_TIME="$2"; shift 2 ;;
        --gpu-mem) GPU_MEM="$2"; shift 2 ;;
        --keep-intermediates) KEEP_INTERMEDIATES=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "[ERROR] Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -n "${REPO_ROOT_OVERRIDE}" ]]; then
    if ! REPO_ROOT="$(cd "${REPO_ROOT_OVERRIDE}" 2>/dev/null && pwd -P)"; then
        echo "[ERROR] --repo-root does not exist or is not readable: ${REPO_ROOT_OVERRIDE}" >&2
        exit 2
    fi
fi
if [[ -z "${REPO_ROOT}" ]] || ! is_repo_root "${REPO_ROOT}"; then
    echo "[ERROR] Could not find the RiboNN_ISM repo root." >&2
    echo "        Run this launcher from the checkout with bash, or pass:" >&2
    echo "        --repo-root /path/to/RiboNN_ISM" >&2
    echo "        Detected script dir: ${SCRIPT_DIR}" >&2
    echo "        SLURM_SUBMIT_DIR : ${SLURM_SUBMIT_DIR:-<unset>}" >&2
    exit 2
fi

case "${SPECIES}" in
    human)
        DEFAULT_GENOME="${REF_ROOT}/human/GRCh38.primary_assembly.genome.fa"
        DEFAULT_GTF="${REF_ROOT}/human/gencode.v44.primary_assembly.annotation.longest_cds_transcripts.gtf.gz"
        DEFAULT_ORF_PREDICTIONS="${REF_ROOT}/human/orfs/gencode.v44.primary_assembly.annotation.longest_cds.transcript_info.orf_predictions.csv.gz"
        ;;
    mouse)
        DEFAULT_GENOME="${REF_ROOT}/mouse/GRCm39.primary_assembly.genome.fa"
        DEFAULT_GTF="${REF_ROOT}/mouse/gencode.vM33.primary_assembly.annotation.longest_cds_transcripts.gtf.gz"
        DEFAULT_ORF_PREDICTIONS="${ANALYSIS_RESULTS_ROOT}/cross_tissue.unmixing.master_table.with_below_tpm_threshold.tsv.gz"
        ;;
    *)
        echo "[ERROR] --species must be human or mouse." >&2
        exit 2
        ;;
esac

if [[ "${TE_COLUMN_SET}" -eq 0 ]]; then
    case "${SPECIES}" in
        human) TE_COLUMN="predicted_TE_normal_brain_tissue" ;;
        mouse) TE_COLUMN="" ;;
    esac
fi
case "${TE_COLUMN}" in
    mean_predicted_TE|mean_all|all|generic|none|"") TE_COLUMN="" ;;
esac

if [[ -n "${INPUT_TABLE}" && ( -n "${TRANSCRIPT_FASTA}" || -n "${GENOME_FASTA}" ) ]]; then
    echo "[ERROR] --input-table cannot be combined with a FASTA source." >&2
    exit 2
fi
if [[ -n "${TRANSCRIPT_FASTA}" && -n "${GENOME_FASTA}" ]]; then
    echo "[ERROR] Choose --transcript-fasta or --genome-fasta, not both." >&2
    exit 2
fi
if [[ "${ALL_UTR5}" -eq 1 && -n "${ORF_PREDICTIONS}" ]]; then
    echo "[ERROR] --all-utr5 cannot be combined with --orf-predictions." >&2
    exit 2
fi

GTF="${GTF:-${DEFAULT_GTF}}"
GENOME_FASTA="${GENOME_FASTA:-${DEFAULT_GENOME}}"
if [[ "${ALL_UTR5}" -eq 0 ]]; then
    ORF_PREDICTIONS="${ORF_PREDICTIONS:-${DEFAULT_ORF_PREDICTIONS}}"
fi
if [[ -n "${TE_COLUMN}" ]]; then
    TE_LABEL="${TE_COLUMN#predicted_TE_}"   # strip leading predicted_TE_ for readability
    TE_DISPLAY="${TE_COLUMN}"
else
    TE_LABEL="mean_predicted_TE"
    TE_DISPLAY="mean_predicted_TE (mean of predicted_TE_* columns)"
fi
RUN_LABEL="orf_start_codon"
DEFAULT_OUTDIR="${REPO_ROOT}/all_utr5_mutagenesis/output/${SPECIES}_orf_starts_${TE_LABEL}"
if [[ "${ALL_UTR5}" -eq 1 ]]; then
    RUN_LABEL="all_utr5"
    DEFAULT_OUTDIR="${REPO_ROOT}/all_utr5_mutagenesis/output/${SPECIES}_all_utr5_${TE_LABEL}"
fi
OUTDIR="${OUTDIR:-${DEFAULT_OUTDIR}}"

if ! [[ "${NUM_SHARDS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] --num-shards must be a positive integer." >&2
    exit 2
fi
if ! [[ "${MAX_CONCURRENT}" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] --max-concurrent must be a positive integer." >&2
    exit 2
fi

if [[ -n "${INPUT_TABLE}" ]]; then
    SOURCE_LABEL="RiboNN table: ${INPUT_TABLE}"
    PREP_SOURCE_ARGS=(--input-table "${INPUT_TABLE}")
elif [[ -n "${TRANSCRIPT_FASTA}" ]]; then
    SOURCE_LABEL="transcript FASTA: ${TRANSCRIPT_FASTA}"
    PREP_SOURCE_ARGS=(--transcript-fasta "${TRANSCRIPT_FASTA}" --gtf "${GTF}")
else
    SOURCE_LABEL="genome FASTA: ${GENOME_FASTA}"
    PREP_SOURCE_ARGS=(--genome-fasta "${GENOME_FASTA}" --gtf "${GTF}")
fi
SCREEN_LABEL="ORF start codons in 5'UTR: ${ORF_PREDICTIONS}"
if [[ "${ALL_UTR5}" -eq 0 ]]; then
    PREP_SOURCE_ARGS+=(--orf-predictions "${ORF_PREDICTIONS}")
else
    SCREEN_LABEL="all retained 5'UTR bases"
fi

echo "=== Targeted 5'UTR mutagenesis ==="
echo "Repo root       : ${REPO_ROOT}"
echo "Species         : ${SPECIES}"
echo "Source          : ${SOURCE_LABEL}"
echo "GTF             : ${GTF}"
echo "Screen          : ${SCREEN_LABEL}"
echo "Run label       : ${RUN_LABEL}"
echo "Output          : ${OUTDIR}"
echo "Array shards    : ${NUM_SHARDS}"
echo "Concurrent GPUs : ${MAX_CONCURRENT}"
echo "GPU partition   : ${GPU_PARTITION}"
echo "Batch / top-k   : ${BATCH_SIZE} / ${TOP_K}"
echo "TE score        : ${TE_DISPLAY}"
echo "Keep large files: ${KEEP_INTERMEDIATES}"

if [[ "${DRY_RUN}" -eq 1 ]]; then
    exit 0
fi

slurm_job_id() {
    # sbatch --parsable may return jobid;cluster on federated systems.
    # Dependencies and local log filenames should use the numeric job ID.
    echo "${1%%;*}"
}

LOG_DIR="${REPO_ROOT}/logs"
mkdir -p "${LOG_DIR}" "${OUTDIR}"
PREP_ARGS_FILE="${OUTDIR}/prepare_source_args.txt"
SUBMISSION_MANIFEST="${OUTDIR}/submission_manifest.tsv"
DEBUG_COMMANDS="${OUTDIR}/slurm_debug_commands.txt"
printf '%s\n' "${PREP_SOURCE_ARGS[@]}" > "${PREP_ARGS_FILE}"

cat > "${SUBMISSION_MANIFEST}" << EOF
key	value
repo_root	${REPO_ROOT}
species	${SPECIES}
source	${SOURCE_LABEL}
gtf	${GTF}
screen	${SCREEN_LABEL}
run_label	${RUN_LABEL}
outdir	${OUTDIR}
num_shards	${NUM_SHARDS}
max_concurrent	${MAX_CONCURRENT}
gpu_partition	${GPU_PARTITION}
gpu_time	${GPU_TIME}
gpu_mem	${GPU_MEM}
batch_size	${BATCH_SIZE}
top_k	${TOP_K}
te_score	${TE_DISPLAY}
keep_intermediates	${KEEP_INTERMEDIATES}
EOF

PREP_JOB_RAW=$(sbatch --parsable \
    --job-name="utr5all_1prep" \
    --output="${LOG_DIR}/utr5all_1_prepare_%j.out" \
    --error="${LOG_DIR}/utr5all_1_prepare_%j.err" \
    --partition=ncpu \
    --cpus-per-task=4 \
    --mem=32G \
    --time=06:00:00 \
    --chdir="${REPO_ROOT}" \
    << EOF
#!/bin/bash
set -eo pipefail
trap 'status=\$?; echo "[ERROR] \${SLURM_JOB_NAME:-prep} failed at line \${LINENO} with exit \${status}" >&2; exit \${status}' ERR
echo "[INFO] Started \$(date -Is)"
echo "[INFO] Host: \$(hostname)"
echo "[INFO] PWD : \$(pwd)"
echo "[INFO] Job : \${SLURM_JOB_ID:-unknown}"

source \$(conda info --base)/etc/profile.d/conda.sh
conda activate ribonn
export LD_LIBRARY_PATH="\${CONDA_PREFIX}/lib:\${LD_LIBRARY_PATH:-}"
export PYTHONWARNINGS="\${PYTHONWARNINGS:+\${PYTHONWARNINGS},}ignore:pkg_resources is deprecated as an API:UserWarning"
echo "[INFO] Python: \$(command -v python)"

mkdir -p "${OUTDIR}" "${OUTDIR}/work" "${OUTDIR}/summaries" "${OUTDIR}/final"
SOURCE_ARGS=()
while IFS= read -r arg; do
    SOURCE_ARGS+=("\${arg}")
done < "${PREP_ARGS_FILE}"

python all_utr5_mutagenesis/prepare_catalog.py \
    --species "${SPECIES}" \
    --outdir "${OUTDIR}" \
    --num-shards "${NUM_SHARDS}" \
    "\${SOURCE_ARGS[@]}"

python run_ribonn_predict.py --species "${SPECIES}" --download-only
touch "${OUTDIR}/prepare_complete"
echo "[INFO] Finished \$(date -Is)"
EOF
)
PREP_JOB="$(slurm_job_id "${PREP_JOB_RAW}")"
echo "Preparation job: ${PREP_JOB_RAW}"

PREP_STATUS_JOB_RAW=$(sbatch --parsable \
    --job-name="utr5all_1prep_status" \
    --output="${LOG_DIR}/utr5all_1_prepare_status_%j.out" \
    --error="${LOG_DIR}/utr5all_1_prepare_status_%j.err" \
    --partition=ncpu \
    --cpus-per-task=1 \
    --mem=1G \
    --time=00:10:00 \
    --dependency="afterany:${PREP_JOB}" \
    --chdir="${REPO_ROOT}" \
    << EOF
#!/bin/bash
set -eo pipefail
REPORT="${OUTDIR}/slurm_prepare_status_${PREP_JOB}.txt"
{
    echo "Preparation job status report"
    echo "generated_at: \$(date -Is)"
    echo "prep_job: ${PREP_JOB_RAW}"
    echo "outdir: ${OUTDIR}"
    echo
    echo "sacct:"
    if command -v sacct >/dev/null 2>&1; then
        sacct -j "${PREP_JOB}" --format=JobID,JobName%24,State,ExitCode,Elapsed,Reason%40 || true
    else
        echo "sacct not available on this node"
    fi
    echo
    if [[ -f "${OUTDIR}/prepare_complete" ]]; then
        echo "prepare_complete marker: present"
    else
        echo "prepare_complete marker: missing"
    fi
    if [[ -f "${OUTDIR}/workflow_manifest.json" ]]; then
        echo "workflow_manifest.json: present"
    else
        echo "workflow_manifest.json: missing"
    fi
    echo
    for log in "${LOG_DIR}/utr5all_1_prepare_${PREP_JOB}.out" "${LOG_DIR}/utr5all_1_prepare_${PREP_JOB}.err"; do
        echo "===== \${log} ====="
        if [[ -s "\${log}" ]]; then
            tail -n 200 "\${log}"
        elif [[ -f "\${log}" ]]; then
            echo "present but empty"
        else
            echo "missing"
        fi
        echo
    done
} > "\${REPORT}" 2>&1
echo "Preparation status report: \${REPORT}"
EOF
)
PREP_STATUS_JOB="$(slurm_job_id "${PREP_STATUS_JOB_RAW}")"
echo "Prep status job: ${PREP_STATUS_JOB_RAW} [afterany:${PREP_JOB}]"

ARRAY_MAX=$((NUM_SHARDS - 1))
ARRAY_JOB_RAW=$(sbatch --parsable \
    --job-name="utr5all_2pred" \
    --output="${LOG_DIR}/utr5all_2_predict_%A_%a.out" \
    --error="${LOG_DIR}/utr5all_2_predict_%A_%a.err" \
    --partition="${GPU_PARTITION}" \
    --gres=gpu:1 \
    --cpus-per-task=4 \
    --mem="${GPU_MEM}" \
    --time="${GPU_TIME}" \
    --array="0-${ARRAY_MAX}%${MAX_CONCURRENT}" \
    --dependency="afterok:${PREP_JOB}" \
    --chdir="${REPO_ROOT}" \
    << EOF
#!/bin/bash
set -eo pipefail
trap 'status=\$?; echo "[ERROR] \${SLURM_JOB_NAME:-predict} task \${SLURM_ARRAY_TASK_ID:-NA} failed at line \${LINENO} with exit \${status}" >&2; exit \${status}' ERR
echo "[INFO] Started \$(date -Is)"
echo "[INFO] Host: \$(hostname)"
echo "[INFO] PWD : \$(pwd)"
echo "[INFO] Job : \${SLURM_ARRAY_JOB_ID:-\${SLURM_JOB_ID:-unknown}} task \${SLURM_ARRAY_TASK_ID:-NA}"

module load CUDA/12.1.1 2>/dev/null || true
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate ribonn
export LD_LIBRARY_PATH="\${CONDA_PREFIX}/lib:\${LD_LIBRARY_PATH:-}"
unset PYTORCH_CUDA_ALLOC_CONF
export PYTHONWARNINGS="\${PYTHONWARNINGS:+\${PYTHONWARNINGS},}ignore:pkg_resources is deprecated as an API:UserWarning"
echo "[INFO] Python: \$(command -v python)"

SHARD_ID=\$(printf '%03d' "\${SLURM_ARRAY_TASK_ID}")
SHARD_CATALOG="${OUTDIR}/catalog/shards/shard_\${SHARD_ID}.transcripts.tsv.gz"
SHARD_WORK="${OUTDIR}/work/shard_\${SHARD_ID}"
SHARD_INPUT="\${SHARD_WORK}/variants.tsv.gz"
SHARD_SCORES="\${SHARD_WORK}/variant_mean_te.tsv.gz"
SHARD_SUMMARY="${OUTDIR}/summaries/shard_\${SHARD_ID}.positions.tsv.gz"
SHARD_ORF_SUMMARY="${OUTDIR}/summaries/shard_\${SHARD_ID}.orfs.tsv.gz"
DONE_FILE="\${SHARD_WORK}/complete.${RUN_LABEL}"

if [[ ! -f "\${SHARD_CATALOG}" ]]; then
    echo "No effective shard \${SHARD_ID}; exiting cleanly."
    exit 0
fi
if [[ -f "\${DONE_FILE}" && -s "\${SHARD_SUMMARY}" ]]; then
    echo "Shard \${SHARD_ID} already complete; skipping."
    exit 0
fi

mkdir -p "\${SHARD_WORK}" "${OUTDIR}/summaries"
python -c "import torch, sys; print('torch:', torch.__version__, 'cuda:', torch.cuda.is_available()); sys.exit('[ERROR] CUDA unavailable') if not torch.cuda.is_available() else print('GPU:', torch.cuda.get_device_name(0))"

PREDICT_TE_ARGS=()
if [[ -n "${TE_COLUMN}" ]]; then
    PREDICT_TE_ARGS=(--te-column "${TE_COLUMN}")
fi

python all_utr5_mutagenesis/make_shard_input.py \
    --catalog "\${SHARD_CATALOG}" \
    --output "\${SHARD_INPUT}"

python all_utr5_mutagenesis/predict_mean_te.py \
    --repo-root "${REPO_ROOT}" \
    --species "${SPECIES}" \
    --input "\${SHARD_INPUT}" \
    --output "\${SHARD_SCORES}" \
    --top-k "${TOP_K}" \
    --batch-size "${BATCH_SIZE}" \
    --num-workers "\${SLURM_CPUS_PER_TASK:-4}" \
    "\${PREDICT_TE_ARGS[@]}"

python all_utr5_mutagenesis/summarize_shard.py \
    --catalog "\${SHARD_CATALOG}" \
    --scores "\${SHARD_SCORES}" \
    --output "\${SHARD_SUMMARY}" \
    --orf-output "\${SHARD_ORF_SUMMARY}"

touch "\${DONE_FILE}"
if [[ "${KEEP_INTERMEDIATES}" -eq 0 ]]; then
    rm -f "\${SHARD_INPUT}" "\${SHARD_SCORES}"
fi
echo "[INFO] Finished \$(date -Is)"
EOF
)
ARRAY_JOB="$(slurm_job_id "${ARRAY_JOB_RAW}")"
echo "GPU array job : ${ARRAY_JOB} [afterok:${PREP_JOB}]"

MERGE_JOB_RAW=$(sbatch --parsable \
    --job-name="utr5all_3merge" \
    --output="${LOG_DIR}/utr5all_3_merge_%j.out" \
    --error="${LOG_DIR}/utr5all_3_merge_%j.err" \
    --partition=ncpu \
    --cpus-per-task=2 \
    --mem=8G \
    --time=02:00:00 \
    --dependency="afterok:${ARRAY_JOB}" \
    --chdir="${REPO_ROOT}" \
    << EOF
#!/bin/bash
set -eo pipefail
trap 'status=\$?; echo "[ERROR] \${SLURM_JOB_NAME:-merge} failed at line \${LINENO} with exit \${status}" >&2; exit \${status}' ERR
echo "[INFO] Started \$(date -Is)"
echo "[INFO] Host: \$(hostname)"
echo "[INFO] PWD : \$(pwd)"
echo "[INFO] Job : \${SLURM_JOB_ID:-unknown}"
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate ribonn
python all_utr5_mutagenesis/merge_summaries.py --outdir "${OUTDIR}"
echo "[INFO] Finished \$(date -Is)"
EOF
)
MERGE_JOB="$(slurm_job_id "${MERGE_JOB_RAW}")"
echo "Merge job     : ${MERGE_JOB} [afterok:${ARRAY_JOB}]"

cat >> "${SUBMISSION_MANIFEST}" << EOF
prep_job	${PREP_JOB_RAW}
prep_status_job	${PREP_STATUS_JOB_RAW}
array_job	${ARRAY_JOB_RAW}
merge_job	${MERGE_JOB_RAW}
prepare_stdout	${LOG_DIR}/utr5all_1_prepare_${PREP_JOB}.out
prepare_stderr	${LOG_DIR}/utr5all_1_prepare_${PREP_JOB}.err
prepare_status_report	${OUTDIR}/slurm_prepare_status_${PREP_JOB}.txt
array_stdout_glob	${LOG_DIR}/utr5all_2_predict_${ARRAY_JOB}_*.out
array_stderr_glob	${LOG_DIR}/utr5all_2_predict_${ARRAY_JOB}_*.err
merge_stdout	${LOG_DIR}/utr5all_3_merge_${MERGE_JOB}.out
merge_stderr	${LOG_DIR}/utr5all_3_merge_${MERGE_JOB}.err
EOF

cat > "${DEBUG_COMMANDS}" << EOF
squeue -j ${PREP_JOB},${PREP_STATUS_JOB},${ARRAY_JOB},${MERGE_JOB}
sacct -j ${PREP_JOB},${PREP_STATUS_JOB},${ARRAY_JOB},${MERGE_JOB} --format=JobID,JobName%24,State,ExitCode,Elapsed,Reason%40
tail -n 200 ${LOG_DIR}/utr5all_1_prepare_${PREP_JOB}.out
tail -n 200 ${LOG_DIR}/utr5all_1_prepare_${PREP_JOB}.err
cat ${OUTDIR}/slurm_prepare_status_${PREP_JOB}.txt
EOF
echo
echo "Monitor:"
echo "  squeue -j ${PREP_JOB},${PREP_STATUS_JOB},${ARRAY_JOB},${MERGE_JOB}"
echo "  sacct  -j ${PREP_JOB},${PREP_STATUS_JOB},${ARRAY_JOB},${MERGE_JOB} --format=JobID,JobName%24,State,ExitCode,Elapsed,Reason%40"
echo "Logs/debug:"
echo "  ${SUBMISSION_MANIFEST}"
echo "  ${DEBUG_COMMANDS}"
echo "  ${LOG_DIR}/utr5all_1_prepare_${PREP_JOB}.out"
echo "  ${LOG_DIR}/utr5all_1_prepare_${PREP_JOB}.err"
echo "  ${OUTDIR}/slurm_prepare_status_${PREP_JOB}.txt"
echo "Final table:"
echo "  ${OUTDIR}/final/all_utr5_position_scores.tsv.gz"
if [[ "${ALL_UTR5}" -eq 0 ]]; then
    echo "  ${OUTDIR}/final/orf_start_codon_scores.tsv.gz"
fi
