#!/bin/bash
# submit_all_ism_scn2a.sh - Submit the full 5'UTR ISM pipeline as chained SLURM jobs.
#
# Despite the historical filename, this launcher now supports:
#   - species: human or mouse
#   - genes: one or more comma-separated gene symbols, default SCN2A,SCN8A
#   - separate data/results/plot outputs for each gene/species pair
#
# Usage:
#   bash submit_all_ism_scn2a.sh [legacy_upstream] [altaug_pos1 altaug_pos2 altaug_pos3]
#   bash submit_all_ism_scn2a.sh --species mouse
#   bash submit_all_ism_scn2a.sh --species human --genes SCN2A,SCN8A --upstream-bases 9999
#   bash submit_all_ism_scn2a.sh --species mouse --fasta /path/genome.fa --gtf /path/annotation.gtf.gz
#
# Reference defaults:
#   human: <ref-root>/human/GRCh38.primary_assembly.genome.fa
#          <ref-root>/human/gencode.v44.primary_assembly.annotation.longest_cds_transcripts.gtf.gz
#   mouse: <ref-root>/mouse/GRCm39.primary_assembly.genome.fa
#          <ref-root>/mouse/gencode.vM33.primary_assembly.annotation.longest_cds_transcripts.gtf.gz
#
# Override defaults with --ref-root, --fasta, --gtf, or:
#   RIBONN_REF_ROOT=/camp/.../ref
#   RIBONN_HUMAN_FASTA=... RIBONN_HUMAN_GTF=...
#   RIBONN_MOUSE_FASTA=... RIBONN_MOUSE_GTF=...

set -eo pipefail

usage() {
    cat << 'EOF'
submit_all_ism_scn2a.sh - Submit the full 5'UTR ISM pipeline as chained SLURM jobs.

Usage:
  bash submit_all_ism_scn2a.sh [legacy_upstream] [altaug_pos1 altaug_pos2 altaug_pos3]
  bash submit_all_ism_scn2a.sh --species mouse
  bash submit_all_ism_scn2a.sh --species human --genes SCN2A,SCN8A --upstream-bases 9999
  bash submit_all_ism_scn2a.sh --species mouse --fasta /path/genome.fa --gtf /path/annotation.gtf.gz

Options:
  --species human|mouse       Model species and default ref/<species> directory.
  --genes SCN2A,SCN8A         Comma-separated genes; default SCN2A,SCN8A.
  --upstream-bases N          Bases before AUG to scan; default 9999, capped by UTR length.
  --max-deletion N            Max growing deletion size; default 15.
  --altaug-positions N ...    Offsets for null test; default -8 -7 -6.
  --ref-root PATH             Reference root containing human/ and mouse/.
  --fasta PATH --gtf PATH     Explicit reference files.
  --skip-null                 Submit prepare/predict/plot jobs only.
EOF
}

# -- repo root ----------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/environment.yml" ]]; then
    REPO_ROOT="${SCRIPT_DIR}"
elif [[ -f "${PWD}/environment.yml" ]]; then
    REPO_ROOT="${PWD}"
else
    echo "[ERROR] Run this script from the RiboNN repo root." >&2
    exit 1
fi

# -- defaults -----------------------------------------------------------------
SPECIES="${RIBONN_SPECIES:-human}"
GENES_CSV="${RIBONN_GENES:-SCN2A,SCN8A}"
UPSTREAM="${RIBONN_UPSTREAM_BASES:-9999}"
MAX_DEL="${RIBONN_MAX_DELETION:-15}"
REF_ROOT="${RIBONN_REF_ROOT:-/camp/lab/ulej/home/shared/oscar_ira_riboloco/ref}"
FASTA_OVERRIDE=""
GTF_OVERRIDE=""
SKIP_NULL=0
ALTAUG_POSITIONS=("-8" "-7" "-6")
POSITIONAL=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --species)
            SPECIES="$2"
            shift 2
            ;;
        --genes|--gene-name)
            GENES_CSV="$2"
            shift 2
            ;;
        --upstream|--upstream-bases)
            UPSTREAM="$2"
            shift 2
            ;;
        --max-deletion|--max-del|max-del)
            MAX_DEL="$2"
            shift 2
            ;;
        --altaug-positions)
            shift
            ALTAUG_POSITIONS=()
            while [[ $# -gt 0 && "$1" != --* ]]; do
                ALTAUG_POSITIONS+=("$1")
                shift
            done
            ;;
        --ref-root)
            REF_ROOT="$2"
            shift 2
            ;;
        --fasta)
            FASTA_OVERRIDE="$2"
            shift 2
            ;;
        --gtf)
            GTF_OVERRIDE="$2"
            shift 2
            ;;
        --skip-null)
            SKIP_NULL=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        human|mouse)
            SPECIES="$1"
            shift
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

# Legacy positional mode: [upstream_bases] [altaug_pos1 altaug_pos2 altaug_pos3]
if [[ ${#POSITIONAL[@]} -gt 0 ]]; then
    UPSTREAM="${POSITIONAL[0]}"
fi
if [[ ${#POSITIONAL[@]} -ge 4 ]]; then
    ALTAUG_POSITIONS=("${POSITIONAL[@]:1:3}")
fi

case "${SPECIES}" in
    human|mouse) ;;
    *)
        echo "[ERROR] --species must be 'human' or 'mouse' (got '${SPECIES}')." >&2
        exit 2
        ;;
esac

if [[ ${#ALTAUG_POSITIONS[@]} -eq 0 ]]; then
    echo "[ERROR] --altaug-positions requires at least one offset." >&2
    exit 2
fi

normalize_gene_symbol() {
    local gene="$1"
    local upper
    local lower

    upper="$(printf '%s' "${gene}" | tr '[:lower:]' '[:upper:]')"
    if [[ "${SPECIES}" == "human" ]]; then
        printf '%s\n' "${upper}"
        return
    fi

    lower="$(printf '%s' "${gene}" | tr '[:upper:]' '[:lower:]')"
    printf '%s%s\n' "$(printf '%s' "${lower:0:1}" | tr '[:lower:]' '[:upper:]')" "${lower:1}"
}

GENES_CSV="${GENES_CSV// /,}"
IFS=',' read -r -a GENES <<< "${GENES_CSV}"
if [[ ${#GENES[@]} -eq 0 ]]; then
    echo "[ERROR] At least one gene must be supplied via --genes." >&2
    exit 2
fi
NORMALIZED_GENES=()
for gene in "${GENES[@]}"; do
    if [[ -n "${gene}" ]]; then
        NORMALIZED_GENES+=("$(normalize_gene_symbol "${gene}")")
    fi
done
GENES=("${NORMALIZED_GENES[@]}")
if [[ ${#GENES[@]} -eq 0 ]]; then
    echo "[ERROR] At least one non-empty gene must be supplied via --genes." >&2
    exit 2
fi

resolve_ref() {
    local kind="$1"
    local species_upper
    local species_env
    local species_override
    local generic_env
    local generic_override
    local default_path
    local species_dir

    species_upper="$(printf '%s' "${SPECIES}" | tr '[:lower:]' '[:upper:]')"
    species_env="RIBONN_${species_upper}_${kind}"
    species_override="${!species_env:-}"
    generic_env="RIBONN_${kind}"
    generic_override="${!generic_env:-}"
    species_dir="${REF_ROOT}/${SPECIES}"

    if [[ "${kind}" == "FASTA" && -n "${FASTA_OVERRIDE}" ]]; then
        printf '%s\n' "${FASTA_OVERRIDE}"
        return
    fi
    if [[ "${kind}" == "GTF" && -n "${GTF_OVERRIDE}" ]]; then
        printf '%s\n' "${GTF_OVERRIDE}"
        return
    fi
    if [[ -n "${species_override}" ]]; then
        printf '%s\n' "${species_override}"
        return
    fi
    if [[ -n "${generic_override}" ]]; then
        printf '%s\n' "${generic_override}"
        return
    fi

    if [[ "${SPECIES}" == "human" && "${kind}" == "FASTA" ]]; then
        default_path="${species_dir}/GRCh38.primary_assembly.genome.fa"
    elif [[ "${SPECIES}" == "human" && "${kind}" == "GTF" ]]; then
        default_path="${species_dir}/gencode.v44.primary_assembly.annotation.longest_cds_transcripts.gtf.gz"
    elif [[ "${SPECIES}" == "mouse" && "${kind}" == "FASTA" ]]; then
        default_path="${species_dir}/GRCm39.primary_assembly.genome.fa"
    else
        default_path="${species_dir}/gencode.vM33.primary_assembly.annotation.longest_cds_transcripts.gtf.gz"
    fi

    if [[ -r "${default_path}" ]]; then
        printf '%s\n' "${default_path}"
        return
    fi

    shopt -s nullglob
    local candidates=()
    if [[ "${kind}" == "FASTA" ]]; then
        candidates=(
            "${species_dir}"/GRC*.primary_assembly.genome.fa
            "${species_dir}"/*.primary_assembly.genome.fa
            "${species_dir}"/*.fa
        )
    else
        candidates=(
            "${species_dir}"/gencode.v*.primary_assembly.annotation.longest_cds_transcripts.gtf.gz
            "${species_dir}"/*longest_cds_transcripts*.gtf.gz
            "${species_dir}"/*.gtf.gz
            "${species_dir}"/*.gtf
        )
    fi
    for candidate in "${candidates[@]}"; do
        if [[ -r "${candidate}" ]]; then
            printf '%s\n' "${candidate}"
            shopt -u nullglob
            return
        fi
    done
    shopt -u nullglob

    printf '%s\n' "${default_path}"
}

FASTA="$(resolve_ref FASTA)"
GTF="$(resolve_ref GTF)"

mkdir -p "${REPO_ROOT}/logs"

echo "=== RiboNN 5'UTR ISM pipeline submission ==="
echo "Repo root       : ${REPO_ROOT}"
echo "Species         : ${SPECIES}"
echo "Genes           : ${GENES[*]}"
echo "Reference FASTA : ${FASTA}"
echo "Reference GTF   : ${GTF}"
echo "Upstream bases  : ${UPSTREAM} (capped to actual UTR length)"
echo "Max deletion    : ${MAX_DEL}"
echo "altAUG positions: ${ALTAUG_POSITIONS[*]}"
if [[ "${SKIP_NULL}" -eq 1 ]]; then
    echo "Null test       : skipped"
fi
echo ""

ALL_JOBS=()

submit_gene_pipeline() {
    local gene="$1"
    local gene_lc
    local run_label
    local input
    local output
    local outdir
    local job_prefix
    local job1
    local job2
    local job3
    local job4

    gene_lc="$(printf '%s' "${gene}" | tr '[:upper:]' '[:lower:]')"
    run_label="${SPECIES}_${gene_lc}"
    input="data/${run_label}_prediction_input.txt"
    output="results/${SPECIES}/${gene_lc}_prediction_output.txt"
    outdir="plots_ism_${run_label}"
    job_prefix="ism_${SPECIES}_${gene_lc}"

    echo "--- Submitting ${gene} (${SPECIES}) ---"
    echo "Input : ${input}"
    echo "Output: ${output}"
    echo "Plots : ${outdir}"

    # Job 1: prepare variant sequences (CPU)
    job1=$(sbatch --parsable \
        --job-name="${job_prefix}_1_prep" \
        --output="${REPO_ROOT}/logs/${job_prefix}_1_prepare_%j.out" \
        --error="${REPO_ROOT}/logs/${job_prefix}_1_prepare_%j.err" \
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
mkdir -p logs data "results/${SPECIES}" "${outdir}"

for ref in "${FASTA}" "${GTF}"; do
    if [[ ! -r "\${ref}" ]]; then
        echo "[ERROR] Reference file is not readable on \$(hostname): \${ref}" >&2
        exit 1
    fi
done

echo "=== Step 1: Generating ${gene} variant sequences (${SPECIES}) ==="
python prepare_scn2a_ism.py \
    --species "${SPECIES}" \
    --gene-name "${gene}" \
    --fasta "${FASTA}" \
    --gtf "${GTF}" \
    --upstream-bases "${UPSTREAM}" \
    --max-deletion "${MAX_DEL}" \
    --truncate-utr3 \
    --output "${input}"
EOF
)
    echo "Job 1 submitted  (prepare)  : ${job1}"

    # Job 2: RiboNN predictions (GPU)
    job2=$(sbatch --parsable \
        --job-name="${job_prefix}_2_pred" \
        --output="${REPO_ROOT}/logs/${job_prefix}_2_predict_%j.out" \
        --error="${REPO_ROOT}/logs/${job_prefix}_2_predict_%j.err" \
        --partition=ga100 \
        --gres=gpu:1 \
        --cpus-per-task=4 \
        --mem=32G \
        --time=06:00:00 \
        --dependency=afterok:${job1} \
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

echo "=== Step 2: Running ${SPECIES} RiboNN predictions for ${gene} ==="
python run_ribonn_predict.py \
    --species "${SPECIES}" \
    --input "${input}" \
    --output "${output}" \
    --top-k 5 \
    --batch-size 1024 \
    --num-workers "\${SLURM_CPUS_PER_TASK:-4}"
EOF
)
    echo "Job 2 submitted  (predict)   : ${job2}  [afterok:${job1}]"

    # Job 3: ISM plots (CPU)
    job3=$(sbatch --parsable \
        --job-name="${job_prefix}_3_plot" \
        --output="${REPO_ROOT}/logs/${job_prefix}_3_plot_%j.out" \
        --error="${REPO_ROOT}/logs/${job_prefix}_3_plot_%j.err" \
        --partition=ncpu \
        --cpus-per-task=2 \
        --mem=8G \
        --time=00:30:00 \
        --dependency=afterok:${job2} \
        --chdir="${REPO_ROOT}" \
        << EOF
#!/bin/bash
set -eo pipefail
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate ribonn
export PYTHONWARNINGS="\${PYTHONWARNINGS:+\${PYTHONWARNINGS},}ignore:pkg_resources is deprecated as an API:UserWarning"

echo "Host: \$(hostname)"

ACTUAL_UPSTREAM=\$(python3 -c "
import re
max_off = 0
with open('${input}') as fh:
    next(fh)
    for line in fh:
        tid = line.split('\t')[0]
        m = re.match(r'^([+-]\d+)_[ACGT]>[ACGT]$', tid)
        if m:
            max_off = max(max_off, abs(int(m.group(1))))
print(max_off if max_off else 15)
")
echo "Detected upstream bases: \${ACTUAL_UPSTREAM}"

echo "=== Step 3: Generating ${gene} ISM plots (${SPECIES}) ==="
python plot_te_changes.py \
    --input "${output}" \
    --outdir "${outdir}" \
    --upstream-bases "\${ACTUAL_UPSTREAM}" \
    --gene-name "${gene}" \
    --species "${SPECIES}"

echo "Plots saved to ${outdir}/"
EOF
)
    echo "Job 3 submitted  (plot)      : ${job3}  [afterok:${job2}]"

    ALL_JOBS+=("${job1}" "${job2}" "${job3}")

    if [[ "${SKIP_NULL}" -eq 0 ]]; then
        # Job 4: altAUG null test (CPU)
        job4=$(sbatch --parsable \
            --job-name="${job_prefix}_4_null" \
            --output="${REPO_ROOT}/logs/${job_prefix}_4_null_%j.out" \
            --error="${REPO_ROOT}/logs/${job_prefix}_4_null_%j.err" \
            --partition=ncpu \
            --cpus-per-task=2 \
            --mem=8G \
            --time=00:30:00 \
            --dependency=afterok:${job2} \
            --chdir="${REPO_ROOT}" \
            << EOF
#!/bin/bash
set -eo pipefail
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate ribonn
export PYTHONWARNINGS="\${PYTHONWARNINGS:+\${PYTHONWARNINGS},}ignore:pkg_resources is deprecated as an API:UserWarning"

echo "Host: \$(hostname)"

echo "=== Step 4: ${gene} altAUG position-matched null test (${SPECIES}) ==="
python analyze_altaug_null.py \
    --input "${output}" \
    --altaug-positions ${ALTAUG_POSITIONS[*]} \
    --outdir "${outdir}" \
    --gene-name "${gene}" \
    --species "${SPECIES}"

echo "Summary: ${outdir}/altaug_null_summary.tsv"
echo "Plot:    ${outdir}/altaug_null_plot.png"
EOF
)
        echo "Job 4 submitted  (null test) : ${job4}  [afterok:${job2}]"
        ALL_JOBS+=("${job4}")
    fi

    echo ""
}

for gene in "${GENES[@]}"; do
    if [[ -z "${gene}" ]]; then
        continue
    fi
    submit_gene_pipeline "${gene}"
done

JOB_CSV="$(IFS=,; echo "${ALL_JOBS[*]}")"

echo "Pipeline submitted. Monitor with:"
echo "  squeue -j ${JOB_CSV}"
echo "  sacct  -j ${JOB_CSV} --format=JobID,JobName,State,ExitCode,Elapsed,Reason"
echo ""
echo "Outputs:"
for gene in "${GENES[@]}"; do
    gene_lc="$(printf '%s' "${gene}" | tr '[:upper:]' '[:lower:]')"
    echo "  ${gene} ${SPECIES}: results/${SPECIES}/${gene_lc}_prediction_output.txt, plots_ism_${SPECIES}_${gene_lc}/"
done
