#!/bin/bash
# submit_riboscanner_master_analysis.sh
#
# Submit the plain-R RiboScanner master-table analysis as a SLURM job.
#
# Usage:
#   bash submit_riboscanner_master_analysis.sh test
#   bash submit_riboscanner_master_analysis.sh full
#   bash submit_riboscanner_master_analysis.sh full --time 12:00:00 --mem 64G
#   bash submit_riboscanner_master_analysis.sh test --rscript /path/to/Rscript
#   bash submit_riboscanner_master_analysis.sh test --riboscanner-env /path/to/envs/RiboScanner
#
# Defaults:
#   test  -> uses scn2a_mouse_centered.fa and writes riboscanner_analysis_scn2a_test/
#   full  -> uses master_table.context_m40_p40.fa.gz and writes riboscanner_analysis/

set -eo pipefail

usage() {
    sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/riboscanner_master_analysis.R" ]]; then
    REPO_ROOT="${SCRIPT_DIR}"
elif [[ -f "${PWD}/riboscanner_master_analysis.R" ]]; then
    REPO_ROOT="${PWD}"
else
    echo "[ERROR] Could not find riboscanner_master_analysis.R. Run from the repo root." >&2
    exit 1
fi

MODE="${1:-test}"
case "${MODE}" in
    test|full)
        shift || true
        ;;
    -h|--help)
        usage
        exit 0
        ;;
    *)
        echo "[ERROR] First argument must be 'test' or 'full' (got '${MODE}')." >&2
        usage >&2
        exit 2
        ;;
esac

PARTITION="${RIBOSCANNER_PARTITION:-ncpu}"
CPUS="${RIBOSCANNER_CPUS:-4}"
MEM="${RIBOSCANNER_MEM:-24G}"
TIME="${RIBOSCANNER_TIME:-08:00:00}"
RIBOSCANNER_DIR="${RIBOSCANNER_DIR:-/camp/lab/ulej/home/users/luscomben/users/iosubi/projects/ag/riboscanner}"
RIBOSCANNER_EXE="${RIBOSCANNER_EXE:-RiboScanner}"
RIBOSCANNER_ENV="${RIBOSCANNER_ENV:-}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
RUN_PREDICTION="${RIBOSCANNER_RUN_PREDICTION:-true}"
ANALYSIS_DIR="riboscanner_analysis"
EXTRA_R_ARGS=()

if [[ "${MODE}" == "test" ]]; then
    TIME="${RIBOSCANNER_TIME:-01:00:00}"
    MEM="${RIBOSCANNER_MEM:-8G}"
    ANALYSIS_DIR="riboscanner_analysis_scn2a_test"
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --partition)
            PARTITION="$2"
            shift 2
            ;;
        --cpus)
            CPUS="$2"
            shift 2
            ;;
        --mem)
            MEM="$2"
            shift 2
            ;;
        --time)
            TIME="$2"
            shift 2
            ;;
        --riboscanner-dir)
            RIBOSCANNER_DIR="$2"
            shift 2
            ;;
        --riboscanner-exe)
            RIBOSCANNER_EXE="$2"
            shift 2
            ;;
        --riboscanner-env)
            RIBOSCANNER_ENV="$2"
            shift 2
            ;;
        --rscript)
            RSCRIPT_BIN="$2"
            shift 2
            ;;
        --analysis-dir)
            ANALYSIS_DIR="$2"
            shift 2
            ;;
        --no-predict)
            RUN_PREDICTION="false"
            shift
            ;;
        --fasta|--test-fasta|--master-table|--prediction-output|--min-sequence-length|--id-col|--prediction-id-col|--class-col|--detected-col|--uorf-col|--gfp-col|--score-col|--uncertainty-col)
            EXTRA_R_ARGS+=("$1" "$2")
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case "${RUN_PREDICTION}" in
    true|TRUE|1|yes|YES)
        NO_PREDICT_ARG=()
        ;;
    false|FALSE|0|no|NO)
        NO_PREDICT_ARG=("--no-predict")
        ;;
    *)
        echo "[ERROR] RIBOSCANNER_RUN_PREDICTION must be true or false (got '${RUN_PREDICTION}')." >&2
        exit 2
        ;;
esac

mkdir -p "${REPO_ROOT}/logs"

echo "=== RiboScanner R analysis submission ==="
echo "Repo root      : ${REPO_ROOT}"
echo "Mode           : ${MODE}"
echo "RiboScanner dir: ${RIBOSCANNER_DIR}"
echo "RiboScanner exe: ${RIBOSCANNER_EXE}"
echo "RiboScanner env: ${RIBOSCANNER_ENV:-auto-detect}"
echo "Rscript        : ${RSCRIPT_BIN}"
echo "Analysis dir   : ${ANALYSIS_DIR}"
echo "Run prediction : $([[ ${#NO_PREDICT_ARG[@]} -eq 0 ]] && echo TRUE || echo FALSE)"
echo "Partition      : ${PARTITION}"
echo "CPUs / mem/time: ${CPUS} / ${MEM} / ${TIME}"
echo ""

JOB=$(sbatch --parsable \
    --job-name="riboscan_${MODE}" \
    --output="${REPO_ROOT}/logs/riboscan_${MODE}_%j.out" \
    --error="${REPO_ROOT}/logs/riboscan_${MODE}_%j.err" \
    --partition="${PARTITION}" \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --time="${TIME}" \
    --chdir="${REPO_ROOT}" \
    << EOF
#!/bin/bash
set -eo pipefail

module load R 2>/dev/null || true

RIBOSCANNER_ENV="${RIBOSCANNER_ENV}"
RIBOSCANNER_EXE_RESOLVED="${RIBOSCANNER_EXE}"
if [[ "\${RIBOSCANNER_EXE_RESOLVED}" != */* ]]; then
    RIBOSCANNER_EXE_RESOLVED="\$(command -v "\${RIBOSCANNER_EXE_RESOLVED}" || true)"
fi
if [[ -z "\${RIBOSCANNER_ENV}" ]] && [[ -n "\${RIBOSCANNER_EXE_RESOLVED}" ]]; then
    RIBOSCANNER_BIN_DIR="\$(cd "\$(dirname "\${RIBOSCANNER_EXE_RESOLVED}")" && pwd)"
    if [[ "\$(basename "\${RIBOSCANNER_BIN_DIR}")" == "bin" ]] && [[ -d "\$(dirname "\${RIBOSCANNER_BIN_DIR}")/lib" ]]; then
        RIBOSCANNER_ENV="\$(dirname "\${RIBOSCANNER_BIN_DIR}")"
    fi
fi
if [[ -n "\${RIBOSCANNER_ENV}" ]]; then
    export PATH="\${RIBOSCANNER_ENV}/bin:\${PATH}"
    export LD_LIBRARY_PATH="\${RIBOSCANNER_ENV}/lib:\${LD_LIBRARY_PATH:-}"
    export CONDA_PREFIX="\${RIBOSCANNER_ENV}"
fi

echo "Host: \$(hostname)"
echo "Working dir: \$(pwd)"
echo "Rscript: ${RSCRIPT_BIN}"
echo "RiboScanner: \$(command -v "${RIBOSCANNER_EXE}" || true)"
echo "RiboScanner env: \${RIBOSCANNER_ENV:-not detected}"
echo "LD_LIBRARY_PATH: \${LD_LIBRARY_PATH:-unset}"

mkdir -p "${RIBOSCANNER_DIR}/${ANALYSIS_DIR}"

if ! "${RSCRIPT_BIN}" --version >/dev/null 2>&1; then
    echo "[ERROR] Rscript is not available or not executable: ${RSCRIPT_BIN}" >&2
    echo "        Try passing --rscript /path/to/Rscript, or load a newer R module." >&2
    exit 1
fi

"${RSCRIPT_BIN}" riboscanner_master_analysis.R \
    --mode "${MODE}" \
    --riboscanner-dir "${RIBOSCANNER_DIR}" \
    --riboscanner-exe "${RIBOSCANNER_EXE}" \
    --riboscanner-env "\${RIBOSCANNER_ENV:-}" \
    --analysis-dir "${ANALYSIS_DIR}" \
    ${NO_PREDICT_ARG[*]} \
    ${EXTRA_R_ARGS[*]}

echo ""
echo "=== Done ==="
echo "Tables: ${RIBOSCANNER_DIR}/${ANALYSIS_DIR}/"
echo "Plots:  ${RIBOSCANNER_DIR}/${ANALYSIS_DIR}/figures/"
EOF
)

echo "Job submitted: ${JOB}"
echo "Monitor with:"
echo "  squeue -j ${JOB}"
echo "  sacct  -j ${JOB} --format=JobID,JobName,State,ExitCode,Elapsed,Reason"
echo ""
echo "Logs:"
echo "  ${REPO_ROOT}/logs/riboscan_${MODE}_${JOB}.out"
echo "  ${REPO_ROOT}/logs/riboscan_${MODE}_${JOB}.err"
