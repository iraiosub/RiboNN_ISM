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
#   bash submit_riboscanner_master_analysis.sh test --cxx-lib-dir /path/to/lib
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
RIBOSCANNER_CXX_LIB_DIR="${RIBOSCANNER_CXX_LIB_DIR:-}"
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
        --cxx-lib-dir)
            RIBOSCANNER_CXX_LIB_DIR="$2"
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
echo "C++ lib dir    : ${RIBOSCANNER_CXX_LIB_DIR:-auto-detect}"
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

RIBOSCANNER_CXX_LIB_DIR="${RIBOSCANNER_CXX_LIB_DIR}"
RIBOSCANNER_CXX_LIB_FILE=""
RIBOSCANNER_CHILD_LD_LIBRARY_PATH=""
RUN_RIBOSCANNER_PREDICTION="${RUN_PREDICTION}"

find_libstdcxx_file() {
    local candidate_dir="\${1}"
    if [[ -e "\${candidate_dir}/libstdc++.so.6" ]]; then
        echo "\${candidate_dir}/libstdc++.so.6"
        return 0
    fi

    local candidate_file
    for candidate_file in "\${candidate_dir}"/libstdc++.so.6*; do
        if [[ -e "\${candidate_file}" ]]; then
            echo "\${candidate_file}"
            return 0
        fi
    done
    return 1
}

libstdcxx_has_cxxabi() {
    local candidate_lib
    candidate_lib="\$(find_libstdcxx_file "\${1}" || true)"
    [[ -e "\${candidate_lib}" ]] || return 1
    command -v strings >/dev/null 2>&1 || return 1
    strings "\${candidate_lib}" | grep -q '^CXXABI_1\\.3\\.11$'
}

report_libstdcxx() {
    local candidate_dir="\${1}"
    local candidate_lib
    candidate_lib="\$(find_libstdcxx_file "\${candidate_dir}" || true)"
    if [[ -e "\${candidate_lib}" ]]; then
        echo "Candidate libstdc++: \$(readlink -f "\${candidate_lib}" 2>/dev/null || echo "\${candidate_lib}")"
        if command -v strings >/dev/null 2>&1; then
            strings "\${candidate_lib}" | grep '^CXXABI_' | tail -n 8 || true
        fi
    else
        echo "Candidate libstdc++ missing: \${candidate_lib}"
    fi
}

diagnose_libstdcxx_install() {
    echo ""
    echo "=== libstdc++ runtime diagnostics ==="
    echo "RIBOSCANNER_ENV=\${RIBOSCANNER_ENV:-unset}"
    echo "RIBOSCANNER_BASE_PREFIX=\${RIBOSCANNER_BASE_PREFIX:-unset}"
    if [[ -n "\${RIBOSCANNER_ENV}" ]] && [[ -x "\${RIBOSCANNER_ENV}/bin/conda" ]]; then
        "\${RIBOSCANNER_ENV}/bin/conda" list -p "\${RIBOSCANNER_ENV}" 'libstdcxx|libgcc|gcc_impl' || true
    elif [[ -n "\${RIBOSCANNER_BASE_PREFIX}" ]] && [[ -x "\${RIBOSCANNER_BASE_PREFIX}/bin/conda" ]]; then
        "\${RIBOSCANNER_BASE_PREFIX}/bin/conda" list -p "\${RIBOSCANNER_ENV}" 'libstdcxx|libgcc|gcc_impl' || true
    elif command -v conda >/dev/null 2>&1; then
        conda list -p "\${RIBOSCANNER_ENV}" 'libstdcxx|libgcc|gcc_impl' || true
    else
        echo "conda command not found in job environment"
    fi

    echo ""
    echo "Searching for libstdc++.so.6* under env/base:"
    for RIBOSCANNER_SEARCH_ROOT in "\${RIBOSCANNER_ENV}" "\${RIBOSCANNER_BASE_PREFIX:-}"; do
        if [[ -n "\${RIBOSCANNER_SEARCH_ROOT}" ]] && [[ -d "\${RIBOSCANNER_SEARCH_ROOT}" ]]; then
            find "\${RIBOSCANNER_SEARCH_ROOT}" -maxdepth 5 -name 'libstdc++.so.6*' -print 2>/dev/null | sed 's/^/  /' || true
        fi
    done
    echo "=== end diagnostics ==="
    echo ""
}

RIBOSCANNER_CXX_CANDIDATES=()
if [[ -n "\${RIBOSCANNER_CXX_LIB_DIR}" ]]; then
    RIBOSCANNER_CXX_CANDIDATES+=("\${RIBOSCANNER_CXX_LIB_DIR}")
elif [[ -n "\${RIBOSCANNER_ENV}" ]]; then
    RIBOSCANNER_CXX_CANDIDATES+=("\${RIBOSCANNER_ENV}/lib")
    RIBOSCANNER_BASE_PREFIX="\$(cd "\${RIBOSCANNER_ENV}/../.." 2>/dev/null && pwd || true)"
    if [[ -n "\${RIBOSCANNER_BASE_PREFIX}" ]]; then
        RIBOSCANNER_CXX_CANDIDATES+=("\${RIBOSCANNER_BASE_PREFIX}/lib")
    fi
fi

if [[ -z "\${RIBOSCANNER_CXX_LIB_DIR}" ]]; then
    for RIBOSCANNER_CANDIDATE_DIR in "\${RIBOSCANNER_CXX_CANDIDATES[@]}"; do
        if libstdcxx_has_cxxabi "\${RIBOSCANNER_CANDIDATE_DIR}"; then
            RIBOSCANNER_CXX_LIB_DIR="\${RIBOSCANNER_CANDIDATE_DIR}"
            RIBOSCANNER_CXX_LIB_FILE="\$(find_libstdcxx_file "\${RIBOSCANNER_CANDIDATE_DIR}")"
            break
        fi
    done
fi

if [[ -z "\${RIBOSCANNER_CXX_LIB_DIR}" ]]; then
    for RIBOSCANNER_CANDIDATE_DIR in "\${RIBOSCANNER_CXX_CANDIDATES[@]}"; do
        RIBOSCANNER_CANDIDATE_FILE="\$(find_libstdcxx_file "\${RIBOSCANNER_CANDIDATE_DIR}" || true)"
        if [[ -n "\${RIBOSCANNER_CANDIDATE_FILE}" ]]; then
            RIBOSCANNER_CXX_LIB_DIR="\${RIBOSCANNER_CANDIDATE_DIR}"
            RIBOSCANNER_CXX_LIB_FILE="\${RIBOSCANNER_CANDIDATE_FILE}"
            echo "[WARN] Could not confirm CXXABI_1.3.11 from strings; using first libstdc++.so.6* candidate and letting torch preflight decide." >&2
            break
        fi
    done
fi

if [[ -n "\${RIBOSCANNER_CXX_LIB_DIR}" ]] && ! libstdcxx_has_cxxabi "\${RIBOSCANNER_CXX_LIB_DIR}"; then
    echo "[WARN] Could not confirm CXXABI_1.3.11 in selected runtime; continuing to torch preflight." >&2
    report_libstdcxx "\${RIBOSCANNER_CXX_LIB_DIR}"
fi

if [[ -n "\${RIBOSCANNER_CXX_LIB_DIR}" ]] && [[ -z "\${RIBOSCANNER_CXX_LIB_FILE}" ]]; then
    RIBOSCANNER_CXX_LIB_FILE="\$(find_libstdcxx_file "\${RIBOSCANNER_CXX_LIB_DIR}" || true)"
fi

if [[ -z "\${RIBOSCANNER_CXX_LIB_DIR}" ]] && [[ "\${RUN_RIBOSCANNER_PREDICTION}" =~ ^(true|TRUE|1|yes|YES)$ ]]; then
    echo "[ERROR] Could not find a libstdc++.so.6 candidate for RiboScanner." >&2
    diagnose_libstdcxx_install
    echo "        Force update with:" >&2
    echo "        \${RIBOSCANNER_BASE_PREFIX:-/path/to/miniconda3}/bin/conda install -y -p \${RIBOSCANNER_ENV} -c conda-forge --override-channels 'libstdcxx-ng>=12' 'libgcc-ng>=12' 'gcc_impl_linux-64>=12'" >&2
    exit 1
fi

if [[ -n "\${RIBOSCANNER_CXX_LIB_DIR}" ]] && [[ -z "\${RIBOSCANNER_CXX_LIB_FILE}" ]]; then
    echo "[ERROR] Could not find a libstdc++.so.6* file in selected C++ lib dir: \${RIBOSCANNER_CXX_LIB_DIR}" >&2
    diagnose_libstdcxx_install
    exit 1
fi

if [[ -n "\${RIBOSCANNER_CXX_LIB_FILE}" ]] && [[ "\$(basename "\${RIBOSCANNER_CXX_LIB_FILE}")" != "libstdc++.so.6" ]]; then
    RIBOSCANNER_RUNTIME_LIB_DIR="${RIBOSCANNER_DIR}/${ANALYSIS_DIR}/runtime_lib"
    mkdir -p "\${RIBOSCANNER_RUNTIME_LIB_DIR}"
    ln -sf "\${RIBOSCANNER_CXX_LIB_FILE}" "\${RIBOSCANNER_RUNTIME_LIB_DIR}/libstdc++.so.6"
    RIBOSCANNER_CXX_LIB_DIR="\${RIBOSCANNER_RUNTIME_LIB_DIR}"
    RIBOSCANNER_CXX_LIB_FILE="\${RIBOSCANNER_RUNTIME_LIB_DIR}/libstdc++.so.6"
fi

RIBOSCANNER_CHILD_LD_LIBRARY_PATH="\${RIBOSCANNER_CXX_LIB_DIR}"
if [[ -n "\${RIBOSCANNER_ENV}" ]]; then
    if [[ -z "\${RIBOSCANNER_CHILD_LD_LIBRARY_PATH}" ]]; then
        RIBOSCANNER_CHILD_LD_LIBRARY_PATH="\${RIBOSCANNER_ENV}/lib"
    elif [[ "\${RIBOSCANNER_ENV}/lib" != "\${RIBOSCANNER_CXX_LIB_DIR}" ]]; then
        RIBOSCANNER_CHILD_LD_LIBRARY_PATH="\${RIBOSCANNER_CHILD_LD_LIBRARY_PATH}:\${RIBOSCANNER_ENV}/lib"
    fi
    for RIBOSCANNER_TORCH_LIB in "\${RIBOSCANNER_ENV}"/lib/python*/site-packages/torch/lib; do
        if [[ -d "\${RIBOSCANNER_TORCH_LIB}" ]]; then
            RIBOSCANNER_CHILD_LD_LIBRARY_PATH="\${RIBOSCANNER_CHILD_LD_LIBRARY_PATH}:\${RIBOSCANNER_TORCH_LIB}"
        fi
    done
fi

echo "Host: \$(hostname)"
echo "Working dir: \$(pwd)"
echo "Rscript: ${RSCRIPT_BIN}"
echo "R analysis script version line:"
grep -n '^analysis_script_version' riboscanner_master_analysis.R || true
echo "R analysis predict command line:"
grep -n 'cmd_args <- c("predict".*"--output"' riboscanner_master_analysis.R || true
echo "RiboScanner: \$(command -v "${RIBOSCANNER_EXE}" || true)"
echo "RiboScanner env: \${RIBOSCANNER_ENV:-not detected}"
echo "Selected C++ lib dir: \${RIBOSCANNER_CXX_LIB_DIR:-not detected}"
echo "Selected C++ lib file: \${RIBOSCANNER_CXX_LIB_FILE:-not detected}"
echo "LD_LIBRARY_PATH: \${LD_LIBRARY_PATH:-unset}"
echo "RiboScanner child LD_LIBRARY_PATH: \${RIBOSCANNER_CHILD_LD_LIBRARY_PATH:-unset}"
echo "RiboScanner child LD_PRELOAD: \${RIBOSCANNER_CXX_LIB_FILE:-unset}"

if [[ -n "\${RIBOSCANNER_ENV}" ]] && [[ "\${RUN_RIBOSCANNER_PREDICTION}" =~ ^(true|TRUE|1|yes|YES)$ ]]; then
    if [[ -n "\${RIBOSCANNER_CXX_LIB_DIR}" ]]; then
        report_libstdcxx "\${RIBOSCANNER_CXX_LIB_DIR}"
    fi

    if [[ -x "\${RIBOSCANNER_ENV}/bin/python" ]]; then
        env PATH="\${RIBOSCANNER_ENV}/bin:\${PATH}" LD_LIBRARY_PATH="\${RIBOSCANNER_CHILD_LD_LIBRARY_PATH}" LD_PRELOAD="\${RIBOSCANNER_CXX_LIB_FILE}" CONDA_PREFIX="\${RIBOSCANNER_ENV}" \
            "\${RIBOSCANNER_ENV}/bin/python" -c "import torch; print('torch import OK', torch.__version__)"
    fi
elif [[ ! "\${RUN_RIBOSCANNER_PREDICTION}" =~ ^(true|TRUE|1|yes|YES)$ ]]; then
    echo "Skipping RiboScanner runtime preflight because prediction is disabled."
fi

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
    --cxx-lib-dir "\${RIBOSCANNER_CXX_LIB_DIR:-}" \
    --analysis-dir "${ANALYSIS_DIR}" \
    --prediction-output "${ANALYSIS_DIR}/riboscanner_predictions.tsv" \
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
