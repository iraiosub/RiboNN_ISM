#!/bin/bash
# setup_hpc.sh – Install RiboNN on CAMP/NEMO HPC (Francis Crick Institute)
#
# Run ONCE from the repo root after cloning:
#   git clone https://github.com/<you>/RiboNN.git
#   cd RiboNN
#   bash setup_hpc.sh
#
# What this does:
#   1. Load mamba (or install miniforge if absent)
#   2. Create/update conda env "ribonn" from environment.yml
#   3. Download model weights from Zenodo (also happens automatically at predict time)
#
# Useful repair mode:
#   bash setup_hpc.sh --recreate

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/environment.yml" ]]; then
    REPO_ROOT="${SCRIPT_DIR}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]] && [[ -f "${SLURM_SUBMIT_DIR}/environment.yml" ]]; then
    REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
elif [[ -f "${PWD}/environment.yml" ]]; then
    REPO_ROOT="${PWD}"
else
    echo "[ERROR] Could not find environment.yml. Run this from the RiboNN repo root, or submit it with sbatch from the repo root." >&2
    exit 1
fi

cd "${REPO_ROOT}"
ENV_NAME="ribonn"
ZENODO_URL="https://zenodo.org/records/17258709/files/weights.zip"
RECREATE_ENV=0

for arg in "$@"; do
    case "${arg}" in
        --recreate|--nuke)
            RECREATE_ENV=1
            ;;
        -h|--help)
            echo "Usage: bash setup_hpc.sh [--recreate]"
            echo "  --recreate  Remove the existing '${ENV_NAME}' conda env before creating it."
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: ${arg}" >&2
            echo "Usage: bash setup_hpc.sh [--recreate]" >&2
            exit 2
            ;;
    esac
done

echo "=== RiboNN HPC Setup ==="
echo "Repo: ${REPO_ROOT}"
echo "Env:  ${ENV_NAME}"
if [[ "${RECREATE_ENV}" -eq 1 ]]; then
    echo "Mode: recreate env"
fi

# ── 1. Make sure mamba/conda is available ────────────────────────────────────
if ! command -v mamba &>/dev/null && ! command -v conda &>/dev/null; then
    echo "Neither mamba nor conda found. Installing Miniforge3 ..."
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    curl -fsSL "${MINIFORGE_URL}" -o /tmp/Miniforge3.sh
    bash /tmp/Miniforge3.sh -b -p "${HOME}/miniforge3"
    rm /tmp/Miniforge3.sh
    source "${HOME}/miniforge3/etc/profile.d/conda.sh"
    conda init bash
    echo "Miniforge installed. Please re-run this script in a new shell."
    exit 0
fi

# Source conda so 'conda activate' works in non-interactive shells
CONDA_BASE="$(conda info --base 2>/dev/null || mamba info --base 2>/dev/null)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Prefer mamba for speed
INSTALLER="conda"
command -v mamba &>/dev/null && INSTALLER="mamba"
echo "Installer: ${INSTALLER}"

# ── 2. Create environment ─────────────────────────────────────────────────────
if [[ "${RECREATE_ENV}" -eq 1 ]] && conda env list | grep -qE "^${ENV_NAME}\s"; then
    if [[ "${CONDA_DEFAULT_ENV:-}" == "${ENV_NAME}" ]]; then
        echo "Deactivating '${ENV_NAME}' before removal ..."
        conda deactivate
    fi
    echo "Removing existing conda env '${ENV_NAME}' ..."
    conda env remove -n "${ENV_NAME}" -y
fi

if conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo "Conda env '${ENV_NAME}' already exists; updating from environment.yml ..."
    ${INSTALLER} env update -n "${ENV_NAME}" -f "${REPO_ROOT}/environment.yml" --prune
else
    echo "Creating conda env '${ENV_NAME}' from environment.yml ..."
    ${INSTALLER} env create -n "${ENV_NAME}" -f "${REPO_ROOT}/environment.yml"
    echo "Environment created."
fi

# Existing environments may contain a stale pip torch build or a too-new MKL
# runtime. Fail loudly; --recreate is the reliable repair path.
if ! TORCH_CHECK_OUTPUT="$(conda run -n "${ENV_NAME}" python -c "
import torch
version = torch.__version__.split('+')[0]
cuda = torch.version.cuda or ''
print(f'  torch build check: torch={torch.__version__} cuda_build={cuda}')
raise SystemExit(0 if version == '1.13.1' and cuda.startswith('11.7') else 1)
" 2>&1)"; then
    echo "${TORCH_CHECK_OUTPUT}" >&2
    if echo "${TORCH_CHECK_OUTPUT}" | grep -q "iJIT_NotifyEvent"; then
        echo "[ERROR] PyTorch import failed because MKL is too new for this PyTorch build." >&2
        echo "        environment.yml now pins mkl=2024.0; run: bash setup_hpc.sh --recreate" >&2
        exit 1
    fi
    echo "[ERROR] The '${ENV_NAME}' env is still importing the wrong PyTorch build." >&2
    echo "        Run: bash setup_hpc.sh --recreate" >&2
    exit 1
fi
echo "${TORCH_CHECK_OUTPUT}"

# ── 3. Download model weights ─────────────────────────────────────────────────
MODELS_DIR="${REPO_ROOT}/models/human"
LEGACY_MODELS_DIR="${REPO_ROOT}/human"
if [[ -f "${MODELS_DIR}/runs.csv" ]] && find "${MODELS_DIR}" -mindepth 2 -maxdepth 2 -name state_dict.pth -print -quit | grep -q .; then
    echo "Model weights already present in ${MODELS_DIR}/"
elif [[ -f "${LEGACY_MODELS_DIR}/runs.csv" ]] && find "${LEGACY_MODELS_DIR}" -mindepth 2 -maxdepth 2 -name state_dict.pth -print -quit | grep -q .; then
    echo "Found model weights in ${LEGACY_MODELS_DIR}/"
    mkdir -p "${REPO_ROOT}/models"
    ln -sfn "../human" "${MODELS_DIR}"
    echo "Linked ${MODELS_DIR} -> ../human"
else
    echo "Downloading model weights from Zenodo ..."
    mkdir -p "${REPO_ROOT}/models" "${REPO_ROOT}/tmp"
    TMP_WEIGHTS_ZIP="${REPO_ROOT}/tmp/weights.zip"
    wget -q --show-progress -O "${TMP_WEIGHTS_ZIP}" "${ZENODO_URL}"
    unzip -oq "${TMP_WEIGHTS_ZIP}" -d "${REPO_ROOT}/models"
    rm -f "${TMP_WEIGHTS_ZIP}"
    echo "Weights extracted to ${MODELS_DIR}/"
fi

# ── 4. Quick smoke test ───────────────────────────────────────────────────────
echo ""
echo "Running quick import test ..."
# LD_LIBRARY_PATH needed on HPC systems where system libstdc++ is older than
# what conda-installed packages (pandas, torch) were compiled against.
CONDA_LIB="$(conda run -n "${ENV_NAME}" python -c 'import sys, os; print(os.path.join(sys.prefix, "lib"))')"
conda run -n "${ENV_NAME}" env LD_LIBRARY_PATH="${CONDA_LIB}:${LD_LIBRARY_PATH}" python -c "
import torch
import importlib.util
import pyfaidx
from src.predict import predict_using_nested_cross_validation_models
if importlib.util.find_spec('pkg_resources') is None:
    raise ModuleNotFoundError('pkg_resources; install setuptools')
print('  torch:', torch.__version__)
print('  torch cuda build:', torch.version.cuda)
print('  pyfaidx:', getattr(pyfaidx, '__version__', 'installed'))
print('  cuda:', torch.cuda.is_available())
print('  RiboNN src import: OK')
"
echo ""
echo "Add to your ~/.bashrc or activate script if needed:"
echo "  export LD_LIBRARY_PATH=\"\${CONDA_PREFIX}/lib:\${LD_LIBRARY_PATH}\""

echo ""
echo "=== Setup complete ==="
echo ""
echo "To run ISM:"
echo "  sbatch submit_ism_scn2a.sh"
echo ""
echo "Or step-by-step:"
echo "  conda activate ${ENV_NAME}"
echo "  python prepare_scn2a_ism.py --truncate-utr3"
echo "  python run_ribonn_predict.py"
echo "  python plot_te_changes.py"
