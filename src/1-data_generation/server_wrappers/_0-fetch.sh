#!/usr/bin/env bash

#SBATCH --job-name=fetch_AM
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%j.out
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8   
#SBATCH --mem=32G
#SBATCH -p ggs

set -Eeuo pipefail

report_error() {
    local status=$?
    local message="[ERROR] Fetch stopped near line ${BASH_LINENO[0]:-unknown} (exit ${status}): ${BASH_COMMAND}"
    trap - ERR
    if [[ -n "${LOG_PATH:-}" ]]; then
        printf '%s\n' "${message}" | tee -a "${LOG_PATH}" >&2 || true
    else
        printf '%s\n' "${message}" >&2
    fi
    exit "${status}"
}
trap report_error ERR

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -f "${SCRIPT_DIR}/__config.sh" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SEARCH_DIR="${SLURM_SUBMIT_DIR}"
    while [[ "${SEARCH_DIR}" != / && ! -f "${SEARCH_DIR}/__config.sh" && ! -f "${SEARCH_DIR}/src/1-data_generation/server_wrappers/__config.sh" ]]; do
        SEARCH_DIR="$(dirname -- "${SEARCH_DIR}")"
    done
    if [[ -f "${SEARCH_DIR}/src/1-data_generation/server_wrappers/__config.sh" ]]; then
        SCRIPT_DIR="${SEARCH_DIR}/src/1-data_generation/server_wrappers"
    elif [[ -f "${SEARCH_DIR}/__config.sh" ]]; then
        SCRIPT_DIR="${SEARCH_DIR}"
    else
        echo "ERROR: Could not locate __config.sh from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}" >&2
        exit 1
    fi
fi
BASE_PATH="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"

# Submit this wrapper with sbatch from the server_wrappers directory.
if [[ -z "${SLURM_JOB_ID:-${SLURM_JOBID:-}}" ]]; then
    echo "ERROR: Run sbatch _0-fetch.sh from the server_wrappers directory." >&2
    exit 1
fi

source "${SCRIPT_DIR}/__config.sh"

RECORD_ID="${ALPHAMISSENSE_RECORD_ID:-10813168}"
RAW_DIR="${BASE_PATH}/data/1-raw/missense_db"
mkdir -p "${RAW_DIR}" "${BASE_PATH}/log"
LOG_PATH="${BASE_PATH}/log/fetch_am_${SLURM_JOB_ID}_$(date +%Y%m%d_%H%M%S).log"

log() {
    printf '%s\n' "$*" | tee -a "${LOG_PATH}"
}

download_table() {
    local file_name="$1"
    local table_path="${RAW_DIR}/${file_name%.gz}"
    local compressed_path="${RAW_DIR}/${file_name}"
    local partial_path="${compressed_path}.part"
    local query_path="${compressed_path}?download=1"
    local url="https://zenodo.org/records/${RECORD_ID}/files/${file_name}"

    if [[ -s "${table_path}" ]]; then
        log "[SKIP] ${table_path} already exists (prepared)."
        return
    fi

    # Recover partial files made by earlier versions of this wrapper.
    if [[ -f "${query_path}" && ! -e "${compressed_path}" ]]; then
        mv -- "${query_path}" "${compressed_path}"
        log "[INFO] Resuming partial download ${file_name}."
    fi
    if [[ -f "${partial_path}" && ! -e "${compressed_path}" ]]; then
        mv -- "${partial_path}" "${compressed_path}"
        log "[INFO] Resuming partial download ${file_name}."
    fi

    log "[DOWNLOAD] ${file_name}"
    local wget_status=0
    if wget --continue --tries=20 --waitretry=5 --timeout=60 --progress=dot:giga \
        --directory-prefix="${RAW_DIR}" "${url}" 2>&1 | tee -a "${LOG_PATH}"; then
        log "[STEP] wget completed: ${file_name}"
    else
        wget_status=$?
        log "[WARN] wget exited with status ${wget_status}; validating the downloaded archive."
        if [[ -s "${compressed_path}" ]] && gzip -t "${compressed_path}" 2>&1 | tee -a "${LOG_PATH}"; then
            log "[STEP] Archive is complete and valid; continuing: ${file_name}"
        else
            log "[ERROR] Download failed or archive is invalid: ${file_name} (wget exit ${wget_status})"
            return "${wget_status}"
        fi
    fi

    if [[ ! -s "${compressed_path}" ]]; then
        log "[ERROR] Download did not produce a non-empty file: ${compressed_path}"
        return 1
    fi
}

log "Fetch output is also being saved to ${LOG_PATH}"
log "[INFO] Project root: ${BASE_PATH}"
log "[INFO] Raw-data folder: ${RAW_DIR}"
log "[INFO] Zenodo record: ${RECORD_ID}"
log "[INFO] Slurm job: ${SLURM_JOB_ID}"

log "[STEP 1/3] Downloading the genome-wide AlphaMissense table."
download_table "AlphaMissense_hg38.tsv.gz"
log "[STEP 2/3] Downloading the amino-acid substitutions table."
download_table "AlphaMissense_aa_substitutions.tsv.gz"

log "[STEP 3/3] Activating operator and preparing both tables."
source "${CONDA_SH}"
conda deactivate || true
conda activate operator

python -u "${SCRIPT_DIR}/../0-fetch_AM.py" "${BASE_PATH}" "${RAW_DIR}" 2>&1 | tee -a "${LOG_PATH}"
log "[DONE] Both AlphaMissense tables were downloaded and prepared."