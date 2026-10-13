#!/usr/bin/env bash

#SBATCH --job-name=missense_parse_AM
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%j.out
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%j.err
#SBATCH --qos=fast
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -p ggs

set -euo pipefail

# Slurm copies this file to a temporary directory. Recover the repository path
# from SLURM_SUBMIT_DIR when that copy can no longer see __config.sh beside it.
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
        echo "ERROR: Could not locate the repository from SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}" >&2
        exit 1
    fi
fi
BASE_PATH="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"

# Submit this wrapper with sbatch from the server_wrappers directory.
if [[ -z "${SLURM_JOB_ID:-${SLURM_JOBID:-}}" ]]; then
    echo "ERROR: Run sbatch _1-parse_AM_data.sh from the server_wrappers directory." >&2
    exit 1
fi

source "${SCRIPT_DIR}/__config.sh"
RAW_DIR="${BASE_PATH}/data/1-raw/missense_db"
for TABLE in AlphaMissense_aa_substitutions.tsv AlphaMissense_hg38.tsv; do
    if [[ ! -f "${RAW_DIR}/${TABLE}" ]]; then
        echo "ERROR: Missing ${RAW_DIR}/${TABLE}; run sbatch _0-fetch.sh and wait for it to finish first." >&2
        exit 1
    fi
done

source "${CONDA_SH}"
conda deactivate || true
conda activate operator
export POLARS_MAX_THREADS="${SLURM_CPUS_PER_TASK:-16}"

srun python -u "${SCRIPT_DIR}/../1-parse_AM_data.py" "${BASE_PATH}" "$@"
echo "Parsing job completed."