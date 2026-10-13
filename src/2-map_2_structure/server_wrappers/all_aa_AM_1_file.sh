#!/bin/bash

# -------------------------
# SLURM configuration (submit from the repository root with sbatch)
# -------------------------
#SBATCH --job-name=am_pathogenicity_2_structure
#SBATCH --output=log/%x_%j.out
#SBATCH --error=log/%x_%j.err

#SBATCH --qos=fast
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -p ggs

# -------------------------
# User configuration
# -------------------------
# Resolve project paths from the original script location.
# Slurm can fall back to the repository root set by SLURM_SUBMIT_DIR.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)" || exit 1
if [[ ! -f "${SCRIPT_DIR}/__config.sh" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    if [[ -f "${SLURM_SUBMIT_DIR}/__config.sh" ]]; then
        SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
    else
        SCRIPT_DIR="${SLURM_SUBMIT_DIR}/src/2-map_2_structure/server_wrappers"
    fi
fi
BASE_PATH="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)" || exit 1
source "${SCRIPT_DIR}/__config.sh" || exit 1
UNIPROT_ID="P10635"

INPUT_TSV="${BASE_PATH}/data/2-processed/missense_jul/aa_substitutions/${UNIPROT_ID}_aa_substitutions.tsv"
INPUT_PDB="${BASE_PATH}/data/2-processed/pdb/${UNIPROT_ID}.pdb"
OUTPUT_PDB="${BASE_PATH}/data/2-processed/stuff/${UNIPROT_ID}_am_pathogenicity_bfactor.pdb"

SCRIPT_PATH="${SCRIPT_DIR}/../all_aa_AM_2_structure.py"

# -------------------------
# Environment setup
# -------------------------
source "${CONDA_SH}" || exit 1
conda deactivate || true
conda activate pymut-env

mkdir -p "$(dirname "$OUTPUT_PDB")"
mkdir -p "${BASE_PATH}/log"

# -------------------------
# Run
# -------------------------
srun python "$SCRIPT_PATH" \
  "$INPUT_TSV" \
  "$INPUT_PDB" \
  "$OUTPUT_PDB"