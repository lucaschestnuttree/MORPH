#!/bin/bash

# -------------------------
# SLURM configuration (submit from the repository root with sbatch)
# -------------------------

#SBATCH --job-name=gnomad_af_2_structure_1file
#SBATCH --output=log/%x_%A_%a.out
#SBATCH --error=log/%x_%A_%a.err

#SBATCH --qos=fast
#SBATCH --time=0:50:00               
#SBATCH --ntasks=1            
#SBATCH --cpus-per-task=1     
#SBATCH --mem=10G             
#SBATCH -p ggs

# -------------------------
# User configuration
# -------------------------

# Define paths
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

# Define directories and paths
INPUT_FILE="${BASE_PATH}/data/2-processed/amissenseXgnomAD/${UNIPROT_ID}_amissenseXgnomAD.tsv"
PDB_FILE="${BASE_PATH}/data/2-processed/pdb/${UNIPROT_ID}.pdb"
OUTPUT_DIR="${BASE_PATH}/data/2-processed/stuff"

SCRIPT_PATH="${SCRIPT_DIR}/../pymut_gnomad_af_2_structure.py"

# -------------------------
# Environment setup
# -------------------------

# Activate conda environment
source "${CONDA_SH}" || exit 1
conda deactivate || true
conda activate pymut-env

# -------------------------
# Run
# -------------------------
srun python "$SCRIPT_PATH" \
  "$INPUT_FILE" \
  "$PDB_FILE" \
  "$OUTPUT_DIR" \
