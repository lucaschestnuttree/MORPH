#!/bin/bash

# -------------------------
# SLURM configuration (submit with sbatch from this directory)
# -------------------------

#SBATCH --job-name=gnomXmiss_1file_clin
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%A_%a.out
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%A_%a.err

#SBATCH --qos=fast
#SBATCH --time=2:00:00               
#SBATCH --ntasks=1            
#SBATCH --cpus-per-task=8     
#SBATCH --mem=128G             
#SBATCH -p ggs

set -euo pipefail

# -------------------------
# User configuration
# -------------------------

# Define paths & parameters

# Work from this wrapper's folder; recover it for Slurm's temporary copy.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
if [[ ! -f __config.sh && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    cd "$SLURM_SUBMIT_DIR" || exit 1
    # Submission may come from the project root or any of its subfolders.
    while [[ ! -d src/1-data_generation/server_wrappers ]]; do
        if [[ "$PWD" == / ]]; then
            echo "ERROR: Cannot find the project from SLURM_SUBMIT_DIR=$SLURM_SUBMIT_DIR" >&2
            exit 1
        fi
        cd .. || exit 1
    done
    cd src/1-data_generation/server_wrappers || exit 1
fi
SCRIPT_DIR="$PWD"
BASE_PATH="$(cd ../../.. && pwd)" || exit 1
source ./__config.sh || exit 1

GNOM_FILE="$BASE_PATH/data/2-processed/missense_jul/gnomAD_AF/P11712_gnomAD.tsv"
OUTPUT_DIR="$BASE_PATH/data/2-processed/amissenseXgnomAD"

HG38_SUBSTITUTIONS_DIR="$BASE_PATH/data/2-processed/missense_jul/hg38"
AA_SUBSTITUTIONS_DIR="$BASE_PATH/data/2-processed/missense_jul/aa_substitutions"

FILTER_SINGLE_BP=1 # Optional flag for filtering single-letter alleles (1 for True, 0 for False)

SCRIPT_PATH="../3-missense_gnomXmiss_1_file.py"

# -------------------------
# Environment setup
# -------------------------

# Load modules
module load apptainer

# Activate conda environment
source "${CONDA_SH}" || exit 1
conda deactivate || true
conda activate hail

# -------------------------
# Run
# -------------------------
srun python "$SCRIPT_PATH" \
  "$GNOM_FILE" \
  "$OUTPUT_DIR" \
  "$VEP_SIF" \
  "$CACHE_DIR" \
  "$HG38_SUBSTITUTIONS_DIR" \
  "$AA_SUBSTITUTIONS_DIR" \
  "$CLINVAR_VCF" \
  "$FILTER_SINGLE_BP"

echo "Job completed!"