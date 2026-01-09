#!/bin/bash

# -------------------------
# SLURM configuration
# -------------------------

#SBATCH --job-name=missense_split_uniprot
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_struct_view/log/%x_%A_%a.out
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_struct_view/log/%x_%A_%a.err

#SBATCH --qos=fast
#SBATCH --time=01:00:00          
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G                
#SBATCH -p ggs


# -------------------------
# User configuration
# -------------------------

# Define paths
BASE_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_struct_view"

TSV_PATH="${BASE_PATH}/data/1-raw/missense_db/AlphaMissense_aa_substitutions.tsv"

SCRIPT_PATH="${BASE_PATH}/src/1-data_generation/2-parse_AM_data.py"

# -------------------------
# Environment setup
# -------------------------
source /pasteur/appa/homes/lchataig/miniconda3/etc/profile.d/conda.sh
conda deactivate || true
conda activate operator

# -------------------------
# Run
# -------------------------

echo "Starting missense split-by-uniprot job"
echo "BASE_PATH: ${BASE_PATH}"
echo "TSV_PATH:  ${TSV_PATH}"
export POLARS_MAX_THREADS="${SLURM_CPUS_PER_TASK:-16}"

srun python "${SCRIPT_PATH}" "${BASE_PATH}" "${TSV_PATH}"

echo "Job completed!"
