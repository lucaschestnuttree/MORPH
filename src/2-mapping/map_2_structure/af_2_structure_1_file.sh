#!/bin/bash

# -------------------------
# SLURM configuration
# -------------------------

#SBATCH --job-name=gnomad_af_2_structure_1file
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%A_%a.out
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%A_%a.err

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
BASE_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view"
UNIPROT_ID="P10635"

# Define directories and paths
INPUT_FILE="${BASE_PATH}/data/2-processed/amissenseXgnomAD/${UNIPROT_ID}_amissenseXgnomAD.tsv"
PDB_FILE="${BASE_PATH}/data/2-processed/pdb/${UNIPROT_ID}.pdb"
OUTPUT_DIR="${BASE_PATH}/data/2-processed/stuff"

SCRIPT_PATH="${BASE_PATH}/src/2-mapping/map_2_structure/pymut_gnomad_af_2_structure.py"

# -------------------------
# Environment setup
# -------------------------

# Activate conda environment
source /pasteur/appa/homes/lchataig/miniconda3/etc/profile.d/conda.sh
conda deactivate || true
conda activate pymut-env

# -------------------------
# Run
# -------------------------
srun python "$SCRIPT_PATH" \
  "$INPUT_FILE" \
  "$PDB_FILE" \
  "$OUTPUT_DIR" \
