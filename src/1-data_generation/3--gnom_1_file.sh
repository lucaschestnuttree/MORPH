#!/bin/bash

# -------------------------
# SLURM configuration
# -------------------------

#SBATCH --job-name=aaa_GnomAD_parse
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_struct_view/log/%x_%A_%a.out
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_struct_view/log/%x_%A_%a.err

#SBATCH --qos=fast
#SBATCH --time=0:05:00               
#SBATCH --ntasks=1            
#SBATCH --cpus-per-task=8   
#SBATCH --mem=5G             
#SBATCH -p ggs

# -------------------------
# User configuration
# -------------------------

# Define paths
BASE_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_struct_view"

INPUT_FILE="${BASE_PATH}/data/2-processed/missense_jul/hg38/Q8WZ42_hg38.tsv"
HT_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/gnomad.joint.v4.1.sites.ht"
OUTPUT_DIR="${BASE_PATH}/data/2-processed/stuff"
LOG_DIR="/pasteur/appa/scratch/lchataig"

SCRIPT_PATH="${BASE_PATH}/src/1-data_generation/3-gnom_1_file.py"

# -------------------------
# Environment setup
# -------------------------
source /pasteur/appa/homes/lchataig/miniconda3/etc/profile.d/conda.sh
conda deactivate || true
conda activate hail

# -------------------------
# Run
# -------------------------

# Monitor memory and CPU usage during the run
echo "Running on file: $tsv_file"
echo "Memory allocated: $SLURM_MEM_PER_NODE"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

# Run the Python script for the specific file
srun --cpu_bind=cores python "${SCRIPT_PATH}" "$INPUT_FILE" "$HT_PATH" "$OUTPUT_DIR" "$LOG_DIR"
