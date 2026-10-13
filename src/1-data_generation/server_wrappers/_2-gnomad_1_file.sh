#!/bin/bash

# -------------------------
# SLURM configuration (submit with sbatch from this directory)
# -------------------------

#SBATCH --job-name=aaa_gnomAD_parse
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%A_%a.out
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%A_%a.err

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

INPUT_FILE="$BASE_PATH/data/2-processed/missense_jul/hg38/Q8WZ42_hg38.tsv"
OUTPUT_DIR="$BASE_PATH/data/2-processed/stuff"

SCRIPT_PATH="../2-gnomad_1_file.py"

# -------------------------
# Environment setup
# -------------------------
source "${CONDA_SH}" || exit 1
conda deactivate || true
conda activate hail

# -------------------------
# Run
# -------------------------

# Monitor memory and CPU usage during the run
echo "Running on file: $INPUT_FILE"
echo "Memory allocated: $SLURM_MEM_PER_NODE"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

# Run the Python script for the specific file
srun --cpu_bind=cores python "${SCRIPT_PATH}" "$INPUT_FILE" "$HT_PATH" "$OUTPUT_DIR" "$LOG_DIR"
