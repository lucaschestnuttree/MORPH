#!/bin/bash

# -------------------------
# SLURM configuration
# -------------------------

#SBATCH --job-name=aaa_GnomAD_parse_array
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view/log/gnomad/%x_%A_%a.out
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view/log/gnomad/%x_%A_%a.err

#SBATCH --qos=fast
#SBATCH --time=0:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH -p ggs

# -------------------------
# User configuration
# -------------------------

# Base project path
BASE_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view"

# Directories and files
TSV_DIR="${BASE_PATH}/data/2-processed/missense_jul/hg38"
OUTPUT_DIR="${BASE_PATH}/data/2-processed/missense_jul/gnomAD_AF"
LOG_DIR="/pasteur/appa/scratch/lchataig"
HT_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/gnomad.joint.v4.1.sites.ht"
SCRIPT_PATH="${BASE_PATH}/src/1-data_generation/3-gnom_1_file.py"

# File that stores the list of UniProt IDs for the array
UNIPROT_IDS_LIST="${BASE_PATH}/src/1-data_generation/gnomad_job_list.txt"

# Max number of concurrent jobs in the array
MAX_CONCURRENT=250

# -------------------------
# Environment setup
# -------------------------

source /pasteur/appa/homes/lchataig/miniconda3/etc/profile.d/conda.sh
conda deactivate || true
conda activate hail

# Ensure output directories exist
mkdir -p "${OUTPUT_DIR}"

# -------------------------
# Job list creation (launcher step)
# -------------------------

# If SLURM_ARRAY_TASK_ID is not set, we are in the launcher step:
# build the list of UniProt IDs that still need processing and resubmit as an array.
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "No SLURM_ARRAY_TASK_ID found – building job list."

    # Collect all UniProt IDs from TSV filenames (e.g., Q8WZ42_hg38.tsv -> Q8WZ42)
    mapfile -t ALL_UNIPROT_IDS < <(find "$TSV_DIR" -name "*.tsv" -printf "%f\n" | cut -d'_' -f1 | sort -u)

    UNIPROT_IDS=()

    for uniprot_id in "${ALL_UNIPROT_IDS[@]}"; do
        output_file="${OUTPUT_DIR}/${uniprot_id}_gnomAD.tsv"
        if [ ! -f "$output_file" ]; then
            UNIPROT_IDS+=("$uniprot_id")
        else
            echo "Output for UniProt ID $uniprot_id already exists. Skipping..."
        fi
    done

    if [ ${#UNIPROT_IDS[@]} -eq 0 ]; then
        echo "No jobs to submit, all outputs already exist."
        exit 0
    fi

    # Save the UniProt IDs list for use by the array tasks
    printf "%s\n" "${UNIPROT_IDS[@]}" > "$UNIPROT_IDS_LIST"
    echo "Saved $((${#UNIPROT_IDS[@]})) UniProt IDs to $UNIPROT_IDS_LIST"

    # Submit the array job
    last_index=$((${#UNIPROT_IDS[@]} - 1))
    echo "Submitting array: 0-${last_index}%${MAX_CONCURRENT}"
    sbatch --array=0-"$last_index"%${MAX_CONCURRENT} "$0"
    exit 0
fi

# -------------------------
# Array task logic
# -------------------------

# Read UniProt IDs from the list file
if [ ! -f "$UNIPROT_IDS_LIST" ]; then
    echo "ERROR: UniProt ID list file not found: $UNIPROT_IDS_LIST"
    exit 1
fi

mapfile -t UNIPROT_IDS < "$UNIPROT_IDS_LIST"

uniprot_id=${UNIPROT_IDS[$SLURM_ARRAY_TASK_ID]}
if [ -z "$uniprot_id" ]; then
    echo "ERROR: No UniProt ID for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Find the corresponding TSV file, assumed pattern: <UNIPROT>_hg38.tsv
input_file=$(find "$TSV_DIR" -name "${uniprot_id}_hg38.tsv" -print -quit)

if [ -f "$input_file" ]; then
    echo "Processing UniProt ID: $uniprot_id"
    echo "Input file: $input_file"
    echo "Memory allocated: $SLURM_MEM_PER_NODE"
    echo "CPUs allocated: $SLURM_CPUS_PER_TASK"

    # Run the Python script for this UniProt ID
    srun --cpu_bind=cores python "${SCRIPT_PATH}" "$input_file" "$HT_PATH" "$OUTPUT_DIR" "$LOG_DIR"
else
    echo "Warning: Input file for UniProt ID $uniprot_id not found in ${TSV_DIR}. Skipping..."
fi
