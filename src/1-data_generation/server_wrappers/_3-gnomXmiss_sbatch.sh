#!/bin/bash

# -------------------------
# SLURM configuration (submit with sbatch from this directory)
# -------------------------

#SBATCH --job-name=missenseXgnomAD
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/gnomXmiss/%x_%A_%a.out
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/gnomXmiss/%x_%A_%a.err

#SBATCH --qos=fast
#SBATCH --time=2:00:00               
#SBATCH --ntasks=1            
#SBATCH --cpus-per-task=8     
#SBATCH --mem=64G             
#SBATCH -p common
# -------------------------
# User configuration
# -------------------------

# Paths and variables
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

GNOM_DIR="$BASE_PATH/data/2-processed/missense_jul/gnomAD_AF"
OUTPUT_DIR="$BASE_PATH/data/2-processed/amissenseXgnomAD"

HG38_SUBSTITUTIONS_DIR="$BASE_PATH/data/2-processed/missense_jul/hg38"
AA_SUBSTITUTIONS_DIR="$BASE_PATH/data/2-processed/missense_jul/aa_substitutions"

FILTER_SINGLE_BP=1

SCRIPT_PATH="../3-missense_gnomXmiss_1_file.py"

# Path for the uniprot IDs list file
UNIPROT_IDS_LIST="../uniprot_ids_list.txt"

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

# If SLURM_ARRAY_TASK_ID is not set, generate the list of jobs to run
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    # Find all UniProt IDs from the GNOM_DIR
    ALL_UNIPROT_IDS=($(find "$GNOM_DIR" -name "*.tsv" | sed 's#.*/##' | cut -d'_' -f1))
    
    # Initialize an empty array for UniProt IDs that need processing
    UNIPROT_IDS=()

    # Loop through the list and filter out those that already have output
    for uniprot_id in "${ALL_UNIPROT_IDS[@]}"; do
        output_file="${OUTPUT_DIR}/${uniprot_id}_amissenseXgnomAD.tsv"
        if [ ! -f "$output_file" ]; then
            UNIPROT_IDS+=("$uniprot_id")
        else
            echo "Output for UniProt ID $uniprot_id already exists. Skipping..."
        fi
    done

    # If there are no UniProt IDs left after filtering, exit
    if [ ${#UNIPROT_IDS[@]} -eq 0 ]; then
        echo "No jobs to submit, all outputs already exist."
        exit 0
    fi

    # Save the UNIPROT_IDS array to the specified file
    echo "${UNIPROT_IDS[@]}" > "$UNIPROT_IDS_LIST"

    # Submit the job array with concurrency limit (n jobs at a time)
    mkdir -p "$BASE_PATH/log/gnomXmiss"
    (cd -- "${BASE_PATH}" && sbatch --chdir="${BASE_PATH}" \
        --array=0-$((${#UNIPROT_IDS[@]}-1))%500 "$SCRIPT_DIR/_3-gnomXmiss_sbatch.sh")
    exit "$?"
fi

# Array job logic: Read UniProt IDs from the saved file
UNIPROT_IDS=($(cat "$UNIPROT_IDS_LIST"))
uniprot_id=${UNIPROT_IDS[$SLURM_ARRAY_TASK_ID]}
gnom_file=$(find "$GNOM_DIR" -name "${uniprot_id}_gnomAD.tsv")

# Check if the input file exists and run the job
if [ -f "$gnom_file" ]; then
    echo "Running for UniProt ID: $uniprot_id"
    echo "GnomAD file: $gnom_file"
    echo "Memory allocated: $SLURM_MEM_PER_NODE"
    echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
    
    # Run the Python script for this UniProt ID
    python "$SCRIPT_PATH" \
      "$gnom_file" \
      "$OUTPUT_DIR" \
      "$VEP_SIF" \
      "$CACHE_DIR" \
      "$HG38_SUBSTITUTIONS_DIR" \
      "$AA_SUBSTITUTIONS_DIR" \
      "$CLINVAR_VCF" \
      "$FILTER_SINGLE_BP"
else
    echo "Warning: GnomAD file for UniProt ID $uniprot_id not found. Skipping..."
fi