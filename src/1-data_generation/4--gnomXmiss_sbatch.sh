#!/bin/bash

# -------------------------
# SLURM configuration
# -------------------------

#SBATCH --job-name=missenseXgnomAD
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view/log/gnomXmiss/%x_%A_%a.out    
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view/log/gnomXmiss/%x_%A_%a.err     

#SBATCH --qos=fast
#SBATCH --time=2:00:00               
#SBATCH --ntasks=1            
#SBATCH --cpus-per-task=8     
#SBATCH --mem=64G             
#SBATCH -p ggs
# -------------------------
# User configuration
# -------------------------

# Paths and variables
BASE_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view"

GNOM_DIR="${BASE_PATH}/data/2-processed/missense_jul/gnomAD_AF"
OUTPUT_DIR="${BASE_PATH}/data/2-processed/amissenseXgnomAD"

VEP_SIF="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/ensemble-VEP/vep.sif"
CACHE_DIR="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/ensemble-VEP/vep_data"
HG38_SUBSTITUTIONS_DIR="${BASE_PATH}/data/2-processed/missense_jul/hg38"
AA_SUBSTITUTIONS_DIR="${BASE_PATH}/data/2-processed/missense_jul/aa_substitutions"
CLINVAR_VCF="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/ensemble-VEP/custom_annotations/clinvar.vcf.gz"

FILTER_SINGLE_BP=1

SCRIPT_PATH="${BASE_PATH}/src/1-data_generation/4-missense_gnomXmiss_1_file.py"

# Path for the uniprot IDs list file
UNIPROT_IDS_LIST="${BASE_PATH}/src/1-data_generation/uniprot_ids_list.txt"

# -------------------------
# Environment setup
# -------------------------

# Load modules
module load apptainer

# Activate conda environment
source /pasteur/appa/homes/lchataig/miniconda3/etc/profile.d/conda.sh
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
    sbatch --array=0-$((${#UNIPROT_IDS[@]}-1))%500 "$0"
    exit 0
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