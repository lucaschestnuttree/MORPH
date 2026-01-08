#!/bin/bash

# -------------------------
# SLURM configuration
# -------------------------

#SBATCH --job-name=gnomXmiss_1file_clin
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view/log/%x_%A_%a.out    
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view/log/%x_%A_%a.err     

#SBATCH --qos=fast
#SBATCH --time=2:00:00               
#SBATCH --ntasks=1            
#SBATCH --cpus-per-task=8     
#SBATCH --mem=128G             
#SBATCH -p ggs

# -------------------------
# User configuration
# -------------------------

# Define paths & parameters

BASE_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view"

GNOM_FILE="${BASE_PATH}/data/2-processed/missense_jul/gnomAD_AF/P15382_gnomAD.tsv"
OUTPUT_DIR="${BASE_PATH}/data/2-processed/amissenseXgnomAD"

VEP_SIF="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/ensemble-VEP/vep.sif"
CACHE_DIR="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/ensemble-VEP/vep_data"
HG38_SUBSTITUTIONS_DIR="${BASE_PATH}/data/2-processed/missense_jul/hg38"
AA_SUBSTITUTIONS_DIR="${BASE_PATH}/data/2-processed/missense_jul/aa_substitutions"
CLINVAR_VCF="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/ensemble-VEP/custom_annotations/clinvar.vcf.gz"

FILTER_SINGLE_BP=1 # Optional flag for filtering single-letter alleles (1 for True, 0 for False)

SCRIPT_PATH="${BASE_PATH}/src/1-data_generation/4-missense_gnomXmiss_1_file.py"

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