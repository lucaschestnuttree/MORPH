#!/bin/bash

# -------------------------
# SLURM configuration
# -------------------------
#SBATCH --job-name=am_pathogenicity_2_structure
#SBATCH --output=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%j.out
#SBATCH --error=/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view/log/%x_%j.err

#SBATCH --qos=fast
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -p ggs

# -------------------------
# User configuration
# -------------------------
BASE_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view"
UNIPROT_ID="P10635"

INPUT_TSV="${BASE_PATH}/data/2-processed/amissenseXgnomAD/${UNIPROT_ID}_amissenseXgnomAD.tsv"
INPUT_PDB="${BASE_PATH}/data/2-processed/pdb/${UNIPROT_ID}.pdb"
OUTPUT_PDB="${BASE_PATH}/data/2-processed/stuff/${UNIPROT_ID}_am_pathogenicity_bfactor.pdb"

SCRIPT_PATH="${BASE_PATH}/src/2-mapping/map_2_structure/your_script.py"

# -------------------------
# Environment setup
# -------------------------
source /pasteur/appa/homes/lchataig/miniconda3/etc/profile.d/conda.sh
conda deactivate || true
conda activate pymut-env

mkdir -p "$(dirname "$OUTPUT_PDB")"
mkdir -p "${BASE_PATH}/log"

# -------------------------
# Run
# -------------------------
srun python "$SCRIPT_PATH" \
  "$INPUT_TSV" \
  "$INPUT_PDB" \
  "$OUTPUT_PDB"