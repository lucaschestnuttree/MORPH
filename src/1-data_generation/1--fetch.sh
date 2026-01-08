#!/bin/bash -l

# -------------------------
# User configuration
# -------------------------

# Define paths
BASE_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view"

SCRIPT_PATH="${BASE_PATH}/src/1-data_generation/1-fetch_data.py"

# -------------------------
# Environment setup
# -------------------------
source /pasteur/appa/homes/lchataig/miniconda3/etc/profile.d/conda.sh
conda deactivate || true
conda activate operator

# -------------------------
# Run
# -------------------------

python "${SCRIPT_PATH}" "${BASE_PATH}"
