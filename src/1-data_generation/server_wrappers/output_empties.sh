#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Paths
# ============================================================
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

IDS="$BASE_PATH/log/missing_amissenseXgnomAD.txt"
GNOMAD_DIR="$BASE_PATH/data/2-processed/missense_jul/gnomAD_AF"

EMPTY_DIR="$BASE_PATH/log/z_empty"
MISSING_DIR="$BASE_PATH/log/z_missing"

# ============================================================
# Setup
# ============================================================
mkdir -p "$EMPTY_DIR" "$MISSING_DIR"

# ============================================================
# Sort files by content
# ============================================================
while read -r uniprot_id; do
  [[ -z "$uniprot_id" ]] && continue

  src="${GNOMAD_DIR}/${uniprot_id}_gnomAD.tsv"

  if [[ ! -f "$src" ]]; then
    echo "WARNING: file not found for ${uniprot_id}" >&2
    continue
  fi

  # Count lines (fast, does not read whole file into memory)
  nlines=$(wc -l < "$src")

  if [[ "$nlines" -le 1 ]]; then
    cp "$src" "$EMPTY_DIR/"
  else
    cp "$src" "$MISSING_DIR/"
  fi
done < "$IDS"

echo "Done."
echo "Empty files â†’   $EMPTY_DIR"
echo "Non-empty files â†’ $MISSING_DIR"
