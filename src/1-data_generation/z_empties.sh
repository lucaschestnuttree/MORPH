#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Paths
# ============================================================
BASE_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view"

IDS="${BASE_PATH}/data/2-processed/missing_amissenseXgnomAD.txt"
GNOMAD_DIR="${BASE_PATH}/data/2-processed/missense_jul/gnomAD_AF"

EMPTY_DIR="${BASE_PATH}/data/2-processed/z_empty"
MISSING_DIR="${BASE_PATH}/data/2-processed/z_missing"

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
echo "Empty files →   $EMPTY_DIR"
echo "Non-empty files → $MISSING_DIR"
