#!/bin/bash

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
FILE="$BASE_PATH/data/1-raw/missense_db/AlphaMissense_hg38.tsv"

UNIPROT_ID="Q4G0S4"

OUT_DIR="$BASE_PATH/log"
mkdir -p "$OUT_DIR"
OUT="${OUT_DIR}/${UNIPROT_ID}_matches.txt"

if [[ ! -f "$FILE" ]]; then
    echo "File not found: $FILE"
    exit 1
fi

echo "Searching for UniProt ID $UNIPROT_ID in $FILE ..."
grep -w "$UNIPROT_ID" "$FILE" | tee "$OUT"

if [[ ! -s "$OUT" ]]; then
    echo "No matches found."
else
    echo "Matches saved to $OUT"
fi
