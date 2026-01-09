#!/bin/bash

BASE_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view"
FILE="${BASE_PATH}/data/1-raw/missense_db/AlphaMissense_hg38.tsv"

UNIPROT_ID="Q4G0S4"

OUT="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/miss_struct_view/src/1-data_generation/${UNIPROT_ID}_matches.txt"

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
