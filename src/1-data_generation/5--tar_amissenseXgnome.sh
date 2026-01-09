#!/usr/bin/env bash
set -euo pipefail

# =========================
# CONFIG
# =========================
BASE_PATH="/pasteur/helix/projects/GGS_JASS/WKD_Lucas/missense_structure_view"
IN_DIR="${BASE_PATH}/data/2-processed/amissenseXgnomAD"
OUT_DIR="${BASE_PATH}/data/2-processed/amissenseXgnomAD_tars_gz"

N=193
TAG="amissenseXgnomAD"
PAD=4
# =========================

mkdir -p "$OUT_DIR"

# Collect files (safe + reproducible ordering)
mapfile -d '' files < <(
  find "$IN_DIR" -maxdepth 1 -type f -name "*_${TAG}.tsv" -print0 | sort -z
)

total="${#files[@]}"
(( total == 0 )) && { echo "No matching files in $IN_DIR"; exit 1; }

parts=$(( (total + N - 1) / N ))
lookup="$OUT_DIR/lookup_${TAG}.tsv"

# UPDATED HEADER
printf "uniprot_id\ttar_file\ttsv_file\n" > "$lookup"

echo "Found $total files"
echo "Chunk size N=$N → $parts archives"

# Use pigz if available (faster), else gzip
if command -v pigz >/dev/null 2>&1; then
  COMPRESSOR="pigz -p ${SLURM_CPUS_PER_TASK:-4}"
else
  COMPRESSOR="gzip -c"
fi

for ((p=0; p<parts; p++)); do
  start=$(( p * N ))
  end=$(( start + N ))
  (( end > total )) && end=$total

  partno=$(printf "%0*d" "$PAD" $((p+1)))
  ofno=$(printf "%0*d" "$PAD" "$parts")
  base="${TAG}_part${partno}of${ofno}"
  targz="$OUT_DIR/${base}.tar.gz"

  listfile="$(mktemp)"
  : > "$listfile"

  for ((i=start; i<end; i++)); do
    bn="$(basename "${files[$i]}")"
    uid="${bn%%_*}"

    printf "%s\n" "$bn" >> "$listfile"
    printf "%s\t%s\t%s\n" "$uid" "$(basename "$targz")" "$bn" >> "$lookup"
  done

  (cd "$IN_DIR" && tar -cf - -T "$listfile") | eval "$COMPRESSOR" > "$targz"
  rm -f "$listfile"

  echo "Wrote $(basename "$targz") with $((end-start)) files"
done

# Sort lookup by UniProt ID (keep header on top)
tmp_lookup="$(mktemp)"
{
  head -n 1 "$lookup"
  tail -n +2 "$lookup" | sort -t $'\t' -k1,1
} > "$tmp_lookup"
mv "$tmp_lookup" "$lookup"

echo "Done."
echo "Lookup table: $lookup"
