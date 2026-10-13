#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Base path
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

# ============================================================
# Inputs
# ============================================================
ALPHA="$BASE_PATH/data/1-raw/missense_db/AlphaMissense_hg38.tsv"

HG38_DIR="$BASE_PATH/data/2-processed/missense_jul/hg38"
AA_DIR="$BASE_PATH/data/2-processed/missense_jul/aa_substitutions"
GNOMAD_DIR="$BASE_PATH/data/2-processed/missense_jul/gnomAD_AF"
AXG_DIR="$BASE_PATH/data/2-processed/amissenseXgnomAD"

# ============================================================
# Outputs
# ============================================================
OUT_INDEX="$BASE_PATH/log/index.tsv"
OUT_COUNTS="$BASE_PATH/log/total_counts.txt"
MISSING_AXG="$BASE_PATH/log/missing_amissenseXgnomAD.txt"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

# ------------------------------------------------------------
# 1) Locate uniprot_id column in AlphaMissense
# ------------------------------------------------------------
col_idx="$(
  awk -F'\t' '
    NR==1{
      for(i=1;i<=NF;i++) if($i=="uniprot_id"){ print i; exit }
      exit 1
    }' "$ALPHA"
)" || {
  echo "ERROR: 'uniprot_id' column not found in $ALPHA" >&2
  exit 1
}

# ------------------------------------------------------------
# 2) Master list of unique uniprot_id
# ------------------------------------------------------------
awk -F'\t' -v c="$col_idx" 'NR>1 && $c!="" {print $c}' "$ALPHA" \
  | sort -u > "$TMPDIR/master.ids"

# ------------------------------------------------------------
# 3) Collect uniprot_id from filenames
# ------------------------------------------------------------

# hg38 files: *_hg38.tsv
find "$HG38_DIR" -type f -name "*_hg38.tsv" \
  -printf "%f\n" \
  | awk -F'_' '{print $1}' \
  | sort -u > "$TMPDIR/hg38.ids"

# aa_substitutions: *_aa_substitutions.tsv
find "$AA_DIR" -type f -name "*_aa_substitutions.tsv" \
  -printf "%f\n" \
  | awk -F'_' '{print $1}' \
  | sort -u > "$TMPDIR/aa.ids"

# gnomAD_AF: *_gnomAD.tsv
find "$GNOMAD_DIR" -type f -name "*_gnomAD.tsv" \
  -printf "%f\n" \
  | awk -F'_' '{print $1}' \
  | sort -u > "$TMPDIR/gnomad.ids"

# amissenseXgnomAD: *_amissenseXgnomAD.tsv
if [[ -d "$AXG_DIR" ]]; then
  find "$AXG_DIR" -type f -name "*_amissenseXgnomAD.tsv" \
    -printf "%f\n" \
    | awk -F'_' '{print $1}' \
    | sort -u > "$TMPDIR/axg.ids"
else
  echo "WARNING: $AXG_DIR does not exist; counting amissenseXgnomAD as 0 files." >&2
  : > "$TMPDIR/axg.ids"
fi

# ------------------------------------------------------------
# 3b) Write IDs missing from amissenseXgnomAD (master - axg)
# ------------------------------------------------------------
mkdir -p "$(dirname "$MISSING_AXG")"
comm -23 \
  <(sort -u "$TMPDIR/master.ids") \
  <(sort -u "$TMPDIR/axg.ids") \
  > "$MISSING_AXG"

# ------------------------------------------------------------
# 4) Write index.tsv
# ------------------------------------------------------------
mkdir -p "$(dirname "$OUT_INDEX")"

{
  printf "uniprot_id\tAlphaMissense_hg38\thg38\taa_substitutions\tgnomAD_AF\tamissenseXgnomAD\n"

  awk -v HG38="$TMPDIR/hg38.ids" \
      -v AA="$TMPDIR/aa.ids" \
      -v GN="$TMPDIR/gnomad.ids" \
      -v AXG="$TMPDIR/axg.ids" '
    BEGIN{
      while((getline x < HG38) > 0) hg38[x]=1
      while((getline x < AA)   > 0) aa[x]=1
      while((getline x < GN)   > 0) gn[x]=1
      while((getline x < AXG)  > 0) axg[x]=1
      close(HG38); close(AA); close(GN); close(AXG)
    }
    {
      id=$0
      printf "%s\t1\t%d\t%d\t%d\t%d\n", id, (id in hg38), (id in aa), (id in gn), (id in axg)
    }
  ' "$TMPDIR/master.ids"
} > "$OUT_INDEX"

# ------------------------------------------------------------
# 5) total_counts.txt
# ------------------------------------------------------------
alpha_rows="$(awk 'END{print NR-1}' "$ALPHA")"
alpha_unique="$(wc -l < "$TMPDIR/master.ids" | tr -d ' ')"
hg38_unique="$(wc -l < "$TMPDIR/hg38.ids" | tr -d ' ')"
aa_unique="$(wc -l < "$TMPDIR/aa.ids" | tr -d ' ')"
gn_unique="$(wc -l < "$TMPDIR/gnomad.ids" | tr -d ' ')"
axg_unique="$(wc -l < "$TMPDIR/axg.ids" | tr -d ' ')"
missing_axg_count="$(wc -l < "$MISSING_AXG" | tr -d ' ')"

hg38_in_alpha="$(awk -F'\t' 'NR>1{s+=$3} END{print s+0}' "$OUT_INDEX")"
aa_in_alpha="$(awk -F'\t'   'NR>1{s+=$4} END{print s+0}' "$OUT_INDEX")"
gn_in_alpha="$(awk -F'\t'   'NR>1{s+=$5} END{print s+0}' "$OUT_INDEX")"
axg_in_alpha="$(awk -F'\t'  'NR>1{s+=$6} END{print s+0}' "$OUT_INDEX")"

mkdir -p "$(dirname "$OUT_COUNTS")"
cat > "$OUT_COUNTS" <<EOF
AlphaMissense_hg38.tsv total rows (excluding header):	${alpha_rows}
AlphaMissense unique uniprot_id:	${alpha_unique}

hg38 folder unique uniprot_id (by filenames):	${hg38_unique}
aa_substitutions folder unique uniprot_id (by filenames):	${aa_unique}
gnomAD_AF folder unique uniprot_id (by filenames):	${gn_unique}
amissenseXgnomAD folder unique uniprot_id (by filenames):	${axg_unique}

Missing from amissenseXgnomAD (in AlphaMissense master list):	${missing_axg_count}

Within AlphaMissense master list (column sums in index.tsv):
  AlphaMissense_hg38:	${alpha_unique}
  hg38:	${hg38_in_alpha}
  aa_substitutions:	${aa_in_alpha}
  gnomAD_AF:	${gn_in_alpha}
  amissenseXgnomAD:	${axg_in_alpha}
EOF

echo "Wrote: $OUT_INDEX"
echo "Wrote: $OUT_COUNTS"
echo "Wrote: $MISSING_AXG"
