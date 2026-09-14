#!/usr/bin/env bash
set -euo pipefail

SRC_ROOT="/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions/100_BACTERIA"
DEST_DIR="/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions/100_BACTERIA_GLOBAL_ANI_RESULTS"

mkdir -p "$DEST_DIR"

LOG_FILE="${DEST_DIR}/copied_global_ani_results_manifest.tsv"

printf "source_file\tcopied_file\n" > "$LOG_FILE"

n_dirs=0
n_files=0

while IFS= read -r -d '' ani_dir; do
    n_dirs=$((n_dirs + 1))

    rel_path="${ani_dir#$SRC_ROOT/}"
    pair_name="${rel_path%%/*}"

    while IFS= read -r -d '' tsv_file; do
        n_files=$((n_files + 1))

        base_name="$(basename "$tsv_file")"
        copied_name="${pair_name}__${base_name}"

        cp "$tsv_file" "${DEST_DIR}/${copied_name}"

        printf "%s\t%s\n" "$tsv_file" "${DEST_DIR}/${copied_name}" >> "$LOG_FILE"

    done < <(find "$ani_dir" -type f -name "*.tsv" -print0)

done < <(find "$SRC_ROOT" -type d -name "GLOBAL_ANI_RESULTS" -print0)

echo "Done."
echo "GLOBAL_ANI_RESULTS folders found: $n_dirs"
echo "TSV files copied: $n_files"
echo "Copied files are in:"
echo "$DEST_DIR"
echo "Manifest written to:"
echo "$LOG_FILE"
