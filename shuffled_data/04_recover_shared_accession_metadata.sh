#!/bin/bash
#SBATCH --job-name=virus_meta_recovery
#SBATCH --account=def-acdoxey
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL

set -euo pipefail
set -x

mkdir -p logs

python -u 04_recover_shared_accession_metadata.py \
  --host-glob "/home/smascar/scratch/VIRUSES/STAT_with_Eukaryota_important_columns_split/**" \
  --pathogen-glob "/home/smascar/scratch/VIRUSES/STAT_with_Viruses_important_columns_split/**" \
  --meta-glob "/home/smascar/projects/def-acdoxey/SRA_DATA/DATA/02_02_2026/SRA_META_DATA-02_02_2026/**" \
  --pairs-csv "virus_host_db_labelled_interactions.csv" \
  --pair-metrics-file "pair_threshold_grid_results/pair_metrics_all_thresholds.parquet" \
  --outdir "/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions/metadata_recovery_virus_hostdb_ht1e-04_pt1e-08" \
  --host-threshold 1e-4 \
  --pathogen-threshold 1e-8 \
  --threads 8 \
  --memory-limit "80GB"

# Optional flags you can add above if needed:
#   --overwrite
#       Deletes the output directory before running.
#
#   --use-pair-metrics-filter
#       Keeps only pairs present in pair_threshold_grid_results/pair_metrics_all_thresholds.parquet
#       at host_threshold=1e-4 and pathogen_threshold=1e-8.
#       Leave this off if that file does not contain this threshold pair.
#
#   --write-large-csv
#       Also writes the large per-accession CSV files.
#       Parquet files are always written; CSV can be huge.
