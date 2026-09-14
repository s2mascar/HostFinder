#!/bin/bash
#SBATCH --job-name=virus_meta_chunked
#SBATCH --account=def-acdoxey
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL

set -euo pipefail
set -x

mkdir -p logs

# Load the same modules you use for DuckDB/Python/Arrow on the cluster.
# If your environment already has these packages, this line can be removed.
module load python gcc arrow || true

python -u 04_recover_shared_accession_metadata_chunked.py \
  --host-glob "/home/smascar/scratch/VIRUSES/STAT_with_Eukaryota_important_columns_split/**" \
  --pathogen-glob "/home/smascar/scratch/VIRUSES/STAT_with_Viruses_important_columns_split/**" \
  --meta-glob "/home/smascar/projects/def-acdoxey/SRA_DATA/DATA/02_02_2026/SRA_META_DATA-02_02_2026/**" \
  --pairs-csv "virus_host_db_labelled_interactions.csv" \
  --pair-metrics-file "pair_threshold_grid_results/pair_metrics_all_thresholds.parquet" \
  --outdir "/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions/metadata_recovery_virus_hostdb_ht1e-04_pt1e-08" \
  --host-threshold 1e-4 \
  --pathogen-threshold 1e-8 \
  --threads 8 \
  --memory-limit "80GB" \
  --n-acc-blocks 512

# Optional flags you can add above if needed:
#   --overwrite-parts              Rebuild chunk outputs from scratch.
#   --write-shared-acc-parts       Save per-accession shared-accession block files. Uses lots of disk.
#   --write-wide-summary           Also make shared_meta_summary_wide.csv/parquet. More expensive.
#   --use-pair-metrics-filter      Keep only pairs present in pair_metrics at this threshold.
