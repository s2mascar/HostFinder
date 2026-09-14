#!/bin/bash
#SBATCH --job-name=pair_acc_finalize
#SBATCH --account=def-acdoxey
#SBATCH --time=00:01:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/pair_acc_finalize_%j.out
#SBATCH --error=logs/pair_acc_finalize_%j.err

set -euo pipefail

mkdir -p logs

PAIRS_CSV="100_host_virus_interactions_for_FDR_trial_2.csv"
#"top_1000_host_virus_predictions.csv"
OUTDIR="/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions/100_host_virus_interactions_for_FDR_trial_2"

#top1000_pair_accessions_array"

TMPDIR_TASK="/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions/duckdb_temp_files/duckdb_tmp_pair_acc_finalize_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR_TASK}"

python finalize_pair_accession_batches.py \
  --pairs-csv "${PAIRS_CSV}" \
  --outdir "${OUTDIR}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --memory-limit "190GB" \
  --temp-dir "${TMPDIR_TASK}" \
  --overwrite

rm -rf "${TMPDIR_TASK}"
