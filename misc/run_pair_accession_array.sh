#!/bin/bash
#SBATCH --job-name=pair_acc_array
#SBATCH --account=def-acdoxey
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --array=1-40%4
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL

set -euo pipefail

mkdir -p logs

PAIRS_CSV="100_Bacteria_05_21.csv"

#1000_HOST_BACTERIA_INTERACTIONS_05_20.csv"

#100_host_virus_interactions_for_FDR_trial_1.csv"

#top_1000_host_virus_predictions.csv"
#"top_1000_host_virus_predictions.csv"
OUTDIR="/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions/100_Bacteria_05_21"

#1000_HOST_BACTERIA_INTERACTIONS_05_20"

#100_host_virus_interactions_for_FDR_trial_1"

#/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions/top1000_pair_accessions_array"

BATCH_SIZE=25

TASK_ID="${SLURM_ARRAY_TASK_ID}"

TMPDIR_TASK="/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions/duckdb_tmp_files/duckdb_tmp_pair_acc_${SLURM_JOB_ID}_${TASK_ID}"
mkdir -p "${TMPDIR_TASK}"

python build_pair_accession_one_batch.py \
  --pairs-csv "${PAIRS_CSV}" \
  --outdir "${OUTDIR}" \
  --batch-id "${TASK_ID}" \
  --batch-size "${BATCH_SIZE}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --memory-limit "30GB" \
  --temp-dir "${TMPDIR_TASK}"

rm -rf "${TMPDIR_TASK}"
