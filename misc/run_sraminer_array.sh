#!/bin/bash
#SBATCH --job-name=sraminer_asm
#SBATCH --account=def-acdoxey
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL


set -euo pipefail

mkdir -p logs

MANIFEST="100_Bacteria_05_21_inputs/sraminer_jobs_manifest.csv"
#1000_HOST_BACTERIA_INTERACTIONS_05_20_inputs/sraminer_jobs_manifest.csv"
#100_host_virus_interactions_for_FDR_trial_1_inputs/sraminer_jobs_manifest.csv"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "PWD: $(pwd)"
echo "Manifest: ${MANIFEST}"
echo

python run_one_sraminer_manifest_job.py \
  --manifest "${MANIFEST}" \
  --array-id "${SLURM_ARRAY_TASK_ID}"
