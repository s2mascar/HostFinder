#!/bin/bash
# Submit the whole shuffled-null pipeline with Slurm dependencies.
# Run this from the directory containing all scripts:
#   bash submit_shuffled_pipeline.sh

set -euo pipefail
mkdir -p logs

prep_job=$(sbatch --parsable 01_prepare_observed_and_baseline.sh)
echo "Submitted prepare job: ${prep_job}"

perm_job=$(sbatch --parsable --dependency=afterok:${prep_job} 02_run_perm_array.sh)
echo "Submitted permutation array job: ${perm_job}"

combine_job=$(sbatch --parsable --dependency=afterok:${perm_job} 03_combine_and_plot.sh)
echo "Submitted combine/plot job: ${combine_job}"

echo "Pipeline submitted. Final outputs will be under: host_virus_abundance_shuffle_null_ht1e-04_pt1e-08/tables"
