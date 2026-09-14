#!/bin/bash
#SBATCH --job-name=hv_combine_plot
#SBATCH --account=def-acdoxey
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user=s2mascar@uwaterloo.ca

set -euo pipefail
mkdir -p logs
module load StdEnv/2023 r/4.5.0 python gcc arrow

python -u 03_combine_perm_hists.py
Rscript 04_plot_observed_vs_null_bins.R
python -u 05_export_hist_bin_counts.py
