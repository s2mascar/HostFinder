#!/bin/bash
#SBATCH --job-name=hv_prepare
#SBATCH --account=def-acdoxey
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL

set -euo pipefail
mkdir -p logs

python -u 01_prepare_observed_and_baseline.py
