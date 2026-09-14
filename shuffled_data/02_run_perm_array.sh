#!/bin/bash
#SBATCH --job-name=hv_perm_array
#SBATCH --account=def-acdoxey
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --array=0-24%20
#SBATCH --mail-user=s2mascar@uwaterloo.ca

set -euo pipefail
mkdir -p logs
module load python gcc arrow

python -u 02_run_perm_array.py
