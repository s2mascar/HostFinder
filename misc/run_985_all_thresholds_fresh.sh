#!/bin/bash
#SBATCH --job-name=vh985_fresh
#SBATCH --account=def-acdoxey
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --array=0-71%8
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL

set -euo pipefail
set -x

mkdir -p logs

export PYTHONUNBUFFERED=1


# ============================================================
# 1,296 threshold combinations
#
# 18 combinations per task
# 1,296 / 18 = 72 tasks
#
# SLURM array indexes:
# 0 through 71
# ============================================================

COMBOS_PER_TASK=18

TASK_ID="${SLURM_ARRAY_TASK_ID}"

COMBO_START=$((TASK_ID * COMBOS_PER_TASK))

COMBO_END=$((COMBO_START + COMBOS_PER_TASK))


# ============================================================
# NODE-LOCAL DUCKDB TEMP DIRECTORY
# ============================================================

TEMP_DIR="${SLURM_TMPDIR:-$PWD}/duckdb_temp_985_${TASK_ID}"

mkdir -p "${TEMP_DIR}"


echo "============================================================"
echo "TASK:        ${TASK_ID}"
echo "COMBO START: ${COMBO_START}"
echo "COMBO END:   ${COMBO_END}"
echo "TEMP DIR:    ${TEMP_DIR}"
echo "============================================================"


python -u rebuild_985_all_thresholds_fresh.py \
    --input-dir virushostdb_985_fresh_build \
    --output-dir virushostdb_985_fresh_build/all_threshold_scores \
    --combo-start "${COMBO_START}" \
    --combo-end "${COMBO_END}" \
    --threads 8 \
    --memory-limit 120GB \
    --temp-dir "${TEMP_DIR}"
