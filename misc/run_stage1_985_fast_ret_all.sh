#!/bin/bash
#SBATCH --job-name=stage1_985_acc
#SBATCH --account=def-acdoxey
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL

set -euo pipefail

mkdir -p logs

export PYTHONUNBUFFERED=1


# ============================================================
# CONFIGURATION
# ============================================================

SCORED_PARQUET="virus_host_db_true_pairs_detected_985_thresholds.parquet"

CHUNK_SIZE=10

REPEATS=100

N_COMBOS=1296


# ============================================================
# NEW SHARED DIRECTORY FOR THE 985-PAIR BENCHMARK
# ============================================================

SHARED_DIR="stage1_985_fast_ret_shared"

GRID_CSV="${SHARED_DIR}/00_threshold_combinations_FULL_GRID.csv"

BASE_UNIVERSE="${SHARED_DIR}/02_fixed_baseline_universe.parquet"


mkdir -p "${SHARED_DIR}"


# ============================================================
# PREPARE THRESHOLD GRID AND FIXED BASELINE UNIVERSE ONCE
# ============================================================

if [[ ! -f "${GRID_CSV}" || ! -f "${BASE_UNIVERSE}" ]]; then

    echo "============================================================"
    echo "PREPARING NEW 985-PAIR STAGE 1 BASELINE"
    echo "============================================================"

    python -u raw_only_threshold_scan_FAST_WITH_RETENTION.py \
      --scored-parquet "${SCORED_PARQUET}" \
      --repeats "${REPEATS}" \
      --raw-score-sources shared_accessions,log_odds_acc \
      --threads 4 \
      --memory-limit 80GB \
      --temp-dir duckdb_temp_stage1_985_prepare \
      --outdir "${SHARED_DIR}" \
      --base-universe-parquet "${BASE_UNIVERSE}" \
      --prepare-only

fi


# ============================================================
# DISPLAY BASELINE UNIVERSE STATS
# ============================================================

echo
echo "============================================================"
echo "BASELINE UNIVERSE STATS"
echo "============================================================"

cat "${SHARED_DIR}/03_fixed_baseline_universe_stats.csv"

echo


# ============================================================
# SCAN ALL 1,296 THRESHOLD COMBINATIONS
# ============================================================

for START in $(seq 0 ${CHUNK_SIZE} $((N_COMBOS - 1))); do

    END=$((START + CHUNK_SIZE))

    CHUNK_ID=$((START / CHUNK_SIZE))


    echo
    echo "============================================================"
    echo "RUNNING THRESHOLD CHUNK ${CHUNK_ID}"
    echo "COMBINATIONS ${START} TO ${END}"
    echo "============================================================"


    python -u raw_only_threshold_scan_FAST_WITH_RETENTION.py \
      --scored-parquet "${SCORED_PARQUET}" \
      --threshold-grid-csv "${GRID_CSV}" \
      --base-universe-parquet "${BASE_UNIVERSE}" \
      --repeats "${REPEATS}" \
      --raw-score-sources shared_accessions,log_odds_acc \
      --combo-start "${START}" \
      --combo-end "${END}" \
      --threads 4 \
      --memory-limit 80GB \
      --temp-dir "duckdb_temp_stage1_985_fast_ret_${CHUNK_ID}" \
      --outdir "stage1_985_fast_ret_chunk_${CHUNK_ID}"

done


echo
echo "============================================================"
echo "STAGE 1 985-PAIR THRESHOLD SCAN COMPLETE"
echo "============================================================"
