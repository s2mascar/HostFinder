#!/bin/bash

#SBATCH --job-name=bioresearch_env_v05
#SBATCH --account=def-acdoxey
#SBATCH --time=01:30:00
#SBATCH --gres=gpu:nvidia_h100_80gb_hbm3_3g.40gb:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48000M
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL


set -euo pipefail


# ============================================================
# ALWAYS RUN FROM THE DIRECTORY CONTAINING THIS SCRIPT
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"


# ============================================================
# PYTHON ENVIRONMENT
# ============================================================

source agent_env/bin/activate


# ============================================================
# OUTPUT DIRECTORIES
# ============================================================

mkdir -p logs
mkdir -p results
mkdir -p data


# ============================================================
# LOG EVERYTHING
# ============================================================

LOG_FILE="logs/bioresearch_env_v05_${SLURM_JOB_ID}.log"
exec > >(tee -a "$LOG_FILE") 2>&1


echo "============================================================"
echo "BIORESEARCHENV v0.5 SLURM JOB"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Working directory: $(pwd)"
echo "Node: $(hostname)"
echo "Start time: $(date)"


echo
echo "============================================================"
echo "GPU"
echo "============================================================"
nvidia-smi


echo
echo "============================================================"
echo "PYTHON"
echo "============================================================"
which python
python --version


# ============================================================
# VIRUS-HOST DB
# ============================================================
# The Biological Context Agent can run without this file, but host-range
# context is much more useful when it is present.
# ============================================================

VIRUSHOST_DB="data/virushostdb.daily.tsv"
VIRUSHOST_URL="https://www.genome.jp/ftp/db/virushostdb/virushostdb.daily.tsv"

if [ ! -s "$VIRUSHOST_DB" ]; then

    echo
    echo "============================================================"
    echo "DOWNLOADING VIRUS-HOST DB"
    echo "============================================================"

    TMP_DB="${VIRUSHOST_DB}.tmp"
    rm -f "$TMP_DB"

    if command -v wget >/dev/null 2>&1; then
        if wget -q -O "$TMP_DB" "$VIRUSHOST_URL"; then
            mv "$TMP_DB" "$VIRUSHOST_DB"
            echo "Virus-Host DB downloaded successfully."
        else
            rm -f "$TMP_DB"
            echo "WARNING: Virus-Host DB download failed."
            echo "The benchmark will continue with NCBI taxonomy context only."
        fi

    elif command -v curl >/dev/null 2>&1; then
        if curl -L --fail -sS "$VIRUSHOST_URL" -o "$TMP_DB"; then
            mv "$TMP_DB" "$VIRUSHOST_DB"
            echo "Virus-Host DB downloaded successfully."
        else
            rm -f "$TMP_DB"
            echo "WARNING: Virus-Host DB download failed."
            echo "The benchmark will continue with NCBI taxonomy context only."
        fi

    else
        echo "WARNING: Neither wget nor curl is available."
        echo "The benchmark will continue with NCBI taxonomy context only."
    fi

else
    echo
    echo "Virus-Host DB found: $VIRUSHOST_DB"
fi


# ============================================================
# SYNTAX CHECK
# ============================================================

echo
echo "============================================================"
echo "PYTHON SYNTAX CHECK"
echo "============================================================"

python -m py_compile \
    bioresearch_env/actions.py \
    bioresearch_env/rewards.py \
    bioresearch_env/relationship_language.py \
    bioresearch_env/biological_context_agent.py \
    bioresearch_env/virus_taxonomy_agent.py \
    bioresearch_env/env.py \
    bioresearch_env/baseline_agent.py \
    search_agent.py \
    evidence_agent.py \
    evaluate_env.py

echo "Syntax check passed."


# ============================================================
# OPTIONAL UPDATED LEGACY BENCHMARK
# ============================================================
# Default is OFF because you already have the historical 6/12 result.
# To run it:
#
# sbatch --export=ALL,RUN_LEGACY=1 run_benchmark.sh
# ============================================================

RUN_LEGACY="${RUN_LEGACY:-0}"

if [ "$RUN_LEGACY" -eq 1 ]; then

    echo
    echo "============================================================"
    echo "STARTING UPDATED LEGACY PIPELINE"
    echo "============================================================"

    python evaluate_pairs.py \
        test_pairs.csv \
        results/legacy_benchmark_v05.csv

    echo
    echo "Updated legacy benchmark complete."

else
    echo
    echo "Skipping updated legacy benchmark."
fi


# ============================================================
# BIORESEARCHENV v0.5 BENCHMARK
# ============================================================

echo
echo "============================================================"
echo "STARTING BIORESEARCHENV v0.5 BENCHMARK"
echo "============================================================"

python evaluate_env.py \
    test_pairs.csv \
    results/bioresearch_env_v05_results.csv


echo
echo "============================================================"
echo "RESULT FILES"
echo "============================================================"

ls -lh results/bioresearch_env_v05_results.csv

if [ -f results/legacy_benchmark_v05.csv ]; then
    ls -lh results/legacy_benchmark_v05.csv
fi


echo
echo "End time: $(date)"
echo "============================================================"
echo "SLURM JOB COMPLETE"
echo "============================================================"
