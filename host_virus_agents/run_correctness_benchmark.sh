#!/bin/bash
#SBATCH --job-name=host-virus-correctness
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G

# Submit from the repository root with the account/GPU options appropriate for
# your allocation. No machine-specific directory or expected pair is embedded.
set -euo pipefail
test -f AGENTS.md
test -f PROJECT_SPEC.md
test -f test_pairs.csv
source agent_env/bin/activate
export MODEL_PATH="${MODEL_PATH:-./models/Qwen3-8B}"
test -d "$MODEL_PATH"
RUN_ID="correctness_${SLURM_JOB_ID:-local}_$(date -u +%Y%m%dT%H%M%SZ)"
RUN_DIR="results/$RUN_ID"
mkdir "$RUN_DIR"
exec > >(tee "$RUN_DIR/run.log") 2>&1
export HOST_ALIAS_CACHE="$RUN_DIR/host_aliases.json"
export VIRUS_TAXONOMY_CACHE="$RUN_DIR/virus_taxonomy.json"
export BIOLOGICAL_CONTEXT_CACHE="$RUN_DIR/biological_context.json"
python -B -m unittest discover -s tests -v
python -B evaluate_env.py test_pairs.csv "$RUN_DIR/env.csv"
python -B evaluate_pairs.py test_pairs.csv "$RUN_DIR/pairs.csv"
python -B validate_live_benchmark.py "$RUN_DIR/env.csv.diagnostics.jsonl" --output "$RUN_DIR/env_acceptance.json"
python -B validate_live_benchmark.py "$RUN_DIR/pairs.csv.diagnostics.jsonl" --output "$RUN_DIR/pairs_acceptance.json"
printf 'Return this complete run directory for evidence review: %s\n' "$RUN_DIR"
