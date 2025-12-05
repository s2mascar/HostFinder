#!/usr/bin/env bash
set -euo pipefail

# Inputs / outputs
INPUT_GLOB='INSERT_NAME/*.parquet'
OUT_DIR="$PWD/alpha_by_acc"

# Resources
THREADS=16
MEM_GB=80
PREFIX_LEN=3

# Use node-local temp if under Slurm; else create a private /tmp dir
if [[ -n "${SLURM_TMPDIR:-}" && -d "$SLURM_TMPDIR" ]]; then
  TMP_DIR="$SLURM_TMPDIR"
  CLEAN_TMP=0
else
  TMP_DIR="$(mktemp -d -p /tmp duckdb-XXXXXX)"
  CLEAN_TMP=1
fi
trap '[[ ${CLEAN_TMP} -eq 1 ]] && rm -rf "$TMP_DIR"' EXIT

mkdir -p "$OUT_DIR"

duckdb -c "
PRAGMA threads=${THREADS};
PRAGMA memory_limit='${MEM_GB}GB';
PRAGMA temp_directory='${TMP_DIR}';

WITH de_duped AS (
  SELECT acc, tax_id, SUM(CAST(total_abundance AS DOUBLE)) AS x
  FROM read_parquet('${INPUT_GLOB}', union_by_name=true)
  WHERE total_abundance > 0
  GROUP BY acc, tax_id
),
base AS (
  SELECT
    acc,
    SUM(x)         AS T,
    SUM(x*ln(x))   AS A,
    SUM(x*x)       AS B,
    MAX(x)         AS M,
    COUNT(*)       AS S
  FROM de_duped
  GROUP BY acc
),
alpha AS (
  SELECT
    acc,
    S                                            AS richness,
    CASE WHEN T > 0 THEN (A/T) - ln(T) END       AS shannon,
    CASE WHEN T > 0 AND B > 0 THEN 1 - (B/(T*T)) END AS simpson_1_minus_D,
    CASE WHEN B > 0 THEN (T*T)/B END             AS inverse_simpson,
    CASE WHEN T > 0 AND S > 1 THEN ((A/T) - ln(T))/ln(S) END AS pielou_evenness,
    CASE WHEN T > 0 THEN M/T END                 AS berger_parker_dominance,
    substr(acc,1,${PREFIX_LEN})                  AS acc_prefix
  FROM base
)
COPY (SELECT * FROM alpha)
TO '${OUT_DIR}'
(WITH (FORMAT PARQUET, COMPRESSION ZSTD, PARTITION_BY (acc_prefix)));
"

echo "Done. Written to: ${OUT_DIR}"

