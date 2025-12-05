#!/bin/bash
#SBATCH --job-name=cooccur_opt
#SBATCH --account=def-acdoxey
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL


# Host-pathogen analysis - threshold-by-threshold processing
# Usage: ./process_flexible_thresholds.sh

set -e  # Exit on error

DUCKDB_CLI="/home/smascar/projects/def-acdoxey/smascar/duckdb"
OUTPUT_DB="host_pathogen_analysis.duckdb"

if [ ! -x "$DUCKDB_CLI" ]; then
  echo "ERROR: duckdb CLI not found or not executable at: $DUCKDB_CLI" >&2
  exit 1
fi

# ============================================
# CONFIGURE YOUR THRESHOLDS HERE (Bash arrays)
# ============================================
# Example subset:
#HOST_THRESHOLDS=(1e-06 1e-05 5e-05 1e-04 5e-04)
#PATHOGEN_THRESHOLDS=(1e-06 1e-05 5e-05 1e-04 5e-04)

# Or uncomment to use all 30:
HOST_THRESHOLDS=(1e-06 1e-05 5e-05 1e-04 5e-04 0.001 0.002 0.003 0.004 0.005 \
                  0.006 0.007 0.008 0.009 0.01 0.02 0.03 0.04 0.05 0.06 \
                  0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 1 5)
PATHOGEN_THRESHOLDS=(1e-06 1e-05 5e-05 1e-04 5e-04 0.001 0.002 0.003 0.004 0.005 \
                      0.006 0.007 0.008 0.009 0.01 0.02 0.03 0.04 0.05 0.06 \
                      0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 1 5)

echo "=========================================="
echo "Host-Pathogen Analysis (threshold-by-threshold)"
echo "=========================================="
echo "Database: $OUTPUT_DB"
echo ""
echo "Host thresholds:      ${HOST_THRESHOLDS[*]}"
echo "Pathogen thresholds:  ${PATHOGEN_THRESHOLDS[*]}"
echo ""

# Convert arrays to comma-separated strings for optional DuckDB tables
HOST_THRESHOLDS_STR=$(IFS=','; echo "${HOST_THRESHOLDS[*]}")
PATHOGEN_THRESHOLDS_STR=$(IFS=','; echo "${PATHOGEN_THRESHOLDS[*]}")

# ============================================
# STEP 1: Setup base tables and threshold tables
# ============================================
echo "[1/6] Setting up base tables and thresholds..."
$DUCKDB_CLI "$OUTPUT_DB" << EOF

SET threads = 16;
SET memory_limit = '180GB';

-- Host-pathogen pairs
CREATE OR REPLACE TABLE host_pathogen_pairs AS
SELECT host_tax_id, path_tax_id, pathogen_root_label, Host_Name, Pathogen_Name
FROM read_csv('manuscript_1_data/100_host_microbe_pairs_11_11.csv');

-- Optional: store thresholds inside DuckDB for reference
CREATE OR REPLACE TABLE host_thresholds AS
SELECT UNNEST([${HOST_THRESHOLDS_STR}]) AS threshold_value;

CREATE OR REPLACE TABLE pathogen_thresholds AS
SELECT UNNEST([${PATHOGEN_THRESHOLDS_STR}]) AS threshold_value;

-- Clean result tables if they exist
DROP TABLE IF EXISTS host_dataset_counts;
DROP TABLE IF EXISTS pathogen_dataset_counts;
DROP TABLE IF EXISTS shared_dataset_counts;

-- New, leaner result tables
CREATE TABLE host_dataset_counts (
    host_tax_id INTEGER,
    host_threshold DOUBLE,
    num_host_datasets INTEGER,
    host_datasets VARCHAR[]
);

CREATE TABLE pathogen_dataset_counts (
    path_tax_id INTEGER,
    pathogen_root_label VARCHAR,
    pathogen_threshold DOUBLE,
    num_pathogen_datasets INTEGER,
    pathogen_datasets VARCHAR[]
);

SELECT 'Loaded ' || COUNT(*) || ' host-pathogen pairs' AS status
FROM host_pathogen_pairs;

SELECT 'Host thresholds: ' || COUNT(*) AS status
FROM host_thresholds;

SELECT 'Pathogen thresholds: ' || COUNT(*) AS status
FROM pathogen_thresholds;

EOF

# ============================================
# STEP 2: Load extracted parquet files as views
# ============================================
echo "[2/6] Creating views over host and pathogen parquet files..."
$DUCKDB_CLI "$OUTPUT_DB" << 'EOF'

CREATE OR REPLACE VIEW host_abundance AS
SELECT 
    tax_id,
    acc,
    total_abundance
FROM read_parquet('host_data_full.parquet');

CREATE OR REPLACE VIEW pathogen_abundance AS
SELECT 
    tax_id,
    acc,
    total_abundance,
    pathogen_type
FROM read_parquet('pathogen_data_full.parquet');

SELECT 'Host records: ' || COUNT(*) AS info FROM host_abundance
UNION ALL
SELECT 'Pathogen records: ' || COUNT(*) AS info FROM pathogen_abundance;

EOF

# ============================================
# STEP 3: Process host dataset counts, one host_threshold at a time
# ============================================
echo "[3/6] Computing host dataset counts per host_threshold..."
for H_THR in "${HOST_THRESHOLDS[@]}"; do
  echo "  -> Host threshold = ${H_THR}"
  $DUCKDB_CLI "$OUTPUT_DB" << EOF
INSERT INTO host_dataset_counts
SELECT
    hp.host_tax_id,
    CAST(${H_THR} AS DOUBLE) AS host_threshold,
    COUNT(DISTINCT ha.acc) AS num_host_datasets,
    LIST(DISTINCT ha.acc) AS host_datasets
FROM (SELECT DISTINCT host_tax_id FROM host_pathogen_pairs) hp
LEFT JOIN host_abundance ha
    ON ha.tax_id = hp.host_tax_id
   AND ha.total_abundance >= ${H_THR}
GROUP BY hp.host_tax_id;
EOF
done

HOST_ROWS=$($DUCKDB_CLI "$OUTPUT_DB" -csv -noheader "SELECT COUNT(*) FROM host_dataset_counts;")
echo "  Completed host_dataset_counts with ${HOST_ROWS} rows."

# ============================================
# STEP 4: Process pathogen dataset counts, one pathogen_threshold at a time
# ============================================
echo "[4/6] Computing pathogen dataset counts per pathogen_threshold..."
for P_THR in "${PATHOGEN_THRESHOLDS[@]}"; do
  echo "  -> Pathogen threshold = ${P_THR}"
  $DUCKDB_CLI "$OUTPUT_DB" << EOF
INSERT INTO pathogen_dataset_counts
SELECT
    hp.path_tax_id,
    hp.pathogen_root_label,
    CAST(${P_THR} AS DOUBLE) AS pathogen_threshold,
    COUNT(DISTINCT pa.acc) AS num_pathogen_datasets,
    LIST(DISTINCT pa.acc) AS pathogen_datasets
FROM (SELECT DISTINCT path_tax_id, pathogen_root_label FROM host_pathogen_pairs) hp
LEFT JOIN pathogen_abundance pa
    ON pa.tax_id = hp.path_tax_id
   AND pa.pathogen_type = hp.pathogen_root_label
   AND pa.total_abundance >= ${P_THR}
GROUP BY hp.path_tax_id, hp.pathogen_root_label;
EOF
done

PATH_ROWS=$($DUCKDB_CLI "$OUTPUT_DB" -csv -noheader "SELECT COUNT(*) FROM pathogen_dataset_counts;")
echo "  Completed pathogen_dataset_counts with ${PATH_ROWS} rows."

# ============================================
# STEP 5: Build all threshold combinations using small aggregated tables
# ============================================
echo "[5/6] Calculating shared dataset counts from aggregated tables..."
$DUCKDB_CLI "$OUTPUT_DB" << 'EOF'

CREATE OR REPLACE TABLE shared_dataset_counts AS
SELECT 
    hp.host_tax_id,
    hp.path_tax_id,
    hdc.host_threshold,
    pdc.pathogen_threshold,
    hdc.num_host_datasets,
    pdc.num_pathogen_datasets,
    COALESCE(
        ARRAY_LENGTH(LIST_INTERSECT(hdc.host_datasets, pdc.pathogen_datasets)),
        0
    ) AS num_shared_datasets
FROM host_pathogen_pairs hp
JOIN host_dataset_counts hdc
    ON hdc.host_tax_id = hp.host_tax_id
JOIN pathogen_dataset_counts pdc
    ON pdc.path_tax_id = hp.path_tax_id
   AND pdc.pathogen_root_label = hp.pathogen_root_label;

SELECT 'Created ' || COUNT(*) || ' shared dataset records' AS status
FROM shared_dataset_counts;

EOF

RESULT_COUNT=$($DUCKDB_CLI "$OUTPUT_DB" -csv -noheader "SELECT COUNT(*) FROM shared_dataset_counts;")
echo "  Shared dataset rows: ${RESULT_COUNT}"

# For reference, expected combinations:
PAIR_COUNT=$($DUCKDB_CLI "$OUTPUT_DB" -csv -noheader "SELECT COUNT(*) FROM host_pathogen_pairs;")
TOTAL_COMBOS=$(( PAIR_COUNT * ${#HOST_THRESHOLDS[@]} * ${#PATHOGEN_THRESHOLDS[@]} ))

# ============================================
# STEP 6: Export results + summary
# ============================================
echo ""
echo "[6/6] Exporting results..."
$DUCKDB_CLI "$OUTPUT_DB" << 'EOF'

COPY shared_dataset_counts
TO 'host_pathogen_analysis_results.csv'
(HEADER, DELIMITER ',');

COPY shared_dataset_counts
TO 'host_pathogen_analysis_results.parquet'
(FORMAT PARQUET);

EOF

echo ""
echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo "Total host-pathogen pairs:     ${PAIR_COUNT}"
echo "Threshold combinations (theoretical): ${TOTAL_COMBOS}"
echo "Total result rows (actual):    ${RESULT_COUNT}"
echo ""
echo "Results saved to:"
echo "  - host_pathogen_analysis_results.csv"
echo "  - host_pathogen_analysis_results.parquet"
echo "  - ${OUTPUT_DB} (persistent database)"
echo ""
echo "Sample results:"
$DUCKDB_CLI "$OUTPUT_DB" << 'EOF'
SELECT * FROM shared_dataset_counts
ORDER BY host_tax_id, path_tax_id, host_threshold, pathogen_threshold
LIMIT 10;
EOF

echo ""
echo "Summary statistics:"
$DUCKDB_CLI "$OUTPUT_DB" << 'EOF'
SELECT 
    MIN(num_host_datasets)          AS min_host,
    MAX(num_host_datasets)          AS max_host,
    ROUND(AVG(num_host_datasets),2) AS avg_host,
    MIN(num_pathogen_datasets)          AS min_pathogen,
    MAX(num_pathogen_datasets)          AS max_pathogen,
    ROUND(AVG(num_pathogen_datasets),2) AS avg_pathogen,
    MIN(num_shared_datasets)          AS min_shared,
    MAX(num_shared_datasets)          AS max_shared,
    ROUND(AVG(num_shared_datasets),2) AS avg_shared
FROM shared_dataset_counts;
EOF

echo ""
echo "Done!"

