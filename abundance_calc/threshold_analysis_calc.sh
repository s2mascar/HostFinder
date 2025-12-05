#!/bin/bash
#SBATCH --job-name=##
#SBATCH --account=##
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --mail-user=##
#SBATCH --mail-type=ALL

DUCKDB="ADD_PATH"
NUM_CHUNKS=100
DB_FILE="threshold_analysis.duckdb"
rm -f "$DB_FILE"

# Initial setup
$DUCKDB "$DB_FILE" <<'EOF'
CREATE OR REPLACE TABLE thresholds(threshold DOUBLE);

INSERT INTO thresholds VALUES
    (1e-06), (1e-05), (5e-05), (1e-04), (5e-04),
    (0.001), (0.002), (0.003), (0.004), (0.005),
    (0.006), (0.007), (0.008), (0.009), (0.01),
    (0.02), (0.03), (0.04), (0.05), (0.06),
    (0.07), (0.08), (0.09), (0.1),
    (0.2), (0.3), (0.4), (0.5), (1.0), (5.0);

CREATE OR REPLACE TABLE host_pathogen_pairs AS
SELECT * FROM read_csv('manuscript_1_data/100_host_microbe_pairs_11_11.csv');

CREATE OR REPLACE TABLE host_list AS
SELECT DISTINCT host_tax_id FROM host_pathogen_pairs;

CREATE OR REPLACE TABLE pathogen_list AS
SELECT DISTINCT path_tax_id FROM host_pathogen_pairs;

CREATE OR REPLACE TABLE host_path_pairs_all AS
SELECT h.host_tax_id, p.path_tax_id
FROM host_list h CROSS JOIN pathogen_list p;

CREATE OR REPLACE TABLE host_data AS
SELECT * FROM read_parquet('host_data_11_20_filtered.parquet');

CREATE OR REPLACE TABLE pathogen_data AS
SELECT * FROM read_parquet('pathogen_data_11_20_filtered.parquet');

CREATE OR REPLACE TABLE host_threshold_counts AS
SELECT
    h.host_tax_id,
    t.threshold AS host_threshold,
    COUNT(DISTINCT CASE WHEN d.total_abundance >= t.threshold THEN d.acc ELSE NULL END) AS host_dataset_count
FROM host_list h
CROSS JOIN thresholds t
LEFT JOIN host_data d ON d.tax_id = h.host_tax_id
GROUP BY h.host_tax_id, t.threshold
ORDER BY h.host_tax_id, t.threshold;

CREATE OR REPLACE TABLE pathogen_threshold_counts AS
SELECT
    p.path_tax_id,
    t.threshold AS path_threshold,
    COUNT(DISTINCT CASE WHEN d.total_abundance >= t.threshold THEN d.acc ELSE NULL END) AS pathogen_dataset_count
FROM pathogen_list p
CROSS JOIN thresholds t
LEFT JOIN pathogen_data d ON d.tax_id = p.path_tax_id
GROUP BY p.path_tax_id, t.threshold
ORDER BY p.path_tax_id, t.threshold;

CREATE OR REPLACE TABLE pair_abundances AS
WITH host_ab AS (
    SELECT acc, tax_id AS host_tax_id, total_abundance AS host_abundance FROM host_data
),
path_ab AS (
    SELECT acc, tax_id AS path_tax_id, total_abundance AS path_abundance FROM pathogen_data
)
SELECT h.host_tax_id, p.path_tax_id, h.acc, h.host_abundance, p.path_abundance
FROM host_ab h
JOIN path_ab p ON p.acc = h.acc
JOIN host_list hl ON hl.host_tax_id = h.host_tax_id
JOIN pathogen_list pl ON pl.path_tax_id = p.path_tax_id;

CREATE OR REPLACE TABLE host_path_pairs_all_chunked AS
SELECT host_tax_id, path_tax_id, NTILE(100) OVER () AS chunk_id
FROM host_path_pairs_all;

CREATE OR REPLACE TABLE pair_threshold_counts (
    host_tax_id BIGINT,
    path_tax_id BIGINT,
    host_threshold DOUBLE,
    path_threshold DOUBLE,
    both_dataset_count BIGINT
);
EOF

# Loop through all chunks
for i in $(seq 1 $NUM_CHUNKS); do
    echo "Processing chunk $i of $NUM_CHUNKS..."
    $DUCKDB "$DB_FILE" <<EOF
INSERT INTO pair_threshold_counts
SELECT
    hp.host_tax_id,
    hp.path_tax_id,
    ht.threshold AS host_threshold,
    pt.threshold AS path_threshold,
    COUNT(DISTINCT CASE
        WHEN pa.host_abundance >= ht.threshold AND pa.path_abundance >= pt.threshold
        THEN pa.acc ELSE NULL
    END) AS both_dataset_count
FROM host_path_pairs_all_chunked hp
CROSS JOIN thresholds ht
CROSS JOIN thresholds pt
LEFT JOIN pair_abundances pa
   ON pa.host_tax_id = hp.host_tax_id AND pa.path_tax_id = hp.path_tax_id
WHERE hp.chunk_id = $i
GROUP BY hp.host_tax_id, hp.path_tax_id, ht.threshold, pt.threshold;
EOF
done

# Final summary
$DUCKDB "$DB_FILE" <<'EOF'
CREATE OR REPLACE TABLE host_pathogen_threshold_summary AS
SELECT
    hp.host_tax_id,
    hp.path_tax_id,
    htc.host_threshold,
    ptc.path_threshold,
    htc.host_dataset_count,
    ptc.pathogen_dataset_count,
    pt.both_dataset_count
FROM host_path_pairs_all hp
JOIN host_threshold_counts htc ON htc.host_tax_id = hp.host_tax_id
JOIN pathogen_threshold_counts ptc ON ptc.path_tax_id = hp.path_tax_id
JOIN pair_threshold_counts pt
    ON pt.host_tax_id = hp.host_tax_id
   AND pt.path_tax_id = hp.path_tax_id
   AND pt.host_threshold = htc.host_threshold
   AND pt.path_threshold = ptc.path_threshold;

COPY host_pathogen_threshold_summary TO 'host_pathogen_threshold_summary.parquet' (FORMAT PARQUET);
EOF

