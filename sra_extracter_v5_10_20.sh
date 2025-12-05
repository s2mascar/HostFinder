#!/bin/bash

# Optimized DuckDB Host-Microbe Abundance Processing
# Usage: ./process_abundances.sh [input_parquet] [stat_path] [output_file]

set -euo pipefail

# Add duckdb to PATH
export PATH="/home/smascar/projects/def-acdoxey/smascar:$PATH"

# Check if duckdb is available
if ! command -v duckdb &> /dev/null; then
    echo "ERROR: duckdb command not found!"
    echo "Please install DuckDB or add it to your PATH"
    exit 1
fi


# Configuration
INPUT_PARQUET="${1:-pathogen_interac_with_sra.parquet}"
STAT_PATH="${2:-/home/smascar/scratch/STAT_sorted_files/*}"
OUTPUT_FILE="${3:-pathogen_all_STAT_data.parquet}"
THREADS="${DUCKDB_THREADS:-4}"
MEMORY="${DUCKDB_MEMORY:-45GB}"

echo "Starting DuckDB processing..."
echo "Input: $INPUT_PARQUET"
echo "STAT files: $STAT_PATH"
echo "Output: $OUTPUT_FILE"
echo "Threads: $THREADS, Memory: $MEMORY"

duckdb -c "
-- Performance tuning
PRAGMA threads=$THREADS;
PRAGMA memory_limit='$MEMORY';
PRAGMA enable_object_cache;
PRAGMA enable_profiling;
PRAGMA preserve_insertion_order=false;

-- Export directly to Parquet for better performance
COPY (
  WITH
  -- 1) Read pairs and explode SRA_datasets
  pairs AS (
    SELECT host_tax_id, microbe_tax_id, SRA_datasets
    FROM read_parquet('$INPUT_PARQUET')
    -- WHERE interaction_type = 'pathogenic'  -- uncomment if needed
  ),
  
  expanded AS (
    SELECT
      p.host_tax_id,
      p.microbe_tax_id,
      UNNEST(p.SRA_datasets) AS SRA_accession
    FROM pairs p
  ),
  
  -- 2) Create filtered lookups (semi-joins for better performance)
  needed_tax_ids AS (
    SELECT DISTINCT host_tax_id AS tax_id FROM pairs
    UNION
    SELECT DISTINCT microbe_tax_id FROM pairs
  ),
  
  needed_accessions AS (
    SELECT DISTINCT SRA_accession FROM expanded
  ),
  
  -- 3) Filtered STAT data with pushed-down predicates
  -- DuckDB will push these filters down into the parquet scan
  stat_filtered AS (
    SELECT s.acc, s.tax_id, s.total_abundance
    FROM read_parquet('$STAT_PATH', 
                      filename=true,
                      hive_partitioning=false) AS s
    WHERE EXISTS (SELECT 1 FROM needed_accessions WHERE SRA_accession = s.acc)
      AND EXISTS (SELECT 1 FROM needed_tax_ids WHERE tax_id = s.tax_id)
  ),
  
  -- 4) Join with conditional logic pushed into join
  abundances AS (
    SELECT
      e.host_tax_id,
      e.microbe_tax_id,
      e.SRA_accession,
      MAX(CASE WHEN s.tax_id = e.host_tax_id 
               THEN s.total_abundance END) AS host_abundance,
      MAX(CASE WHEN s.tax_id = e.microbe_tax_id 
               THEN s.total_abundance END) AS microbe_abundance
    FROM expanded e
    LEFT JOIN stat_filtered s
      ON s.acc = e.SRA_accession
     AND (s.tax_id = e.host_tax_id OR s.tax_id = e.microbe_tax_id)
    GROUP BY e.host_tax_id, e.microbe_tax_id, e.SRA_accession
  )
  
  SELECT 
    host_tax_id,
    microbe_tax_id,
    host_abundance,
    microbe_abundance,
    SRA_accession
  FROM abundances
  WHERE host_abundance IS NOT NULL OR microbe_abundance IS NOT NULL
  -- Optional: remove rows with no abundance data
  
) TO '$OUTPUT_FILE' (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000);

-- Show query profile
-- PRAGMA show_tables_expanded;
" 2>&1 | tee duckdb_processing.log

echo "Processing complete! Output saved to: $OUTPUT_FILE"
echo "Log saved to: duckdb_processing.log"

# Optional: Show basic statistics
duckdb -c "
SELECT 
  COUNT(*) as total_rows,
  COUNT(DISTINCT host_tax_id) as unique_hosts,
  COUNT(DISTINCT microbe_tax_id) as unique_microbes,
  COUNT(DISTINCT SRA_accession) as unique_accessions,
  AVG(host_abundance) as avg_host_abundance,
  AVG(microbe_abundance) as avg_microbe_abundance
FROM read_parquet('$OUTPUT_FILE');
"
