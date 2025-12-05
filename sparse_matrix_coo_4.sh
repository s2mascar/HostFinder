#!/bin/bash
#SBATCH --job-name=cooccur_opt
#SBATCH --account=def-acdoxey
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL


DUCKDB=/home/smascar/projects/def-acdoxey/smascar/duckdb

# Step 1: Create mappings only (very fast)
$DUCKDB << 'EOF'
SET threads TO 32;
SET memory_limit = '180GB';
SET preserve_insertion_order = false;

COPY (
  SELECT DISTINCT acc, 
         ROW_NUMBER() OVER (ORDER BY acc) - 1 AS row_idx
  FROM read_parquet('/home/smascar/scratch/STAT_2025_with_abundance_cleaned.parquet')
) TO 'row_mapping.parquet' (FORMAT PARQUET);

COPY (
  SELECT DISTINCT tax_id,
         ROW_NUMBER() OVER (ORDER BY tax_id) - 1 AS col_idx  
  FROM read_parquet('/home/smascar/scratch/STAT_2025_with_abundance_cleaned.parquet')
) TO 'col_mapping.parquet' (FORMAT PARQUET);
EOF

# Step 2: Use a single-pass approach without materializing tables
$DUCKDB << 'EOF'
SET threads TO 32;
SET memory_limit = '180GB';
SET preserve_insertion_order = false;

-- Direct join without creating intermediate tables
COPY (
  WITH row_map AS (SELECT * FROM 'row_mapping.parquet'),
       col_map AS (SELECT * FROM 'col_mapping.parquet')
  SELECT r.row_idx, c.col_idx, d.total_abundance AS value
  FROM read_parquet('/home/smascar/scratch/STAT_2025_with_abundance_cleaned.parquet') d
  INNER JOIN row_map r ON d.acc = r.acc
  INNER JOIN col_map c ON d.tax_id = c.tax_id
  WHERE d.total_abundance > 0
) TO 'sparse_matrix.parquet' (FORMAT PARQUET, COMPRESSION 'ZSTD');
EOF

echo "Done! Check outputs:"
echo "- row_mapping.parquet"
echo "- col_mapping.parquet"  
echo "- sparse_matrix.parquet"
