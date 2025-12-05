#!/bin/bash

# Add duckdb to PATH
export PATH="ADD_PATH:$PATH"

# Check if duckdb is available
if ! command -v duckdb &> /dev/null; then
    echo "ERROR: duckdb command not found!"
    echo "Please install DuckDB or add it to your PATH"
    exit 1
fi

# Input: your generated parquet with host-microbe pairs and their SRA lists
INPUT_PARQUET="pathogen_interac_with_sra.parquet"

# Path to SRA metadata parquet files
SRA_META_PATH="ADD_PATH/*"

# Intermediate file with unique accessions
UNIQUE_ACC_PARQUET="unique_accessions_pat.parquet"

# Output file with all combined SRA metadata
OUTPUT_PARQUET="pathogen_interac_all_sra_metadata.parquet"

duckdb -cmd "
COPY (
    SELECT DISTINCT UNNEST(SRA_datasets) as acc
    FROM read_parquet('$INPUT_PARQUET')
    WHERE LEN(SRA_datasets) > 0
) TO '$UNIQUE_ACC_PARQUET' (FORMAT PARQUET);

SELECT 'Total unique accessions: ' || COUNT(*) 
FROM read_parquet('$UNIQUE_ACC_PARQUET');
"


duckdb -cmd "
SET memory_limit='8GB';

COPY (
    SELECT s.*
    FROM read_parquet('$SRA_META_PATH') AS s
    WHERE s.acc IN (
        SELECT acc FROM read_parquet('$UNIQUE_ACC_PARQUET')
    )
) TO '$OUTPUT_PARQUET' (FORMAT PARQUET);

SELECT 'Total SRA metadata rows extracted: ' || COUNT(*)
FROM read_parquet('$OUTPUT_PARQUET');
"

#echo ""
#echo "Step 3: (Optional) Adding host-microbe pair information..."
#echo "Creating final output with interaction details..."

#duckdb -cmd "
#SET memory_limit='8GB';

#-- Create expanded interactions table
#CREATE OR REPLACE TEMP TABLE interactions_expanded AS
#SELECT 
#    host_tax_id,
#    microbe_tax_id,
#    interaction_type,
#    UNNEST(SRA_datasets) as acc
#FROM read_parquet('$INPUT_PARQUET')
#WHERE LEN(SRA_datasets) > 0;

#-- Join everything together
#COPY (
#    SELECT 
#        i.host_tax_id,
#        i.microbe_tax_id,
#        i.interaction_type,
#        s.*
#    FROM read_parquet('$OUTPUT_PARQUET') AS s
#    JOIN interactions_expanded AS i ON s.acc = i.acc
#) TO 'pathogen_interac_all_sra_metadata_with_pairs.parquet' (FORMAT PARQUET);
#
#-- Summary statistics
#SELECT 
#    'Total unique host-microbe pairs: ' || COUNT(DISTINCT host_tax_id || '-' || microbe_tax_id) as summary
#FROM read_parquet('pathogen_interac_all_sra_metadata_with_pairs.parquet')
#UNION ALL
#SELECT 
#    'Total SRA datasets: ' || COUNT(DISTINCT acc)
#FROM read_parquet('pathogen_interac_all_sra_metadata_with_pairs.parquet')
#UNION ALL
#SELECT 
#    'Total rows: ' || COUNT(*)
#FROM read_parquet('pathogen_interac_all_sra_metadata_with_pairs.parquet');
#"
echo "- Unique accessions: $UNIQUE_ACC_PARQUET"
echo "- SRA metadata only: $OUTPUT_PARQUET"
#echo "- SRA metadata with pair info: pathogen_interac_all_sra_metadata_with_pairs.parquet"
