#!/bin/bash
#SBATCH --job-name=vh8508_best
#SBATCH --account=def-acdoxey
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL

set -euo pipefail
set -x


# ============================================================
# CONFIGURATION
# ============================================================

DUCKDB_BIN="/home/smascar/scratch/STAT_2025/duckdb"


# ------------------------------------------------------------
# 8,508 DETECTABLE VIRUS-HOST DB POSITIVES
# ------------------------------------------------------------

POSITIVE_PARQUET="/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions/virushostdb_pairs_detectable_1e12.parquet"


# ------------------------------------------------------------
# RAW STAT DATA
# ------------------------------------------------------------

HOST_STAT_GLOB="/home/smascar/scratch/VIRUSES/STAT_with_Eukaryota_important_columns_split/**"

VIRUS_STAT_GLOB="/home/smascar/scratch/VIRUSES/STAT_with_Viruses_important_columns_split/**"


# ------------------------------------------------------------
# OUTPUT
# ------------------------------------------------------------

OUTPUT_DIR="/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions"

OUTPUT_PARQUET="${OUTPUT_DIR}/virushostdb_8508_best_threshold_scored.parquet"


# ------------------------------------------------------------
# SELECTED BEST THRESHOLDS
# ------------------------------------------------------------

HOST_THRESHOLD="0.001"

PATHOGEN_THRESHOLD="1e-8"


# ------------------------------------------------------------
# DUCKDB TEMP DIRECTORY
# ------------------------------------------------------------

TEMP_DIR="${SLURM_TMPDIR:-$PWD}/duckdb_temp_vh8508_best_threshold"


mkdir -p "${OUTPUT_DIR}"

mkdir -p "${TEMP_DIR}"


# ============================================================
# REMOVE OLD OUTPUT
# ============================================================

rm -f "${OUTPUT_PARQUET}"


# ============================================================
# RUN DUCKDB
# ============================================================

"${DUCKDB_BIN}" <<EOF

PRAGMA temp_directory='${TEMP_DIR}';

PRAGMA max_temp_directory_size='2500GB';

PRAGMA enable_progress_bar=true;


SET threads=8;

SET memory_limit='150GB';

SET preserve_insertion_order=false;



-- ============================================================
-- 1. LOAD THE 8,508 KNOWN POSITIVE PAIRS
--
-- Rename virus_taxid -> microbe_taxid to match the schema used
-- in the previous scored benchmark files.
-- ============================================================

CREATE OR REPLACE TEMP TABLE positive_pairs AS

SELECT DISTINCT
    CAST(
        host_taxid AS BIGINT
    ) AS host_taxid,

    CAST(
        virus_taxid AS BIGINT
    ) AS microbe_taxid,

    host_name,

    virus_name

FROM read_parquet(
    '${POSITIVE_PARQUET}'
);



-- ============================================================
-- 2. POSITIVE PAIR QC
-- ============================================================

SELECT
    COUNT(*) AS num_rows,

    COUNT(
        DISTINCT (
            host_taxid,
            microbe_taxid
        )
    ) AS num_pairs,

    COUNT(
        DISTINCT host_taxid
    ) AS num_hosts,

    COUNT(
        DISTINCT microbe_taxid
    ) AS num_viruses

FROM positive_pairs;



-- ============================================================
-- 3. GET THE HOST TAXA USED BY THE 8,508 POSITIVES
-- ============================================================

CREATE OR REPLACE TEMP TABLE benchmark_hosts AS

SELECT DISTINCT
    host_taxid

FROM positive_pairs;



-- ============================================================
-- 4. GET THE VIRUS TAXA USED BY THE 8,508 POSITIVES
-- ============================================================

CREATE OR REPLACE TEMP TABLE benchmark_viruses AS

SELECT DISTINCT
    microbe_taxid

FROM positive_pairs;



-- ============================================================
-- 5. RECALCULATE THE FIXED ACCESSION AND BIOPROJECT UNIVERSES
--
-- These are based on the union of the exact current
-- Eukaryota and Virus STAT Parquet sources.
--
-- Historical current values from the fresh rebuild:
--
--     N_ACC        = 25,796,072
--     N_BIOPROJECT =    592,002
--
-- The query recalculates them instead of blindly hardcoding.
-- ============================================================

CREATE OR REPLACE TEMP TABLE universe_constants AS

WITH combined_stat AS (

    SELECT
        acc,
        bioproject

    FROM read_parquet(
        '${HOST_STAT_GLOB}'
    )


    UNION ALL


    SELECT
        acc,
        bioproject

    FROM read_parquet(
        '${VIRUS_STAT_GLOB}'
    )

)

SELECT
    COUNT(
        DISTINCT acc
    ) AS N_ACC,

    COUNT(
        DISTINCT bioproject
    ) FILTER (

        WHERE bioproject IS NOT NULL

          AND TRIM(
              bioproject
          ) <> ''

    ) AS N_BIOPROJECT

FROM combined_stat;



-- ============================================================
-- 6. DISPLAY UNIVERSE CONSTANTS
-- ============================================================

SELECT *

FROM universe_constants;



-- ============================================================
-- 7. EXTRACT HOST SIGNAL AT HOST THRESHOLD = 0.001
--
-- Restrict immediately to hosts involved in the 8,508 pairs.
--
-- One host can have multiple records for an accession.
-- We keep distinct taxon/accession/BioProject combinations.
-- ============================================================

CREATE OR REPLACE TEMP TABLE host_data AS

SELECT DISTINCT
    CAST(
        p.tax_id AS BIGINT
    ) AS host_taxid,

    p.acc,

    p.bioproject

FROM read_parquet(
    '${HOST_STAT_GLOB}'
) AS p

INNER JOIN benchmark_hosts AS h

    ON CAST(
        p.tax_id AS BIGINT
    ) = h.host_taxid

WHERE p.total_abundance >= ${HOST_THRESHOLD};



-- ============================================================
-- 8. EXTRACT VIRUS SIGNAL AT VIRUS THRESHOLD = 1e-8
-- ============================================================

CREATE OR REPLACE TEMP TABLE virus_data AS

SELECT DISTINCT
    CAST(
        p.tax_id AS BIGINT
    ) AS microbe_taxid,

    p.acc,

    p.bioproject

FROM read_parquet(
    '${VIRUS_STAT_GLOB}'
) AS p

INNER JOIN benchmark_viruses AS v

    ON CAST(
        p.tax_id AS BIGINT
    ) = v.microbe_taxid

WHERE p.total_abundance >= ${PATHOGEN_THRESHOLD};



-- ============================================================
-- 9. HOST COUNTS
-- ============================================================

CREATE OR REPLACE TEMP TABLE host_counts AS

SELECT
    host_taxid,

    COUNT(
        DISTINCT acc
    ) AS num_acc_in_host,

    COUNT(
        DISTINCT bioproject
    ) FILTER (

        WHERE bioproject IS NOT NULL

          AND TRIM(
              bioproject
          ) <> ''

    ) AS num_bioprojects_in_host

FROM host_data

GROUP BY
    host_taxid;



-- ============================================================
-- 10. VIRUS COUNTS
-- ============================================================

CREATE OR REPLACE TEMP TABLE virus_counts AS

SELECT
    microbe_taxid,

    COUNT(
        DISTINCT acc
    ) AS num_acc_in_pathogen,

    COUNT(
        DISTINCT bioproject
    ) FILTER (

        WHERE bioproject IS NOT NULL

          AND TRIM(
              bioproject
          ) <> ''

    ) AS num_bioprojects_in_pathogen

FROM virus_data

GROUP BY
    microbe_taxid;



-- ============================================================
-- 11. IDENTIFY EXACT SHARED ACCESSIONS
--
-- IMPORTANT:
--
-- We only calculate co-occurrence for the 8,508 known
-- Virus-Host DB pairs.
--
-- We do NOT create every possible host x virus combination.
-- ============================================================

CREATE OR REPLACE TEMP TABLE shared_accessions AS

SELECT DISTINCT
    p.host_taxid,

    p.microbe_taxid,

    h.acc

FROM positive_pairs AS p

INNER JOIN host_data AS h

    ON p.host_taxid = h.host_taxid

INNER JOIN virus_data AS v

    ON p.microbe_taxid = v.microbe_taxid

   AND h.acc = v.acc;



-- ============================================================
-- 12. COUNT SHARED ACCESSIONS
-- ============================================================

CREATE OR REPLACE TEMP TABLE shared_accession_counts AS

SELECT
    host_taxid,

    microbe_taxid,

    COUNT(
        DISTINCT acc
    ) AS num_acc_shared

FROM shared_accessions

GROUP BY
    host_taxid,
    microbe_taxid;



-- ============================================================
-- 13. COUNT SHARED BIOPROJECTS
--
-- Definition matches the fresh 985 pipeline:
--
-- BioProjects represented among exact host-virus
-- same-accession co-occurrences.
--
-- We return to host_data for the shared accession and count
-- the distinct BioProjects associated with those accessions.
-- ============================================================

CREATE OR REPLACE TEMP TABLE shared_bioproject_counts AS

SELECT
    s.host_taxid,

    s.microbe_taxid,

    COUNT(
        DISTINCT h.bioproject
    ) FILTER (

        WHERE h.bioproject IS NOT NULL

          AND TRIM(
              h.bioproject
          ) <> ''

    ) AS num_bioprojects_shared

FROM shared_accessions AS s

INNER JOIN host_data AS h

    ON s.host_taxid = h.host_taxid

   AND s.acc = h.acc

GROUP BY
    s.host_taxid,
    s.microbe_taxid;



-- ============================================================
-- 14. BUILD FULL 8,508-PAIR RETENTION-AWARE TABLE
--
-- Every original positive pair is retained.
--
-- If a taxon or pair does not survive the selected threshold,
-- missing counts become zero.
-- ============================================================

CREATE OR REPLACE TEMP TABLE scored_counts AS

SELECT
    p.host_taxid,

    p.microbe_taxid,

    p.host_name,

    p.virus_name,

    CAST(
        ${HOST_THRESHOLD} AS DOUBLE
    ) AS host_threshold,

    CAST(
        ${PATHOGEN_THRESHOLD} AS DOUBLE
    ) AS pathogen_threshold,

    COALESCE(
        h.num_acc_in_host,
        0
    ) AS num_acc_in_host,

    COALESCE(
        v.num_acc_in_pathogen,
        0
    ) AS num_acc_in_pathogen,

    COALESCE(
        s.num_acc_shared,
        0
    ) AS num_acc_shared,

    COALESCE(
        h.num_bioprojects_in_host,
        0
    ) AS num_bioprojects_in_host,

    COALESCE(
        v.num_bioprojects_in_pathogen,
        0
    ) AS num_bioprojects_in_pathogen,

    COALESCE(
        b.num_bioprojects_shared,
        0
    ) AS num_bioprojects_shared,

    CAST(
        1 AS BIGINT
    ) AS Correct_interaction

FROM positive_pairs AS p

LEFT JOIN host_counts AS h

    ON p.host_taxid = h.host_taxid

LEFT JOIN virus_counts AS v

    ON p.microbe_taxid = v.microbe_taxid

LEFT JOIN shared_accession_counts AS s

    ON p.host_taxid = s.host_taxid

   AND p.microbe_taxid = s.microbe_taxid

LEFT JOIN shared_bioproject_counts AS b

    ON p.host_taxid = b.host_taxid

   AND p.microbe_taxid = b.microbe_taxid;



-- ============================================================
-- 15. CALCULATE FREQUENCIES AND LOG-ODDS SCORES
--
-- Accession:
--
-- fh  = host accession frequency
-- fp  = virus accession frequency
-- fhp = shared accession frequency
--
-- log_odds_acc =
--
-- log2(
--
--   (fhp + 1/N)
--   /
--   (fh * fp + 1/N)
--
-- )
--
-- Same equation for BioProjects.
-- ============================================================

COPY (

    WITH constants AS (

        SELECT
            CAST(
                N_ACC AS DOUBLE
            ) AS N_ACC,

            CAST(
                N_BIOPROJECT AS DOUBLE
            ) AS N_BIOPROJECT

        FROM universe_constants

    ),


    frequencies AS (

        SELECT
            c.*,

            u.N_ACC,

            u.N_BIOPROJECT,

            CAST(
                c.num_acc_in_host AS DOUBLE
            ) / u.N_ACC AS fh_acc,

            CAST(
                c.num_acc_in_pathogen AS DOUBLE
            ) / u.N_ACC AS fp_acc,

            CAST(
                c.num_acc_shared AS DOUBLE
            ) / u.N_ACC AS fhp_acc,

            CAST(
                c.num_bioprojects_in_host AS DOUBLE
            ) / u.N_BIOPROJECT AS fh_bioproject,

            CAST(
                c.num_bioprojects_in_pathogen AS DOUBLE
            ) / u.N_BIOPROJECT AS fp_bioproject,

            CAST(
                c.num_bioprojects_shared AS DOUBLE
            ) / u.N_BIOPROJECT AS fhp_bioproject

        FROM scored_counts AS c

        CROSS JOIN constants AS u

    )


    SELECT
        host_taxid,

        microbe_taxid,

        host_name,

        virus_name,

        host_threshold,

        pathogen_threshold,

        num_acc_in_host,

        num_acc_in_pathogen,

        num_acc_shared,

        num_bioprojects_in_host,

        num_bioprojects_in_pathogen,

        num_bioprojects_shared,

        Correct_interaction,

        fh_acc,

        fp_acc,

        fhp_acc,

        LOG2(

            (
                fhp_acc
                +
                (
                    1.0 / N_ACC
                )
            )

            /

            (
                (
                    fh_acc * fp_acc
                )

                +
                (
                    1.0 / N_ACC
                )
            )

        ) AS log_odds_acc,

        fh_bioproject,

        fp_bioproject,

        fhp_bioproject,

        LOG2(

            (
                fhp_bioproject
                +
                (
                    1.0 / N_BIOPROJECT
                )
            )

            /

            (
                (
                    fh_bioproject
                    *
                    fp_bioproject
                )

                +
                (
                    1.0 / N_BIOPROJECT
                )
            )

        ) AS log_odds_bioproject

    FROM frequencies

)

TO '${OUTPUT_PARQUET}'

(
    FORMAT PARQUET,

    COMPRESSION ZSTD,

    ROW_GROUP_SIZE 100000
);



-- ============================================================
-- 16. FINAL OUTPUT QC
-- ============================================================

SELECT
    COUNT(*) AS num_rows,

    COUNT(
        DISTINCT (
            host_taxid,
            microbe_taxid
        )
    ) AS num_pairs,

    COUNT(
        DISTINCT host_taxid
    ) AS num_hosts,

    COUNT(
        DISTINCT microbe_taxid
    ) AS num_viruses,

    COUNT(*) FILTER (
        WHERE Correct_interaction = 1
    ) AS num_positives

FROM read_parquet(
    '${OUTPUT_PARQUET}'
);



-- ============================================================
-- 17. RETENTION QC
--
-- How many of the 8,508 positives still have both taxa
-- represented at the selected thresholds?
-- ============================================================

SELECT
    COUNT(*) AS total_positive_pairs,

    COUNT(*) FILTER (

        WHERE num_acc_in_host > 0

          AND num_acc_in_pathogen > 0

    ) AS pairs_with_both_taxa_detected,

    COUNT(*) FILTER (

        WHERE num_acc_in_host = 0

           OR num_acc_in_pathogen = 0

    ) AS pairs_lost_at_threshold,

    COUNT(*) FILTER (

        WHERE num_acc_shared > 0

    ) AS pairs_with_shared_accessions

FROM read_parquet(
    '${OUTPUT_PARQUET}'
);



-- ============================================================
-- 18. RETENTION PERCENTAGES
-- ============================================================

SELECT
    COUNT(*) AS total_pairs,

    COUNT(*) FILTER (

        WHERE num_acc_in_host > 0

          AND num_acc_in_pathogen > 0

    ) AS retained_pairs,

    ROUND(

        100.0

        *

        COUNT(*) FILTER (

            WHERE num_acc_in_host > 0

              AND num_acc_in_pathogen > 0

        )

        /

        COUNT(*),

        3

    ) AS retained_percent,

    COUNT(*) FILTER (

        WHERE num_acc_shared > 0

    ) AS pairs_with_cooccurrence,

    ROUND(

        100.0

        *

        COUNT(*) FILTER (

            WHERE num_acc_shared > 0

        )

        /

        COUNT(*),

        3

    ) AS cooccurrence_percent

FROM read_parquet(
    '${OUTPUT_PARQUET}'
);



-- ============================================================
-- 19. SCORE QC
-- ============================================================

SELECT
    COUNT(*) FILTER (
        WHERE log_odds_acc IS NULL
    ) AS null_log_odds_acc,

    COUNT(*) FILTER (
        WHERE log_odds_bioproject IS NULL
    ) AS null_log_odds_bioproject,

    COUNT(*) FILTER (

        WHERE num_acc_shared >
              num_acc_in_host

           OR num_acc_shared >
              num_acc_in_pathogen

    ) AS impossible_accession_counts,

    COUNT(*) FILTER (

        WHERE num_bioprojects_shared >
              num_bioprojects_in_host

           OR num_bioprojects_shared >
              num_bioprojects_in_pathogen

    ) AS impossible_bioproject_counts

FROM read_parquet(
    '${OUTPUT_PARQUET}'
);



-- ============================================================
-- 20. SCORE SUMMARY
-- ============================================================

SELECT
    MIN(
        log_odds_acc
    ) AS min_log_odds_acc,

    MEDIAN(
        log_odds_acc
    ) AS median_log_odds_acc,

    AVG(
        log_odds_acc
    ) AS mean_log_odds_acc,

    MAX(
        log_odds_acc
    ) AS max_log_odds_acc,

    MIN(
        log_odds_bioproject
    ) AS min_log_odds_bioproject,

    MEDIAN(
        log_odds_bioproject
    ) AS median_log_odds_bioproject,

    AVG(
        log_odds_bioproject
    ) AS mean_log_odds_bioproject,

    MAX(
        log_odds_bioproject
    ) AS max_log_odds_bioproject

FROM read_parquet(
    '${OUTPUT_PARQUET}'
);



-- ============================================================
-- 21. DISPLAY TOP KNOWN INTERACTIONS
-- ============================================================

SELECT
    host_taxid,

    host_name,

    microbe_taxid,

    virus_name,

    num_acc_shared,

    num_bioprojects_shared,

    log_odds_acc,

    log_odds_bioproject

FROM read_parquet(
    '${OUTPUT_PARQUET}'
)

ORDER BY
    log_odds_acc DESC

LIMIT 30;



-- ============================================================
-- DONE
-- ============================================================

SELECT
    '8,508 VIRUS-HOST DB POSITIVES SCORED AT BEST THRESHOLD'
        AS status;

EOF


# ============================================================
# FINAL OUTPUT
# ============================================================

echo
echo "============================================================"
echo "DONE"
echo "============================================================"

echo
echo "Final scored Parquet:"
echo "${OUTPUT_PARQUET}"

echo
echo "File size:"
ls -lh "${OUTPUT_PARQUET}"

echo
echo "============================================================"
