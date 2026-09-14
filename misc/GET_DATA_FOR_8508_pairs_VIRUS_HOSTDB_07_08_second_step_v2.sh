#!/bin/bash
#SBATCH --job-name=vh8508_cross
#SBATCH --account=def-acdoxey
#SBATCH --time=12:00:00
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

OUTPUT_PARQUET="${OUTPUT_DIR}/virushostdb_8508_all_host_virus_cross_scored.parquet"


# ------------------------------------------------------------
# SELECTED BEST THRESHOLDS
# ------------------------------------------------------------

HOST_THRESHOLD="0.001"

PATHOGEN_THRESHOLD="1e-8"


# ------------------------------------------------------------
# DUCKDB TEMP DIRECTORY
# ------------------------------------------------------------

TEMP_DIR="${SLURM_TMPDIR:-$PWD}/duckdb_temp_vh8508_cross"


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
-- These pairs define the TRUE POSITIVE interaction keys.
--
-- virus_taxid is renamed to microbe_taxid.
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
-- 2. CREATE UNIQUE POSITIVE PAIR KEYS
--
-- This table is used ONLY to label the known interactions.
--
-- Correct_interaction = 1
-- ============================================================

CREATE OR REPLACE TEMP TABLE positive_keys AS

SELECT DISTINCT
    host_taxid,
    microbe_taxid

FROM positive_pairs;



-- ============================================================
-- 3. POSITIVE PAIR QC
-- ============================================================

SELECT
    COUNT(*) AS num_positive_rows,

    COUNT(
        DISTINCT (
            host_taxid,
            microbe_taxid
        )
    ) AS num_positive_pairs,

    COUNT(
        DISTINCT host_taxid
    ) AS num_hosts,

    COUNT(
        DISTINCT microbe_taxid
    ) AS num_viruses

FROM positive_pairs;



-- ============================================================
-- 4. GET ALL UNIQUE BENCHMARK HOST TAXA
--
-- These are ALL host taxids represented among the
-- 8,508 Virus-Host DB positive pairs.
-- ============================================================

CREATE OR REPLACE TEMP TABLE benchmark_hosts AS

SELECT
    host_taxid,

    MIN(
        host_name
    ) AS host_name

FROM positive_pairs

GROUP BY
    host_taxid;



-- ============================================================
-- 5. GET ALL UNIQUE BENCHMARK VIRUS TAXA
--
-- These are ALL virus taxids represented among the
-- 8,508 Virus-Host DB positive pairs.
-- ============================================================

CREATE OR REPLACE TEMP TABLE benchmark_viruses AS

SELECT
    microbe_taxid,

    MIN(
        virus_name
    ) AS virus_name

FROM positive_pairs

GROUP BY
    microbe_taxid;



-- ============================================================
-- 6. CROSS EVERY HOST BY EVERY VIRUS
--
-- THIS IS THE MAIN CHANGE.
--
-- Example:
--
--     Host 1 x Virus A
--     Host 1 x Virus B
--     Host 1 x Virus C
--     Host 2 x Virus A
--     Host 2 x Virus B
--     Host 2 x Virus C
--
-- The original 8,508 pairs are labelled 1.
--
-- Every other host-virus combination is labelled 0.
-- ============================================================

CREATE OR REPLACE TEMP TABLE candidate_pairs AS

SELECT
    h.host_taxid,

    v.microbe_taxid,

    h.host_name,

    v.virus_name,

    CAST(

        CASE

            WHEN p.host_taxid IS NOT NULL
            THEN 1

            ELSE 0

        END

        AS BIGINT

    ) AS Correct_interaction

FROM benchmark_hosts AS h

CROSS JOIN benchmark_viruses AS v

LEFT JOIN positive_keys AS p

    ON h.host_taxid = p.host_taxid

   AND v.microbe_taxid = p.microbe_taxid;



-- ============================================================
-- 7. CANDIDATE CROSS QC
--
-- Expected:
--
-- total_pairs =
--
--     number of hosts
--     x
--     number of viruses
--
-- num_positives should equal exactly 8,508.
--
-- Everything else is unlabelled.
-- ============================================================

SELECT
    COUNT(*) AS total_pairs,

    COUNT(
        DISTINCT host_taxid
    ) AS num_hosts,

    COUNT(
        DISTINCT microbe_taxid
    ) AS num_viruses,

    COUNT(*) FILTER (

        WHERE Correct_interaction = 1

    ) AS num_positives,

    COUNT(*) FILTER (

        WHERE Correct_interaction = 0

    ) AS num_unlabelled

FROM candidate_pairs;



-- ============================================================
-- 8. RECALCULATE FIXED ACCESSION AND BIOPROJECT UNIVERSES
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
-- 9. DISPLAY UNIVERSE CONSTANTS
-- ============================================================

SELECT *

FROM universe_constants;



-- ============================================================
-- 10. EXTRACT HOST SIGNAL AT HOST THRESHOLD = 0.001
--
-- Only hosts represented in the 8,508 benchmark pairs
-- are needed.
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
-- 11. EXTRACT VIRUS SIGNAL AT VIRUS THRESHOLD = 1e-8
--
-- Only viruses represented in the 8,508 benchmark pairs
-- are needed.
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
-- 12. HOST COUNTS
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
-- 13. VIRUS COUNTS
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
-- 14. IDENTIFY ALL HOST-VIRUS SHARED ACCESSIONS
--
-- IMPORTANT CHANGE:
--
-- The old code started from positive_pairs.
--
-- Therefore, only known positive pairs could ever
-- have shared accession counts.
--
-- Here we join ALL benchmark hosts and ALL benchmark viruses
-- on accession.
--
-- Any host-virus combination detected in the same accession
-- is recovered.
-- ============================================================

CREATE OR REPLACE TEMP TABLE shared_accessions AS

SELECT DISTINCT
    h.host_taxid,

    v.microbe_taxid,

    h.acc

FROM host_data AS h

INNER JOIN virus_data AS v

    ON h.acc = v.acc;



-- ============================================================
-- 15. COUNT SHARED ACCESSIONS
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
-- 16. COUNT SHARED BIOPROJECTS
--
-- BioProjects represented among exact host-virus
-- same-accession co-occurrences.
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
-- 17. BUILD FULL HOST x VIRUS RETENTION-AWARE TABLE
--
-- EVERY host-virus combination is retained.
--
-- Known Virus-Host DB pair:
--
--     Correct_interaction = 1
--
-- All other combinations:
--
--     Correct_interaction = 0
--
-- Missing threshold counts become zero.
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

    p.Correct_interaction

FROM candidate_pairs AS p

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
-- 18. CALCULATE FREQUENCIES AND LOG-ODDS SCORES
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
-- 19. FINAL OUTPUT QC
--
-- Critical check:
--
-- num_positives MUST = 8,508
--
-- num_unlabelled = every other crossed combination
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

    ) AS num_positives,

    COUNT(*) FILTER (

        WHERE Correct_interaction = 0

    ) AS num_unlabelled

FROM read_parquet(
    '${OUTPUT_PARQUET}'
);



-- ============================================================
-- 20. VERIFY ALL 8,508 POSITIVE PAIRS WERE RETAINED
-- ============================================================

SELECT
    COUNT(*) AS expected_positive_pairs,

    COUNT(*) FILTER (

        WHERE o.Correct_interaction = 1

    ) AS positive_pairs_in_output

FROM positive_keys AS p

LEFT JOIN read_parquet(
    '${OUTPUT_PARQUET}'
) AS o

    ON p.host_taxid = o.host_taxid

   AND p.microbe_taxid = o.microbe_taxid;



-- ============================================================
-- 21. LABEL DISTRIBUTION
-- ============================================================

SELECT
    Correct_interaction,

    COUNT(*) AS num_pairs,

    ROUND(

        100.0
        *
        COUNT(*)
        /
        SUM(
            COUNT(*)
        ) OVER (),

        6

    ) AS percent_of_pairs

FROM read_parquet(
    '${OUTPUT_PARQUET}'
)

GROUP BY
    Correct_interaction

ORDER BY
    Correct_interaction DESC;



-- ============================================================
-- 22. POSITIVE RETENTION QC
--
-- Only evaluate known positive pairs here.
-- ============================================================

SELECT
    COUNT(*) AS total_positive_pairs,

    COUNT(*) FILTER (

        WHERE num_acc_in_host > 0

          AND num_acc_in_pathogen > 0

    ) AS positive_pairs_with_both_taxa_detected,

    COUNT(*) FILTER (

        WHERE num_acc_in_host = 0

           OR num_acc_in_pathogen = 0

    ) AS positive_pairs_lost_at_threshold,

    COUNT(*) FILTER (

        WHERE num_acc_shared > 0

    ) AS positive_pairs_with_shared_accessions

FROM read_parquet(
    '${OUTPUT_PARQUET}'
)

WHERE Correct_interaction = 1;



-- ============================================================
-- 23. POSITIVE VS UNLABELLED SCORE SUMMARY
-- ============================================================

SELECT
    Correct_interaction,

    COUNT(*) AS num_pairs,

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
)

GROUP BY
    Correct_interaction

ORDER BY
    Correct_interaction DESC;



-- ============================================================
-- 24. CO-OCCURRENCE SUMMARY BY LABEL
-- ============================================================

SELECT
    Correct_interaction,

    COUNT(*) AS total_pairs,

    COUNT(*) FILTER (

        WHERE num_acc_in_host > 0

          AND num_acc_in_pathogen > 0

    ) AS both_taxa_detected,

    COUNT(*) FILTER (

        WHERE num_acc_shared > 0

    ) AS pairs_with_shared_accessions,

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
)

GROUP BY
    Correct_interaction

ORDER BY
    Correct_interaction DESC;



-- ============================================================
-- 25. SCORE QC
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
-- 26. DISPLAY TOP KNOWN POSITIVE INTERACTIONS
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

WHERE Correct_interaction = 1

ORDER BY
    log_odds_acc DESC

LIMIT 30;



-- ============================================================
-- 27. DISPLAY TOP UNLABELLED INTERACTIONS
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

WHERE Correct_interaction = 0

ORDER BY
    log_odds_acc DESC

LIMIT 30;



-- ============================================================
-- DONE
-- ============================================================

SELECT
    'ALL BENCHMARK HOST x VIRUS PAIRS SCORED; 8,508 POSITIVES RETAINED'
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
