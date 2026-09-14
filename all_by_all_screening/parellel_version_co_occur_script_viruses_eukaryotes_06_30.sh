#!/bin/bash
#SBATCH --job-name=cooccur_virus_array
#SBATCH --account=def-acdoxey
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --array=0-9
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL

set -euo pipefail
set -x

# ---------------- cfg ----------------
DUCKDB_BIN="/home/smascar/scratch/STAT_2025/duckdb"

# INPUTS (host list is the same for every array task)
HOST_ID_FILE="eukaryota_tax_ids.txt"

# virus split files live here (created beforehand via split -n l/10 ...)
BACT_SPLIT_DIR="virus_splits"

# Parquet directories (partitioned folders are OK; we use /** globs)
parquet_path_host="/home/smascar/scratch/VIRUSES/STAT_with_Eukaryota_important_columns_split/"
parquet_path_pathogen="/home/smascar/scratch/VIRUSES/STAT_with_Viruses_important_columns_split/"

# THRESHOLDS (fixed per side)
HOST_T="5e-03"
PATH_T="1e-06"

# BATCH SIZES (within each array task)
HOST_BATCH_SIZE=1000
PATH_BATCH_SIZE=1000

# OUTPUT (task-specific)
OUT_DIR="cooccur_virus_outputs"
mkdir -p "$OUT_DIR"
TASK_ID="$(printf '%02d' "${SLURM_ARRAY_TASK_ID:-0}")"
output_file="${OUT_DIR}/host_virus_cooccur_part_${TASK_ID}.parquet"
batch_output_dir="${OUT_DIR}/batch_outputs_virus_part_${TASK_ID}_${SLURM_JOB_ID:-$$}"

# Use node-local scratch for DuckDB temp files
TMPDIR_DUCK="${SLURM_TMPDIR:-$PWD/.ducktmp_virus.$$}"
mkdir -p "$TMPDIR_DUCK"
mkdir -p "$batch_output_dir"
CLEANUP_TMP=$([[ -z "${SLURM_TMPDIR:-}" ]] && echo 1 || echo 0)
trap '[[ "$CLEANUP_TMP" = 1 ]] && rm -rf "$TMPDIR_DUCK" 2>/dev/null || true' EXIT

# ---------------- resources -----------
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
  MEM_MB=$SLURM_MEM_PER_NODE
elif [[ -n "${SLURM_MEM_PER_CPU:-}" && -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
  MEM_MB=$(( SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK ))
else
  MEM_MB=70000
fi
MEM_LIMIT_GB=$(( (MEM_MB * 9) / (10 * 1024) ))
(( MEM_LIMIT_GB < 1 )) && MEM_LIMIT_GB=1

MAX_THREADS="${DUCKDB_MAX_THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
(( THREADS > MAX_THREADS )) && THREADS="$MAX_THREADS"

# ---------------- checks -------------
[[ -s "$HOST_ID_FILE" ]] || { echo "Missing/empty $HOST_ID_FILE"; exit 1; }

# Resolve virus split file for this array task (supports with/without .txt extension)
PATH_ID_FILE_BASE="${BACT_SPLIT_DIR}/virus_tax_ids_${TASK_ID}"
if [[ -f "${PATH_ID_FILE_BASE}.txt" ]]; then
  PATH_ID_FILE="${PATH_ID_FILE_BASE}.txt"
elif [[ -f "${PATH_ID_FILE_BASE}" ]]; then
  PATH_ID_FILE="${PATH_ID_FILE_BASE}"
else
  echo "Missing virus split for task ${TASK_ID}: ${PATH_ID_FILE_BASE}[.txt]"
  exit 1
fi
[[ -s "$PATH_ID_FILE" ]] || { echo "Missing/empty $PATH_ID_FILE"; exit 1; }

echo "[$(date +'%T')] Array task: ${TASK_ID}"
echo "[$(date +'%T')] Host IDs: $HOST_ID_FILE"
echo "[$(date +'%T')] virus IDs (this task): $PATH_ID_FILE"
echo "[$(date +'%T')] Using DuckDB temp dir: $TMPDIR_DUCK"
df -h "$TMPDIR_DUCK" || true

# Normalize parquet dir to glob patterns (recursive)
parquet_glob_host="${parquet_path_host%/}/**"
parquet_glob_pathogen="${parquet_path_pathogen%/}/**"

# ---------------- batching -----------
echo "[$(date +'%T')] Creating batch files in: $batch_output_dir"

split -l "$HOST_BATCH_SIZE" --numeric-suffixes=1 --suffix-length=3 \
  "$HOST_ID_FILE" "${batch_output_dir}/host_batch_"

split -l "$PATH_BATCH_SIZE" --numeric-suffixes=1 --suffix-length=3 \
  "$PATH_ID_FILE" "${batch_output_dir}/path_batch_"

host_batches=($(ls "${batch_output_dir}"/host_batch_* | sort -V))
path_batches=($(ls "${batch_output_dir}"/path_batch_* | sort -V))

echo "[$(date +'%T')] Created ${#host_batches[@]} host batches and ${#path_batches[@]} pathogen batches"

batch_count=0
total_batches=$((${#host_batches[@]} * ${#path_batches[@]}))

# ---------------- main loop ----------
for host_batch in "${host_batches[@]}"; do
  for path_batch in "${path_batches[@]}"; do
    batch_count=$((batch_count + 1))
    echo "[$(date +'%T')] Processing batch ${batch_count}/${total_batches}: $(basename "$host_batch") x $(basename "$path_batch")"

    batch_output="${batch_output_dir}/batch_${batch_count}_output.parquet"

    "$DUCKDB_BIN" <<EOF
PRAGMA temp_directory='${TMPDIR_DUCK}';
PRAGMA max_temp_directory_size='2500GB';
PRAGMA enable_progress_bar=true;
SET threads=${THREADS};
SET memory_limit='${MEM_LIMIT_GB}GB';
SET preserve_insertion_order=false;

-- ID lists for this batch
CREATE TEMP TABLE hosts AS
SELECT DISTINCT CAST(TRIM(col) AS INTEGER) AS host_taxid
FROM read_csv('${host_batch}', header=false, columns={'col':'VARCHAR'})
WHERE TRIM(col) ~ '^[0-9]+\$';

CREATE TEMP TABLE pathogens AS
SELECT DISTINCT CAST(TRIM(col) AS INTEGER) AS pathogen_taxid
FROM read_csv('${path_batch}', header=false, columns={'col':'VARCHAR'})
WHERE TRIM(col) ~ '^[0-9]+\$';

-- Filter and de-dup host data with abundance threshold (INCLUDING BIOPROJECT)
CREATE TEMP TABLE host_data AS
SELECT
  tax_id,
  acc,
  bioproject,
  MAX(total_abundance) AS total_abundance
FROM read_parquet('${parquet_glob_host}')
WHERE tax_id IN (SELECT host_taxid FROM hosts)
  AND total_abundance >= ${HOST_T}
GROUP BY 1,2,3;

-- Filter and de-dup pathogen data with abundance threshold (INCLUDING BIOPROJECT)
CREATE TEMP TABLE path_data AS
SELECT
  tax_id,
  acc,
  bioproject,
  MAX(total_abundance) AS total_abundance
FROM read_parquet('${parquet_glob_pathogen}')
WHERE tax_id IN (SELECT pathogen_taxid FROM pathogens)
  AND total_abundance >= ${PATH_T}
GROUP BY 1,2,3;

ANALYZE host_data;
ANALYZE path_data;

-- Per-side counts (already filtered by thresholds)
CREATE TEMP TABLE host_counts AS
SELECT tax_id AS host_taxid, COUNT(*) AS num_acc_in_host
FROM host_data
GROUP BY 1;

CREATE TEMP TABLE pathogen_counts AS
SELECT tax_id AS pathogen_taxid, COUNT(*) AS num_acc_in_pathogen
FROM path_data
GROUP BY 1;

-- Per-side BIOPROJECT diversity (unique bioproject IDs after filtering)
CREATE TEMP TABLE host_bioprojects AS
SELECT
  tax_id AS host_taxid,
  COUNT(DISTINCT bioproject) AS num_bioprojects_in_host
FROM host_data
WHERE bioproject IS NOT NULL
  AND bioproject <> ''
GROUP BY 1;

CREATE TEMP TABLE pathogen_bioprojects AS
SELECT
  tax_id AS pathogen_taxid,
  COUNT(DISTINCT bioproject) AS num_bioprojects_in_pathogen
FROM path_data
WHERE bioproject IS NOT NULL
  AND bioproject <> ''
GROUP BY 1;

-- Shared accs (both sides already meet thresholds)
CREATE TEMP TABLE shared_counts AS
SELECT
  h.tax_id AS host_taxid,
  p.tax_id AS pathogen_taxid,
  COUNT(*) AS num_acc_shared
FROM host_data h
JOIN path_data p USING (acc)
GROUP BY 1,2;

-- Shared BIOPROJECTS (unique bioproject IDs among shared accs)
CREATE TEMP TABLE shared_bioprojects AS
SELECT
  h.tax_id AS host_taxid,
  p.tax_id AS pathogen_taxid,
  COUNT(DISTINCT h.bioproject) AS num_bioprojects_shared
FROM host_data h
JOIN path_data p USING (acc)
WHERE h.bioproject IS NOT NULL
  AND h.bioproject <> ''
GROUP BY 1,2;

-- Final rows for this batch
COPY (
  SELECT
    h.host_taxid,
    p.pathogen_taxid,
    CAST(${HOST_T} AS DOUBLE) AS host_threshold,
    CAST(${PATH_T} AS DOUBLE) AS pathogen_threshold,

    COALESCE(hc.num_acc_in_host, 0)             AS num_acc_in_host,
    COALESCE(pc.num_acc_in_pathogen, 0)         AS num_acc_in_pathogen,
    COALESCE(sc.num_acc_shared, 0)              AS num_acc_shared,

    COALESCE(hb.num_bioprojects_in_host, 0)     AS num_bioprojects_in_host,
    COALESCE(pb.num_bioprojects_in_pathogen, 0) AS num_bioprojects_in_pathogen,
    COALESCE(sb.num_bioprojects_shared, 0)      AS num_bioprojects_shared

  FROM hosts h
  CROSS JOIN pathogens p

  LEFT JOIN host_counts          hc ON hc.host_taxid      = h.host_taxid
  LEFT JOIN pathogen_counts      pc ON pc.pathogen_taxid  = p.pathogen_taxid
  LEFT JOIN shared_counts        sc ON sc.host_taxid      = h.host_taxid
                                   AND sc.pathogen_taxid = p.pathogen_taxid

  LEFT JOIN host_bioprojects     hb ON hb.host_taxid      = h.host_taxid
  LEFT JOIN pathogen_bioprojects pb ON pb.pathogen_taxid  = p.pathogen_taxid
  LEFT JOIN shared_bioprojects   sb ON sb.host_taxid      = h.host_taxid
                                   AND sb.pathogen_taxid = p.pathogen_taxid

  ORDER BY h.host_taxid, p.pathogen_taxid
) TO '${batch_output}' (FORMAT PARQUET);
EOF

    echo "[$(date +'%T')] Batch ${batch_count}/${total_batches} completed"
  done
done

# Cleanup ID batch text files (keep batch parquet outputs until after final combine)
rm -f "${batch_output_dir}"/host_batch_* "${batch_output_dir}"/path_batch_*

# ---------------- combine ------------
echo "[$(date +'%T')] Combining batch parquet files into: ${output_file}"

"$DUCKDB_BIN" <<EOF
PRAGMA temp_directory='${TMPDIR_DUCK}';
SET threads=${THREADS};
SET memory_limit='${MEM_LIMIT_GB}GB';

COPY (
  SELECT * FROM read_parquet('${batch_output_dir}/batch_*_output.parquet')
) TO '${output_file}' (FORMAT PARQUET);
EOF

# Cleanup per-batch parquet outputs to save space
rm -f "${batch_output_dir}"/batch_*_output.parquet
rmdir "${batch_output_dir}" 2>/dev/null || true

echo "[$(date +'%T')] Done. Output for task ${TASK_ID}: ${output_file}"
ls -lh "${output_file}" || true

# Parquet row count isn't available via wc; this is just a harmless check
wc -l "${output_file}" 2>/dev/null || echo "Row count not available for parquet file"


