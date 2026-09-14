#!/usr/bin/env python3
"""
Recover SRA accession metadata for labelled virus-host pairs.

Purpose
-------
This script takes a labelled host-virus pair CSV, identifies SRA accessions where
both the host and virus/pathogen are detected above fixed abundance thresholds,
joins those shared accessions to SRA metadata, and summarizes metadata fields
such as assay_type and librarysource.

Main outputs
------------
- shared_acc_selected.parquet
- shared_acc_selected_counts.csv/parquet
- shared_acc_selected_with_meta.parquet
- shared_meta_summary_long.csv/parquet
- shared_meta_summary_wide.csv/parquet
- metadata_composition_pair_counts.csv/parquet
- background_metadata_counts.csv/parquet
- metadata_enrichment_by_label.csv/parquet

Notes
-----
The script is virus-specific by default, but column names use the generic term
"pathogen" internally. It accepts labelled CSVs with either microbe_taxid or
pathogen_taxid columns.
"""

from __future__ import annotations

import argparse
import glob
import shutil
from pathlib import Path
from typing import Iterable

import duckdb
import pandas as pd


# =========================================================
# DEFAULT CONFIGURATION
# =========================================================

DEFAULT_HOST_GLOB = "STAT_with_Eukaryota_important_columns_split/**"
DEFAULT_PATHOGEN_GLOB = "STAT_with_Viruses_important_columns_split/**"
DEFAULT_META_GLOB = "/home/smascar/projects/def-acdoxey/SRA_DATA/DATA/02_02_2026/SRA_META_DATA-02_02_2026/*"
DEFAULT_PAIRS_CSV = "virus_host_db_labelled_interactions.csv"

# Optional. This is not required unless --use-pair-metrics-filter is supplied.
DEFAULT_PAIR_METRICS_FILE = "pair_threshold_grid_results/pair_metrics_all_thresholds.parquet"

DEFAULT_OUTDIR = "/home/smascar/projects/def-acdoxey/smascar/novel_host_virus_predictions/metadata_recovery_virus_hostdb_ht1e-04_pt1e-08"

DEFAULT_HOST_THRESHOLD = 1e-4
DEFAULT_PATHOGEN_THRESHOLD = 1e-8

ACC_COL = "acc"
TAX_COL = "tax_id"
ABUND_COL = "total_abundance"

DEFAULT_THREADS = 8
DEFAULT_MEMORY_LIMIT = "250GB"

# Only fields found in the SRA metadata parquet will be used.
DEFAULT_METADATA_FIELDS = [
    "assay_type",
    "instrument",
    "librarylayout",
    "libraryselection",
    "librarysource",
    "platform",
    "organism",
    "geo_loc_name_country_calc",
]


# =========================================================
# HELPERS
# =========================================================

def resolve_parquet_files(root_or_glob: str) -> list[str]:
    """
    Resolve parquet inputs robustly.

    This intentionally supports both conventional files ending in .parquet and
    extensionless parquet part files, because some SRA/STAT exports are valid
    parquet files but do not have a .parquet suffix.

    Accepted inputs:
      - /path/to/dir
      - /path/to/dir/*
      - /path/to/dir/**
      - /path/to/dir/**/*.parquet
      - explicit glob patterns
    """
    pattern = root_or_glob.strip()

    def good_file(x: str) -> bool:
        path = Path(x)
        if not path.is_file():
            return False
        name = path.name
        if name.startswith(".") or name.startswith("_"):
            return False
        if name.endswith(('.crc', '.json', '.txt', '.csv', '.tsv', '.log', '.md')):
            return False
        return True

    def unique_sorted(xs):
        return sorted({str(Path(x)) for x in xs if good_file(x)})

    # Build candidate glob patterns. Search for .parquet first.
    parquet_patterns = []
    all_file_patterns = []

    if pattern.endswith("/**"):
        base = pattern[:-3].rstrip("/")
        parquet_patterns.extend([
            f"{base}/*.parquet",
            f"{base}/**/*.parquet",
        ])
        all_file_patterns.extend([
            f"{base}/*",
            f"{base}/**/*",
        ])
    elif "*" in pattern:
        parquet_patterns.append(pattern)
        all_file_patterns.append(pattern)
        # If the user supplied a directory wildcard like /dir/*, also allow nested.
        if pattern.endswith("/*"):
            base = pattern[:-2].rstrip("/")
            parquet_patterns.append(f"{base}/**/*.parquet")
            all_file_patterns.append(f"{base}/**/*")
    else:
        pth = Path(pattern)
        if pth.exists() and pth.is_dir():
            parquet_patterns.extend([
                str(pth / "*.parquet"),
                str(pth / "**" / "*.parquet"),
            ])
            all_file_patterns.extend([
                str(pth / "*"),
                str(pth / "**" / "*"),
            ])
        else:
            parquet_patterns.append(pattern)
            all_file_patterns.append(pattern)

    parquet_files = []
    for pat in parquet_patterns:
        parquet_files.extend(glob.glob(pat, recursive=True))
    parquet_files = unique_sorted(parquet_files)
    if parquet_files:
        return parquet_files

    # Fallback: no .parquet suffixes were found, so try extensionless parquet
    # part files. DuckDB will validate them when read_parquet() is called.
    all_files = []
    for pat in all_file_patterns:
        all_files.extend(glob.glob(pat, recursive=True))
    all_files = unique_sorted(all_files)
    if all_files:
        print(
            f"  Warning: no files ending in .parquet found for {root_or_glob}. "
            f"Using {len(all_files):,} regular non-hidden file(s) instead; "
            "these must be valid parquet files."
        )
        return all_files

    clean_pattern = pattern.replace("/**", "")
    raise FileNotFoundError(
        f"No parquet or candidate extensionless parquet files found for: {root_or_glob}\n"
        "Check the directory path and run something like:\n"
        f"  find {clean_pattern} -maxdepth 3 -type f | head\n"
        f"  find {clean_pattern} -type f -name '*.parquet' | wc -l"
    )


def path_string(path_obj: Path | str) -> str:
    return Path(path_obj).as_posix()


def sql_quote_identifier(col: str) -> str:
    return '"' + col.replace('"', '""') + '"'


def sql_quote_literal(x: str) -> str:
    return "'" + x.replace("'", "''") + "'"


def standardize_pair_df(df: pd.DataFrame, require_label: bool = True) -> pd.DataFrame:
    """
    Standardize pair CSV columns.

    Accepted host columns:
      host_taxid, host_tax_id, host

    Accepted pathogen columns:
      pathogen_taxid, pathogen_tax_id, microbe_taxid, microbe_tax_id,
      virus_taxid, virus_tax_id, pathogen, microbe, virus

    Accepted label columns:
      Correct_interaction, correct_interaction, label

    Output columns:
      host_taxid, pathogen_taxid, Correct_interaction
    """
    host_col = None
    pathogen_col = None
    label_col = None

    for c in df.columns:
        cl = c.lower().strip()

        if cl in {"host_taxid", "host_tax_id", "host"}:
            host_col = c

        if cl in {
            "pathogen_taxid",
            "pathogen_tax_id",
            "microbe_taxid",
            "microbe_tax_id",
            "virus_taxid",
            "virus_tax_id",
            "pathogen",
            "microbe",
            "virus",
        }:
            pathogen_col = c

        if cl in {"correct_interaction", "label"}:
            label_col = c

    if host_col is None or pathogen_col is None:
        raise ValueError(
            "Could not find host/pathogen taxid columns in PAIRS_CSV. "
            "Expected host_taxid plus pathogen_taxid/microbe_taxid/virus_taxid."
        )

    if require_label and label_col is None:
        raise ValueError(
            "Could not find Correct_interaction or label column in PAIRS_CSV. "
            "This script needs labels so metadata can be stratified by positive/control."
        )

    out = pd.DataFrame({
        "host_taxid": pd.to_numeric(df[host_col], errors="coerce"),
        "pathogen_taxid": pd.to_numeric(df[pathogen_col], errors="coerce"),
    })

    if label_col is not None:
        out["Correct_interaction"] = pd.to_numeric(df[label_col], errors="coerce")
    else:
        out["Correct_interaction"] = 1

    out = out.dropna(subset=["host_taxid", "pathogen_taxid"]).copy()
    out["host_taxid"] = out["host_taxid"].astype("int64")
    out["pathogen_taxid"] = out["pathogen_taxid"].astype("int64")

    if require_label:
        out = out.dropna(subset=["Correct_interaction"]).copy()
        out["Correct_interaction"] = out["Correct_interaction"].astype("int64")
        valid_labels = set(out["Correct_interaction"].unique())
        if not valid_labels.issubset({0, 1}):
            raise ValueError(f"Correct_interaction must contain only 0/1. Found: {valid_labels}")

    out = out.drop_duplicates(subset=["host_taxid", "pathogen_taxid", "Correct_interaction"])
    return out


def get_table_columns(con: duckdb.DuckDBPyConnection, relation_sql: str) -> list[str]:
    """Return columns for a SQL relation/view."""
    return [row[0] for row in con.execute(f"DESCRIBE SELECT * FROM {relation_sql} LIMIT 0").fetchall()]


def get_available_metadata_fields(
    con: duckdb.DuckDBPyConnection,
    requested_fields: Iterable[str],
) -> list[str]:
    meta_cols = set(get_table_columns(con, "meta_src"))
    if ACC_COL not in meta_cols:
        raise ValueError(f"Metadata source does not contain required accession column: {ACC_COL}")

    available = [c for c in requested_fields if c in meta_cols]
    missing = [c for c in requested_fields if c not in meta_cols]

    if missing:
        print("WARNING: These requested metadata fields were not found and will be skipped:")
        for c in missing:
            print(f"  - {c}")

    if not available:
        raise ValueError("None of the requested metadata fields were found in the metadata source.")

    print("Metadata fields that will be summarized:")
    for c in available:
        print(f"  - {c}")

    return available


def normalize_pair_metrics_if_requested(
    con: duckdb.DuckDBPyConnection,
    pair_metrics_file: Path,
    host_threshold: float,
    pathogen_threshold: float,
) -> bool:
    """
    Normalize a pair metrics file into a temp view with columns:
      host_taxid, pathogen_taxid, host_threshold, pathogen_threshold

    Returns True if a pair metrics filter was created.
    """
    if not pair_metrics_file.exists():
        raise FileNotFoundError(f"Missing pair metrics file: {pair_metrics_file}")

    con.execute(f"""
    CREATE OR REPLACE VIEW pair_metrics_raw AS
    SELECT *
    FROM read_parquet('{path_string(pair_metrics_file)}');
    """)

    cols = set(get_table_columns(con, "pair_metrics_raw"))

    required_host = "host_taxid"
    if required_host not in cols:
        raise ValueError("pair_metrics file is missing host_taxid")

    pathogen_col = None
    for c in ["pathogen_taxid", "microbe_taxid", "virus_taxid"]:
        if c in cols:
            pathogen_col = c
            break
    if pathogen_col is None:
        raise ValueError("pair_metrics file is missing pathogen_taxid/microbe_taxid/virus_taxid")

    host_threshold_col = "host_threshold" if "host_threshold" in cols else None
    pathogen_threshold_col = None
    for c in ["pathogen_threshold", "microbe_threshold", "virus_threshold"]:
        if c in cols:
            pathogen_threshold_col = c
            break

    if host_threshold_col is None or pathogen_threshold_col is None:
        raise ValueError(
            "pair_metrics file must contain host_threshold and pathogen_threshold/microbe_threshold."
        )

    eps = 1e-15
    con.execute(f"""
    CREATE OR REPLACE TEMP TABLE pair_metrics_filter AS
    SELECT DISTINCT
        CAST(host_taxid AS BIGINT) AS host_taxid,
        CAST({sql_quote_identifier(pathogen_col)} AS BIGINT) AS pathogen_taxid
    FROM pair_metrics_raw
    WHERE ABS(CAST({sql_quote_identifier(host_threshold_col)} AS DOUBLE) - {host_threshold}) <= {eps}
      AND ABS(CAST({sql_quote_identifier(pathogen_threshold_col)} AS DOUBLE) - {pathogen_threshold}) <= {eps};
    """)

    n = con.execute("SELECT COUNT(*) FROM pair_metrics_filter").fetchone()[0]
    print(f"Pair metrics filter rows at selected threshold: {n:,}")
    if n == 0:
        raise ValueError(
            "The pair metrics file has zero rows at the selected thresholds. "
            "Either use the correct threshold pair or run without --use-pair-metrics-filter."
        )

    return True


def write_df_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


# =========================================================
# ARGUMENTS
# =========================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Recover SRA metadata for labelled host-virus pairs."
    )

    parser.add_argument("--host-glob", default=DEFAULT_HOST_GLOB)
    parser.add_argument("--pathogen-glob", default=DEFAULT_PATHOGEN_GLOB)
    parser.add_argument("--meta-glob", default=DEFAULT_META_GLOB)
    parser.add_argument("--pairs-csv", default=DEFAULT_PAIRS_CSV)
    parser.add_argument("--pair-metrics-file", default=DEFAULT_PAIR_METRICS_FILE)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)

    parser.add_argument("--host-threshold", type=float, default=DEFAULT_HOST_THRESHOLD)
    parser.add_argument("--pathogen-threshold", type=float, default=DEFAULT_PATHOGEN_THRESHOLD)

    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--memory-limit", default=DEFAULT_MEMORY_LIMIT)

    parser.add_argument(
        "--metadata-fields",
        default=",".join(DEFAULT_METADATA_FIELDS),
        help="Comma-separated SRA metadata fields to summarize.",
    )

    parser.add_argument(
        "--use-pair-metrics-filter",
        action="store_true",
        help=(
            "Optional: keep only pairs that are present in pair_metrics_all_thresholds.parquet "
            "at the selected thresholds. Not required for metadata recovery."
        ),
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete the output directory before running.",
    )

    parser.add_argument(
        "--write-large-csv",
        action="store_true",
        help=(
            "Also write large per-accession CSV files. Parquet is always written. "
            "Leave this off for large runs."
        ),
    )

    return parser.parse_args()


# =========================================================
# MAIN
# =========================================================

def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir)
    pairs_csv = Path(args.pairs_csv)
    pair_metrics_file = Path(args.pair_metrics_file)
    metadata_fields = [x.strip() for x in args.metadata_fields.split(",") if x.strip()]

    if outdir.exists() and args.overwrite:
        print(f"Deleting old output directory: {outdir}")
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not pairs_csv.exists():
        raise FileNotFoundError(f"Missing labelled pair CSV: {pairs_csv}")

    print("\n========== SRA metadata recovery for labelled host-virus pairs ==========")
    print(f"Host STAT glob:        {args.host_glob}")
    print(f"Pathogen STAT glob:    {args.pathogen_glob}")
    print(f"Metadata glob:         {args.meta_glob}")
    print(f"Pairs CSV:             {pairs_csv}")
    print(f"Pair metrics file:     {pair_metrics_file}")
    print(f"Use pair metrics filt: {args.use_pair_metrics_filter}")
    print(f"Host threshold:        {args.host_threshold}")
    print(f"Pathogen threshold:    {args.pathogen_threshold}")
    print(f"Output directory:      {outdir}")
    print(f"Threads:               {args.threads}")
    print(f"DuckDB memory limit:   {args.memory_limit}\n")

    print("Resolving parquet files...")
    host_files = resolve_parquet_files(args.host_glob)
    pathogen_files = resolve_parquet_files(args.pathogen_glob)
    meta_files = resolve_parquet_files(args.meta_glob)
    print(f"  Host parquet files:     {len(host_files):,}")
    print(f"  Pathogen parquet files: {len(pathogen_files):,}")
    print(f"  Metadata parquet files: {len(meta_files):,}")

    con = duckdb.connect()
    con.execute(f"PRAGMA threads={args.threads};")
    con.execute(f"PRAGMA memory_limit='{args.memory_limit}';")
    con.execute("PRAGMA preserve_insertion_order=false;")

    con.read_parquet(host_files).create_view("host_src")
    con.read_parquet(pathogen_files).create_view("pathogen_src")
    con.read_parquet(meta_files).create_view("meta_src")

    # -----------------------------------------------------
    # 1. Load threshold and labelled pairs
    # -----------------------------------------------------
    selected_thresholds = pd.DataFrame({
        "host_threshold": [float(args.host_threshold)],
        "pathogen_threshold": [float(args.pathogen_threshold)],
    })
    write_df_csv(selected_thresholds, outdir / "selected_thresholds.csv")
    con.register("selected_thresholds_df", selected_thresholds)

    selected_pairs = standardize_pair_df(pd.read_csv(pairs_csv), require_label=True)
    if selected_pairs.empty:
        raise ValueError("No valid labelled pairs found in pairs CSV.")

    write_df_csv(selected_pairs, outdir / "selected_pairs_with_correct_interaction.csv")
    print(f"\nSelected labelled pairs: {len(selected_pairs):,}")
    print("Label counts:")
    print(selected_pairs["Correct_interaction"].value_counts().sort_index())
    con.register("selected_pairs_raw_df", selected_pairs)

    # -----------------------------------------------------
    # 2. Optional pair-metrics filtering
    # -----------------------------------------------------
    if args.use_pair_metrics_filter:
        normalize_pair_metrics_if_requested(
            con=con,
            pair_metrics_file=pair_metrics_file,
            host_threshold=args.host_threshold,
            pathogen_threshold=args.pathogen_threshold,
        )
        selected_pairs_in_metrics = con.execute("""
            SELECT DISTINCT
                s.host_taxid,
                s.pathogen_taxid,
                s.Correct_interaction
            FROM selected_pairs_raw_df s
            JOIN pair_metrics_filter p
              ON s.host_taxid = p.host_taxid
             AND s.pathogen_taxid = p.pathogen_taxid;
        """).df()
        if selected_pairs_in_metrics.empty:
            raise ValueError("No selected pairs remained after pair metrics filtering.")
        print(f"Selected pairs after pair-metrics filter: {len(selected_pairs_in_metrics):,}")
    else:
        selected_pairs_in_metrics = selected_pairs.copy()
        print("Skipping pair-metrics filter. Metadata recovery will use all pairs in PAIRS_CSV.")

    selected_pairs_in_metrics["host_taxid"] = selected_pairs_in_metrics["host_taxid"].astype("int64")
    selected_pairs_in_metrics["pathogen_taxid"] = selected_pairs_in_metrics["pathogen_taxid"].astype("int64")
    selected_pairs_in_metrics["Correct_interaction"] = selected_pairs_in_metrics["Correct_interaction"].astype("int64")

    write_df_csv(
        selected_pairs_in_metrics,
        outdir / "selected_pairs_used_with_correct_interaction.csv",
    )
    con.unregister("selected_pairs_raw_df")
    con.register("selected_pairs_df", selected_pairs_in_metrics)

    # -----------------------------------------------------
    # 3. Selected host/pathogen taxid lists
    # -----------------------------------------------------
    con.execute("""
    CREATE OR REPLACE TEMP TABLE selected_hosts AS
    SELECT DISTINCT host_taxid
    FROM selected_pairs_df;
    """)

    con.execute("""
    CREATE OR REPLACE TEMP TABLE selected_pathogens AS
    SELECT DISTINCT pathogen_taxid
    FROM selected_pairs_df;
    """)

    n_hosts = con.execute("SELECT COUNT(*) FROM selected_hosts").fetchone()[0]
    n_pathogens = con.execute("SELECT COUNT(*) FROM selected_pathogens").fetchone()[0]
    print(f"Selected host taxids:     {n_hosts:,}")
    print(f"Selected pathogen taxids: {n_pathogens:,}")

    # -----------------------------------------------------
    # 4. Aggregate max abundance for selected taxa only
    # -----------------------------------------------------
    print("\nBuilding aggregated host abundance table for selected hosts...")
    con.execute(f"""
    COPY (
        SELECT
            src.{ACC_COL} AS acc,
            CAST(src.{TAX_COL} AS BIGINT) AS host_taxid,
            MAX(CAST(src.{ABUND_COL} AS DOUBLE)) AS max_abundance
        FROM host_src src
        JOIN selected_hosts h
          ON CAST(src.{TAX_COL} AS BIGINT) = h.host_taxid
        WHERE src.{ACC_COL} IS NOT NULL
        GROUP BY 1, 2
    )
    TO '{path_string(outdir / "host_agg_selected.parquet")}'
    (FORMAT PARQUET);
    """)

    print("Building aggregated pathogen abundance table for selected pathogens...")
    con.execute(f"""
    COPY (
        SELECT
            src.{ACC_COL} AS acc,
            CAST(src.{TAX_COL} AS BIGINT) AS pathogen_taxid,
            MAX(CAST(src.{ABUND_COL} AS DOUBLE)) AS max_abundance
        FROM pathogen_src src
        JOIN selected_pathogens p
          ON CAST(src.{TAX_COL} AS BIGINT) = p.pathogen_taxid
        WHERE src.{ACC_COL} IS NOT NULL
        GROUP BY 1, 2
    )
    TO '{path_string(outdir / "pathogen_agg_selected.parquet")}'
    (FORMAT PARQUET);
    """)

    # -----------------------------------------------------
    # 5. Build thresholded detections
    # -----------------------------------------------------
    print("Building host detections for selected threshold...")
    con.execute(f"""
    COPY (
        SELECT
            a.acc,
            a.host_taxid,
            {args.host_threshold}::DOUBLE AS host_threshold
        FROM read_parquet('{path_string(outdir / "host_agg_selected.parquet")}') a
        WHERE a.max_abundance >= {args.host_threshold}
    )
    TO '{path_string(outdir / "host_detects_selected_threshold.parquet")}'
    (FORMAT PARQUET);
    """)

    print("Building pathogen detections for selected threshold...")
    con.execute(f"""
    COPY (
        SELECT
            a.acc,
            a.pathogen_taxid,
            {args.pathogen_threshold}::DOUBLE AS pathogen_threshold
        FROM read_parquet('{path_string(outdir / "pathogen_agg_selected.parquet")}') a
        WHERE a.max_abundance >= {args.pathogen_threshold}
    )
    TO '{path_string(outdir / "pathogen_detects_selected_threshold.parquet")}'
    (FORMAT PARQUET);
    """)

    host_det_n = con.execute(f"SELECT COUNT(*) FROM read_parquet('{path_string(outdir / 'host_detects_selected_threshold.parquet')}')").fetchone()[0]
    path_det_n = con.execute(f"SELECT COUNT(*) FROM read_parquet('{path_string(outdir / 'pathogen_detects_selected_threshold.parquet')}')").fetchone()[0]
    print(f"Host detection rows:     {host_det_n:,}")
    print(f"Pathogen detection rows: {path_det_n:,}")

    # -----------------------------------------------------
    # 6. Recover shared accessions for selected pairs
    # -----------------------------------------------------
    print("Recovering shared accessions for selected pairs...")
    con.execute(f"""
    COPY (
        SELECT
            p.host_taxid,
            p.pathogen_taxid,
            p.Correct_interaction,
            {args.host_threshold}::DOUBLE AS host_threshold,
            {args.pathogen_threshold}::DOUBLE AS pathogen_threshold,
            h.acc
        FROM selected_pairs_df p
        JOIN read_parquet('{path_string(outdir / "host_detects_selected_threshold.parquet")}') h
          ON p.host_taxid = h.host_taxid
        JOIN read_parquet('{path_string(outdir / "pathogen_detects_selected_threshold.parquet")}') v
          ON p.pathogen_taxid = v.pathogen_taxid
         AND h.acc = v.acc
    )
    TO '{path_string(outdir / "shared_acc_selected.parquet")}'
    (FORMAT PARQUET);
    """)

    if args.write_large_csv:
        con.execute(f"""
        COPY (
            SELECT *
            FROM read_parquet('{path_string(outdir / "shared_acc_selected.parquet")}')
        )
        TO '{path_string(outdir / "shared_acc_selected.csv")}'
        (HEADER, DELIMITER ',');
        """)

    # -----------------------------------------------------
    # 7. Summarize shared accession counts per pair
    # -----------------------------------------------------
    print("Summarizing shared accession counts per pair...")
    con.execute(f"""
    COPY (
        SELECT
            host_taxid,
            pathogen_taxid,
            Correct_interaction,
            host_threshold,
            pathogen_threshold,
            COUNT(*)::BIGINT AS n_shared_acc
        FROM read_parquet('{path_string(outdir / "shared_acc_selected.parquet")}')
        GROUP BY 1, 2, 3, 4, 5
        ORDER BY Correct_interaction DESC, n_shared_acc DESC
    )
    TO '{path_string(outdir / "shared_acc_selected_counts.parquet")}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT *
        FROM read_parquet('{path_string(outdir / "shared_acc_selected_counts.parquet")}')
    )
    TO '{path_string(outdir / "shared_acc_selected_counts.csv")}'
    (HEADER, DELIMITER ',');
    """)

    # -----------------------------------------------------
    # 8. Build metadata by accession and background metadata counts
    # -----------------------------------------------------
    print("Building accession-level metadata table...")
    available_fields = get_available_metadata_fields(con, metadata_fields)

    meta_select_parts = []
    for c in available_fields:
        q = sql_quote_identifier(c)
        meta_select_parts.append(f"ANY_VALUE({q}) AS {q}")
    meta_select = ",\n            ".join(meta_select_parts)

    con.execute(f"""
    COPY (
        SELECT
            {ACC_COL} AS acc,
            {meta_select}
        FROM meta_src
        WHERE {ACC_COL} IS NOT NULL
        GROUP BY {ACC_COL}
    )
    TO '{path_string(outdir / "meta_by_acc.parquet")}'
    (FORMAT PARQUET);
    """)

    print("Building SRA background metadata composition...")
    background_parts = []
    for field in available_fields:
        q = sql_quote_identifier(field)
        background_parts.append(f"""
        SELECT
            {sql_quote_literal(field)} AS field_name,
            COALESCE(NULLIF(CAST({q} AS VARCHAR), ''), 'Unknown') AS field_value,
            COUNT(*)::BIGINT AS n_sra_accessions
        FROM read_parquet('{path_string(outdir / "meta_by_acc.parquet")}')
        GROUP BY 1, 2
        """)
    background_sql = "\nUNION ALL\n".join(background_parts)

    con.execute(f"""
    COPY (
        WITH bg AS (
            {background_sql}
        ), totals AS (
            SELECT field_name, SUM(n_sra_accessions)::BIGINT AS total_sra_accessions
            FROM bg
            GROUP BY 1
        )
        SELECT
            bg.field_name,
            bg.field_value,
            bg.n_sra_accessions,
            totals.total_sra_accessions,
            CASE
                WHEN totals.total_sra_accessions > 0
                THEN bg.n_sra_accessions * 1.0 / totals.total_sra_accessions
                ELSE NULL
            END AS fraction_sra_accessions
        FROM bg
        JOIN totals USING (field_name)
        ORDER BY field_name, n_sra_accessions DESC
    )
    TO '{path_string(outdir / "background_metadata_counts.parquet")}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT * FROM read_parquet('{path_string(outdir / "background_metadata_counts.parquet")}')
    )
    TO '{path_string(outdir / "background_metadata_counts.csv")}'
    (HEADER, DELIMITER ',');
    """)

    # -----------------------------------------------------
    # 9. Join shared accessions to metadata
    # -----------------------------------------------------
    print("Joining shared accessions to metadata...")
    meta_cols = ",\n            ".join(
        [f"m.{sql_quote_identifier(c)} AS {sql_quote_identifier(c)}" for c in available_fields]
    )

    con.execute(f"""
    COPY (
        SELECT
            s.host_taxid,
            s.pathogen_taxid,
            s.Correct_interaction,
            s.host_threshold,
            s.pathogen_threshold,
            s.acc,
            {meta_cols}
        FROM read_parquet('{path_string(outdir / "shared_acc_selected.parquet")}') s
        LEFT JOIN read_parquet('{path_string(outdir / "meta_by_acc.parquet")}') m
          USING (acc)
    )
    TO '{path_string(outdir / "shared_acc_selected_with_meta.parquet")}'
    (FORMAT PARQUET);
    """)

    if args.write_large_csv:
        con.execute(f"""
        COPY (
            SELECT *
            FROM read_parquet('{path_string(outdir / "shared_acc_selected_with_meta.parquet")}')
        )
        TO '{path_string(outdir / "shared_acc_selected_with_meta.csv")}'
        (HEADER, DELIMITER ',');
        """)

    # -----------------------------------------------------
    # 10. Wide and long metadata summaries
    # -----------------------------------------------------
    print("Building wide metadata summary by pair and threshold...")
    group_cols = ",\n            ".join([sql_quote_identifier(c) for c in available_fields])

    con.execute(f"""
    COPY (
        SELECT
            host_taxid,
            pathogen_taxid,
            Correct_interaction,
            host_threshold,
            pathogen_threshold,
            {group_cols},
            COUNT(*)::BIGINT AS n_shared_acc
        FROM read_parquet('{path_string(outdir / "shared_acc_selected_with_meta.parquet")}')
        GROUP BY
            host_taxid,
            pathogen_taxid,
            Correct_interaction,
            host_threshold,
            pathogen_threshold,
            {group_cols}
    )
    TO '{path_string(outdir / "shared_meta_summary_wide.parquet")}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT *
        FROM read_parquet('{path_string(outdir / "shared_meta_summary_wide.parquet")}')
    )
    TO '{path_string(outdir / "shared_meta_summary_wide.csv")}'
    (HEADER, DELIMITER ',');
    """)

    print("Building long metadata summary by field and value...")
    union_parts = []
    for field in available_fields:
        q = sql_quote_identifier(field)
        union_parts.append(f"""
        SELECT
            host_taxid,
            pathogen_taxid,
            Correct_interaction,
            host_threshold,
            pathogen_threshold,
            {sql_quote_literal(field)} AS field_name,
            COALESCE(NULLIF(CAST({q} AS VARCHAR), ''), 'Unknown') AS field_value,
            COUNT(*)::BIGINT AS n_shared_acc
        FROM read_parquet('{path_string(outdir / "shared_acc_selected_with_meta.parquet")}')
        GROUP BY 1, 2, 3, 4, 5, 6, 7
        """)
    long_sql = "\nUNION ALL\n".join(union_parts)

    con.execute(f"""
    COPY (
        {long_sql}
    )
    TO '{path_string(outdir / "shared_meta_summary_long.parquet")}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT *
        FROM read_parquet('{path_string(outdir / "shared_meta_summary_long.parquet")}')
    )
    TO '{path_string(outdir / "shared_meta_summary_long.csv")}'
    (HEADER, DELIMITER ',');
    """)

    # -----------------------------------------------------
    # 11. Metadata composition by label
    # -----------------------------------------------------
    print("Building metadata composition table by label...")
    con.execute(f"""
    COPY (
        SELECT
            Correct_interaction,
            host_threshold,
            pathogen_threshold,
            field_name,
            field_value,
            COUNT(DISTINCT CAST(host_taxid AS VARCHAR) || '__' || CAST(pathogen_taxid AS VARCHAR))::BIGINT AS n_unique_pairs,
            SUM(n_shared_acc)::BIGINT AS n_shared_accessions
        FROM read_parquet('{path_string(outdir / "shared_meta_summary_long.parquet")}')
        GROUP BY 1, 2, 3, 4, 5
        ORDER BY Correct_interaction DESC, field_name, n_shared_accessions DESC
    )
    TO '{path_string(outdir / "metadata_composition_pair_counts.parquet")}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT *
        FROM read_parquet('{path_string(outdir / "metadata_composition_pair_counts.parquet")}')
    )
    TO '{path_string(outdir / "metadata_composition_pair_counts.csv")}'
    (HEADER, DELIMITER ',');
    """)

    print("Building metadata enrichment table compared with SRA background...")
    con.execute(f"""
    COPY (
        WITH comp AS (
            SELECT *
            FROM read_parquet('{path_string(outdir / "metadata_composition_pair_counts.parquet")}')
        ), totals AS (
            SELECT
                Correct_interaction,
                host_threshold,
                pathogen_threshold,
                field_name,
                SUM(n_shared_accessions)::BIGINT AS total_shared_accessions_for_field
            FROM comp
            GROUP BY 1, 2, 3, 4
        ), labelled AS (
            SELECT
                c.Correct_interaction,
                c.host_threshold,
                c.pathogen_threshold,
                c.field_name,
                c.field_value,
                c.n_unique_pairs,
                c.n_shared_accessions,
                t.total_shared_accessions_for_field,
                CASE
                    WHEN t.total_shared_accessions_for_field > 0
                    THEN c.n_shared_accessions * 1.0 / t.total_shared_accessions_for_field
                    ELSE NULL
                END AS fraction_shared_accessions
            FROM comp c
            JOIN totals t
              ON c.Correct_interaction = t.Correct_interaction
             AND c.host_threshold = t.host_threshold
             AND c.pathogen_threshold = t.pathogen_threshold
             AND c.field_name = t.field_name
        ), bg AS (
            SELECT *
            FROM read_parquet('{path_string(outdir / "background_metadata_counts.parquet")}')
        )
        SELECT
            l.Correct_interaction,
            l.host_threshold,
            l.pathogen_threshold,
            l.field_name,
            l.field_value,
            l.n_unique_pairs,
            l.n_shared_accessions,
            l.total_shared_accessions_for_field,
            l.fraction_shared_accessions,
            bg.n_sra_accessions,
            bg.total_sra_accessions,
            bg.fraction_sra_accessions,
            CASE
                WHEN bg.fraction_sra_accessions > 0
                THEN l.fraction_shared_accessions / bg.fraction_sra_accessions
                ELSE NULL
            END AS enrichment_vs_sra_background,
            CASE
                WHEN bg.fraction_sra_accessions > 0 AND l.fraction_shared_accessions > 0
                THEN LN(l.fraction_shared_accessions / bg.fraction_sra_accessions) / LN(2.0)
                ELSE NULL
            END AS log2_enrichment_vs_sra_background
        FROM labelled l
        LEFT JOIN bg
          ON l.field_name = bg.field_name
         AND l.field_value = bg.field_value
        ORDER BY Correct_interaction DESC, field_name, enrichment_vs_sra_background DESC NULLS LAST
    )
    TO '{path_string(outdir / "metadata_enrichment_by_label.parquet")}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT *
        FROM read_parquet('{path_string(outdir / "metadata_enrichment_by_label.parquet")}')
    )
    TO '{path_string(outdir / "metadata_enrichment_by_label.csv")}'
    (HEADER, DELIMITER ',');
    """)

    # -----------------------------------------------------
    # 12. Diagnostics
    # -----------------------------------------------------
    print("\nDiagnostics:")

    diagnostics = con.execute(f"""
    WITH pair_counts AS (
        SELECT
            Correct_interaction,
            COUNT(DISTINCT CAST(host_taxid AS VARCHAR) || '__' || CAST(pathogen_taxid AS VARCHAR)) AS n_pairs_with_shared_acc,
            COUNT(*) AS n_shared_acc_rows
        FROM read_parquet('{path_string(outdir / "shared_acc_selected.parquet")}')
        GROUP BY 1
    ), all_pairs AS (
        SELECT
            Correct_interaction,
            COUNT(*) AS n_pairs_input
        FROM selected_pairs_df
        GROUP BY 1
    )
    SELECT
        a.Correct_interaction,
        a.n_pairs_input,
        COALESCE(p.n_pairs_with_shared_acc, 0) AS n_pairs_with_shared_acc,
        COALESCE(p.n_shared_acc_rows, 0) AS n_shared_acc_rows
    FROM all_pairs a
    LEFT JOIN pair_counts p USING (Correct_interaction)
    ORDER BY a.Correct_interaction;
    """).df()

    print("\nPairs with at least one shared accession by label:")
    print(diagnostics)
    write_df_csv(diagnostics, outdir / "diagnostics_pair_recovery_by_label.csv")

    meta_counts = con.execute(f"""
        SELECT
            Correct_interaction,
            field_name,
            COUNT(DISTINCT field_value) AS n_distinct_values,
            SUM(n_unique_pairs) AS sum_pair_field_values,
            SUM(n_shared_accessions) AS sum_shared_accessions
        FROM read_parquet('{path_string(outdir / "metadata_composition_pair_counts.parquet")}')
        GROUP BY 1, 2
        ORDER BY 1, 2
    """).df()

    print("\nMetadata field summary by label:")
    print(meta_counts)
    write_df_csv(meta_counts, outdir / "diagnostics_metadata_field_summary_by_label.csv")

    con.close()

    print("\nDone.")
    print(f"Outputs written to: {outdir.resolve()}")
    print("\nMain outputs:")
    print("  shared_acc_selected_counts.csv")
    print("  shared_meta_summary_long.csv")
    print("  metadata_composition_pair_counts.csv")
    print("  background_metadata_counts.csv")
    print("  metadata_enrichment_by_label.csv")


if __name__ == "__main__":
    main()
