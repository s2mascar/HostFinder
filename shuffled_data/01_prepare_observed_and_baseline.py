from pathlib import Path
import glob
import shutil

import duckdb

from shuffle_null_config import (
    HOST_GLOB,
    PATHOGEN_GLOB,
    ORIGINAL_PAIR_PARQUET,
    BASE,
    TABLES,
    OVERWRITE_OUTDIR,
    ACC_COL,
    TAX_COL,
    ABUND_COL,
    PAIR_HOST_COL,
    PAIR_PATHOGEN_COL,
    HOST_THRESHOLD,
    PATHOGEN_THRESHOLD,
    N_TOTAL_OVERRIDE,
    BIN_WIDTH,
    PREPARE_THREADS,
    PREPARE_MEMORY_LIMIT,
)


def p(path_obj: Path) -> str:
    return path_obj.as_posix()


def resolve_parquet_files(root_or_glob: str) -> list[str]:
    """
    Accepts:
      - some_dir/**
      - some_dir/**/*.parquet
      - some_dir
      - explicit glob patterns
    and returns all parquet files underneath.
    """
    pattern = root_or_glob.strip()

    if pattern.endswith("/**"):
        pattern1 = pattern.rstrip("/") + "/*.parquet"
        files = sorted(f for f in glob.glob(pattern1, recursive=True) if f.endswith(".parquet"))
        if files:
            return files

        pattern2 = pattern.rstrip("/") + "/**/*.parquet"
        files = sorted(f for f in glob.glob(pattern2, recursive=True) if f.endswith(".parquet"))
        if files:
            return files

    if ".parquet" in pattern or "*" in pattern:
        files = sorted(f for f in glob.glob(pattern, recursive=True) if f.endswith(".parquet"))
        if files:
            return files

    path_obj = Path(root_or_glob)
    if path_obj.exists() and path_obj.is_dir():
        files = sorted(str(x) for x in path_obj.rglob("*.parquet"))
        if files:
            return files

    fallback = str(Path(root_or_glob) / "**" / "*.parquet")
    files = sorted(f for f in glob.glob(fallback, recursive=True) if f.endswith(".parquet"))
    if not files:
        raise FileNotFoundError(f"No parquet files found for: {root_or_glob}")
    return files


def binned_log_odds_sql(a_expr: str, n_host_expr: str, n_pathogen_expr: str, n_total: int) -> str:
    """
    Returns a DuckDB SQL expression that calculates log2 observed/expected and
    bins it using BIN_WIDTH. A pseudocount of 1/N is used, matching the earlier
    thesis log-odds calculation.
    """
    return f"""
        ROUND(
            FLOOR(
                (
                    LN(
                        (
                            (({a_expr}) * 1.0 / {n_total}) + (1.0 / {n_total})
                        ) /
                        (
                            ((({n_host_expr}) * 1.0 / {n_total}) * (({n_pathogen_expr}) * 1.0 / {n_total})) + (1.0 / {n_total})
                        )
                    ) / LN(2.0)
                ) / {BIN_WIDTH}
            ) * {BIN_WIDTH},
            10
        )
    """


def main() -> None:
    if BASE.exists() and OVERWRITE_OUTDIR:
        print(f"Deleting old output directory: {BASE}")
        shutil.rmtree(BASE)

    TABLES.mkdir(parents=True, exist_ok=True)

    print("\n========== Host-virus shuffled-null prepare step ==========")
    print(f"Host threshold:     {HOST_THRESHOLD:g}")
    print(f"Pathogen threshold: {PATHOGEN_THRESHOLD:g}")
    print(f"Output directory:   {BASE}")

    print("\nResolving parquet files...")
    host_files = resolve_parquet_files(HOST_GLOB)
    pathogen_files = resolve_parquet_files(PATHOGEN_GLOB)

    print(f"  Host parquet files found:     {len(host_files):,}")
    print(f"  Pathogen parquet files found: {len(pathogen_files):,}")

    original_pair_path = Path(ORIGINAL_PAIR_PARQUET)
    if not original_pair_path.exists():
        raise FileNotFoundError(f"Missing original pair parquet: {ORIGINAL_PAIR_PARQUET}")

    con = duckdb.connect()
    con.execute(f"PRAGMA threads={PREPARE_THREADS};")
    con.execute(f"PRAGMA memory_limit='{PREPARE_MEMORY_LIMIT}';")
    con.execute("PRAGMA preserve_insertion_order=false;")

    con.read_parquet(host_files).create_view("host_src")
    con.read_parquet(pathogen_files).create_view("pathogen_src")
    con.read_parquet(str(original_pair_path)).create_view("pair_src")

    acc_universe_path = TABLES / "acc_universe.parquet"
    host_universe_path = TABLES / "host_universe.parquet"
    pathogen_universe_path = TABLES / "pathogen_universe.parquet"
    host_detects_path = TABLES / "host_detects_thresholded.parquet"
    pathogen_detects_path = TABLES / "pathogen_detects_thresholded.parquet"
    host_counts_path = TABLES / "host_counts.parquet"
    pathogen_counts_path = TABLES / "pathogen_counts.parquet"
    observed_host_detects_path = TABLES / "observed_host_detects.parquet"
    observed_pathogen_detects_path = TABLES / "observed_pathogen_detects.parquet"
    observed_nonzero_pairs_path = TABLES / "observed_nonzero_pairs.parquet"
    observed_scored_nonzero_pairs_path = TABLES / "observed_scored_nonzero_pairs.parquet"
    zero_bin_baseline_path = TABLES / "zero_bin_baseline.parquet"
    observed_hist_fixed_path = TABLES / "observed_hist_FIXED.parquet"

    print("\n[1/7] Building accession universe...")
    con.execute(f"""
    COPY (
        SELECT
            acc,
            ROW_NUMBER() OVER (ORDER BY acc) - 1 AS acc_idx
        FROM (
            SELECT DISTINCT {ACC_COL} AS acc
            FROM host_src
            WHERE {ACC_COL} IS NOT NULL

            UNION

            SELECT DISTINCT {ACC_COL} AS acc
            FROM pathogen_src
            WHERE {ACC_COL} IS NOT NULL
        )
    )
    TO '{p(acc_universe_path)}'
    (FORMAT PARQUET);
    """)

    n_total_acc_computed = con.execute(
        f"SELECT COUNT(*) FROM read_parquet('{p(acc_universe_path)}')"
    ).fetchone()[0]
    n_total_acc = N_TOTAL_OVERRIDE if N_TOTAL_OVERRIDE is not None else n_total_acc_computed

    print(f"    Total accession universe from raw STAT files: {n_total_acc_computed:,}")
    if N_TOTAL_OVERRIDE is not None:
        print(f"    Using N_TOTAL_OVERRIDE for log-odds/permutations: {n_total_acc:,}")

    print("\n[2/7] Building host and pathogen universe from original pair parquet...")
    con.execute(f"""
    COPY (
        SELECT DISTINCT CAST({PAIR_HOST_COL} AS BIGINT) AS host_taxid
        FROM pair_src
        WHERE {PAIR_HOST_COL} IS NOT NULL
        ORDER BY 1
    )
    TO '{p(host_universe_path)}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT DISTINCT CAST({PAIR_PATHOGEN_COL} AS BIGINT) AS pathogen_taxid
        FROM pair_src
        WHERE {PAIR_PATHOGEN_COL} IS NOT NULL
        ORDER BY 1
    )
    TO '{p(pathogen_universe_path)}'
    (FORMAT PARQUET);
    """)

    n_hosts_universe = con.execute(
        f"SELECT COUNT(*) FROM read_parquet('{p(host_universe_path)}')"
    ).fetchone()[0]
    n_pathogens_universe = con.execute(
        f"SELECT COUNT(*) FROM read_parquet('{p(pathogen_universe_path)}')"
    ).fetchone()[0]
    full_pair_universe = int(n_hosts_universe) * int(n_pathogens_universe)

    print(f"    Host universe size:     {n_hosts_universe:,}")
    print(f"    Pathogen universe size: {n_pathogens_universe:,}")
    print(f"    Full pair universe:     {full_pair_universe:,}")

    print("\n[3/7] Extracting above-threshold host and pathogen detections...")
    con.execute(f"""
    COPY (
        WITH host_max AS (
            SELECT
                h.{ACC_COL} AS acc,
                CAST(h.{TAX_COL} AS BIGINT) AS host_taxid,
                MAX(CAST(h.{ABUND_COL} AS DOUBLE)) AS max_abundance
            FROM host_src h
            JOIN read_parquet('{p(host_universe_path)}') u
              ON CAST(h.{TAX_COL} AS BIGINT) = u.host_taxid
            WHERE h.{ACC_COL} IS NOT NULL
            GROUP BY 1, 2
        )
        SELECT DISTINCT
            u.acc_idx,
            h.host_taxid
        FROM host_max h
        JOIN read_parquet('{p(acc_universe_path)}') u
          ON h.acc = u.acc
        WHERE h.max_abundance >= {HOST_THRESHOLD}
    )
    TO '{p(host_detects_path)}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        WITH pathogen_max AS (
            SELECT
                psrc.{ACC_COL} AS acc,
                CAST(psrc.{TAX_COL} AS BIGINT) AS pathogen_taxid,
                MAX(CAST(psrc.{ABUND_COL} AS DOUBLE)) AS max_abundance
            FROM pathogen_src psrc
            JOIN read_parquet('{p(pathogen_universe_path)}') u
              ON CAST(psrc.{TAX_COL} AS BIGINT) = u.pathogen_taxid
            WHERE psrc.{ACC_COL} IS NOT NULL
            GROUP BY 1, 2
        )
        SELECT DISTINCT
            u.acc_idx,
            psrc.pathogen_taxid
        FROM pathogen_max psrc
        JOIN read_parquet('{p(acc_universe_path)}') u
          ON psrc.acc = u.acc
        WHERE psrc.max_abundance >= {PATHOGEN_THRESHOLD}
    )
    TO '{p(pathogen_detects_path)}'
    (FORMAT PARQUET);
    """)

    n_host_detects = con.execute(
        f"SELECT COUNT(*) FROM read_parquet('{p(host_detects_path)}')"
    ).fetchone()[0]
    n_pathogen_detects = con.execute(
        f"SELECT COUNT(*) FROM read_parquet('{p(pathogen_detects_path)}')"
    ).fetchone()[0]
    print(f"    Host detections:     {n_host_detects:,}")
    print(f"    Pathogen detections: {n_pathogen_detects:,}")

    print("\n[4/7] Building count tables over the full original universe...")
    con.execute(f"""
    COPY (
        WITH signal_counts AS (
            SELECT
                host_taxid,
                COUNT(*)::BIGINT AS n_host_acc
            FROM read_parquet('{p(host_detects_path)}')
            GROUP BY 1
        )
        SELECT
            u.host_taxid,
            COALESCE(s.n_host_acc, 0)::BIGINT AS n_host_acc
        FROM read_parquet('{p(host_universe_path)}') u
        LEFT JOIN signal_counts s USING (host_taxid)
        ORDER BY 1
    )
    TO '{p(host_counts_path)}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        WITH signal_counts AS (
            SELECT
                pathogen_taxid,
                COUNT(*)::BIGINT AS n_pathogen_acc
            FROM read_parquet('{p(pathogen_detects_path)}')
            GROUP BY 1
        )
        SELECT
            u.pathogen_taxid,
            COALESCE(s.n_pathogen_acc, 0)::BIGINT AS n_pathogen_acc
        FROM read_parquet('{p(pathogen_universe_path)}') u
        LEFT JOIN signal_counts s USING (pathogen_taxid)
        ORDER BY 1
    )
    TO '{p(pathogen_counts_path)}'
    (FORMAT PARQUET);
    """)

    print("\n[5/7] Computing observed overlaps and observed pair-level log-odds...")
    con.execute(f"""
    COPY (
        SELECT DISTINCT acc_idx, host_taxid
        FROM read_parquet('{p(host_detects_path)}')
    )
    TO '{p(observed_host_detects_path)}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT DISTINCT acc_idx, pathogen_taxid
        FROM read_parquet('{p(pathogen_detects_path)}')
    )
    TO '{p(observed_pathogen_detects_path)}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT
            h.host_taxid,
            psrc.pathogen_taxid,
            COUNT(*)::BIGINT AS a_obs
        FROM read_parquet('{p(observed_host_detects_path)}') h
        JOIN read_parquet('{p(observed_pathogen_detects_path)}') psrc
          USING (acc_idx)
        GROUP BY 1, 2
    )
    TO '{p(observed_nonzero_pairs_path)}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT
            o.host_taxid,
            o.pathogen_taxid,
            o.a_obs,
            h.n_host_acc,
            psrc.n_pathogen_acc,
            LN(
                (
                    (o.a_obs * 1.0 / {n_total_acc}) + (1.0 / {n_total_acc})
                ) /
                (
                    ((h.n_host_acc * 1.0 / {n_total_acc}) * (psrc.n_pathogen_acc * 1.0 / {n_total_acc})) + (1.0 / {n_total_acc})
                )
            ) / LN(2.0) AS thesis_log_odds
        FROM read_parquet('{p(observed_nonzero_pairs_path)}') o
        JOIN read_parquet('{p(host_counts_path)}') h USING (host_taxid)
        JOIN read_parquet('{p(pathogen_counts_path)}') psrc USING (pathogen_taxid)
    )
    TO '{p(observed_scored_nonzero_pairs_path)}'
    (FORMAT PARQUET);
    """)

    n_obs_pairs = con.execute(
        f"SELECT COUNT(*) FROM read_parquet('{p(observed_nonzero_pairs_path)}')"
    ).fetchone()[0]
    print(f"    Observed nonzero host-pathogen pairs: {n_obs_pairs:,}")

    print("\n[6/7] Building zero-overlap baseline over the full host x pathogen universe...")
    zero_bin_expr = binned_log_odds_sql("0.0", "hf.n_host_acc", "pf.n_pathogen_acc", n_total_acc)
    con.execute(f"""
    COPY (
        WITH host_count_freq AS (
            SELECT
                n_host_acc,
                COUNT(*)::BIGINT AS n_host_taxa
            FROM read_parquet('{p(host_counts_path)}')
            GROUP BY 1
        ),
        pathogen_count_freq AS (
            SELECT
                n_pathogen_acc,
                COUNT(*)::BIGINT AS n_pathogen_taxa
            FROM read_parquet('{p(pathogen_counts_path)}')
            GROUP BY 1
        ),
        pair_bins AS (
            SELECT
                {zero_bin_expr} AS bin_left,
                (hf.n_host_taxa * pf.n_pathogen_taxa)::BIGINT AS n_pairs
            FROM host_count_freq hf
            CROSS JOIN pathogen_count_freq pf
        )
        SELECT
            bin_left,
            SUM(n_pairs)::BIGINT AS n_pairs
        FROM pair_bins
        GROUP BY 1
        ORDER BY 1
    )
    TO '{p(zero_bin_baseline_path)}'
    (FORMAT PARQUET);
    """)

    print("\n[7/7] Building fixed observed histogram over the full host x pathogen universe...")
    zero_bin_expr = binned_log_odds_sql("0.0", "n_host_acc", "n_pathogen_acc", n_total_acc)
    actual_bin_expr = binned_log_odds_sql("a_obs", "n_host_acc", "n_pathogen_acc", n_total_acc)
    con.execute(f"""
    COPY (
        WITH x AS (
            SELECT
                o.host_taxid,
                o.pathogen_taxid,
                o.a_obs,
                h.n_host_acc,
                psrc.n_pathogen_acc
            FROM read_parquet('{p(observed_nonzero_pairs_path)}') o
            JOIN read_parquet('{p(host_counts_path)}') h USING (host_taxid)
            JOIN read_parquet('{p(pathogen_counts_path)}') psrc USING (pathogen_taxid)
        ),
        scored AS (
            SELECT
                {zero_bin_expr} AS zero_bin,
                {actual_bin_expr} AS actual_bin
            FROM x
        ),
        remove_zero AS (
            SELECT
                zero_bin AS bin_left,
                -COUNT(*)::BIGINT AS delta
            FROM scored
            GROUP BY 1
        ),
        add_actual AS (
            SELECT
                actual_bin AS bin_left,
                COUNT(*)::BIGINT AS delta
            FROM scored
            GROUP BY 1
        ),
        delta_sum AS (
            SELECT
                bin_left,
                SUM(delta)::BIGINT AS delta
            FROM (
                SELECT * FROM remove_zero
                UNION ALL
                SELECT * FROM add_actual
            )
            GROUP BY 1
        ),
        all_bins AS (
            SELECT bin_left FROM read_parquet('{p(zero_bin_baseline_path)}')
            UNION
            SELECT bin_left FROM delta_sum
        )
        SELECT
            a.bin_left,
            (COALESCE(z.n_pairs, 0) + COALESCE(d.delta, 0))::BIGINT AS n_pairs
        FROM all_bins a
        LEFT JOIN read_parquet('{p(zero_bin_baseline_path)}') z USING (bin_left)
        LEFT JOIN delta_sum d USING (bin_left)
        ORDER BY 1
    )
    TO '{p(observed_hist_fixed_path)}'
    (FORMAT PARQUET);
    """)

    observed_total_pairs = con.execute(
        f"SELECT SUM(n_pairs) FROM read_parquet('{p(observed_hist_fixed_path)}')"
    ).fetchone()[0]
    if int(observed_total_pairs) != full_pair_universe:
        raise RuntimeError(
            f"Observed histogram total {observed_total_pairs:,} does not equal full pair universe {full_pair_universe:,}."
        )

    max_obs_bin = con.execute(
        f"SELECT MAX(bin_left) FROM read_parquet('{p(observed_hist_fixed_path)}')"
    ).fetchone()[0]

    run_info_path = TABLES / "run_info.txt"
    run_info_path.write_text(
        "Host-virus shuffled-null run\n"
        f"HOST_THRESHOLD={HOST_THRESHOLD:g}\n"
        f"PATHOGEN_THRESHOLD={PATHOGEN_THRESHOLD:g}\n"
        f"BIN_WIDTH={BIN_WIDTH:g}\n"
        f"N_TOTAL_ACC={n_total_acc}\n"
        f"N_HOSTS={n_hosts_universe}\n"
        f"N_PATHOGENS={n_pathogens_universe}\n"
        f"FULL_PAIR_UNIVERSE={full_pair_universe}\n"
        f"OBSERVED_NONZERO_PAIRS={n_obs_pairs}\n"
        f"MAX_OBSERVED_BIN={max_obs_bin}\n"
    )

    print("\nDone.")
    print(f"    Total accessions used in score: {n_total_acc:,}")
    print(f"    Host universe size:             {n_hosts_universe:,}")
    print(f"    Pathogen universe size:         {n_pathogens_universe:,}")
    print(f"    Full pair universe:             {full_pair_universe:,}")
    print(f"    Observed histogram total:       {int(observed_total_pairs):,}")
    print(f"    Max observed histogram bin:     {max_obs_bin}")
    print(f"    Outputs written to:             {BASE}")

    con.close()


if __name__ == "__main__":
    main()
