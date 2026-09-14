from pathlib import Path
import glob

import duckdb

from shuffle_null_config import (
    TABLES,
    PERM_HISTS,
    N_PERM_TOTAL,
    COMBINE_THREADS,
    COMBINE_MEMORY_LIMIT,
)


ALL_PERM_FIXED = TABLES / "all_perm_hist_FIXED.parquet"
NULL_SUMMARY_FIXED = TABLES / "null_hist_summary_FIXED.parquet"
NULL_SUMMARY_QUANTILES_FIXED = TABLES / "null_hist_summary_quantiles_FIXED.parquet"


def p(path_obj: Path) -> str:
    return path_obj.as_posix()


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)

    hist_files = sorted(glob.glob(str(PERM_HISTS / "perm_*_hist_FIXED.parquet")))
    if not hist_files:
        raise FileNotFoundError(
            f"No perm_*_hist_FIXED.parquet files found in {PERM_HISTS}. Run the permutation array first."
        )

    print("\n========== Combining permutation histograms ==========")
    print(f"Found {len(hist_files)} permutation histograms.")
    if len(hist_files) != N_PERM_TOTAL:
        print(f"WARNING: expected {N_PERM_TOTAL} permutation histograms, found {len(hist_files)}.")

    con = duckdb.connect()
    con.execute(f"PRAGMA threads={COMBINE_THREADS};")
    con.execute(f"PRAGMA memory_limit='{COMBINE_MEMORY_LIMIT}';")
    con.execute("PRAGMA preserve_insertion_order=false;")

    glob_path = p(PERM_HISTS / "perm_*_hist_FIXED.parquet")

    con.execute(f"""
    COPY (
        SELECT
            CAST(perm AS INTEGER) AS perm,
            CAST(bin_left AS DOUBLE) AS bin_left,
            CAST(n_pairs AS BIGINT) AS n_pairs
        FROM read_parquet('{glob_path}')
        ORDER BY perm, bin_left
    )
    TO '{p(ALL_PERM_FIXED)}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT
            bin_left,
            AVG(n_pairs) AS mean_n_pairs,
            STDDEV_SAMP(n_pairs) AS sd_n_pairs,
            MIN(n_pairs) AS min_n_pairs,
            MAX(n_pairs) AS max_n_pairs
        FROM read_parquet('{glob_path}')
        GROUP BY 1
        ORDER BY 1
    )
    TO '{p(NULL_SUMMARY_FIXED)}'
    (FORMAT PARQUET);
    """)

    con.execute(f"""
    COPY (
        SELECT
            bin_left,
            AVG(n_pairs) AS mean_n_pairs,
            STDDEV_SAMP(n_pairs) AS sd_n_pairs,
            MIN(n_pairs) AS min_n_pairs,
            QUANTILE_CONT(n_pairs, 0.05) AS q05_n_pairs,
            QUANTILE_CONT(n_pairs, 0.50) AS median_n_pairs,
            QUANTILE_CONT(n_pairs, 0.95) AS q95_n_pairs,
            MAX(n_pairs) AS max_n_pairs
        FROM read_parquet('{glob_path}')
        GROUP BY 1
        ORDER BY 1
    )
    TO '{p(NULL_SUMMARY_QUANTILES_FIXED)}'
    (FORMAT PARQUET);
    """)

    n_perm = con.execute(
        f"SELECT COUNT(DISTINCT perm) FROM read_parquet('{glob_path}')"
    ).fetchone()[0]
    con.close()

    print("Done.")
    print(f"Combined permutations: {n_perm}")
    print(f"Saved: {ALL_PERM_FIXED}")
    print(f"Saved: {NULL_SUMMARY_FIXED}")
    print(f"Saved: {NULL_SUMMARY_QUANTILES_FIXED}")


if __name__ == "__main__":
    main()
