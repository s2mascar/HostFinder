from pathlib import Path

import duckdb
import pandas as pd

from shuffle_null_config import TABLES, EXPORT_THREADS, EXPORT_MEMORY_LIMIT


OBS_HIST = TABLES / "observed_hist_FIXED.parquet"
ALL_PERM_HIST = TABLES / "all_perm_hist_FIXED.parquet"

OUT_LONG = TABLES / "observed_and_null_integer_logodds_bins_long.csv"
OUT_WIDE = TABLES / "observed_and_null_integer_logodds_bins_wide.csv"

# Export integer-width bins for quick spreadsheet inspection.
# Values below LOWER_EDGE are grouped as < LOWER_EDGE.
# Values >= UPPER_EDGE are grouped as >= UPPER_EDGE.
LOWER_EDGE = -10
UPPER_EDGE = 10


def p(path_obj: Path) -> str:
    return path_obj.as_posix()


def label_bin(bin_left: float) -> tuple[int, str, object, object]:
    if bin_left < LOWER_EDGE:
        return 0, f"< {LOWER_EDGE}", None, LOWER_EDGE
    if bin_left >= UPPER_EDGE:
        return UPPER_EDGE - LOWER_EDGE + 1, f">= {UPPER_EDGE}", UPPER_EDGE, None
    start = int(bin_left // 1)
    end = start + 1
    return start - LOWER_EDGE + 1, f"[{start},{end})", start, end


def main() -> None:
    for req in [OBS_HIST, ALL_PERM_HIST]:
        if not req.exists():
            raise FileNotFoundError(
                f"Missing required file: {req}\nRun prepare, permutation, and combine steps first."
            )

    con = duckdb.connect()
    con.execute(f"PRAGMA threads={EXPORT_THREADS};")
    con.execute(f"PRAGMA memory_limit='{EXPORT_MEMORY_LIMIT}';")

    df = con.execute(f"""
        WITH observed_src AS (
            SELECT
                'observed' AS kind,
                NULL::INTEGER AS perm,
                'observed' AS source_id,
                CAST(bin_left AS DOUBLE) AS bin_left,
                CAST(n_pairs AS BIGINT) AS n_pairs
            FROM read_parquet('{p(OBS_HIST)}')
        ),
        perm_src AS (
            SELECT
                'null' AS kind,
                CAST(perm AS INTEGER) AS perm,
                'null_' || LPAD(CAST(perm AS VARCHAR), 3, '0') AS source_id,
                CAST(bin_left AS DOUBLE) AS bin_left,
                CAST(n_pairs AS BIGINT) AS n_pairs
            FROM read_parquet('{p(ALL_PERM_HIST)}')
        )
        SELECT * FROM observed_src
        UNION ALL
        SELECT * FROM perm_src
    """).df()
    con.close()

    bin_info = df["bin_left"].apply(label_bin)
    df["row_order"] = [x[0] for x in bin_info]
    df["bin_label"] = [x[1] for x in bin_info]
    df["bin_start"] = [x[2] for x in bin_info]
    df["bin_end"] = [x[3] for x in bin_info]

    totals = (
        df.groupby(["kind", "perm", "source_id"], dropna=False)["n_pairs"]
        .sum()
        .reset_index()
        .rename(columns={"n_pairs": "total_pairs"})
    )

    long_bins = (
        df.groupby(
            ["kind", "perm", "source_id", "row_order", "bin_label", "bin_start", "bin_end"],
            dropna=False,
        )["n_pairs"]
        .sum()
        .reset_index()
        .merge(totals, on=["kind", "perm", "source_id"], how="left")
    )
    long_bins["fraction_of_total"] = long_bins["n_pairs"] / long_bins["total_pairs"]
    long_bins = long_bins.sort_values(["kind", "perm", "row_order"], na_position="first")

    long_bins.to_csv(OUT_LONG, index=False)

    wide_input = long_bins.copy()
    # pivot_table can drop rows where an index value is NA, so use a sentinel
    # for the observed row and convert it back after pivoting.
    wide_input["perm_for_pivot"] = wide_input["perm"].where(wide_input["perm"].notna(), -1).astype(int)
    wide = wide_input.pivot_table(
        index=["kind", "perm_for_pivot", "source_id", "total_pairs"],
        columns="bin_label",
        values="n_pairs",
        aggfunc="sum",
        fill_value=0,
    ).reset_index()
    wide = wide.rename(columns={"perm_for_pivot": "perm"})
    wide.loc[wide["kind"] == "observed", "perm"] = pd.NA
    wide.columns.name = None
    wide.to_csv(OUT_WIDE, index=False)

    check = long_bins.groupby(["kind", "perm", "source_id"], dropna=False)["n_pairs"].sum().reset_index()
    check = check.merge(totals, on=["kind", "perm", "source_id"], how="left")
    bad = check.loc[check["n_pairs"] != check["total_pairs"]]

    print(f"Saved long CSV: {OUT_LONG}")
    print(f"Saved wide CSV: {OUT_WIDE}")
    if bad.empty:
        print("Sanity check passed: exported bins sum to total_pairs for every observed/null source.")
    else:
        print("WARNING: some exported rows do not sum to total_pairs.")
        print(bad)


if __name__ == "__main__":
    main()
