from pathlib import Path
import math
import os
import shutil

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from shuffle_null_config import (
    TABLES,
    PERM_HISTS,
    PERM_PAIRS,
    TEMP,
    PERM_THREADS,
    PERM_MEMORY_LIMIT,
    SEED,
    BIN_WIDTH,
    N_TOTAL_OVERRIDE,
    N_PERM_TOTAL,
    N_ARRAY_TASKS,
    HOST_TAXA_PER_CHUNK,
    PATHOGEN_TAXA_PER_CHUNK,
    SAVE_PERM_NONZERO_PAIRS,
    KEEP_TEMP_DETECT_DIRS,
    OVERWRITE_EXISTING_PERMS,
)


ACC_FILE = TABLES / "acc_universe.parquet"
HOST_COUNTS_FILE = TABLES / "host_counts.parquet"
PATHOGEN_COUNTS_FILE = TABLES / "pathogen_counts.parquet"
ZERO_BASE_FILE = TABLES / "zero_bin_baseline.parquet"

PERM_HISTS.mkdir(parents=True, exist_ok=True)
PERM_PAIRS.mkdir(parents=True, exist_ok=True)
TEMP.mkdir(parents=True, exist_ok=True)


def p(path_obj: Path) -> str:
    return path_obj.as_posix()


def task_id_from_env() -> int:
    """Return 0-based Slurm array task ID. Defaults to 0 for local testing."""
    return int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))


def perms_for_task(task_id: int, n_perm_total: int, n_tasks: int) -> list[int]:
    """Split 1..N_PERM_TOTAL across 0..N_ARRAY_TASKS-1."""
    perms_per_task = math.ceil(n_perm_total / n_tasks)
    start = task_id * perms_per_task + 1
    end = min((task_id + 1) * perms_per_task, n_perm_total)
    if start > n_perm_total:
        return []
    return list(range(start, end + 1))


def write_chunked_detects_from_counts(
    counts_df: pd.DataFrame,
    tax_col: str,
    count_col: str,
    n_total_acc: int,
    rng: np.random.Generator,
    out_dir: Path,
    taxa_per_chunk: int,
) -> None:
    """
    For each taxon, keep its detection count fixed but randomly choose which
    acc_idx rows contain that taxon in this permutation.
    """
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for start in range(0, len(counts_df), taxa_per_chunk):
        sub = counts_df.iloc[start:start + taxa_per_chunk]

        acc_parts = []
        tax_parts = []

        for taxid, k in zip(sub[tax_col].to_numpy(), sub[count_col].to_numpy()):
            k = int(k)
            if k <= 0:
                continue
            if k > n_total_acc:
                raise ValueError(
                    f"Taxon {taxid} has count k={k:,}, which is larger than N={n_total_acc:,}."
                )

            acc_idx_perm = rng.choice(n_total_acc, size=k, replace=False)
            acc_parts.append(acc_idx_perm.astype(np.int64))
            tax_parts.append(np.full(k, taxid, dtype=np.int64))

        if not acc_parts:
            continue

        table = pa.table({
            "acc_idx": np.concatenate(acc_parts),
            tax_col: np.concatenate(tax_parts),
        })
        pq.write_table(table, out_dir / f"chunk_{start:06d}.parquet")


def binned_log_odds_sql(a_expr: str, n_host_expr: str, n_pathogen_expr: str, n_total: int) -> str:
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
    for req in [ACC_FILE, HOST_COUNTS_FILE, PATHOGEN_COUNTS_FILE, ZERO_BASE_FILE]:
        if not req.exists():
            raise FileNotFoundError(
                f"Missing required file: {req}\nRun 01_prepare_observed_and_baseline.py first."
            )

    task_id = task_id_from_env()
    assigned_perms = perms_for_task(task_id, N_PERM_TOTAL, N_ARRAY_TASKS)

    print("\n========== Host-virus shuffled-null permutation step ==========")
    print(f"SLURM_ARRAY_TASK_ID = {task_id}")
    print(f"Assigned permutations: {assigned_perms}")

    if not assigned_perms:
        print("No permutations assigned to this task. Exiting.")
        return

    print("\nLoading count tables...")
    host_counts_df = pd.read_parquet(HOST_COUNTS_FILE).sort_values("host_taxid").reset_index(drop=True)
    pathogen_counts_df = pd.read_parquet(PATHOGEN_COUNTS_FILE).sort_values("pathogen_taxid").reset_index(drop=True)

    # Zero-count taxa remain represented in zero_bin_baseline. We skip them here
    # because there is nothing to randomly assign for those taxa.
    host_counts_nonzero = host_counts_df.loc[host_counts_df["n_host_acc"] > 0].reset_index(drop=True)
    pathogen_counts_nonzero = pathogen_counts_df.loc[pathogen_counts_df["n_pathogen_acc"] > 0].reset_index(drop=True)

    n_total_acc_computed = len(pd.read_parquet(ACC_FILE, columns=["acc_idx"]))
    n_total_acc = N_TOTAL_OVERRIDE if N_TOTAL_OVERRIDE is not None else n_total_acc_computed
    print(f"Total accession universe used in score/permutations: {n_total_acc:,}")
    print(f"Hosts with nonzero detections:     {len(host_counts_nonzero):,}")
    print(f"Pathogens with nonzero detections: {len(pathogen_counts_nonzero):,}")

    zero_base_total = int(
        duckdb.sql(f"SELECT SUM(n_pairs) FROM read_parquet('{p(ZERO_BASE_FILE)}')").fetchone()[0]
    )

    con = duckdb.connect()
    con.execute(f"PRAGMA threads={PERM_THREADS};")
    con.execute(f"PRAGMA memory_limit='{PERM_MEMORY_LIMIT}';")
    con.execute("PRAGMA preserve_insertion_order=false;")

    for perm in assigned_perms:
        perm_hist_file = PERM_HISTS / f"perm_{perm:03d}_hist_FIXED.parquet"
        perm_pairs_file = PERM_PAIRS / f"perm_{perm:03d}_nonzero_pairs.parquet"

        if perm_hist_file.exists() and not OVERWRITE_EXISTING_PERMS:
            print(f"Permutation {perm}: histogram already exists, skipping.")
            continue

        print(f"\nRunning permutation {perm}...")
        rng = np.random.default_rng(SEED + perm)

        host_temp_dir = TEMP / f"host_perm_detects_{perm:03d}"
        pathogen_temp_dir = TEMP / f"pathogen_perm_detects_{perm:03d}"

        write_chunked_detects_from_counts(
            counts_df=host_counts_nonzero,
            tax_col="host_taxid",
            count_col="n_host_acc",
            n_total_acc=n_total_acc,
            rng=rng,
            out_dir=host_temp_dir,
            taxa_per_chunk=HOST_TAXA_PER_CHUNK,
        )

        write_chunked_detects_from_counts(
            counts_df=pathogen_counts_nonzero,
            tax_col="pathogen_taxid",
            count_col="n_pathogen_acc",
            n_total_acc=n_total_acc,
            rng=rng,
            out_dir=pathogen_temp_dir,
            taxa_per_chunk=PATHOGEN_TAXA_PER_CHUNK,
        )

        host_glob = host_temp_dir / "*.parquet"
        pathogen_glob = pathogen_temp_dir / "*.parquet"

        perm_nonzero_query = f"""
            SELECT
                h.host_taxid,
                psrc.pathogen_taxid,
                COUNT(*)::BIGINT AS a_perm
            FROM read_parquet('{p(host_glob)}') h
            JOIN read_parquet('{p(pathogen_glob)}') psrc
              USING (acc_idx)
            GROUP BY 1, 2
        """

        if SAVE_PERM_NONZERO_PAIRS:
            con.execute(f"""
            COPY (
                {perm_nonzero_query}
            )
            TO '{p(perm_pairs_file)}'
            (FORMAT PARQUET);
            """)
            perm_nonzero_source = f"read_parquet('{p(perm_pairs_file)}')"
        else:
            perm_nonzero_source = f"({perm_nonzero_query})"

        zero_bin_expr = binned_log_odds_sql("0.0", "n_host_acc", "n_pathogen_acc", n_total_acc)
        actual_bin_expr = binned_log_odds_sql("a_perm", "n_host_acc", "n_pathogen_acc", n_total_acc)

        con.execute(f"""
        COPY (
            WITH perm_nonzero_pairs AS (
                SELECT * FROM {perm_nonzero_source}
            ),
            x AS (
                SELECT
                    p.host_taxid,
                    p.pathogen_taxid,
                    p.a_perm,
                    h.n_host_acc,
                    pc.n_pathogen_acc
                FROM perm_nonzero_pairs p
                JOIN read_parquet('{p(HOST_COUNTS_FILE)}') h USING (host_taxid)
                JOIN read_parquet('{p(PATHOGEN_COUNTS_FILE)}') pc USING (pathogen_taxid)
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
                SELECT bin_left FROM read_parquet('{p(ZERO_BASE_FILE)}')
                UNION
                SELECT bin_left FROM delta_sum
            )
            SELECT
                {perm}::INTEGER AS perm,
                a.bin_left,
                (COALESCE(z.n_pairs, 0) + COALESCE(d.delta, 0))::BIGINT AS n_pairs
            FROM all_bins a
            LEFT JOIN read_parquet('{p(ZERO_BASE_FILE)}') z USING (bin_left)
            LEFT JOIN delta_sum d USING (bin_left)
            ORDER BY 2
        )
        TO '{p(perm_hist_file)}'
        (FORMAT PARQUET);
        """)

        perm_total = int(con.execute(
            f"SELECT SUM(n_pairs) FROM read_parquet('{p(perm_hist_file)}')"
        ).fetchone()[0])
        if perm_total != zero_base_total:
            raise RuntimeError(
                f"Permutation {perm} histogram total {perm_total:,} does not equal full pair universe {zero_base_total:,}."
            )

        if not KEEP_TEMP_DETECT_DIRS:
            shutil.rmtree(host_temp_dir, ignore_errors=True)
            shutil.rmtree(pathogen_temp_dir, ignore_errors=True)

        print(f"Finished permutation {perm}: {perm_hist_file.name}")

    con.close()
    print("\nTask complete.")


if __name__ == "__main__":
    main()
