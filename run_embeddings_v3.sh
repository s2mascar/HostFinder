#!/bin/bash
# run_jl_embeddings_true.sh
set -euo pipefail

# ---- Tunables ----
EMBEDDING_DIM=${EMBEDDING_DIM:-512}
DATA_FILE=${DATA_FILE:-"/home/smascar/scratch/STAT_2025_with_abundance_cleaned.parquet"}
MEMORY=${MEMORY:-"120GB"}     # bump if you can
THREADS=${THREADS:-48}        # bump if you can
OUTDIR=${OUTDIR:-"./jl_out"}  # where results go (not a temp dir)
SHARDS=${SHARDS:-64}          # for dataset embeddings; 64–256 is typical
PARQUET_COMPRESSION=${PARQUET_COMPRESSION:-"ZSTD"}

mkdir -p "$OUTDIR"

echo "=== Config ==="
echo "EMBEDDING_DIM=$EMBEDDING_DIM  DATA_FILE=$DATA_FILE"
echo "MEMORY=$MEMORY  THREADS=$THREADS  SHARDS=$SHARDS"
echo "OUTDIR=$OUTDIR"
echo

# A tiny helper to run embedded Python with env vars expanded
py() { python3 - "$@" <<'PYCODE'
import os, sys, math, duckdb, pyarrow as pa, pyarrow.parquet as pq, numpy as np
EMBEDDING_DIM = int(os.environ['EMBEDDING_DIM'])
DATA_FILE     = os.environ['DATA_FILE']
MEMORY        = os.environ['MEMORY']
THREADS       = int(os.environ['THREADS'])
OUTDIR        = os.environ['OUTDIR']
SHARDS        = int(os.environ['SHARDS'])
COMPRESS      = os.environ['PARQUET_COMPRESSION']

# Scale for JL (±1 / sqrt(k))
SCALE = 1.0 / math.sqrt(EMBEDDING_DIM)

# Connect once
con = duckdb.connect()
con.execute(f"SET memory_limit='{MEMORY}'")
con.execute(f"SET threads={THREADS}")

# 1) Build random sign matrix for taxa: (tax_id × k) with entries in {+SCALE, -SCALE}
print("Step 1: Create random signs per tax_id (dense JL)...", flush=True)
taxa = con.execute(f"""
  SELECT DISTINCT tax_id
  FROM read_parquet('{DATA_FILE}')
  ORDER BY tax_id
""").fetchnumpy()['tax_id']
n_taxa = len(taxa)
print(f"  found {n_taxa:,} taxa", flush=True)

rng = np.random.default_rng(42)
# int8 signs then scale to float32 on write to save memory while building
signs = rng.choice(np.int8([-1, 1]), size=(n_taxa, EMBEDDING_DIM), replace=True)
# Write as Parquet (float32) to avoid a giant join materialization in RAM
# We keep it columnar for fast joins
cols = {'tax_id': taxa}
# store in chunks to keep peak memory lower
chunk = 64
for i0 in range(0, EMBEDDING_DIM, chunk):
    i1 = min(i0+chunk, EMBEDDING_DIM)
    for j in range(i0, i1):
        cols[f"sign_{j}"] = (signs[:, j].astype('float32') * SCALE)
    table = pa.table(cols)
    path = f"{OUTDIR}/random_signs_{i0:04d}_{i1-1:04d}.parquet"
    pq.write_table(table, path, compression=COMPRESS)
    # free sign columns to keep memory predictable
    for j in range(i0, i1):
        del cols[f"sign_{j}"]
# Also write a tiny manifest of ids (handy later)
pq.write_table(pa.table({'tax_id': taxa}), f"{OUTDIR}/taxa_ids.parquet", compression=COMPRESS)
print("  ✓ random_signs_* written", flush=True)

# 2) Species embeddings: emb_species[tax_id, i] = Σ_acc total_abundance * sign_i(acc) / √k
#    Implemented as 512 grouped SUMs over tax_id. No temp dir needed.
print("Step 2: Species embeddings via JL over acc...", flush=True)

# Build the 512 SUM(...) expressions that depend on acc only (no big joins).
exprs = []
for i in range(EMBEDDING_DIM):
    # stable per-dimension sign from acc using hash of concatenated string
    # CASE WHEN (hash(acc || '#i') & 1)=0 THEN -1 ELSE +1 END
    sign_expr = f"(CASE WHEN (hash(d.acc || '#{i}') & 1) = 0 THEN -1.0 ELSE 1.0 END) * {SCALE}"
    exprs.append(f"SUM(d.total_abundance * {sign_expr}) AS emb_{i}")

species_sql = f"""
COPY (
  SELECT d.tax_id,
         {", ".join(exprs)}
  FROM read_parquet('{DATA_FILE}') d
  GROUP BY d.tax_id
  ORDER BY d.tax_id
) TO '{OUTDIR}/species_embeddings.parquet' (FORMAT PARQUET, COMPRESSION '{COMPRESS}');
"""
con.execute(species_sql)
print("  ✓ species_embeddings.parquet", flush=True)

# 3) Dataset embeddings at scale:
#    y_acc[i] = Σ_tax total_abundance * sign_i(tax_id) / √k
#    We avoid an all-acc giant GROUP BY by physically partitioning once by acc-hash.
print("Step 3a: One-pass partition by acc hash...", flush=True)
con.execute(f"""
COPY (
  SELECT *, (hash(acc) % {SHARDS})::INT AS shard
  FROM read_parquet('{DATA_FILE}')
) TO '{OUTDIR}/acc_sharded' (FORMAT PARQUET, PARTITION_BY (shard), COMPRESSION '{COMPRESS}');
""")
print("  ✓ acc_sharded/ created (partitioned by shard)", flush=True)

# 3b) For each shard, aggregate with JL via joins to random_signs_* in manageable chunks of dims.
print("Step 3b: Compute dataset embeddings per shard...", flush=True)

# Prepare schema string for output wide table
emb_cols = ", ".join([f"emb_{i} FLOAT" for i in range(EMBEDDING_DIM)])
con.execute(f"CREATE TABLE IF NOT EXISTS dataset_embeddings_all(acc VARCHAR, {emb_cols});")
con.execute("DELETE FROM dataset_embeddings_all;")  # clean slate

# Process shards sequentially; within each shard, do dim-blocks to keep query complexity sane
DIM_BLOCK = 64  # 64-dim blocks; adjust if you want fewer/larger SQL expressions
for shard in range(SHARDS):
    shard_path = f"{OUTDIR}/acc_sharded/shard={shard}"
    # sanity: skip empty shards
    try:
        _ = pq.ParquetFile(f"{shard_path}/part-0.parquet")
    except Exception:
        continue

    print(f"  - shard {shard}/{SHARDS-1}", flush=True)
    # We will build the wide 512D result by iterating blocks and LEFT JOINing on acc
    # Start with an "acc only" temp table for this shard
    con.execute("DROP TABLE IF EXISTS _acc_tmp")
    con.execute(f"""
      CREATE TEMP TABLE _acc_tmp AS
      SELECT DISTINCT acc FROM read_parquet('{shard_path}/*.parquet')
      ORDER BY acc
    """)

    # Working table for this shard
    con.execute("DROP TABLE IF EXISTS _acc_emb")
    con.execute("CREATE TEMP TABLE _acc_emb AS SELECT acc FROM _acc_tmp")

    for i0 in range(0, EMBEDDING_DIM, DIM_BLOCK):
        i1 = min(i0 + DIM_BLOCK, EMBEDDING_DIM)
        # Union the sign chunks for tax_id we wrote earlier
        sign_inputs = " , ".join([f"read_parquet('{OUTDIR}/random_signs_{j:04d}_{min(j+DIM_BLOCK,EMBEDDING_DIM)-1:04d}.parquet')"
                                  for j in range(i0, i1, DIM_BLOCK)])
        # Join shard with signs and aggregate per acc
        agg_cols = ",\n       ".join([f"SUM(d.total_abundance * s.sign_{i}) AS emb_{i}" for i in range(i0, i1)])
        con.execute(f"""
          CREATE TEMP TABLE _block AS
          SELECT d.acc, {agg_cols}
          FROM read_parquet('{shard_path}/*.parquet') d
          JOIN ({sign_inputs}) s USING (tax_id)
          GROUP BY d.acc
        """)
        # Merge this block into working table
        sel_cols = "e.acc" + "".join([f", b.emb_{i}" for i in range(i0, i1)])
        set_cols = ", ".join([f"emb_{i} = COALESCE(b.emb_{i}, emb_{i})" for i in range(i0, i1)])
        # Extend _acc_emb to have these columns (if first block), otherwise left join-update
        if i0 == 0:
            # Build initial wide table with NULLs and then COALESCE assign
            con.execute(f"ALTER TABLE _acc_emb ADD COLUMN {', ADD COLUMN '.join([f'emb_{i} FLOAT' for i in range(i0, i1)])}")
            con.execute(f"""
              UPDATE _acc_emb AS e
              SET {set_cols}
              FROM _block b
              WHERE e.acc = b.acc
            """)
        else:
            con.execute(f"ALTER TABLE _acc_emb ADD COLUMN {', ADD COLUMN '.join([f'emb_{i} FLOAT' for i in range(i0, i1)])}")
            con.execute(f"""
              UPDATE _acc_emb AS e
              SET {set_cols}
              FROM _block b
              WHERE e.acc = b.acc
            """)
        con.execute("DROP TABLE _block")

    # Append this shard’s rows into the final table
    con.execute("INSERT INTO dataset_embeddings_all SELECT * FROM _acc_emb")
    con.execute("DROP TABLE _acc_emb")
    con.execute("DROP TABLE _acc_tmp")

# 3c) Write final dataset embeddings to Parquet
con.execute(f"""
COPY dataset_embeddings_all
TO '{OUTDIR}/dataset_embeddings.parquet' (FORMAT PARQUET, COMPRESSION '{COMPRESS}');
""")
print("  ✓ dataset_embeddings.parquet", flush=True)

con.close()
print("=== Done ===", flush=True)
PYCODE
}

# -------- run --------
py
echo
echo "Outputs:"
ls -lh "$OUTDIR" | sed 's/^/  /'

