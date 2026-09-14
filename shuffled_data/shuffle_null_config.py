from pathlib import Path

# =========================================================
# Shuffled host-virus null model configuration
# =========================================================
# Edit this file if you want to change paths, thresholds, number of
# permutations, or compute resources used by DuckDB.

# Raw STAT parquet collections
HOST_GLOB = "/home/smascar/scratch/VIRUSES/STAT_with_Eukaryota_important_columns_split/**"
PATHOGEN_GLOB = "/home/smascar/scratch/VIRUSES/STAT_with_Viruses_important_columns_split/**"

# Original all-pairs parquet used only to define the taxon universe.
# It should contain host_taxid and pathogen_taxid columns.
ORIGINAL_PAIR_PARQUET = "/home/smascar/scratch/VIRUSES/host_viruses_co_occurences_fixed_thresholds_batched_v1_labeled_with_logodds__with_correct_interaction.parquet"

# Output directory for this one-threshold shuffled-null run.
# I changed the name from host_bacteria... to host_virus... so it is clear
# that this run is for viruses/pathogens, not bacteria.
BASE = Path("host_virus_abundance_shuffle_null_ht1e-03_pt1e-08")
TABLES = BASE / "tables"
PERM_HISTS = BASE / "perm_hists"
PERM_PAIRS = BASE / "perm_nonzero_pairs"
TEMP = BASE / "temp_perm_detects"

# If True, the prepare script deletes BASE before rebuilding the observed data.
# Keep True when changing thresholds so old outputs do not get mixed in.
OVERWRITE_OUTDIR = True

# Raw STAT column names
ACC_COL = "acc"
TAX_COL = "tax_id"
ABUND_COL = "total_abundance"

# Original pair parquet column names
PAIR_HOST_COL = "host_taxid"
PAIR_PATHOGEN_COL = "pathogen_taxid"

# Fixed abundance thresholds for this run
HOST_THRESHOLD = 1e-3
PATHOGEN_THRESHOLD = 1e-8

# If you want to force the exact accession universe used in an older log-odds
# calculation, set this to an integer. Otherwise the code uses the accession
# union from the host and pathogen STAT parquet files.
N_TOTAL_OVERRIDE = None

# Histogram resolution for stored log-odds histograms.
# Keep 0.1. You can make wider display bins later in the R plotting script.
BIN_WIDTH = 0.1

# DuckDB settings
PREPARE_THREADS = 8
PREPARE_MEMORY_LIMIT = "250GB"
PERM_THREADS = 8
PERM_MEMORY_LIMIT = "250GB"
COMBINE_THREADS = 8
COMBINE_MEMORY_LIMIT = "50GB"
EXPORT_THREADS = 8
EXPORT_MEMORY_LIMIT = "50GB"

# Permutation settings
SEED = 123
N_PERM_TOTAL = 100
N_ARRAY_TASKS = 25

# Chunk sizes for temporary shuffled-detection parquet chunks
HOST_TAXA_PER_CHUNK = 250
PATHOGEN_TAXA_PER_CHUNK = 500

# Optional outputs from permutation step
SAVE_PERM_NONZERO_PAIRS = False
KEEP_TEMP_DETECT_DIRS = False
OVERWRITE_EXISTING_PERMS = False
