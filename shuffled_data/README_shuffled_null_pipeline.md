# Host-virus shuffled-null pipeline

This is the cleaned shuffled-null pipeline only. It does not run metadata recovery.

## Fixed thresholds

The run is configured in `shuffle_null_config.py` with:

```python
HOST_THRESHOLD = 1e-4
PATHOGEN_THRESHOLD = 1e-8
```

The output directory is:

```text
host_virus_abundance_shuffle_null_ht1e-04_pt1e-08/
```

## What the pipeline does

1. `01_prepare_observed_and_baseline.py`
   - Reads host STAT parquet files and virus/pathogen STAT parquet files.
   - Applies the fixed host/pathogen thresholds.
   - Builds the accession universe.
   - Counts how many accessions contain each host and each pathogen.
   - Computes observed shared accessions for real host-pathogen pairs.
   - Builds the full observed host × pathogen log-odds histogram using a zero-overlap baseline.

2. `02_run_perm_array.py`
   - Runs 100 total permutations split across 25 Slurm array tasks.
   - For each taxon, keeps its number of detected accessions fixed.
   - Randomly reassigns that taxon to accession rows.
   - Recomputes shuffled host-pathogen overlaps and saves one histogram per permutation.

3. `03_combine_perm_hists.py`
   - Combines the 100 permutation histograms.
   - Writes the full null histogram table and null summary files.

4. `04_plot_observed_vs_null_bins.R`
   - Makes non-cumulative log-odds bin plots like `0–0.5`, `0.5–1`, `1–1.5`, etc.
   - Also writes cumulative tail-count outputs.

5. `05_export_hist_bin_counts.py`
   - Exports integer-width histogram bins to long and wide CSV files for quick checking.

## How to run

The easiest way is:

```bash
bash submit_shuffled_pipeline.sh
```

Or manually:

```bash
mkdir -p logs
sbatch 01_prepare_observed_and_baseline.sh
sbatch 02_run_perm_array.sh
sbatch 03_combine_and_plot.sh
```

If you run manually, wait until the prepare job finishes before starting the permutation array, and wait until the permutation array finishes before running the combine/plot step.

## Main outputs

```text
host_virus_abundance_shuffle_null_ht1e-04_pt1e-08/tables/observed_hist_FIXED.parquet
host_virus_abundance_shuffle_null_ht1e-04_pt1e-08/tables/all_perm_hist_FIXED.parquet
host_virus_abundance_shuffle_null_ht1e-04_pt1e-08/tables/null_hist_summary_FIXED.parquet
host_virus_abundance_shuffle_null_ht1e-04_pt1e-08/tables/observed_vs_null_logodds_bins.csv
host_virus_abundance_shuffle_null_ht1e-04_pt1e-08/tables/observed_vs_null_logodds_bins_log.png
host_virus_abundance_shuffle_null_ht1e-04_pt1e-08/tables/observed_vs_null_tail_counts.csv
```

## How to adjust the displayed histogram bins

Edit this part of `04_plot_observed_vs_null_bins.R`:

```r
breaks <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, Inf)
labels <- c("0–0.5", "0.5–1", "1–1.5", "1.5–2", "2–2.5", "2.5–3", "3+")
```

Changing the plot bins does not require rerunning the permutations. Only rerun:

```bash
Rscript 04_plot_observed_vs_null_bins.R
```
