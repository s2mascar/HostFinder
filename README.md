# HostFinder
HostFinder is a bioinformatics framework for identifying and prioritizing potential **host–microbe interactions** from large-scale public sequencing data in the NCBI Sequence Read Archive (SRA).

The framework uses taxonomic co-occurrence patterns across millions of sequencing datasets to determine whether a host and microorganism are observed together more frequently than expected by chance. HostFinder was developed as part of my MSc research in Bioinformatics at the University of Waterloo.

## Overview

Public sequencing repositories contain biological information beyond the original purpose of individual sequencing experiments. HostFinder treats the SRA as a large observational dataset and searches for recurring associations between host and microbial taxa.

The general hypothesis is:

> If two organisms repeatedly co-occur across independent sequencing datasets more often than expected from their individual prevalence, that co-occurrence may provide evidence of a biological interaction.

HostFinder processes SRA taxonomic abundance data, calculates host–microbe co-occurrence across different abundance thresholds, engineers statistical association features, and evaluates whether these signals can distinguish known biological interactions from control pairs.

## Workflow

The current HostFinder workflow consists of five major stages:

```text
NCBI SRA taxonomic profiles
            │
            ▼
      Pre-processing
            │
            ▼
   Abundance calculation
            │
            ▼
 Alpha-diversity analysis
            │
            ▼
    Feature engineering
            │
            ▼
Host–microbe co-occurrence
      and evaluation
```

The corresponding components of the repository are:

```text
HostFinder/
│
├── pre_process/
│   └── Scripts for preparing and filtering SRA-derived taxonomic data
│
├── abundance_calc/
│   └── Scripts for calculating taxonomic abundance and occurrence
│
├── alpha_diversity_calc/
│   └── Scripts for calculating within-sample taxonomic diversity
│
├── feature_eng/
│   └── Scripts for generating features used to characterize interactions
│
├── co_occurence_stats.R
│   └── Co-occurrence scoring, statistical analysis, and evaluation
│
├── LICENSE
└── README.md
```

## Co-occurrence Analysis

For each candidate host–microbe pair, HostFinder tracks the number of SRA datasets containing:

* the host
* the microorganism
* both the host and microorganism

These counts are calculated across different taxonomic abundance thresholds.

A simplified representation is:

```text
N       = total number of SRA datasets
N_h     = datasets containing the host
N_m     = datasets containing the microorganism
N_hm    = datasets containing both
```

These values are used to determine whether the observed host–microbe overlap is greater than expected from their individual prevalence.

## Association Scoring

The main analysis script, `co_occurence_stats.R`, calculates several measures of association between host and microbial taxa.

These include:

* Log-odds score
* Hypergeometric probability
* Jaccard similarity
* Sørensen–Dice coefficient
* Ochiai coefficient
* Phi coefficient
* Mutual information
* Fisher-based statistics

### Log-Odds

One of the primary HostFinder scores compares observed host–microbe co-occurrence with the co-occurrence expected from their individual frequencies.

For host \(h\) and microbe \(m\):

```text
f_h  = frequency of the host
f_m  = frequency of the microorganism
f_hm = frequency of host–microbe co-occurrence
```

The general form of the score is:

$$
\text{log-odds} =
\log_2
\left(
\frac{f_{hm} + k}
{f_h f_m + k}
\right)
$$

where \(k\) is a small smoothing constant.

A positive score indicates that the host and microorganism occur together more frequently than expected from their individual prevalence.

## Abundance Thresholds

Taxonomic detections in sequencing datasets can vary substantially in abundance. HostFinder therefore evaluates interactions across combinations of host and microbial abundance thresholds rather than relying on a single cutoff.

For each threshold combination, the pipeline recalculates:

```text
Host prevalence
Microbe prevalence
Shared prevalence
Association score
```

This makes it possible to determine how the strength and predictive performance of an interaction changes as detection criteria become more or less stringent.

## Benchmarking Known Interactions

To determine whether co-occurrence provides meaningful biological information, HostFinder compares known host–microbe interactions with negative and background interaction sets.

The current evaluation framework includes interaction categories such as:

```text
Positive Control
Negative Control
PHI-Base
Commensal
Random
```

PHI-Base provides experimentally supported pathogen–host relationships that can be used to test whether known biological associations receive stronger HostFinder scores than unrelated organism pairs.

## Model Evaluation

HostFinder evaluates the ability of co-occurrence scores to distinguish interaction classes across abundance thresholds.

The analysis includes:

* ROC curves
* ROC AUC
* threshold-specific performance
* optimal score cutoffs
* comparisons between interaction classes
* AUC heatmaps across host and microbial abundance thresholds
* score-distribution diagnostics

For example, interaction classes can be evaluated using comparisons such as:

```text
Positive Control vs Negative Control
Positive Control vs Random
Commensal vs Negative Control
Commensal vs Random
PHI-Base vs Negative Control
PHI-Base vs Random
```

This allows the pipeline to identify abundance thresholds under which co-occurrence provides the strongest signal of a biological interaction.

## Input Data

The main co-occurrence analysis currently expects two CSV files:

```text
host_pathogen_threshold_summary.csv
100_host_microbe_pairs_11_11_corrected.csv
```

### `host_pathogen_threshold_summary.csv`

Contains host–microbe occurrence counts across abundance threshold combinations.

Core fields include information equivalent to:

```text
host_tax_id
path_tax_id
host_threshold
path_threshold
host_dataset_count
pathogen_dataset_count
both_dataset_count
```

### `100_host_microbe_pairs_11_11_corrected.csv`

Contains metadata describing the benchmark interaction pairs, including organism identities and interaction labels.

These files are generated from earlier stages of the HostFinder pipeline and are not necessarily included in the repository because the underlying SRA-derived datasets can be very large.

## Running the Co-occurrence Analysis

The statistical analysis can be run from R using:

```bash
Rscript co_occurence_stats.R
```

The script expects the required input CSV files to be available in the working directory.

The R analysis currently uses packages including:

```r
dplyr
tidyr
ggplot2
pROC
tibble
```

These can be installed using:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "pROC",
  "tibble"
))
```



## License

This project is distributed under the **MIT License**. See [`LICENSE`](LICENSE) for details.
