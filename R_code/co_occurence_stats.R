#!/usr/bin/env Rscript

###############################################################################
# Host–microbe interaction scoring and evaluation pipeline
#
# This script:
#   1. Loads host–microbe co-occurrence counts and interaction metadata.
#   2. Computes hypergeometric P-values and log-odds scores.
#   3. Adds co-occurrence metrics (Jaccard, Sørensen–Dice, Ochiai, phi, MI, Fisher).
#   4. Samples a balanced subset of pairs per interaction type.
#   5. Evaluates performance by abundance thresholds (ROC curves, AUC heatmaps).
#   6. Produces diagnostic plots of log-odds distributions.
#
# Inputs (CSV):
#   - host_pathogen_threshold_summary.csv
#   - 100_host_microbe_pairs_11_11_corrected.csv
#
# Key outputs (in-memory objects):
#   - df_with_hypergeo : full dataset with scores
#   - df_subset        : sampled pairs with metrics
#   - all_roc_curves   : ROC curve points per threshold
#   - auc_labels       : AUC summary per threshold (PHI-Base vs Negative Control)
#   - auc_grid         : AUC summary by comparison and thresholds
#
# Configure file paths as needed.
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
library(tibble)

set.seed(1)  # for reproducible sampling

DATA_FILE     <- "host_pathogen_threshold_summary.csv"
METADATA_FILE <- "100_host_microbe_pairs_11_11_corrected.csv"

if (!file.exists(DATA_FILE)) {
  stop("Data file not found: ", DATA_FILE)
}
if (!file.exists(METADATA_FILE)) {
  stop("Metadata file not found: ", METADATA_FILE)
}

###############################################################################
# Load data
###############################################################################

df_data  <- read.csv(DATA_FILE, check.names = FALSE)
metadata <- read.csv(METADATA_FILE, check.names = FALSE)

df_data <- as.data.frame(df_data)
metadata <- as.data.frame(metadata)

###############################################################################
# Clean and standardize columns
###############################################################################

# Main data (host–microbe counts per threshold)
df_data <- df_data %>%
  rename(
    host_tax_id                        = host_tax_id,
    microbe_tax_id                     = path_tax_id,
    host_threshold                     = host_threshold,
    microbe_threshold                  = path_threshold,
    num_host_datasets_at_threshold     = host_dataset_count,
    num_microbe_datasets_at_threshold  = pathogen_dataset_count,
    num_shared_datasets_both_threshold = both_dataset_count
  )

# Metadata (interaction labels and names)
metadata <- metadata %>%
  rename(
    microbe_tax_id     = path_tax_id,
    Microbe_Name       = Pathogen_Name,
    microbe_root_label = pathogen_root_label
  )

###############################################################################
# Sanity checks on host / microbe coverage
###############################################################################

host_meta    <- unique(metadata$host_tax_id)
microbe_meta <- unique(metadata$microbe_tax_id)
host_data    <- unique(df_data$host_tax_id)
microbe_data <- unique(df_data$microbe_tax_id)

# These calls are mainly for interactive inspection
setdiff(host_meta, host_data)     # in metadata but not in data
setdiff(host_data, host_meta)     # in data but not in metadata
setdiff(microbe_meta, microbe_data)
setdiff(microbe_data, microbe_meta)

length(host_meta)
length(microbe_meta)
length(host_data)
length(microbe_data)

###############################################################################
# Hypergeometric score and log-odds
###############################################################################

# Total number of datasets (update as appropriate)
N_total_datasets <- 29289124L

# Hypergeometric P(X >= observed overlap)
df_with_hypergeo <- df_data %>%
  mutate(
    Hyp_Geo_Score = phyper(
      q          = num_shared_datasets_both_threshold - 1,
      m          = num_host_datasets_at_threshold,
      n          = N_total_datasets - num_host_datasets_at_threshold,
      k          = num_microbe_datasets_at_threshold,
      lower.tail = FALSE
    )
  ) %>%
  mutate(
    Hyp_Geo_Score = if_else(Hyp_Geo_Score == 0, 1e-323, Hyp_Geo_Score),
    neg_log_p     = -log10(Hyp_Geo_Score)
  )

# Log-odds (global background)
k_pseudo <- 1 / N_total_datasets

df_with_hypergeo <- df_with_hypergeo %>%
  mutate(
    fh      = num_host_datasets_at_threshold    / N_total_datasets,
    fp      = num_microbe_datasets_at_threshold / N_total_datasets,
    fhp     = num_shared_datasets_both_threshold / N_total_datasets,
    log_odds = log2((fhp + k_pseudo) / (fh * fp + k_pseudo))
  )

# Alternative log-odds variant (local background)
df_with_hypergeo <- df_with_hypergeo %>%
  mutate(
    fh_2 = num_host_datasets_at_threshold /
      (num_host_datasets_at_threshold + num_microbe_datasets_at_threshold),
    fp_2 = num_microbe_datasets_at_threshold /
      (num_host_datasets_at_threshold + num_microbe_datasets_at_threshold),
    fhp_2 = num_shared_datasets_both_threshold /
      (num_host_datasets_at_threshold + num_microbe_datasets_at_threshold),
    log_odds_2 = log2(
      (fhp_2 + 1 / (num_host_datasets_at_threshold +
                      num_microbe_datasets_at_threshold)) /
        (fh_2 * fp_2 + 1 / (num_host_datasets_at_threshold +
                              num_microbe_datasets_at_threshold))
    )
  )

###############################################################################
# Resolve duplicate host–microbe pairs in metadata
###############################################################################

# Inspect duplicates per pair and interaction type
dups <- metadata %>%
  count(host_tax_id, microbe_tax_id, Interaction_type, name = "n_pairs") %>%
  filter(n_pairs > 1)

dups

# Priority order for resolving duplicates (edit to match actual labels)
priority_order <- c("Positive", "Negative", "Commensal", "Phi_Base")

metadata_prioritized <- metadata %>%
  mutate(
    Interaction_type_prio = factor(Interaction_type, levels = priority_order)
  ) %>%
  group_by(host_tax_id, microbe_tax_id) %>%
  arrange(Interaction_type_prio, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(host_tax_id, microbe_tax_id, Interaction_type)

###############################################################################
# Join interaction labels and label non-matches as Random
###############################################################################

df_with_labels <- df_with_hypergeo %>%
  left_join(
    metadata_prioritized,
    by = c("host_tax_id", "microbe_tax_id")
  ) %>%
  mutate(
    Interaction_type = tidyr::replace_na(Interaction_type, "Random")
  )

df_with_labels %>%
  count(Interaction_type)

###############################################################################
# Sample balanced subset of 100 pairs per interaction type (if available)
###############################################################################

sampled_pairs <- df_with_labels %>%
  distinct(host_tax_id, microbe_tax_id, Interaction_type) %>%
  group_by(Interaction_type) %>%
  slice_sample(n = min(100L, n())) %>%  # avoids error if < 100 pairs
  ungroup()

df_subset <- df_with_labels %>%
  inner_join(
    sampled_pairs,
    by = c("host_tax_id", "microbe_tax_id", "Interaction_type")
  )

df_subset %>%
  distinct(host_tax_id, microbe_tax_id, Interaction_type) %>%
  count(Interaction_type)

###############################################################################
# Co-occurrence metrics (Jaccard, Sørensen–Dice, Ochiai, phi, MI)
###############################################################################

# Universe size for 2x2 table
N_universe <- 34724088L

# Base counts and Jaccard / Sørensen–Dice
df_subset <- df_subset %>%
  mutate(
    a = as.numeric(num_shared_datasets_both_threshold),
    b = as.numeric(num_host_datasets_at_threshold)    - a,
    c = as.numeric(num_microbe_datasets_at_threshold) - a
  ) %>%
  mutate(
    b = pmax(b, 0),
    c = pmax(c, 0),
    Jaccard = if_else(a + b + c > 0, a / (a + b + c), NA_real_),
    Sorensen_Dice = if_else(2 * a + b + c > 0,
                            2 * a / (2 * a + b + c),
                            NA_real_)
  )

# Ochiai / cosine, phi, and C-score
df_subset <- df_subset %>%
  mutate(
    d = as.numeric(N_universe) - (a + b + c),
    ochiai = if_else(
      (a + b) > 0 & (a + c) > 0,
      a / sqrt((a + b) * (a + c)),
      NA_real_
    ),
    phi = if_else(
      (a + b) > 0 & (a + c) > 0 & (b + d) > 0 & (c + d) > 0,
      (a * d - b * c) / sqrt((a + b) * (a + c) * (b + d) * (c + d)),
      NA_real_
    ),
    c_score = (num_host_datasets_at_threshold - a) *
      (num_microbe_datasets_at_threshold - a)
  )

# Mutual information (2x2)
df_subset <- df_subset %>%
  mutate(
    p_a = a / N_universe,
    p_b = b / N_universe,
    p_c = c / N_universe,
    p_d = d / N_universe,
    p_row1 = p_a + p_b,
    p_row2 = p_c + p_d,
    p_col1 = p_a + p_c,
    p_col2 = p_b + p_d,
    mutual_info = (
      ifelse(p_a > 0, p_a * log2(p_a / (p_row1 * p_col1)), 0) +
        ifelse(p_b > 0, p_b * log2(p_b / (p_row1 * p_col2)), 0) +
        ifelse(p_c > 0, p_c * log2(p_c / (p_row2 * p_col1)), 0) +
        ifelse(p_d > 0, p_d * log2(p_d / (p_row2 * p_col2)), 0)
    )
  )

###############################################################################
# ROC curves by (host_threshold, microbe_threshold)
#   Comparison: PHI-Base (positive) vs Negative Control (negative)
###############################################################################

threshold_combos <- df_subset %>%
  distinct(host_threshold, microbe_threshold)

roc_curve_list <- list()
auc_label_list <- list()

for (i in seq_len(nrow(threshold_combos))) {
  host_thresh    <- threshold_combos$host_threshold[i]
  microbe_thresh <- threshold_combos$microbe_threshold[i]
  
  # Data for this threshold combo
  long_labeled <- df_subset %>%
    filter(
      host_threshold    == host_thresh,
      microbe_threshold == microbe_thresh,
      Interaction_type %in% c("PHI-Base", "Negative Control")
    ) %>%
    mutate(
      label = ifelse(Interaction_type == "PHI-Base", 1L, 0L)
    ) %>%
    filter(!is.na(log_odds))
  
  n_obs <- nrow(long_labeled)
  
  # Need at least one positive and one negative
  if (n_obs == 0 || length(unique(long_labeled$label)) < 2) {
    auc_label_list[[length(auc_label_list) + 1]] <- tibble(
      host_threshold    = host_thresh,
      microbe_threshold = microbe_thresh,
      AUC       = NA_real_,
      N_total   = n_obs,
      N_correct = NA_integer_
    )
    next
  }
  
  roc_obj <- tryCatch({
    roc(
      response  = long_labeled$label,
      predictor = long_labeled$log_odds,
      quiet     = TRUE
    )
  }, error = function(e) NULL)
  
  if (is.null(roc_obj)) {
    auc_label_list[[length(auc_label_list) + 1]] <- tibble(
      host_threshold    = host_thresh,
      microbe_threshold = microbe_thresh,
      AUC       = NA_real_,
      N_total   = n_obs,
      N_correct = NA_integer_
    )
    next
  }
  
  auc_val <- as.numeric(auc(roc_obj))
  
  # Best cutoff (Youden index)
  best_coords <- coords(
    roc_obj,
    x   = "best",
    ret = c("threshold", "tp", "tn", "fp", "fn", "accuracy"),
    transpose = FALSE
  )
  
  N_correct <- best_coords$tp + best_coords$tn
  
  # Store ROC curve points
  roc_curve_list[[length(roc_curve_list) + 1]] <- tibble(
    specificity       = rev(roc_obj$specificities),
    sensitivity       = rev(roc_obj$sensitivities),
    host_threshold    = host_thresh,
    microbe_threshold = microbe_thresh,
    AUC               = auc_val
  )
  
  # Store AUC summary
  auc_label_list[[length(auc_label_list) + 1]] <- tibble(
    host_threshold    = host_thresh,
    microbe_threshold = microbe_thresh,
    AUC       = auc_val,
    N_total   = n_obs,
    N_correct = as.integer(N_correct)
  )
}

all_roc_curves <- bind_rows(roc_curve_list)
auc_labels     <- bind_rows(auc_label_list)

# ROC plot faceted by thresholds
ggplot(all_roc_curves, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_grid(
    rows = vars(microbe_threshold),
    cols = vars(host_threshold),
    labeller = labeller(
      host_threshold    = function(x) paste0("host=", x),
      microbe_threshold = function(x) paste0("microbe=", x)
    )
  ) +
  geom_text(
    data = auc_labels,
    aes(
      x = 0.6,
      y = 0.2,
      label = paste0(
        "AUC = ", round(AUC, 2),
        "\nN = ", N_correct
      )
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  labs(
    title = "ROC curves by host and microbe thresholds",
    x = "False positive rate (1 - specificity)",
    y = "True positive rate (sensitivity)"
  ) +
  theme_minimal()

###############################################################################
# AUC heatmaps for multiple interaction-type comparisons
###############################################################################

comparisons <- tribble(
  ~comparison,      ~pos_class,          ~neg_class,
  "Pos vs Neg",     "Positive Control",  "Negative Control",
  "Pos vs Random",  "Positive Control",  "Random",
  "Comm vs Neg",    "Commensal",         "Negative Control",
  "Comm vs Random", "Commensal",         "Random",
  "PHI vs Neg",     "PHI-Base",          "Negative Control",
  "PHI vs Random",  "PHI-Base",          "Random"
)

threshold_combos <- df_subset %>%
  distinct(host_threshold, microbe_threshold)

auc_results_list <- list()

for (j in seq_len(nrow(comparisons))) {
  comp_name <- comparisons$comparison[j]
  pos_type  <- comparisons$pos_class[j]
  neg_type  <- comparisons$neg_class[j]
  
  for (i in seq_len(nrow(threshold_combos))) {
    host_thresh    <- threshold_combos$host_threshold[i]
    microbe_thresh <- threshold_combos$microbe_threshold[i]
    
    long_labeled <- df_subset %>%
      filter(
        host_threshold    == host_thresh,
        microbe_threshold == microbe_thresh,
        Interaction_type %in% c(pos_type, neg_type)
      ) %>%
      mutate(
        label = ifelse(Interaction_type == pos_type, 1L, 0L)
      ) %>%
      filter(!is.na(log_odds))
    
    n_obs <- nrow(long_labeled)
    
    if (n_obs == 0 || length(unique(long_labeled$label)) < 2) {
      auc_results_list[[length(auc_results_list) + 1]] <- tibble(
        comparison        = comp_name,
        host_threshold    = host_thresh,
        microbe_threshold = microbe_thresh,
        AUC       = NA_real_,
        N_total   = n_obs,
        N_correct = NA_integer_
      )
      next
    }
    
    roc_obj <- tryCatch({
      roc(
        response  = long_labeled$label,
        predictor = long_labeled$log_odds,
        quiet     = TRUE
      )
    }, error = function(e) NULL)
    
    if (is.null(roc_obj)) {
      auc_results_list[[length(auc_results_list) + 1]] <- tibble(
        comparison        = comp_name,
        host_threshold    = host_thresh,
        microbe_threshold = microbe_thresh,
        AUC       = NA_real_,
        N_total   = n_obs,
        N_correct = NA_integer_
      )
      next
    }
    
    auc_val <- as.numeric(auc(roc_obj))
    
    best_coords <- coords(
      roc_obj,
      x   = "best",
      ret = c("tp", "tn", "fp", "fn", "accuracy"),
      transpose = FALSE
    )
    
    N_correct <- best_coords$tp + best_coords$tn
    
    auc_results_list[[length(auc_results_list) + 1]] <- tibble(
      comparison        = comp_name,
      host_threshold    = host_thresh,
      microbe_threshold = microbe_thresh,
      AUC       = auc_val,
      N_total   = n_obs,
      N_correct = as.integer(N_correct)
    )
  }
}

auc_grid <- bind_rows(auc_results_list)

auc_grid_plot <- auc_grid %>%
  mutate(
    host_threshold    = factor(host_threshold,    levels = sort(unique(host_threshold))),
    microbe_threshold = factor(microbe_threshold, levels = sort(unique(microbe_threshold)))
  )

ggplot(auc_grid_plot,
       aes(x = host_threshold, y = microbe_threshold, fill = AUC)) +
  geom_tile() +
  facet_wrap(~ comparison) +
  scale_fill_viridis_c(option = "magma", na.value = "grey90") +
  labs(
    title = "AUC heatmaps by interaction-type comparison",
    x = "Host abundance threshold",
    y = "Microbe abundance threshold",
    fill = "AUC"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )

###############################################################################
# Log-odds distribution by interaction type
###############################################################################

df_subset %>%
  filter(!is.na(log_odds)) %>%
  ggplot(aes(x = log_odds)) +
  geom_histogram(bins = 50) +
  facet_wrap(~ Interaction_type, scales = "free_y") +
  labs(
    title = "Distribution of log-odds by interaction type",
    x = "log-odds",
    y = "Count"
  ) +
  theme_minimal()

###############################################################################
# Fisher's exact test for each pair
###############################################################################

df_subset <- df_subset %>%
  rowwise() %>%
  mutate(
    a = num_shared_datasets_both_threshold,
    b = num_host_datasets_at_threshold - a,
    c = num_microbe_datasets_at_threshold - a,
    d = b + c,  # minimal assumption to close the table
    fisher_p = fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
  ) %>%
  ungroup() %>%
  select(-a, -b, -c, -d)

