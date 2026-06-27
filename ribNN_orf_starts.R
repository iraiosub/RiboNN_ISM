#!/usr/bin/env/Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mgcv))

orf_definition_colors <- c(
  "Not detected (expression-filtered)" = "#6B6B6B",
  "Not detected (low expression)"      = "#C7C7C7",
  "uORF"                               = "#76baa6",
  "uoORF"                              = "#74669d",
  "Canonical CDS"                      = "#b3cde3"
)

orf_definition_levels <- c(
  "Not detected (low expression)",
  "Not detected (expression-filtered)",
  "uORF",
  "uoORF"
)

# Sub-stratified palette: theoreticals split by structural type (uORF vs uoORF)
orf_definition_sub_colors <- c(
  "Not det. (low expr.) · uORF"    = "#C7DDD8",   # light teal-grey
  "Not det. (low expr.) · uoORF"   = "#C7C2D8",   # light purple-grey
  "Not det. (filt.) · uORF"        = "#6B8B85",   # mid teal-grey
  "Not det. (filt.) · uoORF"       = "#6B6580",   # mid purple-grey
  "uORF"                            = "#76baa6",
  "uoORF"                           = "#74669d"
)

orf_definition_sub_levels <- c(
  "Not det. (low expr.) · uORF",
  "Not det. (low expr.) · uoORF",
  "Not det. (filt.) · uORF",
  "Not det. (filt.) · uoORF",
  "uORF",
  "uoORF"
)

direction_colors <- c(
  "increase"  = "#d73027",
  "unchanged" = "#fee090",
  "decrease"  = "#4575b4"
)

codon_pos_labels <- c("1 (A)", "2 (T/U)", "3 (G)")


# ============================================================
# Load data
# ============================================================

# Paste a run output folder here if you want the script to ignore the defaults.
# Use the run folder that contains summaries/, not final/ and not a .tsv.gz file.
INPUT_RESULTS_DIR <- "/Volumes/lab-ulej/home/shared/oscar_ira_riboloco/RiboNN_ISM/all_utr5_mutagenesis/output/mouse_whole_atg_deletions_mean_predicted_TE"
# Example:
# INPUT_RESULTS_DIR <- "/nemo/lab/ulej/home/shared/oscar_ira_riboloco/RiboNN_ISM/all_utr5_mutagenesis/output/mouse_whole_atg_deletions_mean_predicted_TE"

# Usually leave this blank. Set only if you want a non-default master table.
INPUT_MASTER_TABLE <- ""

default_shared_root <- if (dir.exists("/camp/lab/ulej/home/shared/oscar_ira_riboloco")) {
  "/camp/lab/ulej/home/shared/oscar_ira_riboloco"
} else {
  "/Volumes/lab-ulej/home/shared/oscar_ira_riboloco"
}
shared_root <- Sys.getenv("RIBONN_SHARED_ROOT", default_shared_root)
analysis_results.dir <- Sys.getenv(
  "RIBONN_ANALYSIS_RESULTS_DIR",
  file.path(shared_root, "analysis_results")
)
output_root.dir <- Sys.getenv(
  "RIBONN_UTR5_OUTPUT_ROOT",
  file.path(shared_root, "RiboNN_ISM", "all_utr5_mutagenesis", "output")
)

input_results_dir <- trimws(INPUT_RESULTS_DIR)
input_master_table <- trimws(INPUT_MASTER_TABLE)

infer_species_from_path <- function(path) {
  normalized_path <- normalizePath(path, mustWork = FALSE)
  folder <- basename(normalized_path)
  if (folder %in% c("summaries", "final")) {
    folder <- basename(dirname(normalized_path))
  }
  if (str_detect(folder, "^mouse_")) {
    return("mouse")
  }
  if (str_detect(folder, "^human_")) {
    return("human")
  }
  Sys.getenv("RIBONN_ORF_SPECIES", "human")
}

species <- if (nzchar(input_results_dir)) {
  infer_species_from_path(input_results_dir)
} else {
  Sys.getenv("RIBONN_ORF_SPECIES", "human")
}
if (!species %in% c("human", "mouse")) {
  stop("Species must be 'human' or 'mouse'. Set RIBONN_ORF_SPECIES if the path is ambiguous.")
}
default_te_label <- if (species == "mouse") "mean_predicted_TE" else "normal_brain_tissue"
te_label <- Sys.getenv("RIBONN_ORF_TE_LABEL", default_te_label)
orf_screen <- Sys.getenv("RIBONN_ORF_SCREEN", "orf_starts")
if (orf_screen %in% c("whole-atg-deletion", "whole_atg_deletions")) {
  orf_screen <- "whole_atg_deletion"
}
if (!orf_screen %in% c("orf_starts", "whole_atg_deletion")) {
  stop("RIBONN_ORF_SCREEN must be 'orf_starts' or 'whole_atg_deletion'")
}
default_results_name <- if (orf_screen == "whole_atg_deletion") {
  paste0(species, "_whole_atg_deletions_", te_label)
} else {
  paste0(species, "_orf_starts_", te_label)
}
ism_results.dir <- Sys.getenv(
  "RIBONN_ORF_RESULTS_DIR",
  file.path(output_root.dir, default_results_name)
)
if (nzchar(input_results_dir)) {
  ism_results.dir <- input_results_dir
}
results_source <- if (nzchar(input_results_dir)) {
  "INPUT_RESULTS_DIR"
} else if (nzchar(Sys.getenv("RIBONN_ORF_RESULTS_DIR"))) {
  "RIBONN_ORF_RESULTS_DIR"
} else {
  "default"
}
ism_summary.dir <- if (dir.exists(file.path(ism_results.dir, "summaries"))) {
  file.path(ism_results.dir, "summaries")
} else {
  ism_results.dir
}

default_master_tables <- c(
  human = file.path(
    analysis_results.dir,
    "human_brain.unmixing.master_table.with_below_tpm_threshold.tsv.gz"
  ),
  mouse = file.path(
    analysis_results.dir,
    "cross_tissue.unmixing.master_table.with_below_tpm_threshold.tsv.gz"
  )
)
master_table <- if (nzchar(input_master_table)) {
  input_master_table
} else {
  Sys.getenv("RIBONN_ORF_MASTER_TABLE", default_master_tables[[species]])
}

cat("Species:", species, "\n")
cat("Results source:", results_source, "\n")
cat("ISM summaries:", ism_summary.dir, "\n")
cat("Master table:", master_table, "\n")

ism.ls <- list.files(ism_summary.dir, pattern = "\\.orfs\\.tsv\\.gz$", full.names = TRUE)
if (length(ism.ls) == 0) {
  stop("No shard ORF summaries found in ", ism_summary.dir)
}
ism.df <- rbindlist(lapply(ism.ls, fread))

pos.ls <- list.files(ism_summary.dir, pattern = "\\.positions\\.tsv\\.gz$", full.names = TRUE)
if (length(pos.ls) == 0) {
  stop("No shard position summaries found in ", ism_summary.dir)
}
pos.df <- rbindlist(lapply(pos.ls, fread))

whole_atg_mode <- "whole_atg_deletion_te_change" %in% names(ism.df)
cat("Screen mode:", if_else(whole_atg_mode, "whole_atg_deletion", "orf_start_codon"), "\n")

master.df <- fread(master_table)
required_master_cols <- c("transcript_id", "orf_id", "orf_label", "tpm_threshold_group", "orf_rel_to_cds")
missing_master_cols <- setdiff(required_master_cols, names(master.df))
if (length(missing_master_cols) > 0) {
  stop("Master table is missing required columns: ", paste(missing_master_cols, collapse = ", "))
}

master.df <- master.df %>%
  mutate(orf_definition = case_when(
    orf_label == "Upstream ORF" ~ "uORF",
    orf_label == "Upstream overlapping ORF" ~ "uoORF",
    orf_label == "5'UTR theoretical (AUG + in-frame stop)" & tpm_threshold_group == "above_threshold" ~ "Not detected (expression-filtered)",
    orf_label == "5'UTR theoretical (AUG + in-frame stop)" & tpm_threshold_group == "below_threshold" ~ "Not detected (low expression)",
    TRUE ~ NA_character_
  ))

master_orf_def.df <- master.df %>%
  filter(!is.na(orf_definition)) %>%
  distinct(orf_id, orf_definition)

# Sub-stratified definition: theoreticals split by orf_rel_to_cds (uORF vs uoORF structure)
master_orf_def_sub.df <- master.df %>%
  filter(!is.na(orf_definition)) %>%
  mutate(orf_definition_sub = case_when(
    orf_definition == "uORF"  ~ "uORF",
    orf_definition == "uoORF" ~ "uoORF",
    orf_definition == "Not detected (low expression)"      & orf_rel_to_cds == "uORF"  ~
      "Not det. (low expr.) · uORF",
    orf_definition == "Not detected (low expression)"      & orf_rel_to_cds == "uoORF" ~
      "Not det. (low expr.) · uoORF",
    orf_definition == "Not detected (expression-filtered)" & orf_rel_to_cds == "uORF"  ~
      "Not det. (filt.) · uORF",
    orf_definition == "Not detected (expression-filtered)" & orf_rel_to_cds == "uoORF" ~
      "Not det. (filt.) · uoORF",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(orf_definition_sub)) %>%
  distinct(orf_id, orf_definition_sub)


# ============================================================
# Prepare ORF-level table
# ============================================================

ism.df <- ism.df %>%
  semi_join(master.df, by = "transcript_id") %>%
  tidyr::unite("orf_id", transcript_id, orf_start_1based, orf_stop_1based, orf_frame,
               sep = "_", remove = TRUE) %>%
  left_join(master_orf_def.df, by = "orf_id") %>%
  filter(!is.na(orf_definition)) %>%
  mutate(orf_definition = factor(orf_definition, levels = orf_definition_levels))

cat("ORF-level rows after joining master table:", nrow(ism.df), "\n")


# ============================================================
# Prepare position-level table
# ============================================================

pos.df <- pos.df %>%
  tidyr::unite("orf_id", transcript_id, orf_start_1based, orf_stop_1based, orf_frame,
               sep = "_", remove = FALSE) %>%
  left_join(master_orf_def.df, by = "orf_id") %>%
  filter(!is.na(orf_definition)) %>%
  mutate(orf_definition = factor(orf_definition, levels = orf_definition_levels))

if (!whole_atg_mode) {
  # In orf_start_codon mode every row is one of the 3 ATG positions
  pos.df <- pos.df %>%
    mutate(
      pos_label = factor(codon_pos_labels[orf_position_in_start_codon],
                         levels = codon_pos_labels)
    )
}

cat("Position-level rows after joining master table:", nrow(pos.df), "\n")

if (whole_atg_mode) {


# ============================================================
# Whole-ATG deletion mode
# ============================================================

suppressPackageStartupMessages(library(cowplot))

SHOW_LOW_EXPR <- FALSE
WHOLE_ATG_BOOTSTRAP_REPS <- as.integer(Sys.getenv(
  "RIBONN_WHOLE_ATG_BOOTSTRAP_REPS",
  "0"
))
WHOLE_ATG_BOOTSTRAP_SEED <- as.integer(Sys.getenv(
  "RIBONN_WHOLE_ATG_BOOTSTRAP_SEED",
  "1"
))
if (is.na(WHOLE_ATG_BOOTSTRAP_REPS) || WHOLE_ATG_BOOTSTRAP_REPS < 0) {
  stop("RIBONN_WHOLE_ATG_BOOTSTRAP_REPS must be an integer >= 0.")
}
if (is.na(WHOLE_ATG_BOOTSTRAP_SEED)) {
  stop("RIBONN_WHOLE_ATG_BOOTSTRAP_SEED must be an integer.")
}
set.seed(WHOLE_ATG_BOOTSTRAP_SEED)
cat("Whole-ATG HL bootstrap reps:", WHOLE_ATG_BOOTSTRAP_REPS, "\n")

find_phylop_col <- function(df) {
  exact_match <- names(df)[tolower(names(df)) == "phylop"]
  if (length(exact_match) > 0) {
    return(exact_match[[1]])
  }

  loose_match <- names(df)[str_detect(tolower(names(df)), "phylop")]
  if (length(loose_match) > 0) {
    return(loose_match[[1]])
  }

  NA_character_
}

find_gene_label_col <- function(df) {
  candidates <- c(
    "gene_name",
    "gene_symbol",
    "symbol",
    "external_gene_name",
    "mgi_symbol",
    "hgnc_symbol"
  )
  for (candidate in candidates) {
    exact_match <- names(df)[tolower(names(df)) == tolower(candidate)]
    if (length(exact_match) > 0) {
      return(exact_match[[1]])
    }
  }

  NA_character_
}

phylop_col <- find_phylop_col(master.df)
if (is.na(phylop_col)) {
  stop("Master table is missing a phyloP column needed for whole-ATG deletion scatter plots.")
}
cat("phyloP column:", phylop_col, "\n")

PUB_SUB_LEVELS <- if (SHOW_LOW_EXPR) orf_definition_sub_levels else
  grep("low expr", orf_definition_sub_levels, value = TRUE, invert = TRUE)

PUB_SUB_COLORS <- orf_definition_sub_colors[PUB_SUB_LEVELS]

FACET_GROUP <- c(
  "uORF"                         = "uORF",
  "uoORF"                        = "uoORF",
  "Not det. (filt.) · uORF"     = "uORF",
  "Not det. (filt.) · uoORF"    = "uoORF",
  "Not det. (low expr.) · uORF"  = "uORF",
  "Not det. (low expr.) · uoORF" = "uoORF"
)

add_facet_group <- function(df) {
  df %>%
    mutate(facet_group = factor(
      FACET_GROUP[as.character(orf_definition_sub)],
      levels = c("uORF", "uoORF")
    ))
}

add_dist_bin <- function(df, offset_col) {
  df %>% mutate(
    distance_bin = factor(
      case_when(
        .data[[offset_col]] >= -30  ~ "≤30 nt",
        .data[[offset_col]] >= -60  ~ "30–60 nt",
        .data[[offset_col]] >= -120 ~ "60–120 nt",
        TRUE                         ~ ">120 nt"
      ),
      levels = c("≤30 nt", "30–60 nt", "60–120 nt", ">120 nt")
    )
  )
}

pub_base_theme <- theme_classic(base_size = 11) +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(4, "mm"),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9),
    plot.title       = element_text(size = 11, face = "bold"),
    plot.subtitle    = element_text(size = 8,  colour = "grey40"),
    strip.background = element_blank(),
    strip.text       = element_text(size = 9, face = "bold")
  )

pub_box_theme <- pub_base_theme +
  theme(
    axis.text.x        = element_text(angle = 30, hjust = 1, size = 8),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)
  )

orf_start_dist.df <- pos.df %>%
  mutate(
    orf_start_1based_num = as.numeric(orf_start_1based),
    orf_stop_1based_num  = as.numeric(orf_stop_1based),
    orf_length_nt = abs(orf_stop_1based_num - orf_start_1based_num) + 1
  ) %>%
  dplyr::select(orf_id, orf_start_offset = offset_from_cds_start, orf_length_nt)

master_phylop.df <- master.df %>%
  filter(!is.na(orf_definition)) %>%
  mutate(phylop = suppressWarnings(as.numeric(.data[[phylop_col]]))) %>%
  group_by(orf_id) %>%
  summarise(
    phylop = {
      finite_phylop <- phylop[is.finite(phylop)]
      if (length(finite_phylop) == 0) NA_real_ else mean(finite_phylop)
    },
    .groups = "drop"
  )

ism_whole.df <- ism.df %>%
  left_join(master_orf_def_sub.df, by = "orf_id") %>%
  left_join(master_phylop.df, by = "orf_id") %>%
  filter(!is.na(orf_definition_sub)) %>%
  { if (!SHOW_LOW_EXPR) filter(., !grepl("low expr", orf_definition_sub)) else . } %>%
  mutate(
    orf_definition_sub = factor(orf_definition_sub, levels = PUB_SUB_LEVELS),
    detected_group     = factor(
      if_else(orf_definition_sub %in% c("uORF", "uoORF"), "Detected", "Not detected"),
      levels = c("Detected", "Not detected")
    ),
    whole_atg_delta = as.numeric(whole_atg_deletion_te_change)
  ) %>%
  add_facet_group() %>%
  left_join(orf_start_dist.df, by = "orf_id") %>%
  add_dist_bin("orf_start_offset") %>%
  filter(is.finite(whole_atg_delta))

gene_label_col <- find_gene_label_col(ism_whole.df)
if (!is.na(gene_label_col)) {
  ism_whole.df <- ism_whole.df %>%
    mutate(plot_gene_name = as.character(.data[[gene_label_col]]))
  cat("Gene label column:", gene_label_col, "\n")
} else {
  master_gene_label_col <- find_gene_label_col(master.df)
  if (!is.na(master_gene_label_col)) {
    master_gene_labels.df <- master.df %>%
      filter(!is.na(orf_definition)) %>%
      transmute(
        orf_id,
        plot_gene_name = as.character(.data[[master_gene_label_col]])
      ) %>%
      distinct(orf_id, plot_gene_name)
    ism_whole.df <- ism_whole.df %>%
      left_join(master_gene_labels.df, by = "orf_id")
    cat("Gene label column:", master_gene_label_col, "(master table)\n")
  } else {
    ism_whole.df <- ism_whole.df %>%
      mutate(plot_gene_name = NA_character_)
    warning("No gene-name column found; scatter labels will use generic top-ORF labels.")
  }
}

WHOLE_ATG_PHYLOP_COR_METHOD <- "spearman"
WHOLE_ATG_PHYLOP_COR_LABEL <- if_else(
  WHOLE_ATG_PHYLOP_COR_METHOD == "pearson",
  "Pearson r",
  "Spearman rho"
)

format_corr_p <- function(p) {
  ifelse(is.na(p), "NA", format.pval(p, digits = 2, eps = 1e-3))
}

whole_atg_phylop_cor_one_group <- function(df) {
  complete_df <- df %>%
    filter(is.finite(phylop), is.finite(whole_atg_delta))

  n_complete <- nrow(complete_df)
  base_row <- tibble(
    method = WHOLE_ATG_PHYLOP_COR_METHOD,
    n_complete = n_complete,
    correlation = NA_real_,
    p_value = NA_real_
  )

  if (
    n_complete < 3 ||
    n_distinct(complete_df$phylop) < 2 ||
    n_distinct(complete_df$whole_atg_delta) < 2
  ) {
    return(base_row %>%
      mutate(label = paste0(WHOLE_ATG_PHYLOP_COR_LABEL, "=NA\nn=", n_complete)))
  }

  cor_test <- tryCatch(
    suppressWarnings(cor.test(
      complete_df$phylop,
      complete_df$whole_atg_delta,
      method = WHOLE_ATG_PHYLOP_COR_METHOD,
      exact = FALSE
    )),
    error = function(err) err
  )

  if (inherits(cor_test, "error")) {
    return(base_row %>%
      mutate(label = paste0(WHOLE_ATG_PHYLOP_COR_LABEL, "=NA\nn=", n_complete)))
  }

  base_row %>%
    mutate(
      correlation = unname(cor_test$estimate),
      p_value = cor_test$p.value,
      label = paste0(
        WHOLE_ATG_PHYLOP_COR_LABEL, "=", sprintf("%.2f", correlation),
        "\nP=", format_corr_p(p_value),
        "\nn=", n_complete
      )
    )
}

whole_atg_phylop_cor_labels <- function(df, group_cols = character()) {
  if (length(group_cols) == 0) {
    return(whole_atg_phylop_cor_one_group(df))
  }

  df %>%
    group_by(across(all_of(group_cols))) %>%
    group_modify(~whole_atg_phylop_cor_one_group(.x)) %>%
    ungroup()
}

whole_atg_phylop.df <- ism_whole.df %>%
  filter(
    is.finite(phylop),
    is.finite(whole_atg_delta),
    as.character(orf_definition_sub) %in% PUB_SUB_LEVELS,
    detected_group %in% c("Detected", "Not detected")
  ) %>%
  mutate(
    orf_definition = fct_drop(orf_definition),
    orf_definition_sub = fct_drop(orf_definition_sub),
    detected_group = fct_drop(detected_group),
    facet_group = fct_drop(facet_group)
  )

cat("Whole-ATG phyloP scatter complete rows:", nrow(whole_atg_phylop.df), "\n")

whole_atg_phylop_plot.df <- whole_atg_phylop.df %>%
  arrange(detected_group == "Detected")

whole_atg_phylop_overall_corr.df <- whole_atg_phylop_cor_labels(whole_atg_phylop.df)
whole_atg_phylop_orf_corr.df <- whole_atg_phylop_cor_labels(
  whole_atg_phylop.df,
  "orf_definition"
)
whole_atg_phylop_sub_corr.df <- whole_atg_phylop_cor_labels(
  whole_atg_phylop.df,
  c("facet_group", "detected_group")
)
whole_atg_phylop_distance_corr.df <- whole_atg_phylop_cor_labels(
  whole_atg_phylop.df %>% filter(!is.na(distance_bin)),
  "distance_bin"
)
whole_atg_phylop_distance_structure_corr.df <- whole_atg_phylop_cor_labels(
  whole_atg_phylop.df %>% filter(!is.na(distance_bin), !is.na(facet_group)),
  c("facet_group", "distance_bin")
)

whole_atg_phylop_corr_stats.df <- bind_rows(
  whole_atg_phylop_overall_corr.df %>%
    mutate(stratification = "pooled"),
  whole_atg_phylop_distance_corr.df %>%
    mutate(stratification = "distance_from_cds"),
  whole_atg_phylop_distance_structure_corr.df %>%
    mutate(stratification = "distance_from_cds_by_uORF_uoORF"),
  whole_atg_phylop_orf_corr.df %>%
    mutate(stratification = "orf_definition"),
  whole_atg_phylop_sub_corr.df %>%
    mutate(stratification = "orf_definition_by_uORF_uoORF")
) %>%
  relocate(stratification)

print(whole_atg_phylop_corr_stats.df)

phylop_scatter_theme <- pub_base_theme +
  theme(
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
    legend.position = "right",
    plot.margin = margin(5.5, 12, 5.5, 5.5)
  )

phylop_target_genes <- c("scn2a", "scn8a")

label_orf_name <- function(df) {
  case_when(
    !is.na(df$plot_gene_name) & nzchar(df$plot_gene_name) ~ df$plot_gene_name,
    TRUE ~ "max |ΔTE| ORF"
  )
}

grouped_slice_max_abs_delta <- function(df, group_cols) {
  if (length(group_cols) == 0) {
    return(df %>%
      slice_max(order_by = abs(whole_atg_delta), n = 1, with_ties = FALSE))
  }

  df %>%
    group_by(across(all_of(group_cols))) %>%
    slice_max(order_by = abs(whole_atg_delta), n = 1, with_ties = FALSE) %>%
    ungroup()
}

grouped_slice_max_delta <- function(df, group_cols) {
  if (length(group_cols) == 0) {
    return(df %>%
      slice_max(order_by = whole_atg_delta, n = 1, with_ties = FALSE))
  }

  df %>%
    group_by(across(all_of(group_cols))) %>%
    slice_max(order_by = whole_atg_delta, n = 1, with_ties = FALSE) %>%
    ungroup()
}

target_gene_phylop_labels <- function(df, group_cols = character()) {
  target_df <- df %>%
    filter(str_to_lower(plot_gene_name) %in% phylop_target_genes)

  if (length(group_cols) == 0) {
    target_df <- target_df %>%
      group_by(plot_gene_name) %>%
      slice_max(order_by = abs(whole_atg_delta), n = 1, with_ties = FALSE) %>%
      ungroup()
  } else {
    target_df <- target_df %>%
      group_by(across(all_of(c(group_cols, "plot_gene_name")))) %>%
      slice_max(order_by = abs(whole_atg_delta), n = 1, with_ties = FALSE) %>%
      ungroup()
  }

  target_df %>%
    mutate(
      plot_label = plot_gene_name,
      label_priority = 1L
    )
}

top_delta_phylop_labels <- function(df, group_cols = character()) {
  grouped_slice_max_abs_delta(df, group_cols) %>%
    mutate(
      plot_label = paste0(label_orf_name(.), "\nmax |ΔTE|"),
      label_priority = 2L
    )
}

top_detected_delta_phylop_labels <- function(df, group_cols = character()) {
  df %>%
    filter(detected_group == "Detected") %>%
    grouped_slice_max_delta(group_cols) %>%
    mutate(
      plot_label = paste0(label_orf_name(.), "\nhighest detected ΔTE"),
      label_priority = 3L
    )
}

phylop_label_df <- function(df, group_cols = character(), include_top = TRUE) {
  label_df <- target_gene_phylop_labels(df, group_cols)
  if (include_top) {
    label_df <- bind_rows(
      label_df,
      top_delta_phylop_labels(df, group_cols),
      top_detected_delta_phylop_labels(df, group_cols)
    )
  }

  label_df %>%
    arrange(label_priority) %>%
    distinct(across(all_of(c(group_cols, "orf_id"))), .keep_all = TRUE)
}

phylop_label_layer <- function(label_df, size = 2.2) {
  geom_label(
    data = label_df,
    aes(x = phylop, y = whole_atg_delta, label = plot_label),
    inherit.aes = FALSE,
    size = size,
    label.size = 0.12,
    label.padding = unit(0.08, "lines"),
    fill = "white",
    alpha = 0.75,
    colour = "grey15",
    show.legend = FALSE
  )
}

whole_atg_phylop_overall_label.df <- phylop_label_df(whole_atg_phylop.df)
whole_atg_phylop_orf_label.df <- phylop_label_df(
  whole_atg_phylop.df,
  "orf_definition"
)
whole_atg_phylop_distance_label.df <- phylop_label_df(
  whole_atg_phylop.df %>% filter(!is.na(distance_bin)),
  "distance_bin"
)
whole_atg_phylop_distance_structure_label.df <- phylop_label_df(
  whole_atg_phylop.df %>% filter(!is.na(distance_bin), !is.na(facet_group)),
  c("facet_group", "distance_bin")
)
whole_atg_phylop_sub_label.df <- phylop_label_df(
  whole_atg_phylop.df,
  c("facet_group", "detected_group")
)

whole_atg_phylop_scatter.gg <- ggplot(
  whole_atg_phylop_plot.df,
  aes(x = phylop, y = whole_atg_delta)
) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed",
             linewidth = 0.35) +
  geom_point(colour = "#3A6EA5", alpha = 0.35, size = 0.75, stroke = 0) +
  phylop_label_layer(whole_atg_phylop_overall_label.df, size = 2.3) +
  geom_text(
    data = whole_atg_phylop_overall_corr.df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.15,
    size = 3.0,
    colour = "grey20"
  ) +
  labs(
    x = "phyloP",
    y = "ΔTE (whole ATG deletion)",
    title = "phyloP vs whole-ATG deletion ΔTE"
  ) +
  phylop_scatter_theme

whole_atg_phylop_orf_type_scatter.gg <- ggplot(
  whole_atg_phylop_plot.df,
  aes(x = phylop, y = whole_atg_delta, colour = orf_definition)
) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed",
             linewidth = 0.35) +
  geom_point(alpha = 0.35, size = 0.65, stroke = 0) +
  phylop_label_layer(whole_atg_phylop_orf_label.df, size = 2.2) +
  geom_text(
    data = whole_atg_phylop_orf_corr.df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.15,
    size = 2.6,
    colour = "grey20"
  ) +
  facet_wrap(~orf_definition, ncol = 2) +
  scale_colour_manual(values = orf_definition_colors, name = "ORF class") +
  labs(
    x = "phyloP",
    y = "ΔTE (whole ATG deletion)",
    title = "phyloP vs whole-ATG deletion ΔTE by ORF class"
  ) +
  phylop_scatter_theme +
  guides(colour = "none")

whole_atg_phylop_distance_scatter.gg <- whole_atg_phylop_plot.df %>%
  filter(!is.na(distance_bin)) %>%
  ggplot(aes(x = phylop, y = whole_atg_delta)) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed",
             linewidth = 0.35) +
  geom_point(colour = "#3A6EA5", alpha = 0.35, size = 0.6, stroke = 0) +
  phylop_label_layer(whole_atg_phylop_distance_label.df, size = 2.1) +
  geom_text(
    data = whole_atg_phylop_distance_corr.df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.15,
    size = 2.4,
    colour = "grey20"
  ) +
  facet_wrap(~distance_bin, nrow = 1) +
  labs(
    x = "phyloP",
    y = "ΔTE (whole ATG deletion)",
    title = "phyloP vs whole-ATG deletion ΔTE by ORF start distance from CDS"
  ) +
  phylop_scatter_theme

distance_structure_colors <- PUB_SUB_COLORS
distance_structure_colors[grepl("^Not det\\.", names(distance_structure_colors))] <- "#AFAFAF"

whole_atg_phylop_distance_structure_scatter.gg <- whole_atg_phylop_plot.df %>%
  filter(!is.na(distance_bin), !is.na(facet_group)) %>%
  ggplot(aes(x = phylop, y = whole_atg_delta, colour = orf_definition_sub)) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed",
             linewidth = 0.35) +
  geom_point(alpha = 0.45, size = 0.85, stroke = 0) +
  phylop_label_layer(whole_atg_phylop_distance_structure_label.df, size = 1.9) +
  geom_text(
    data = whole_atg_phylop_distance_structure_corr.df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.15,
    size = 2.1,
    colour = "grey20"
  ) +
  facet_grid(facet_group ~ distance_bin) +
  scale_colour_manual(values = distance_structure_colors, name = "ORF class") +
  labs(
    x = "phyloP",
    y = "ΔTE (whole ATG deletion)",
    title = "phyloP vs whole-ATG deletion ΔTE by distance and uORF/uoORF structure"
  ) +
  phylop_scatter_theme +
  theme(legend.position = "bottom")

whole_atg_phylop_orf_type_structure_scatter.gg <- ggplot(
  whole_atg_phylop_plot.df,
  aes(x = phylop, y = whole_atg_delta, colour = orf_definition_sub)
) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed",
             linewidth = 0.35) +
  geom_point(alpha = 0.35, size = 0.6, stroke = 0) +
  phylop_label_layer(whole_atg_phylop_sub_label.df, size = 2.1) +
  geom_text(
    data = whole_atg_phylop_sub_corr.df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.15,
    size = 2.4,
    colour = "grey20"
  ) +
  facet_grid(facet_group ~ detected_group) +
  scale_colour_manual(values = PUB_SUB_COLORS, name = "ORF class") +
  labs(
    x = "phyloP",
    y = "ΔTE (whole ATG deletion)",
    title = "phyloP vs whole-ATG deletion ΔTE by ORF class and uORF/uoORF structure"
  ) +
  phylop_scatter_theme

pub_whole_atg_phylop.gg <- cowplot::plot_grid(
  whole_atg_phylop_scatter.gg,
  whole_atg_phylop_distance_scatter.gg,
  whole_atg_phylop_distance_structure_scatter.gg,
  whole_atg_phylop_orf_type_scatter.gg,
  whole_atg_phylop_orf_type_structure_scatter.gg,
  ncol = 1,
  rel_heights = c(1.0, 1.0, 1.25, 1.35, 1.45),
  align = "v",
  axis = "l"
)

format_mw_p <- function(p) {
  ifelse(is.na(p), "NA", format.pval(p, digits = 2, eps = 1e-3))
}

safe_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  median(x, na.rm = TRUE)
}

pairwise_detected_minus_undetected <- function(det_values, undet_values) {
  as.numeric(outer(det_values, undet_values, "-"))
}

bootstrap_hl_ci <- function(det_values, undet_values, n_boot) {
  boot_hl <- replicate(n_boot, {
    det_boot <- sample(det_values, length(det_values), replace = TRUE)
    undet_boot <- sample(undet_values, length(undet_values), replace = TRUE)
    safe_median(pairwise_detected_minus_undetected(det_boot, undet_boot))
  })
  stats::quantile(boot_hl, probs = c(0.025, 0.975),
                  na.rm = TRUE, names = FALSE)
}

whole_atg_mw_one_bin <- function(df) {
  det_values <- df$whole_atg_delta[df$detected_group == "Detected"]
  undet_values <- df$whole_atg_delta[df$detected_group == "Not detected"]
  n_detected <- length(det_values)
  n_undetected <- length(undet_values)
  median_detected <- safe_median(det_values)
  median_undetected <- safe_median(undet_values)

  base_row <- tibble(
    status = "ok",
    message = NA_character_,
    test = "Mann-Whitney/Wilcoxon rank-sum",
    n_detected = n_detected,
    n_undetected = n_undetected,
    median_detected = median_detected,
    median_undetected = median_undetected,
    median_difference = median_detected - median_undetected,
    hodges_lehmann_difference = NA_real_,
    hodges_lehmann_ci_low = NA_real_,
    hodges_lehmann_ci_high = NA_real_,
    bootstrap_reps = WHOLE_ATG_BOOTSTRAP_REPS,
    common_language_p_detected_gt_not_detected = NA_real_,
    wilcox_p = NA_real_
  )

  if (n_detected < 3 || n_undetected < 3) {
    return(base_row %>%
      mutate(
        status = "skipped",
        message = "Need at least 3 detected and 3 undetected ORFs."
      ))
  }

  wt <- tryCatch(
    wilcox.test(
      whole_atg_delta ~ detected_group,
      data = df,
      exact = FALSE,
      alternative = "two.sided"
    ),
    error = function(err) err
  )
  if (inherits(wt, "error")) {
    return(base_row %>%
      mutate(status = "error", message = conditionMessage(wt)))
  }

  if (WHOLE_ATG_BOOTSTRAP_REPS == 0) {
    return(base_row %>%
      mutate(
        message = "Hodges-Lehmann/bootstrap skipped; set RIBONN_WHOLE_ATG_BOOTSTRAP_REPS > 0 to enable.",
        wilcox_p = wt$p.value
      ))
  }

  pairwise_diffs <- pairwise_detected_minus_undetected(det_values, undet_values)
  hl_ci <- bootstrap_hl_ci(det_values, undet_values, WHOLE_ATG_BOOTSTRAP_REPS)

  base_row %>%
    mutate(
      hodges_lehmann_difference = safe_median(pairwise_diffs),
      hodges_lehmann_ci_low = hl_ci[[1]],
      hodges_lehmann_ci_high = hl_ci[[2]],
      common_language_p_detected_gt_not_detected = mean(pairwise_diffs > 0),
      wilcox_p = wt$p.value
    )
}

whole_atg_bin_mw_stats.df <- ism_whole.df %>%
  filter(!is.na(distance_bin)) %>%
  group_by(facet_group, distance_bin) %>%
  group_modify(~whole_atg_mw_one_bin(.x)) %>%
  ungroup() %>%
  mutate(
    p_adjust_method = "BH across tested whole-ATG distance bins",
    n_tests_bh = sum(!is.na(wilcox_p)),
    wilcox_p_adj_bh = p.adjust(wilcox_p, method = "BH"),
    mw_label = paste0("MW FDR=", format_mw_p(wilcox_p_adj_bh))
  )

print(whole_atg_bin_mw_stats.df)

whole_atg_bin_mw_labels.df <- whole_atg_bin_mw_stats.df %>%
  filter(status == "ok")

whole_atg_deletion_box.gg <- ism_whole.df %>%
  ggplot(aes(x = detected_group, y = whole_atg_delta,
             fill = orf_definition_sub, colour = orf_definition_sub)) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed",
             linewidth = 0.35) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5,
               colour = "grey30", linewidth = 0.35) +
  geom_jitter(width = 0.14, alpha = 0.35, size = 0.7, stroke = 0,
              show.legend = FALSE) +
  facet_wrap(~facet_group, ncol = 2) +
  scale_fill_manual(values   = PUB_SUB_COLORS, name = "ORF class") +
  scale_colour_manual(values = PUB_SUB_COLORS, name = "ORF class") +
  labs(x = NULL, y = "ΔTE (whole ATG deletion)") +
  pub_box_theme +
  guides(fill = "none", colour = "none")

whole_atg_deletion_dist_bins.gg <- ism_whole.df %>%
  ggplot(aes(x = detected_group, y = whole_atg_delta,
             fill = orf_definition_sub, colour = orf_definition_sub)) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed",
             linewidth = 0.35) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5,
               colour = "grey30", linewidth = 0.35) +
  geom_jitter(width = 0.12, alpha = 0.35, size = 0.55, stroke = 0,
              show.legend = FALSE) +
  stat_summary(
    fun.data = function(x) data.frame(y = Inf, label = paste0("n=", length(x))),
    geom = "text",
    vjust = 1.5, size = 2.0, colour = "grey30", show.legend = FALSE
  ) +
  geom_text(
    data = whole_atg_bin_mw_labels.df,
    aes(x = 1.5, y = Inf, label = mw_label),
    inherit.aes = FALSE,
    vjust = 3.3,
    size = 2.0,
    colour = "grey25"
  ) +
  facet_grid(facet_group ~ distance_bin) +
  scale_fill_manual(values = PUB_SUB_COLORS, name = "ORF class") +
  scale_colour_manual(values = PUB_SUB_COLORS, name = "ORF class") +
  labs(x = NULL, y = "ΔTE (whole ATG deletion)") +
  pub_box_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "top")

make_whole_atg_distance_plot <- function(structural_type) {
  structural_levels <- PUB_SUB_LEVELS[
    FACET_GROUP[PUB_SUB_LEVELS] == structural_type
  ]
  structural_colors <- PUB_SUB_COLORS[structural_levels]
  plot_df <- ism_whole.df %>%
    filter(facet_group == structural_type, !is.na(distance_bin)) %>%
    mutate(orf_definition_sub = fct_drop(orf_definition_sub))
  stat_df <- whole_atg_bin_mw_labels.df %>%
    filter(facet_group == structural_type)

  ggplot(plot_df, aes(x = detected_group, y = whole_atg_delta,
                      fill = orf_definition_sub, colour = orf_definition_sub)) +
    geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed",
               linewidth = 0.35) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5,
                 colour = "grey30", linewidth = 0.35) +
    geom_jitter(width = 0.12, alpha = 0.35, size = 0.55, stroke = 0,
                show.legend = FALSE) +
    stat_summary(
      fun.data = function(x) data.frame(y = Inf, label = paste0("n=", length(x))),
      geom = "text",
      vjust = 1.5, size = 2.0, colour = "grey30", show.legend = FALSE
    ) +
    geom_text(
      data = stat_df,
      aes(x = 1.5, y = Inf, label = mw_label),
      inherit.aes = FALSE,
      vjust = 3.3,
      size = 2.0,
      colour = "grey25"
    ) +
    facet_wrap(~distance_bin, nrow = 1) +
    scale_fill_manual(values = structural_colors, name = structural_type) +
    scale_colour_manual(values = structural_colors, name = structural_type) +
    labs(x = NULL, y = "ΔTE (whole ATG deletion)",
         title = structural_type) +
    pub_box_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "top")
}

pub_dist_bins_uorf.gg <- make_whole_atg_distance_plot("uORF")
pub_dist_bins_uo_orf.gg <- make_whole_atg_distance_plot("uoORF")
pub_dist_bins.gg <- cowplot::plot_grid(pub_dist_bins_uorf.gg,
                                       pub_dist_bins_uo_orf.gg,
                                       ncol = 2, align = "hv")

cat("Whole-ATG deletion mode: skipping reference-base mutation heatmaps.\n")
hm_a.gg <- NULL
box_a.gg <- NULL
pub_orf_class.gg <- NULL
pub_hm_sub.gg <- NULL
pub_hm_del.gg <- NULL
heatmap.gg <- NULL
pub_line.gg <- NULL
within_aug_staircase.gg <- NULL
within_aug_contrast.gg <- NULL
within_aug_interaction_coef.gg <- NULL

pub_whole_atg_deletion.gg <- cowplot::plot_grid(
  whole_atg_deletion_box.gg,
  pub_dist_bins.gg,
  ncol = 1,
  rel_heights = c(1.0, 1.2),
  align = "v",
  axis = "l"
)

} else {


# ============================================================
# Violin + boxplot (original)
# ============================================================

substitution.gg <- ggplot(ism.df, aes(x = orf_definition, y = orf_substitution_te_change_mean_3nt)) +
  geom_violin(aes(fill = orf_definition)) +
  geom_boxplot(fill = "white", alpha = 0.4, width = 0.1) +
  scale_fill_manual(values = orf_definition_colors) +
  theme_classic()

deletion.gg <- ggplot(ism.df, aes(x = orf_definition, y = orf_deletion_te_change_mean_3nt)) +
  geom_violin(aes(fill = orf_definition)) +
  geom_boxplot(fill = "white", alpha = 0.4, width = 0.1) +
  scale_fill_manual(values = orf_definition_colors) +
  theme_classic()


# ============================================================
# Stacked bar charts: direction of effect per ORF class
# ============================================================

direction_bar <- function(df, col, title) {
  df %>%
    count(orf_definition, direction = .data[[col]]) %>%
    group_by(orf_definition) %>%
    mutate(
      pct       = n / sum(n),
      direction = factor(direction, levels = c("increase", "unchanged", "decrease"))
    ) %>%
    ungroup() %>%
    ggplot(aes(x = orf_definition, y = pct, fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = direction_colors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(title = title, x = NULL, y = "Fraction of ORFs", fill = "Direction") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

sub_dir_bar.gg <- direction_bar(
  ism.df, "orf_substitution_direction_vs_wt",
  "Substitution effect on TE by ORF class"
)

del_dir_bar.gg <- direction_bar(
  ism.df, "orf_deletion_direction_vs_wt",
  "Deletion effect on TE by ORF class"
)


# ============================================================
# Scatter: substitution vs deletion TE change per ORF
# ============================================================

scatter.gg <- ggplot(ism.df,
                     aes(x = orf_substitution_te_change_mean_3nt,
                         y = orf_deletion_te_change_mean_3nt,
                         colour = orf_definition)) +
  geom_point(alpha = 0.25, size = 0.8, stroke = 0) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40",
              linewidth = 0.4) +
  scale_colour_manual(values = orf_definition_colors) +
  facet_wrap(~orf_definition, ncol = 2) +
  labs(
    x      = "Mean substitution ΔTE (ATG avg)",
    y      = "Mean deletion ΔTE (ATG avg)",
    title  = "Substitution vs deletion TE change per ORF start codon"
  ) +
  theme_classic() +
  guides(colour = "none")


# ============================================================
# Codon-position summary heatmaps (4 groups × 3 positions)
# ============================================================

summarise_by_codon_pos <- function(df, value_col) {
  df %>%
    group_by(orf_definition, pos_label) %>%
    summarise(mean_change = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(orf_definition = fct_rev(orf_definition))
}

codon_heatmap_theme <- list(
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                       midpoint = 0, name = "Mean ΔTE"),
  theme_classic(),
  theme(
    axis.line        = element_blank(),
    panel.border     = element_rect(fill = NA, colour = "grey80"),
    legend.position  = "right"
  )
)



# # ============================================================
# # Individual ORF meta-heatmap: one row per ORF, 3 codon positions
# # ============================================================
#
# orf_wide.df <- pos.df %>%
#   dplyr::select(orf_id, orf_definition, orf_position_in_start_codon,
#                 substitution_te_change_mean, deletion_te_change) %>%
#   pivot_wider(
#     names_from  = orf_position_in_start_codon,
#     values_from = c(substitution_te_change_mean, deletion_te_change),
#     names_prefix = "pos"
#   ) %>%
#   filter(complete.cases(.)) %>%
#   arrange(orf_definition, substitution_te_change_mean_pos1) %>%
#   mutate(orf_rank = row_number())
#
# clamp_lim <- quantile(
#   abs(c(orf_wide.df$substitution_te_change_mean_pos1,
#         orf_wide.df$substitution_te_change_mean_pos2,
#         orf_wide.df$substitution_te_change_mean_pos3)),
#   0.99, na.rm = TRUE
# )


# ============================================================
# Meta-position heatmap: ΔTE vs distance of ORF start from CDS
# Use pos-1 (A of ATG) offset as the representative ORF distance
# ============================================================

orf_dist.df <- pos.df %>%
  filter(orf_position_in_start_codon == 1) %>%           # A of ATG = ORF start
  mutate(dist_bin = floor(offset_from_cds_start / 25) * 25)   # 25-nt bins

meta_dist_sub <- orf_dist.df %>%
  group_by(orf_definition, dist_bin) %>%
  summarise(mean_change = mean(substitution_te_change_mean, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 5)

meta_dist_del <- orf_dist.df %>%
  group_by(orf_definition, dist_bin) %>%
  summarise(mean_change = mean(deletion_te_change, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 5)

meta_dist_theme <- list(
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                       midpoint = 0, name = "Mean ΔTE"),
  scale_x_continuous(expand = c(0, 0)),
  theme_classic(),
  theme(
    axis.line    = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey80")
  )
)



# ============================================================
# Trend line: mean ΔTE vs ORF start distance (smoothed)
# ============================================================

meta_trend_sub.gg <- orf_dist.df %>%
  ggplot(aes(x = offset_from_cds_start, y = substitution_te_change_mean,
             colour = orf_definition, fill = orf_definition)) +
  geom_smooth(method = "loess", span = 0.3, linewidth = 0.8, alpha = 0.15) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.4) +
  scale_colour_manual(values = orf_definition_colors) +
  scale_fill_manual(values = orf_definition_colors) +
  labs(
    x      = "ORF start offset from CDS (nt)",
    y      = "Mean substitution ΔTE",
    colour = "ORF class", fill = "ORF class",
    title  = "Substitution sensitivity across 5’UTR (by ORF start position)"
  ) +
  theme_classic()

meta_trend_del.gg <- orf_dist.df %>%
  ggplot(aes(x = offset_from_cds_start, y = deletion_te_change,
             colour = orf_definition, fill = orf_definition)) +
  geom_smooth(method = "loess", span = 0.3, linewidth = 0.8, alpha = 0.15) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.4) +
  scale_colour_manual(values = orf_definition_colors) +
  scale_fill_manual(values = orf_definition_colors) +
  labs(
    x      = "ORF start offset from CDS (nt)",
    y      = "Mean deletion ΔTE",
    colour = "ORF class", fill = "ORF class",
    title  = "Deletion sensitivity across 5’UTR (by ORF start position)"
  ) +
  theme_classic()


# ============================================================
# Sub-stratified heatmaps: theoreticals split by uORF vs uoORF structure
# ============================================================

# Join sub-stratified definition onto position data
pos_sub.df <- pos.df %>%
  dplyr::select(-orf_definition) %>%
  left_join(master_orf_def_sub.df, by = "orf_id") %>%
  filter(!is.na(orf_definition_sub)) %>%
  mutate(orf_definition_sub = factor(orf_definition_sub, levels = orf_definition_sub_levels))

# -- Codon-position heatmaps (6 groups × 3 ATG positions) --

summarise_by_codon_pos_sub <- function(df, value_col) {
  df %>%
    group_by(orf_definition_sub, pos_label) %>%
    summarise(mean_change = mean(.data[[value_col]], na.rm = TRUE),
              n           = n(),
              .groups = "drop") %>%
    mutate(orf_definition_sub = fct_rev(orf_definition_sub))
}

codon_sub_heatmap_theme <- list(
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                       midpoint = 0, name = "Mean ΔTE"),
  theme_classic(),
  theme(
    axis.line    = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey80")
  )
)

codon_sub_sub.gg <- summarise_by_codon_pos_sub(pos_sub.df, "substitution_te_change_mean") %>%
  ggplot(aes(x = pos_label, y = orf_definition_sub, fill = mean_change)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.4f", mean_change)), size = 3) +
  labs(x = "ATG codon position", y = NULL,
       title = "Mean substitution ΔTE per codon position\n(theoreticals sub-stratified by structure)") +
  codon_sub_heatmap_theme

codon_del_sub.gg <- summarise_by_codon_pos_sub(pos_sub.df, "deletion_te_change") %>%
  ggplot(aes(x = pos_label, y = orf_definition_sub, fill = mean_change)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.4f", mean_change)), size = 3) +
  labs(x = "ATG codon position", y = NULL,
       title = "Mean deletion ΔTE per codon position\n(theoreticals sub-stratified by structure)") +
  codon_sub_heatmap_theme

# -- Meta-distance heatmaps (6 groups × distance bins) --

orf_dist_sub.df <- pos_sub.df %>%
  filter(orf_position_in_start_codon == 1) %>%
  mutate(dist_bin = floor(offset_from_cds_start / 25) * 25)

meta_dist_sub_grid_sub <- orf_dist_sub.df %>%
  group_by(orf_definition_sub, dist_bin) %>%
  summarise(mean_change = mean(substitution_te_change_mean, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 5)

meta_dist_sub_grid_del <- orf_dist_sub.df %>%
  group_by(orf_definition_sub, dist_bin) %>%
  summarise(mean_change = mean(deletion_te_change, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 5)

meta_dist_sub_theme <- list(
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                       midpoint = 0, name = "Mean ΔTE"),
  scale_x_continuous(expand = c(0, 0)),
  theme_classic(),
  theme(
    axis.line    = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey80")
  )
)

meta_dist_sub_sub.gg <- ggplot(meta_dist_sub_grid_sub,
                               aes(x = dist_bin, y = fct_rev(orf_definition_sub), fill = mean_change)) +
  geom_tile() +
  labs(
    x     = "Distance of ORF start from CDS (nt, negative = upstream)",
    y     = NULL,
    title = "Mean substitution ΔTE vs ORF start distance\n(theoreticals sub-stratified by structure)"
  ) +
  meta_dist_sub_theme

meta_dist_sub_del.gg <- ggplot(meta_dist_sub_grid_del,
                               aes(x = dist_bin, y = fct_rev(orf_definition_sub), fill = mean_change)) +
  geom_tile() +
  labs(
    x     = "Distance of ORF start from CDS (nt, negative = upstream)",
    y     = NULL,
    title = "Mean deletion ΔTE vs ORF start distance\n(theoreticals sub-stratified by structure)"
  ) +
  meta_dist_sub_theme


# ============================================================
# Final publication plots
# ============================================================

suppressPackageStartupMessages(library(cowplot))

SHOW_LOW_EXPR <- FALSE

PUB_SUB_LEVELS <- if (SHOW_LOW_EXPR) orf_definition_sub_levels else
  grep("low expr", orf_definition_sub_levels, value = TRUE, invert = TRUE)

PUB_SUB_COLORS <- orf_definition_sub_colors[PUB_SUB_LEVELS]

DIST_CROP <- -500
BIN_WIDTH <- 50
MIN_N_BIN <- 15

FACET_GROUP <- c(
  "uORF"                         = "uORF",
  "uoORF"                        = "uoORF",
  "Not det. (filt.) · uORF"     = "uORF",
  "Not det. (filt.) · uoORF"    = "uoORF",
  "Not det. (low expr.) · uORF"  = "uORF",
  "Not det. (low expr.) · uoORF" = "uoORF"
)

add_facet_group <- function(df) {
  df %>%
    mutate(facet_group = factor(
      FACET_GROUP[as.character(orf_definition_sub)],
      levels = c("uORF", "uoORF")
    ))
}

pub_base_theme <- theme_classic(base_size = 11) +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(4, "mm"),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9),
    plot.title       = element_text(size = 11, face = "bold"),
    plot.subtitle    = element_text(size = 8,  colour = "grey40"),
    strip.background = element_blank(),
    strip.text       = element_text(size = 9, face = "bold")
  )


prepare_hm_sub <- function(df) {
  df %>%
    filter(dist_bin >= DIST_CROP) %>%
    { if (!SHOW_LOW_EXPR) filter(., !grepl("low expr", orf_definition_sub)) else . } %>%
    add_facet_group() %>%
    mutate(orf_definition_sub = fct_rev(
      factor(orf_definition_sub, levels = PUB_SUB_LEVELS)
    ))
}

pub_fill_scale <- scale_fill_gradient2(
  low = "#74669d", mid = "grey92", high = "#76baa6", midpoint = 0,
  name = "Mean ΔTE",
  guide = guide_colourbar(barheight = unit(22, "mm"), barwidth = unit(3, "mm"))
)

pub_hm_theme <- theme_classic(base_size = 11) +
  theme(
    axis.line        = element_blank(),
    panel.border     = element_rect(fill = NA, colour = "grey70", linewidth = 0.5),
    panel.spacing.y  = unit(6, "pt"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 9),
    axis.title       = element_text(size = 10),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    plot.title       = element_text(size = 11, face = "bold"),
    strip.background = element_blank(),
    strip.text.y     = element_text(size = 9, face = "bold", angle = 0, hjust = 0)
  )

hm_x_scale <- scale_x_continuous(
  breaks = seq(DIST_CROP, 0, by = 100),
  expand = expansion(add = c(12.5, 12.5))
)

pub_hm_sub.gg <- prepare_hm_sub(meta_dist_sub_grid_sub) %>%
  ggplot(aes(x = dist_bin, y = orf_definition_sub, fill = mean_change)) +
  geom_tile(height = 0.85, colour = "white", linewidth = 0.25) +
  facet_grid(facet_group ~ ., scales = "free_y", space = "free_y") +
  hm_x_scale +
  scale_y_discrete(expand = expansion(add = 0.4)) +
  labs(x = NULL, y = NULL, title = "Substitution  —  mean ΔTE") +
  pub_fill_scale + pub_hm_theme

pub_hm_del.gg <- prepare_hm_sub(meta_dist_sub_grid_del) %>%
  ggplot(aes(x = dist_bin, y = orf_definition_sub, fill = mean_change)) +
  geom_tile(height = 0.85, colour = "white", linewidth = 0.25) +
  facet_grid(facet_group ~ ., scales = "free_y", space = "free_y") +
  hm_x_scale +
  scale_y_discrete(expand = expansion(add = 0.4)) +
  labs(x = "ORF start offset from CDS (nt)", y = NULL,
       title = "Deletion  —  mean ΔTE") +
  pub_fill_scale + pub_hm_theme
heatmap.gg <- cowplot::plot_grid(pub_hm_sub.gg, pub_hm_del.gg, ncol = 1, align = "v")

bin_summary_sub <- function(df, value_col) {
  df %>%
    filter(offset_from_cds_start >= DIST_CROP) %>%
    { if (!SHOW_LOW_EXPR) filter(., !grepl("low expr", orf_definition_sub)) else . } %>%
    add_facet_group() %>%
    mutate(
      orf_definition_sub = factor(orf_definition_sub, levels = PUB_SUB_LEVELS),
      dist_bin = floor(offset_from_cds_start / BIN_WIDTH) * BIN_WIDTH + BIN_WIDTH / 2
    ) %>%
    group_by(facet_group, orf_definition_sub, dist_bin) %>%
    summarise(
      mean_val = mean(.data[[value_col]], na.rm = TRUE),
      se_val   = sd(.data[[value_col]], na.rm = TRUE) / sqrt(n()),
      n        = n(),
      .groups  = "drop"
    ) %>%
    filter(n >= MIN_N_BIN)
}

sub_line_sub.df <- bin_summary_sub(orf_dist_sub.df, "substitution_te_change_mean")
del_line_sub.df <- bin_summary_sub(orf_dist_sub.df, "deletion_te_change")

pub_line_layers <- list(
  geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed",
             linewidth = 0.35),
  geom_ribbon(aes(ymin = mean_val - 1.96 * se_val,
                  ymax = mean_val + 1.96 * se_val,
                  fill = orf_definition_sub),
              alpha = 0.15, colour = NA),
  geom_line(aes(colour = orf_definition_sub), linewidth = 0.75),
  facet_wrap(~facet_group, ncol = 2),
  scale_colour_manual(values = PUB_SUB_COLORS, name = "ORF class"),
  scale_fill_manual(values   = PUB_SUB_COLORS, name = "ORF class"),
  scale_x_continuous(
    breaks = seq(DIST_CROP, 0, by = 100),
    expand = expansion(mult = c(0.01, 0.02))
  ),
  coord_cartesian(xlim = c(DIST_CROP, 0)),
  pub_base_theme
)

pub_line_sub.gg <- ggplot(sub_line_sub.df, aes(x = dist_bin, y = mean_val)) +
  pub_line_layers +
  labs(x = NULL, y = "Mean substitution ΔTE", title = "Substitution")

pub_line_del.gg <- ggplot(del_line_sub.df, aes(x = dist_bin, y = mean_val)) +
  pub_line_layers +
  labs(x = "ORF start offset from CDS (nt)", y = "Mean deletion ΔTE",
       title = "Deletion")

pub_line.gg <- cowplot::plot_grid(pub_line_sub.gg, pub_line_del.gg, ncol = 1, align = "v")


# ============================================================
# ATG codon sensitivity summary plots
# ============================================================

ism_sub.df <- ism.df %>%
  left_join(master_orf_def_sub.df, by = "orf_id") %>%
  filter(!is.na(orf_definition_sub)) %>%
  { if (!SHOW_LOW_EXPR) filter(., !grepl("low expr", orf_definition_sub)) else . } %>%
  mutate(
    orf_definition_sub = factor(orf_definition_sub, levels = PUB_SUB_LEVELS),
    detected_group     = factor(
      if_else(orf_definition_sub %in% c("uORF", "uoORF"), "Detected", "Not detected"),
      levels = c("Detected", "Not detected")
    )
  ) %>%
  add_facet_group()

orf_start_dist.df <- pos.df %>%
  filter(orf_position_in_start_codon == 1) %>%
  mutate(
    orf_start_1based_num = as.numeric(orf_start_1based),
    orf_stop_1based_num  = as.numeric(orf_stop_1based),
    orf_length_nt = abs(orf_stop_1based_num - orf_start_1based_num) + 1
  ) %>%
  dplyr::select(orf_id, orf_start_offset = offset_from_cds_start, orf_length_nt)

add_dist_bin <- function(df, offset_col) {
  df %>% mutate(
    distance_bin = factor(
      case_when(
        .data[[offset_col]] >= -30  ~ "≤30 nt",
        .data[[offset_col]] >= -60  ~ "30–60 nt",
        .data[[offset_col]] >= -120 ~ "60–120 nt",
        TRUE                         ~ ">120 nt"
      ),
      levels = c("≤30 nt", "30–60 nt", "60–120 nt", ">120 nt")
    )
  )
}

ism_sub.df <- ism_sub.df %>%
  left_join(orf_start_dist.df, by = "orf_id") %>%
  add_dist_bin("orf_start_offset")

pos_sub_aug.df <- pos_sub.df %>%
  { if (!SHOW_LOW_EXPR) filter(., !grepl("low expr", orf_definition_sub)) else . } %>%
  mutate(
    orf_definition_sub = factor(orf_definition_sub, levels = PUB_SUB_LEVELS),
    detected_group     = factor(
      if_else(orf_definition_sub %in% c("uORF", "uoORF"), "Detected", "Not detected"),
      levels = c("Detected", "Not detected")
    )
  ) %>%
  add_facet_group() %>%
  left_join(orf_start_dist.df, by = "orf_id") %>%
  add_dist_bin("orf_start_offset")

pub_box_theme <- pub_base_theme +
  theme(
    axis.text.x        = element_text(angle = 30, hjust = 1, size = 8),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)
  )

pos_hm_theme <- theme_classic(base_size = 10) +
  theme(
    axis.line        = element_blank(),
    panel.border     = element_rect(fill = NA, colour = "grey70", linewidth = 0.4),
    panel.spacing    = unit(4, "pt"),
    axis.text        = element_text(size = 8),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    axis.title       = element_text(size = 9),
    strip.background = element_blank(),
    strip.text       = element_text(size = 8, face = "bold"),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7)
  )

codon_hm_fill <- scale_fill_gradient2(
  low = "#74669d", mid = "grey92", high = "#76baa6", midpoint = 0,
  name = "Mean ΔTE",
  guide = guide_colourbar(barheight = unit(16, "mm"), barwidth = unit(3, "mm"))
)

codon_pos_long_summary <- function(df, group_cols) {
  df %>%
    group_by(across(all_of(c(group_cols, "pos_label")))) %>%
    summarise(
      mean_sub = mean(substitution_te_change_mean, na.rm = TRUE),
      mean_del = mean(deletion_te_change,          na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    pivot_longer(cols = c(mean_sub, mean_del),
                 names_to = "metric", values_to = "mean_val") %>%
    mutate(metric = recode(metric,
                           mean_sub = "Substitution",
                           mean_del = "Deletion"))
}

mutation_matrix_long <- function(df) {
  required_cols <- c(
    "substitution_A_te_change",
    "substitution_C_te_change",
    "substitution_G_te_change",
    "substitution_T_te_change",
    "deletion_te_change"
  )
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Position summaries are missing per-alt-base ΔTE columns: ",
      paste(missing_cols, collapse = ", "),
      ". Rerun the mutagenesis summarization with the updated ",
      "summarize_shard.py before making mutation-matrix heatmaps."
    )
  }
  to_numeric <- function(x) suppressWarnings(as.numeric(x))

  df %>%
    mutate(
      ref_label = recode(ref_base, T = "U"),
      substitution_A_te_change = if_else(
        ref_base == "A", 0, to_numeric(substitution_A_te_change)
      ),
      substitution_C_te_change = if_else(
        ref_base == "C", 0, to_numeric(substitution_C_te_change)
      ),
      substitution_G_te_change = if_else(
        ref_base == "G", 0, to_numeric(substitution_G_te_change)
      ),
      substitution_T_te_change = if_else(
        ref_base == "T", 0, to_numeric(substitution_T_te_change)
      ),
      deletion_te_change = to_numeric(deletion_te_change)
    ) %>%
    pivot_longer(
      cols = all_of(required_cols),
      names_to = "mutation",
      values_to = "delta_te"
    ) %>%
    mutate(
      mutation = recode(
        mutation,
        substitution_A_te_change = "A",
        substitution_C_te_change = "C",
        substitution_G_te_change = "G",
        substitution_T_te_change = "U",
        deletion_te_change = "del"
      ),
      mutation = factor(mutation, levels = c("A", "C", "G", "U", "del")),
      ref_label = factor(ref_label, levels = c("A", "U", "G"))
    ) %>%
    group_by(distance_bin, detected_group, ref_label, mutation) %>%
    summarise(
      mean_val = mean(delta_te, na.rm = TRUE),
      n = sum(!is.na(delta_te)),
      .groups = "drop"
    )
}

orf_effect_long <- function(df) {
  df %>%
    pivot_longer(
      cols = c(orf_substitution_te_change_mean_3nt,
               orf_deletion_te_change_mean_3nt),
      names_to = "effect_type",
      values_to = "delta_te"
    ) %>%
    mutate(
      effect_type = recode(
        effect_type,
        orf_substitution_te_change_mean_3nt = "Substitution",
        orf_deletion_te_change_mean_3nt = "Deletion"
      )
    )
}


n_a.df <- count(ism_sub.df, facet_group, detected_group, orf_definition_sub,
                name = "n_orfs")


make_distance_bin_plot <- function(structural_type, heatmap_low, heatmap_high) {
  structural_levels <- PUB_SUB_LEVELS[
    FACET_GROUP[PUB_SUB_LEVELS] == structural_type
  ]
  structural_colors <- PUB_SUB_COLORS[structural_levels]
  box_df <- ism_sub.df %>%
    filter(facet_group == structural_type, !is.na(distance_bin)) %>%
    mutate(orf_definition_sub = fct_drop(orf_definition_sub))
  hm_df <- pos_sub_aug.df %>%
    filter(facet_group == structural_type, !is.na(distance_bin)) %>%
    mutate(orf_definition_sub = fct_drop(orf_definition_sub))

  box_plot <- box_df %>%
    orf_effect_long() %>%
    ggplot(aes(x = detected_group, y = delta_te,
               fill = orf_definition_sub, colour = orf_definition_sub)) +
    geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed",
               linewidth = 0.35) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5,
                 colour = "grey30", linewidth = 0.35) +
    geom_jitter(width = 0.12, alpha = 0.35, size = 0.55, stroke = 0,
                show.legend = FALSE) +
    stat_summary(
      fun.data = function(x) data.frame(y = Inf, label = paste0("n=", length(x))),
      geom = "text",
      vjust = 1.5, size = 2.0, colour = "grey30", show.legend = FALSE
    ) +
    facet_grid(effect_type ~ distance_bin) +
    scale_fill_manual(values = structural_colors, name = structural_type) +
    scale_colour_manual(values = structural_colors, name = structural_type) +
    labs(x = NULL, y = "ΔTE (ATG codon avg)",
         title = structural_type) +
    pub_box_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "top")

  heatmap_fill <- scale_fill_gradient2(
    low = heatmap_low, mid = "grey92", high = heatmap_high, midpoint = 0,
    name = "Mean ΔTE",
    guide = guide_colourbar(barheight = unit(14, "mm"), barwidth = unit(3, "mm"))
  )

  heatmap_plot <- mutation_matrix_long(hm_df) %>%
    ggplot(aes(x = mutation, y = fct_rev(ref_label), fill = mean_val)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", mean_val)), size = 2.2) +
    facet_grid(distance_bin ~ detected_group) +
    heatmap_fill +
    scale_x_discrete(expand = expansion(add = 0)) +
    scale_y_discrete(expand = expansion(add = 0.3)) +
    labs(x = NULL, y = "Reference\nbase") +
    pos_hm_theme

  cowplot::plot_grid(box_plot, heatmap_plot, ncol = 1,
                     rel_heights = c(3.0, 1.6), align = "v")
}

pub_dist_bins_uorf.gg <- make_distance_bin_plot(
  "uORF",
  heatmap_low = "#dceee9",
  heatmap_high = "#76baa6"
)
pub_dist_bins_uo_orf.gg <- make_distance_bin_plot(
  "uoORF",
  heatmap_low = "#e4e0ef",
  heatmap_high = "#74669d"
)
pub_dist_bins.gg <- cowplot::plot_grid(pub_dist_bins_uorf.gg,
                                       pub_dist_bins_uo_orf.gg,
                                       ncol = 2, align = "hv")

}



#################

if (whole_atg_mode) {
  whole_atg_plot_dir_env <- trimws(Sys.getenv("RIBONN_ORF_PLOT_DIR", ""))
  whole_atg_plot.dir <- if (nzchar(whole_atg_plot_dir_env)) {
    whole_atg_plot_dir_env
  } else {
    file.path(ism_results.dir, "plots")
  }
  dir.create(whole_atg_plot.dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(
    whole_atg_phylop_corr_stats.df,
    file.path(whole_atg_plot.dir, "whole_atg_phylop_correlation_stats.tsv"),
    sep = "\t"
  )

  ggsave(
    file.path(whole_atg_plot.dir, "whole_atg_phylop_scatter.png"),
    whole_atg_phylop_scatter.gg,
    width = 5.5, height = 4, dpi = 300
  )
  ggsave(
    file.path(whole_atg_plot.dir, "whole_atg_phylop_by_cds_distance.png"),
    whole_atg_phylop_distance_scatter.gg,
    width = 8, height = 3.5, dpi = 300
  )
  ggsave(
    file.path(whole_atg_plot.dir, "whole_atg_phylop_by_cds_distance_and_orf_type.png"),
    whole_atg_phylop_distance_structure_scatter.gg,
    width = 9, height = 5.5, dpi = 300
  )
  ggsave(
    file.path(whole_atg_plot.dir, "whole_atg_phylop_by_orf_type.png"),
    whole_atg_phylop_orf_type_scatter.gg,
    width = 7, height = 5.5, dpi = 300
  )
  ggsave(
    file.path(whole_atg_plot.dir, "whole_atg_phylop_by_orf_type_and_structure.png"),
    whole_atg_phylop_orf_type_structure_scatter.gg,
    width = 7, height = 5.5, dpi = 300
  )
  ggsave(
    file.path(whole_atg_plot.dir, "whole_atg_phylop_all_panels.png"),
    pub_whole_atg_phylop.gg,
    width = 8, height = 17, dpi = 300
  )

  cat("Saved whole-ATG phyloP plots to:", whole_atg_plot.dir, "\n")

  # print(whole_atg_phylop_scatter.gg)
  print(whole_atg_phylop_distance_scatter.gg)
  print(whole_atg_phylop_distance_structure_scatter.gg)
  print(whole_atg_phylop_orf_type_scatter.gg)
  print(whole_atg_phylop_orf_type_structure_scatter.gg)

}
