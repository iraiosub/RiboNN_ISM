#!/usr/bin/env/Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

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

species <- Sys.getenv("RIBONN_ORF_SPECIES", "human")
if (!species %in% c("human", "mouse")) {
  stop("RIBONN_ORF_SPECIES must be 'human' or 'mouse'")
}
default_te_label <- if (species == "mouse") "mean_predicted_TE" else "normal_brain_tissue"
te_label <- Sys.getenv("RIBONN_ORF_TE_LABEL", default_te_label)
ism_results.dir <- Sys.getenv(
  "RIBONN_ORF_RESULTS_DIR",
  file.path(output_root.dir, paste0(species, "_orf_starts_", te_label))
)
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
master_table <- Sys.getenv("RIBONN_ORF_MASTER_TABLE", default_master_tables[[species]])

cat("Species:", species, "\n")
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

# In orf_start_codon mode every row is one of the 3 ATG positions
pos.df <- pos.df %>%
  tidyr::unite("orf_id", transcript_id, orf_start_1based, orf_stop_1based, orf_frame,
               sep = "_", remove = FALSE) %>%
  left_join(master_orf_def.df, by = "orf_id") %>%
  filter(!is.na(orf_definition)) %>%
  mutate(
    orf_definition   = factor(orf_definition, levels = orf_definition_levels),
    pos_label        = factor(codon_pos_labels[orf_position_in_start_codon],
                              levels = codon_pos_labels)
  )

cat("Position-level rows after joining master table:", nrow(pos.df), "\n")


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

codon_sub.gg <- summarise_by_codon_pos(pos.df, "substitution_te_change_mean") %>%
  ggplot(aes(x = pos_label, y = orf_definition, fill = mean_change)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.4f", mean_change)), size = 3) +
  labs(x = "ATG codon position", y = NULL,
       title = "Mean substitution ΔTE per codon position and ORF class") +
  codon_heatmap_theme

codon_del.gg <- summarise_by_codon_pos(pos.df, "deletion_te_change") %>%
  ggplot(aes(x = pos_label, y = orf_definition, fill = mean_change)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.4f", mean_change)), size = 3) +
  labs(x = "ATG codon position", y = NULL,
       title = "Mean deletion ΔTE per codon position and ORF class") +
  codon_heatmap_theme


# ============================================================
# Codon-position violin: per-position distribution per group
# ============================================================

codon_pos_violin_sub.gg <- ggplot(pos.df,
                                  aes(x = pos_label, y = substitution_te_change_mean, fill = orf_definition)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  geom_boxplot(fill = "white", alpha = 0.5, width = 0.08, outlier.shape = NA) +
  scale_fill_manual(values = orf_definition_colors) +
  facet_wrap(~orf_definition, ncol = 2) +
  labs(x = "ATG codon position", y = "Substitution ΔTE",
       title = "Substitution TE change per codon position") +
  theme_classic() +
  guides(fill = "none")

codon_pos_violin_del.gg <- ggplot(pos.df,
                                  aes(x = pos_label, y = deletion_te_change, fill = orf_definition)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  geom_boxplot(fill = "white", alpha = 0.5, width = 0.08, outlier.shape = NA) +
  scale_fill_manual(values = orf_definition_colors) +
  facet_wrap(~orf_definition, ncol = 2) +
  labs(x = "ATG codon position", y = "Deletion ΔTE",
       title = "Deletion TE change per codon position") +
  theme_classic() +
  guides(fill = "none")


# ============================================================
# Individual ORF meta-heatmap: one row per ORF, 3 codon positions
# ============================================================

orf_wide.df <- pos.df %>%
  dplyr::select(orf_id, orf_definition, orf_position_in_start_codon,
                substitution_te_change_mean, deletion_te_change) %>%
  pivot_wider(
    names_from  = orf_position_in_start_codon,
    values_from = c(substitution_te_change_mean, deletion_te_change),
    names_prefix = "pos"
  ) %>%
  filter(complete.cases(.)) %>%
  arrange(orf_definition, substitution_te_change_mean_pos1) %>%
  mutate(orf_rank = row_number())

clamp_lim <- quantile(
  abs(c(orf_wide.df$substitution_te_change_mean_pos1,
        orf_wide.df$substitution_te_change_mean_pos2,
        orf_wide.df$substitution_te_change_mean_pos3)),
  0.99, na.rm = TRUE
)

orf_long.df <- orf_wide.df %>%
  pivot_longer(
    cols         = matches("^(substitution_te_change_mean|deletion_te_change)_pos\\d"),
    names_to     = c("metric", "position"),
    names_pattern = "(substitution_te_change_mean|deletion_te_change)_pos(\\d+)"
  ) %>%
  mutate(
    metric        = recode(metric,
                           "substitution_te_change_mean" = "Substitution",
                           "deletion_te_change"          = "Deletion"),
    pos_label     = factor(codon_pos_labels[as.integer(position)],
                           levels = codon_pos_labels),
    value_clamped = pmin(pmax(value, -clamp_lim), clamp_lim)
  )

meta_heatmap.gg <- ggplot(orf_long.df,
                          aes(x = pos_label, y = orf_rank, fill = value_clamped)) +
  geom_raster() +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                       midpoint = 0, name = "ΔTE",
                       limits = c(-clamp_lim, clamp_lim)) +
  facet_grid(orf_definition ~ metric, scales = "free_y", space = "free_y") +
  labs(
    x     = "ATG codon position",
    y     = "ORF (ranked by pos-1 substitution effect)",
    title = "Per-ORF ΔTE across ATG codon positions"
  ) +
  theme_classic() +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    strip.background = element_blank(),
    panel.border     = element_rect(fill = NA, colour = "grey80"),
    panel.spacing.x  = unit(4, "pt")
  )


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

meta_dist_sub.gg <- ggplot(meta_dist_sub,
                           aes(x = dist_bin, y = fct_rev(orf_definition), fill = mean_change)) +
  geom_tile() +
  labs(
    x     = "Distance of ORF start from CDS (nt, negative = upstream)",
    y     = NULL,
    title = "Mean substitution ΔTE vs ORF start distance from CDS"
  ) +
  meta_dist_theme

meta_dist_del.gg <- ggplot(meta_dist_del,
                           aes(x = dist_bin, y = fct_rev(orf_definition), fill = mean_change)) +
  geom_tile() +
  labs(
    x     = "Distance of ORF start from CDS (nt, negative = upstream)",
    y     = NULL,
    title = "Mean deletion ΔTE vs ORF start distance from CDS"
  ) +
  meta_dist_theme


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
