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


# Load ISM results, ORF-level
ism_results.dir <- "/Volumes/lab-ulej/home/shared/oscar_ira_riboloco/RiboNN_ISM/all_utr5_mutagenesis/output/human_orf_starts/summaries"
ism.ls <- list.files(ism_results.dir, pattern = ".orfs.tsv.gz", full.names = T)
ism.df <- rbindlist(lapply(ism.ls, fread))


# Positions also available: shard_051.positions.tsv.gz

# Load human brain mastertable of ORFs
human.master.df <- fread("/Volumes/lab-ulej/home/shared/oscar_ira_riboloco/analysis_results/human_brain.unmixing.master_table.with_below_tpm_threshold.tsv.gz")


human.master.df <- human.master.df %>%
  mutate(orf_definition = case_when(orf_label == "Upstream ORF" ~ "uORF",
                                    orf_label == "Upstream overlapping ORF" ~ "uoORF",
                                    orf_label == "5'UTR theoretical (AUG + in-frame stop)" & tpm_threshold_group == "above_threshold" ~ "Not detected (expression-filtered)",
                                    orf_label == "5'UTR theoretical (AUG + in-frame stop)" & tpm_threshold_group == "below_threshold" ~ "Not detected (low expression)",
                                    TRUE ~ "suspicious"))

# Only keep tx in the master-table
ism.df <- semi_join(ism.df, human.master.df, by = "transcript_id")
nrow(ism.df)

# Derive orf_id for the ISM table
ism.df <- ism.df %>%
  tidyr::unite("orf_id", transcript_id, orf_start_1based, orf_stop_1based, orf_frame, sep = "_")

ism.df <- left_join(ism.df, dplyr::select(human.master.df, orf_id, orf_definition), by = "orf_id")

ism.df <- ism.df %>%
  dplyr::filter(!is.na(orf_definition))

head(ism.df)

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
















