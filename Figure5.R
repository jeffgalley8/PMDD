# ============================================================
# PICRUSt2 / MetaCyc CYTOKINE HEATMAPS — STANDALONE FIGURE SCRIPT
# ============================================================
# Main-effects heatmap: IL-1beta, IL-6, IL-8 (no TNFalpha main effect).
# Interaction heatmap: IL-1beta, IL-6, IL-8, TNFalpha x PMDD.
# Significant associations are colored; robust associations have asterisks.

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ------------------------------------------------------------
# 1. Read existing ANCOM-BC2 results
# ------------------------------------------------------------
main_file <- "PICRUSt_cytokine_main_ALL.csv"
interaction_file <- "PICRUSt_cytokine_PMDD_interaction_ALL.csv"
metacyc_file <- "MetaCyc-pathway.txt"

needed_files <- c(main_file, interaction_file, metacyc_file)
if (any(!file.exists(needed_files))) {
  stop("Missing required file(s): ",
       paste(needed_files[!file.exists(needed_files)], collapse = ", "))
}

func_cyto_main_results <- read_csv(main_file, show_col_types = FALSE)
func_cyto_int_results <- read_csv(interaction_file, show_col_types = FALSE)

required_columns <- c("pathway", "cytokine", "lfc", "significant", "robust")
stopifnot(all(required_columns %in% names(func_cyto_main_results)),
          all(required_columns %in% names(func_cyto_int_results)))

# ------------------------------------------------------------
# 2. Combine all FDR-significant main and interaction hits
# ------------------------------------------------------------
func_main_sig <- func_cyto_main_results |>
  filter(significant == TRUE) |>
  mutate(analysis = "Main", effect = cytokine)

func_int_sig <- func_cyto_int_results |>
  filter(significant == TRUE) |>
  mutate(analysis = "Interaction", effect = paste0(cytokine, " × PMDD"))

func_sig <- bind_rows(func_main_sig, func_int_sig)

# Exact final order: omit TNFalpha MAIN effect, retain its interaction.
effect_order <- c(
  "IL-1β", "IL-6", "IL-8",
  "IL-1β × PMDD", "IL-6 × PMDD", "IL-8 × PMDD", "TNFα × PMDD"
)

# Full pathway x effect grid; nonsignificant cells remain NA.
heat_func_sig <- tidyr::expand_grid(
  pathway = unique(func_sig$pathway),
  effect = effect_order
) |>
  left_join(
    func_sig |>
      filter(effect %in% effect_order) |>
      select(pathway, effect, lfc, robust),
    by = c("pathway", "effect")
  )

if (anyDuplicated(heat_func_sig[c("pathway", "effect")])) {
  stop("Duplicate pathway/effect combinations after joining results. Check input CSVs.")
}

# ------------------------------------------------------------
# 3. MetaCyc ID -> pathway-name annotations
# ------------------------------------------------------------
# The original mapping file is a tab-delimited, two-column file, no header.
metacyc_map <- read.delim(
  metacyc_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)
if (ncol(metacyc_map) != 2L) {
  stop("Expected exactly two tab-separated columns in ", metacyc_file)
}
colnames(metacyc_map) <- c("pathway", "pathway_name")
metacyc_map <- metacyc_map |>
  mutate(across(c(pathway, pathway_name), trimws)) |>
  filter(!is.na(pathway), pathway != "") |>
  distinct(pathway, .keep_all = TRUE)

sig_pathways <- unique(func_sig$pathway)
unmapped <- setdiff(sig_pathways, metacyc_map$pathway)
cat("Significant unique pathways:", length(sig_pathways), "\n")
cat("Annotated in MetaCyc map:", length(sig_pathways) - length(unmapped), "\n")
if (length(unmapped) > 0) {
  warning("Unmapped MetaCyc IDs (displaying ID as fallback): ",
          paste(unmapped, collapse = ", "))
}

heat_func_sig <- heat_func_sig |>
  left_join(metacyc_map, by = "pathway") |>
  mutate(
    pathway_label = ifelse(
      is.na(pathway_name) | pathway_name == "",
      pathway,
      paste0(pathway_name, " (", pathway, ")")
    )
  )

# Final label-cleaning function from the figure workflow.
clean_pathway_label <- function(x) {
  x |>
    stringr::str_remove("\\s*\\([^()]*PWY[^()]*\\)\\s*$") |>
    stringr::str_remove(
      stringr::regex("^superpathway of\\s+", ignore_case = TRUE)
    ) |>
    stringr::str_wrap(width = 34)
}

label_lookup <- heat_func_sig |>
  distinct(pathway, pathway_label) |>
  mutate(pathway_label_short = clean_pathway_label(pathway_label))

# ------------------------------------------------------------
# 4. Separate main and interaction pathways
# ------------------------------------------------------------
heat_func_split <- heat_func_sig |>
  mutate(
    effect_type = ifelse(
      grepl("PMDD", effect, ignore.case = TRUE),
      "Interaction", "Main"
    )
  )

main_hit_ids <- heat_func_split |>
  filter(effect_type == "Main", !is.na(lfc)) |>
  distinct(pathway) |>
  pull(pathway)

interaction_hit_ids <- heat_func_split |>
  filter(effect_type == "Interaction", !is.na(lfc)) |>
  distinct(pathway) |>
  pull(pathway)

heat_main <- heat_func_split |>
  filter(effect_type == "Main", pathway %in% main_hit_ids)

heat_interaction <- heat_func_split |>
  filter(effect_type == "Interaction", pathway %in% interaction_hit_ids)

# ------------------------------------------------------------
# 5. LFC matrices and independent Ward.D2 row clustering
# ------------------------------------------------------------
main_mat <- heat_main |>
  select(pathway, effect, lfc) |>
  pivot_wider(names_from = effect, values_from = lfc) |>
  column_to_rownames("pathway") |>
  as.matrix()
main_mat <- main_mat[, effect_order[1:3], drop = FALSE]

main_cluster <- main_mat
main_cluster[is.na(main_cluster)] <- 0
hc_main <- hclust(dist(main_cluster), method = "ward.D2")

int_mat <- heat_interaction |>
  select(pathway, effect, lfc) |>
  pivot_wider(names_from = effect, values_from = lfc) |>
  column_to_rownames("pathway") |>
  as.matrix()
int_mat <- int_mat[, effect_order[4:7], drop = FALSE]

int_cluster <- int_mat
int_cluster[is.na(int_cluster)] <- 0
hc_int <- hclust(dist(int_cluster), method = "ward.D2")

# ------------------------------------------------------------
# 6. Display labels and robust-asterisk matrices
# ------------------------------------------------------------
main_labels <- label_lookup$pathway_label_short[
  match(rownames(main_mat), label_lookup$pathway)
]

robust_main <- heat_main |>
  select(pathway, effect, robust) |>
  distinct() |>
  pivot_wider(names_from = effect, values_from = robust, values_fill = FALSE) |>
  column_to_rownames("pathway") |>
  as.matrix()
robust_main <- robust_main[rownames(main_mat), colnames(main_mat), drop = FALSE]

int_labels <- label_lookup$pathway_label_short[
  match(rownames(int_mat), label_lookup$pathway)
]

robust_int <- heat_interaction |>
  select(pathway, effect, robust) |>
  distinct() |>
  pivot_wider(names_from = effect, values_from = robust, values_fill = FALSE) |>
  column_to_rownames("pathway") |>
  as.matrix()
robust_int <- robust_int[rownames(int_mat), colnames(int_mat), drop = FALSE]

# ------------------------------------------------------------
# 7. Final ComplexHeatmap definitions
# ------------------------------------------------------------
col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

ht_main_effects <- Heatmap(
  main_mat,
  name = "LFC", col = col_fun, na_col = "grey95",
  cluster_rows = as.dendrogram(hc_main), cluster_columns = FALSE,
  show_row_dend = TRUE, row_dend_side = "left",
  row_dend_width = unit(18, "mm"),
  row_labels = main_labels, row_names_gp = gpar(fontsize = 7),
  column_title = "Cytokine Main Effects",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  column_names_rot = 45, column_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "grey90", lwd = 0.4),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(robust_main[i, j]) && robust_main[i, j]) {
      grid.points(x, y, pch = 8, size = unit(3.5, "mm"))
    }
  }
)

int_col_labels <- gsub(" × PMDD", "", colnames(int_mat), fixed = TRUE)

ht_interactions <- Heatmap(
  int_mat,
  name = "LFC", col = col_fun, na_col = "grey95",
  cluster_rows = as.dendrogram(hc_int), cluster_columns = FALSE,
  show_row_dend = TRUE, row_dend_side = "left",
  row_dend_width = unit(18, "mm"),
  row_labels = int_labels, row_names_gp = gpar(fontsize = 7),
  column_title = "PMDD Interaction",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  column_labels = int_col_labels,
  column_names_rot = 45, column_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "grey90", lwd = 0.4),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(robust_int[i, j]) && robust_int[i, j]) {
      grid.points(x, y, pch = 8, size = unit(3.5, "mm"))
    }
  }
)

# ------------------------------------------------------------
# 8. Final exports — exact filenames, sizes, resolution
# ------------------------------------------------------------
pdf("PICRUSt_Cytokine_Main_Effects.pdf", width = 8.5, height = 13)
draw(ht_main_effects)
dev.off()

png("PICRUSt_Cytokine_Main_Effects.png",
    width = 8.5, height = 13, units = "in", res = 600)
draw(ht_main_effects)
dev.off()

pdf("PICRUSt_PMDD_Interactions.pdf", width = 8.5, height = 12)
draw(ht_interactions)
dev.off()

png("PICRUSt_PMDD_Interactions.png",
    width = 8.5, height = 12, units = "in", res = 600)
draw(ht_interactions)
dev.off()