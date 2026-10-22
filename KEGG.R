############################################################
# PMDD PICRUSt2 KO-derived KEGG pathway analysis
# Cytokine x PMDD interactions (luteal samples only)
#
# Inputs:
#   picrust_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz
#   metadata.tsv
############################################################

# ==========================================================
# 1. PACKAGES AND PATHS
# ==========================================================
library(ggpicrust2)
library(ANCOMBC)
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)

ko_file <- "picrust_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"
metadata_file <- "metadata.tsv"
outdir <- "KEGG_cytokine_interactions"

# ==========================================================
# 2. METADATA AND KO-TO-KEGG PATHWAY CONVERSION
# ==========================================================
md <- as.data.frame(readr::read_tsv(metadata_file, show_col_types = FALSE))
required <- c("sample-id", "PMDD", "Phase", "blood_il1", "blood_il6",
              "blood_il8", "blood_tnfa")
missing_cols <- setdiff(required, names(md))
if (length(missing_cols)) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
if (anyNA(md$`sample-id`) || anyDuplicated(md$`sample-id`)) {
  stop("Metadata sample IDs must be nonmissing and unique.")
}
rownames(md) <- md$`sample-id`

# Reference group = Control (C); comparison group = PMDD (P).
md$PMDD <- factor(as.character(md$PMDD), levels = c("C", "P"))
md$Phase <- factor(as.character(md$Phase), levels = c("F", "L"))
if (anyNA(md$PMDD) || anyNA(md$Phase)) {
  stop("PMDD must use C/P and Phase must use F/L in metadata.tsv.")
}

kegg_abundance <- as.data.frame(ggpicrust2::ko2kegg_abundance(ko_file))
if (is.null(rownames(kegg_abundance)) || anyDuplicated(rownames(kegg_abundance))) {
  stop("KEGG pathway IDs must be present and unique.")
}
missing_samples <- setdiff(rownames(md), colnames(kegg_abundance))
if (length(missing_samples)) {
  stop("Samples absent from KO-derived KEGG matrix: ", paste(missing_samples, collapse = ", "))
}
kegg_pmdd <- kegg_abundance[, rownames(md), drop = FALSE]
stopifnot(identical(colnames(kegg_pmdd), rownames(md)))

# ==========================================================
# 3. LUTEAL CYTOKINES AND PHYLOSEQ OBJECT
# ==========================================================
md_luteal <- md[md$Phase == "L", , drop = FALSE]
if (!nrow(md_luteal)) stop("No luteal samples found.")
cytokines <- c(il1_z = "blood_il1", il6_z = "blood_il6",
               il8_z = "blood_il8", tnfa_z = "blood_tnfa")
cytokine_labels <- c(il1_z = "IL-1β", il6_z = "IL-6",
                     il8_z = "IL-8", tnfa_z = "TNFα")
for (cyt in names(cytokines)) {
  values <- suppressWarnings(as.numeric(as.character(md_luteal[[cytokines[[cyt]]]])))
  sd_value <- sd(values, na.rm = TRUE)
  if (!is.finite(sd_value) || sd_value == 0) {
    stop("Cannot standardize cytokine ", cytokines[[cyt]], ": insufficient variation.")
  }
  md_luteal[[cyt]] <- as.numeric(scale(values))
}

kegg_luteal <- as.matrix(kegg_pmdd[, rownames(md_luteal), drop = FALSE])
storage.mode(kegg_luteal) <- "numeric"
if (anyNA(kegg_luteal) || any(!is.finite(kegg_luteal)) || any(kegg_luteal < 0)) {
  stop("KEGG abundance matrix contains missing, non-finite, or negative values.")
}
ps_kegg <- phyloseq(
  otu_table(kegg_luteal, taxa_are_rows = TRUE),
  sample_data(md_luteal)
)

# ==========================================================
# 4. KEGG PATHWAY NAMES AND BIOLOGICAL CLASSES
# ==========================================================
# ggpicrust2 may query KEGG for annotations; internet access may be needed.
annotation_input <- tibble(feature = rownames(kegg_luteal), p_adjust = 0)
annotation_raw <- ggpicrust2::pathway_annotation(
  pathway = "KO", daa_results_df = annotation_input, ko_to_kegg = TRUE
)
annotation <- annotation_raw |>
  dplyr::select(pathway = feature, pathway_name, pathway_class) |>
  dplyr::distinct()
if (anyDuplicated(annotation$pathway)) {
  stop("KEGG annotation has multiple rows for a pathway; resolve before joining.")
}
write.csv(annotation, file.path(outdir, "KEGG_pathway_annotations.csv"), row.names = FALSE)

# ==========================================================
# 5. CYTOKINE x PMDD ANCOM-BC2 INTERACTION MODELS
# ==========================================================
# PMDD reference = C (Control). No random effect: one luteal sample/participant.
run_interaction <- function(cyt) {
  md_ps <- as.data.frame(sample_data(ps_kegg))
  keep <- rownames(md_ps)[!is.na(md_ps[[cyt]])]
  ps_sub <- prune_samples(keep, ps_kegg)
  groups <- table(as.data.frame(sample_data(ps_sub))$PMDD)
  if (length(groups) != 2 || any(groups == 0)) {
    stop("Both PMDD groups must be represented for ", cyt)
  }
  message("Running ", cyt, " interaction model (n=", nsamples(ps_sub), ")")
  fit <- ANCOMBC::ancombc2(
    data = ps_sub,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = paste0(cyt, " * PMDD"),
    rand_formula = NULL,
    p_adj_method = "BH",
    prv_cut = 0.10,
    lib_cut = 0,
    s0_perc = 0.05,
    group = "PMDD",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 1,
    verbose = TRUE,
    global = FALSE
  )
  result <- list(fit = fit, n = nsamples(ps_sub))
  saveRDS(result, file.path(outdir, paste0("fit_interaction_", cyt, ".rds")))
  result
}
interaction_fits <- setNames(lapply(names(cytokines), run_interaction), names(cytokines))

# ==========================================================
# 6. EXTRACT INTERACTION RESULTS AND GROUP-SPECIFIC SLOPES
# ==========================================================
extract_interaction <- function(obj, cyt) {
  res <- obj$fit$res
  candidates <- grep(paste0("^lfc_", cyt, ":PMDD"), names(res), value = TRUE)
  if (length(candidates) != 1) {
    stop("Expected one interaction coefficient for ", cyt,
         "; found: ", paste(candidates, collapse = ", "))
  }
  term <- sub("^lfc_", "", candidates)
  required_cols <- c("taxon", paste0(c("lfc_", "se_", "W_", "p_", "q_",
                                     "diff_", "diff_robust_"), term),
                     paste0("lfc_", cyt))
  missing <- setdiff(required_cols, names(res))
  if (length(missing)) {
    stop("Missing ANCOM-BC2 result columns for ", cyt, ": ",
         paste(missing, collapse = ", "))
  }
  control_slope <- res[[paste0("lfc_", cyt)]]
  interaction_lfc <- res[[paste0("lfc_", term)]]
  pmdd_slope <- control_slope + interaction_lfc
  tibble(
    cytokine = cytokine_labels[[cyt]],
    n = obj$n,
    pathway = res$taxon,
    lfc = interaction_lfc,
    se = res[[paste0("se_", term)]],
    W = res[[paste0("W_", term)]],
    p_value = res[[paste0("p_", term)]],
    q_value = res[[paste0("q_", term)]],
    significant = res[[paste0("diff_", term)]],
    robust = res[[paste0("diff_robust_", term)]],
    control_slope = control_slope,
    pmdd_slope = pmdd_slope
  )
}
interaction_results <- dplyr::bind_rows(lapply(names(cytokines), function(cyt) {
  extract_interaction(interaction_fits[[cyt]], cyt)
})) |>
  dplyr::left_join(annotation, by = "pathway")

# ==========================================================
# 7. SIGNIFICANT INTERACTIONS AND SUMMARY TABLES
# ==========================================================
significant_interactions <- interaction_results |>
  dplyr::filter(!is.na(significant) & significant) |>
  dplyr::mutate(
    pattern = dplyr::case_when(
      is.na(control_slope) | is.na(pmdd_slope) ~ NA_character_,
      sign(control_slope) != sign(pmdd_slope) ~ "Opposite direction",
      abs(pmdd_slope) > abs(control_slope) ~ "Stronger in PMDD",
      abs(pmdd_slope) < abs(control_slope) ~ "Weaker in PMDD",
      TRUE ~ "No magnitude change"
    ),
    control_direction = dplyr::case_when(
      control_slope > 0 ~ "Positive", control_slope < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    pmdd_direction = dplyr::case_when(
      pmdd_slope > 0 ~ "Positive", pmdd_slope < 0 ~ "Negative",
      TRUE ~ "Zero"
    )
  ) |>
  dplyr::arrange(cytokine, q_value)

hit_counts <- interaction_results |>
  dplyr::group_by(cytokine) |>
  dplyr::summarise(n_interactions = sum(significant, na.rm = TRUE), .groups = "drop")
pattern_summary <- significant_interactions |>
  dplyr::count(cytokine, pattern, name = "n") |>
  dplyr::group_by(cytokine) |>
  dplyr::mutate(total = sum(n), percent = round(100 * n / total, 1)) |>
  dplyr::ungroup()
class_summary <- significant_interactions |>
  dplyr::count(cytokine, pathway_class, sort = TRUE, name = "n")

# ==========================================================
# 8. EXPORT RESULTS AND SESSION INFORMATION
# ==========================================================
write.csv(interaction_results,
          file.path(outdir, "KEGG_cytokine_PMDD_interactions_ALL.csv"), row.names = FALSE)
write.csv(significant_interactions,
          file.path(outdir, "KEGG_cytokine_PMDD_interactions_SIGNIFICANT.csv"), row.names = FALSE)
saveRDS(interaction_fits, file.path(outdir, "KEGG_interaction_fits.rds"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))