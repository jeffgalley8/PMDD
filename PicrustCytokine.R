
# ============================================================
# PICRUSt2 CYTOKINE ANALYSIS
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(phyloseq)
library(ANCOMBC)

# ------------------------------------------------------------
# 1. Import PICRUSt2 pathways and metadata
# ------------------------------------------------------------

PICRUST_TSV <- "pathways_out/path_abun_unstrat.tsv"
META_TSV <- "metadata.tsv"

pic <- read_tsv(PICRUST_TSV, show_col_types = FALSE) |>
  rename(function_id = 1)

func_mat <- pic |>
  column_to_rownames("function_id") |>
  as.matrix()

storage.mode(func_mat) <- "numeric"

md_pic <- read_tsv(META_TSV, show_col_types = FALSE)

md_pic <- md_pic |>
  mutate(
    PMDD = case_when(
      PMDD %in% c(0, "0", "C", "Control") ~ "Control",
      PMDD %in% c(2, "2", "P", "PMDD") ~ "PMDD",
      TRUE ~ NA_character_
    ),
    PMDD = factor(PMDD, levels = c("Control", "PMDD")),

    Phase = case_when(
      timepoint %in% c(1, "1", "F") ~ "F",
      timepoint %in% c(2, "2", "L") ~ "L",
      TRUE ~ as.character(timepoint)
    ),
    Phase = factor(Phase, levels = c("F", "L")),

    enrollment_id = factor(enrollment_id)
  )

# ------------------------------------------------------------
# 2. Construct functional phyloseq object
# ------------------------------------------------------------

build_func_phyloseq <- function(func_mat, md) {

  md$`sample-id` <- trimws(as.character(md$`sample-id`))
  colnames(func_mat) <- trimws(colnames(func_mat))

  samps <- intersect(
    colnames(func_mat),
    md$`sample-id`
  )

  cat("Overlapping samples:", length(samps), "\n")

  func2 <- func_mat[, samps, drop = FALSE]

  md2 <- md |>
    filter(`sample-id` %in% samps) |>
    arrange(match(`sample-id`, samps)) |>
    as.data.frame()

  rownames(md2) <- md2$`sample-id`

  tax <- data.frame(
    Function = rownames(func2)
  )
  rownames(tax) <- rownames(func2)

  phyloseq(
    otu_table(func2, taxa_are_rows = TRUE),
    sample_data(md2),
    tax_table(as.matrix(tax))
  )
}

ps_func <- build_func_phyloseq(
  func_mat,
  md_pic
)

ps_func

# ------------------------------------------------------------
# 3. Restrict to luteal phase and standardize cytokines
# ------------------------------------------------------------

cytokines <- c(
  "blood_il1",
  "blood_il6",
  "blood_tnfa",
  "blood_il8"
)

ps_func_L <- subset_samples(
  ps_func,
  Phase == "L"
)

sample_df <- as(
  sample_data(ps_func_L),
  "data.frame"
)

sample_df <- sample_df |>
  mutate(
    il1_z  = as.numeric(scale(blood_il1)),
    il6_z  = as.numeric(scale(blood_il6)),
    tnfa_z = as.numeric(scale(blood_tnfa)),
    il8_z  = as.numeric(scale(blood_il8))
  )

sample_data(ps_func_L) <- sample_data(sample_df)

cytokines_z <- c(
  "il1_z",
  "il6_z",
  "tnfa_z",
  "il8_z"
)

# ------------------------------------------------------------
# 4. Main cytokine models
# ------------------------------------------------------------

run_func_cyto_main <- function(ps, cytokine) {

  meta <- as.data.frame(
    sample_data(ps)
  )

  keep_ids <- rownames(meta)[
    !is.na(meta[[cytokine]]) &
    !is.na(meta$PMDD)
  ]

  ps_sub <- prune_samples(
    keep_ids,
    ps
  )

  fit <- ancombc2(
    data = ps_sub,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = cytokine,
    rand_formula = NULL,
    p_adj_method = "fdr",
    prv_cut = 0.10,
    lib_cut = 0,
    group = "PMDD",
    struc_zero = TRUE,
    neg_lb = TRUE,
    global = FALSE
  )

  list(
    fit = fit,
    n = nsamples(ps_sub)
  )
}

func_cyto_main <- lapply(
  cytokines_z,
  function(cyt) {
    run_func_cyto_main(
      ps_func_L,
      cyt
    )
  }
)

names(func_cyto_main) <- cytokines_z

# ------------------------------------------------------------
# 5. Cytokine x PMDD interaction models
# ------------------------------------------------------------

run_func_cyto_int <- function(ps, cytokine) {

  meta <- as.data.frame(
    sample_data(ps)
  )

  keep_ids <- rownames(meta)[
    !is.na(meta[[cytokine]]) &
    !is.na(meta$PMDD)
  ]

  ps_sub <- prune_samples(
    keep_ids,
    ps
  )

  fit <- ancombc2(
    data = ps_sub,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = paste0(
      cytokine,
      " * PMDD"
    ),
    rand_formula = NULL,
    p_adj_method = "fdr",
    prv_cut = 0.10,
    lib_cut = 0,
    group = "PMDD",
    struc_zero = TRUE,
    neg_lb = TRUE,
    global = FALSE
  )

  list(
    fit = fit,
    n = nsamples(ps_sub)
  )
}

func_cyto_int <- lapply(
  cytokines_z,
  function(cyt) {
    run_func_cyto_int(
      ps_func_L,
      cyt
    )
  }
)

names(func_cyto_int) <- cytokines_z

# ------------------------------------------------------------
# 6. Extract main and interaction results
# ------------------------------------------------------------

cyto_labels <- c(
  il1_z = "IL-1β",
  il6_z = "IL-6",
  il8_z = "IL-8",
  tnfa_z = "TNFα"
)

extract_func_main <- function(obj, cyt) {

  res <- obj$fit$res

  data.frame(
    cytokine = cyto_labels[[cyt]],
    n = obj$n,
    pathway = res$taxon,
    lfc = res[[paste0("lfc_", cyt)]],
    p_value = res[[paste0("p_", cyt)]],
    q_value = res[[paste0("q_", cyt)]],
    significant = res[[paste0("diff_", cyt)]],
    robust = res[[paste0("diff_robust_", cyt)]]
  )
}

extract_func_int <- function(obj, cyt) {

  res <- obj$fit$res

  term <- paste0(cyt, ":PMDDPMDD")

  data.frame(
    cytokine = cyto_labels[[cyt]],
    n = obj$n,
    pathway = res$taxon,
    lfc = res[[paste0("lfc_", term)]],
    p_value = res[[paste0("p_", term)]],
    q_value = res[[paste0("q_", term)]],
    significant = res[[paste0("diff_", term)]],
    robust = res[[paste0("diff_robust_", term)]]
  )
}

func_cyto_main_results <- bind_rows(
  lapply(
    cytokines_z,
    function(cyt) {
      extract_func_main(
        func_cyto_main[[cyt]],
        cyt
      )
    }
  )
)

func_cyto_int_results <- bind_rows(
  lapply(
    cytokines_z,
    function(cyt) {
      extract_func_int(
        func_cyto_int[[cyt]],
        cyt
      )
    }
  )
)

# ------------------------------------------------------------
# 7. Export results
# ------------------------------------------------------------

write.csv(
  func_cyto_main_results,
  "PICRUSt_cytokine_main_ALL.csv",
  row.names = FALSE
)

write.csv(
  func_cyto_main_results |>
    filter(significant == TRUE),
  "PICRUSt_cytokine_main_SIG.csv",
  row.names = FALSE
)

write.csv(
  func_cyto_int_results,
  "PICRUSt_cytokine_PMDD_interaction_ALL.csv",
  row.names = FALSE
)

write.csv(
  func_cyto_int_results |>
    filter(significant == TRUE),
  "PICRUSt_cytokine_PMDD_interaction_SIG.csv",
  row.names = FALSE
)
