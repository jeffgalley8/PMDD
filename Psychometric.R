
# ================================================================
# PSYCHOMETRIC ASSOCIATIONS WITH ALPHA, BETA, AND TAXA
# ================================================================
#
# ================================================================

library(readr)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(vegan)
library(phyloseq)
library(ANCOMBC)

OUTDIR <- "psychometric_analysis"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

psych_vars <- c(
  "drsp_tot", "drsp_mood", "ace_tot", "tleq_tot",
  "qids_tot", "stais_tot", "pss_tot", "pirs_tot",
  "fcqs_tot", "phq9_tot", "gad7_tot"
)
alpha_vars <- c("shannon", "observed", "faith", "evenness")

# ------------------------------------------------
# 1. Metadata and psychometric standardization
# ------------------------------------------------
md <- read_tsv("metadata.tsv", show_col_types = FALSE)

required <- c("sample-id", "enrollment_id", "PMDD", psych_vars, alpha_vars)
missing_cols <- setdiff(required, names(md))
if (length(missing_cols)) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}
if (anyDuplicated(md$`sample-id`)) stop("Duplicate sample IDs in metadata.tsv")

md <- md |>
  mutate(enrollment_id = factor(enrollment_id))

# The PMDD reference category for interaction models is Control.
pmdd_raw <- as.character(md$PMDD)
pmdd_raw[pmdd_raw %in% c("C")] <- "Control"
pmdd_raw[pmdd_raw %in% c("P")] <- "PMDD"
if (any(!is.na(pmdd_raw) & !pmdd_raw %in% c("Control", "PMDD"))) {
  stop("Unexpected PMDD values; check FinalMap.tsv")
}
md$PMDD <- factor(pmdd_raw, levels = c("Control", "PMDD"))

md_psych <- md |>
  mutate(across(all_of(psych_vars), ~ as.numeric(scale(as.numeric(.x)))))

# ------------------------------------------------
# 2. Alpha diversity: one psychometric per LME
#    alpha ~ psychometric + (1 | enrollment_id)
# ------------------------------------------------
run_psych_alpha <- function(data, alpha_col, psych) {
  form <- as.formula(paste0(alpha_col, " ~ ", psych,
                            " + (1 | enrollment_id)"))
  mod <- lmer(form, data = data, REML = FALSE)
  s <- summary(mod)$coefficients
  mf <- model.frame(mod)
  data.frame(
    metric = alpha_col,
    psychometric = psych,
    estimate = s[psych, "Estimate"],
    std_error = s[psych, "Std. Error"],
    df = s[psych, "df"],
    t_value = s[psych, "t value"],
    p_value = s[psych, "Pr(>|t|)"],
    n_obs = nobs(mod),
    n_subjects = n_distinct(mf$enrollment_id)
  )
}

psych_alpha_results <- bind_rows(lapply(alpha_vars, function(a) {
  bind_rows(lapply(psych_vars, function(p) run_psych_alpha(md_psych, a, p)))
})) |>
  group_by(metric) |>
  mutate(p_adj_fdr = p.adjust(p_value, method = "fdr")) |>
  ungroup() |>
  arrange(metric, p_adj_fdr)

write.csv(psych_alpha_results,
          file.path(OUTDIR, "psychometric_vs_alpha_timevarying_LME.csv"),
          row.names = FALSE)

# ------------------------------------------------
# 3. Beta diversity: one psychometric per PERMANOVA
#    distance ~ psychometric
#    9,999 permutations, restricted within subjects
# ------------------------------------------------
read_qiime_dist <- function(path) {
  m <- read.delim(path, check.names = FALSE, row.names = 1)
  m <- as.matrix(m)
  storage.mode(m) <- "numeric"
  as.dist(m)
}

dist_bray <- read_qiime_dist("braycurtis/distance-matrix.tsv")
dist_uwuf <- read_qiime_dist("unweightunifrac/distance-matrix.tsv")

run_beta_psych <- function(dist_obj, metadata, psych, nperm = 9999) {
  md_filt <- metadata |>
    filter(!is.na(.data[[psych]]), !is.na(enrollment_id))
  samps <- intersect(labels(dist_obj), md_filt$`sample-id`)
  if (length(samps) < 3) stop("Too few overlapping samples for ", psych)
  dist_sub <- as.dist(as.matrix(dist_obj)[samps, samps, drop = FALSE])
  md_sub <- md_filt |>
    filter(`sample-id` %in% samps) |>
    arrange(match(`sample-id`, samps))
  stopifnot(identical(as.character(md_sub$`sample-id`), samps))

  form <- as.formula(paste0("dist_sub ~ ", psych))
  ad <- adonis2(form, data = md_sub, permutations = nperm,
               strata = md_sub$enrollment_id)
  out <- as.data.frame(ad)[1, , drop = FALSE]
  out$psychometric <- psych
  out$n_obs <- length(samps)
  out$n_subjects <- n_distinct(md_sub$enrollment_id)
  out
}

results_bray_psych <- bind_rows(lapply(psych_vars, function(p) {
  run_beta_psych(dist_bray, md_psych, p)
})) |>
  mutate(p_adj_fdr = p.adjust(`Pr(>F)`, method = "fdr")) |>
  arrange(p_adj_fdr)

results_uwuf_psych <- bind_rows(lapply(psych_vars, function(p) {
  run_beta_psych(dist_uwuf, md_psych, p)
})) |>
  mutate(p_adj_fdr = p.adjust(`Pr(>F)`, method = "fdr")) |>
  arrange(p_adj_fdr)

write.csv(results_bray_psych,
          file.path(OUTDIR, "beta_psychometric_bray_timevarying_PERMANOVA.csv"),
          row.names = FALSE)
write.csv(results_uwuf_psych,
          file.path(OUTDIR, "beta_psychometric_unweighted_unifrac_timevarying_PERMANOVA.csv"),
          row.names = FALSE)

# ------------------------------------------------
# 4. Taxa: align phyloseq with the same metadata
# ------------------------------------------------
ps_psych <- prune_samples(sample_names(ps_all) %in% md_psych$`sample-id`,
                          ps_all)
md_taxa <- md_psych |>
  filter(`sample-id` %in% sample_names(ps_psych)) |>
  arrange(match(`sample-id`, sample_names(ps_psych))) |>
  as.data.frame()
rownames(md_taxa) <- md_taxa$`sample-id`
stopifnot(identical(rownames(md_taxa), sample_names(ps_psych)))
sample_data(ps_psych) <- sample_data(md_taxa)

TAXADIR <- file.path(OUTDIR, "psychometric_taxa_ancombc")
dir.create(TAXADIR, showWarnings = FALSE, recursive = TRUE)

# Keep all ANCOM-BC2 settings the same for both model sets.
run_psych_taxa <- function(ps, var, interaction = FALSE) {
  formula_text <- if (interaction) paste0(var, " * PMDD") else var
  ancombc2(
    data = ps,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = formula_text,
    rand_formula = "(1 | enrollment_id)",
    p_adj_method = "fdr",
    prv_cut = 0.10,
    lib_cut = 0,
    group = "PMDD",
    struc_zero = TRUE,
    neg_lb = TRUE,
    global = FALSE
  )
}

# Extract one coefficient's full, nominal, and FDR-significant taxa.
# In the interaction models, extract the interaction coefficient only.
extract_psych_taxa <- function(fit, psych, interaction = FALSE) {
  res <- as.data.frame(fit$res)
  if (interaction) {
    p_candidates <- grep("^p_", names(res), value = TRUE)
    p_candidates <- p_candidates[
      grepl(psych, p_candidates, fixed = TRUE) &
      grepl("PMDD", p_candidates, fixed = TRUE) &
      grepl(":", p_candidates, fixed = TRUE)
    ]
    if (length(p_candidates) != 1) {
      stop("Cannot identify unique interaction term for ", psych)
    }
    term <- sub("^p_", "", p_candidates)
  } else {
    term <- psych
  }

  p_col <- paste0("p_", term)
  q_col <- paste0("q_", term)
  lfc_col <- paste0("lfc_", term)
  diff_col <- paste0("diff_", term)
  required_cols <- c(p_col, q_col, lfc_col)
  if (!all(required_cols %in% names(res))) {
    stop("Missing ANCOM-BC2 result columns for ", psych, ": ", term)
  }

  # Retain the full original result columns and add readable identifiers.
  res$psychometric <- psych
  res$tested_term <- term
  res$tested_lfc <- res[[lfc_col]]
  res$tested_p <- res[[p_col]]
  res$tested_q <- res[[q_col]]
  if (diff_col %in% names(res)) {
    res$tested_diff <- res[[diff_col]]
  } else {
    res$tested_diff <- !is.na(res$tested_q) & res$tested_q < 0.05
  }
  res <- res |> arrange(tested_p)
  list(
    all = res,
    nominal = res |> filter(!is.na(tested_p), tested_p < 0.05),
    fdr = res |> filter(!is.na(tested_diff), tested_diff)
  )
}

summary_rows <- list()
for (pv in psych_vars) {
  fit_main <- run_psych_taxa(ps_psych, pv, interaction = FALSE)
  saveRDS(fit_main, file.path(TAXADIR, paste0("ancombc2_psych_", pv, ".rds")))
  main <- extract_psych_taxa(fit_main, pv, interaction = FALSE)
  write.csv(main$all, file.path(TAXADIR, paste0("all_taxa_psych_", pv, ".csv")), row.names = FALSE)
  write.csv(main$nominal, file.path(TAXADIR, paste0("p05_taxa_psych_", pv, ".csv")), row.names = FALSE)
  write.csv(main$fdr, file.path(TAXADIR, paste0("FDR_taxa_psych_", pv, ".csv")), row.names = FALSE)
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    model = "main", psychometric = pv,
    n_p05 = nrow(main$nominal), n_FDR = nrow(main$fdr)
  )

  fit_int <- run_psych_taxa(ps_psych, pv, interaction = TRUE)
  saveRDS(fit_int, file.path(TAXADIR, paste0("ancombc2_psych_PMDD_", pv, ".rds")))
  int <- extract_psych_taxa(fit_int, pv, interaction = TRUE)
  write.csv(int$all, file.path(TAXADIR, paste0("all_taxa_psych_PMDD_", pv, ".csv")), row.names = FALSE)
  write.csv(int$nominal, file.path(TAXADIR, paste0("p05_taxa_psych_PMDD_", pv, ".csv")), row.names = FALSE)
  write.csv(int$fdr, file.path(TAXADIR, paste0("FDR_taxa_psych_PMDD_", pv, ".csv")), row.names = FALSE)
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    model = "interaction", psychometric = pv,
    n_p05 = nrow(int$nominal), n_FDR = nrow(int$fdr)
  )
}
