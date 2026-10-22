# Steroid-microbiome analysis
# Inputs: metadata.tsv, QIIME 2 distance matrices, and ps_all.

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(lme4)
library(lmerTest)
library(vegan)
library(ANCOMBC)

# ============================================================
# 1. METADATA
# ============================================================
md <- read_tsv('metadata.tsv', show_col_types = FALSE)

# Set factor levels
md <- md |>
  mutate(
    PMDD = case_when(
      as.character(PMDD) %in% c('C') ~ 'Control',
      as.character(PMDD) %in% c('P') ~ 'PMDD',
      TRUE ~ NA_character_
    ),
    PMDD = factor(PMDD, levels = c('Control', 'PMDD')),
    Phase = factor(as.character(Phase), levels = c('F', 'L')),
    enrollment_id = factor(enrollment_id)
  )
stopifnot(!anyNA(md$PMDD), !anyNA(md$Phase))

raw_steroids <- c('Allo', 'PA', 'Epi', 'IsoAllo')

stopifnot(all(vapply(md[raw_steroids], function(x)
  all(is.na(x) | (is.finite(x) & x > 0)), logical(1))))

# Log-transform and standardize steroid concentrations
md <- md |>
  mutate(
    Allo_z = as.numeric(scale(log(Allo))),
    PA_z = as.numeric(scale(log(PA))),
    Epi_z = as.numeric(scale(log(Epi))),
    IsoAllo_z = as.numeric(scale(log(IsoAllo)))
  )
steroids <- c('Allo_z', 'PA_z', 'Epi_z', 'IsoAllo_z')
metrics <- c('shannon', 'observed', 'faith', 'evenness')
print(sapply(md[steroids], function(x) c(
  n = sum(!is.na(x)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))))

# ============================================================
# 2. ALPHA DIVERSITY
# ============================================================
run_steroid_alpha_models <- function(data, outcome, steroid) {
  f0 <- as.formula(paste0(outcome, ' ~ ', steroid,
                          ' + PMDD + Phase + (1 | enrollment_id)'))
  f1 <- as.formula(paste0(outcome, ' ~ ', steroid,
                          ' + PMDD + Phase + demos_age + bmi + (1 | enrollment_id)'))
  f2 <- as.formula(paste0(outcome, ' ~ ', steroid,
                          ' + PMDD + Phase + demos_age + bmi + Kcal + Fibe + omega',
                          ' + (1 | enrollment_id)'))
  list(M0 = lmer(f0, data = data, REML = FALSE),
       M1 = lmer(f1, data = data, REML = FALSE),
       M2 = lmer(f2, data = data, REML = FALSE))
}

extract_steroid_result <- function(model, outcome, steroid, model_name) {
  cf <- summary(model)$coefficients[steroid, ]
  data.frame(metric = outcome, steroid = steroid, model = model_name,
             estimate = unname(cf['Estimate']),
             SE = unname(cf['Std. Error']),
             df = unname(cf['df']),
             t_value = unname(cf['t value']),
             p_value = unname(cf['Pr(>|t|)']),
             n = nobs(model),
             n_subjects = length(unique(model.frame(model)$enrollment_id)))
}

results_steroid_models <- bind_rows(lapply(metrics, function(met) {
  bind_rows(lapply(steroids, function(st) {
    mods <- run_steroid_alpha_models(md, met, st)
    bind_rows(extract_steroid_result(mods$M0, met, st, 'M0'),
              extract_steroid_result(mods$M1, met, st, 'M1'),
              extract_steroid_result(mods$M2, met, st, 'M2'))
  }))
})) |>
  mutate(steroid = recode(steroid, Allo_z = 'Allo', PA_z = 'PA',
                          Epi_z = 'Epi', IsoAllo_z = 'IsoAllo')) |>
  group_by(model, metric) |>
  mutate(p_adj_fdr_metric = p.adjust(p_value, method = 'fdr')) |>
  ungroup() |>
  group_by(model) |>
  mutate(p_adj_fdr_global = p.adjust(p_value, method = 'fdr')) |>
  ungroup() |>
  arrange(model, metric, p_value)

print(as.data.frame(results_steroid_models))
print(as.data.frame(results_steroid_models |>
  filter(p_value < 0.05 | p_adj_fdr_metric < 0.05 | p_adj_fdr_global < 0.05)))
write.csv(results_steroid_models,
          'alpha_steroid_M0_M1_M2.csv', row.names = FALSE)

# ============================================================
# 3. BETA DIVERSITY
# ============================================================
read_qiime_dist <- function(path) {
  m <- read.delim(path, check.names = FALSE, row.names = 1)
  m <- as.matrix(m)
  storage.mode(m) <- 'numeric'
  as.dist(m)
}
dist_bray <- read_qiime_dist('braycurtis/distance-matrix.tsv')
dist_uwuf <- read_qiime_dist('unweightunifrac/distance-matrix.tsv')

align_dist_and_md <- function(dist_obj, metadata) {
  samps <- intersect(attr(dist_obj, 'Labels'), metadata$`sample-id`)
  if (length(samps) < 3) stop('Too few overlapping samples')
  dist_sub <- as.dist(as.matrix(dist_obj)[samps, samps])
  md_sub <- metadata |>
    filter(`sample-id` %in% samps) |>
    arrange(match(`sample-id`, samps))
  stopifnot(identical(md_sub$`sample-id`, attr(dist_sub, 'Labels')))
  list(dist = dist_sub, md = md_sub)
}

run_steroid_beta <- function(dist_obj, metadata, steroid,
                                      distance_name, interaction = FALSE,
                                      nperm = 9999) {
  md_filt <- metadata |>
    filter(!is.na(.data[[steroid]]), !is.na(PMDD),
           !is.na(Phase), !is.na(enrollment_id))
  al <- align_dist_and_md(dist_obj, md_filt)
  d <- al$dist
  md_sub <- al$md
  rhs <- if (interaction) paste0(steroid, ' * PMDD + Phase') else
    paste0(steroid, ' + PMDD + Phase')
  form <- as.formula(paste('d ~', rhs))
  fit <- adonis2(form, data = md_sub, permutations = nperm,
                strata = md_sub$enrollment_id, by = 'terms')
  fit_df <- as.data.frame(fit) |> rownames_to_column('term')
  if (interaction) {
    fit_df <- fit_df |>
      filter(grepl(steroid, term, fixed = TRUE),
             grepl('PMDD', term, fixed = TRUE), grepl(':', term, fixed = TRUE))
  } else {
    fit_df <- fit_df |> filter(term == steroid)
  }
  if (nrow(fit_df) != 1) stop('Cannot identify beta term: ', steroid)
  fit_df |>
    transmute(distance = distance_name, steroid = steroid,
              term = term, n_obs = nrow(md_sub),
              n_subjects = n_distinct(md_sub$enrollment_id),
              Df = Df, SumOfSqs = SumOfSqs, R2 = R2,
              F = F, p_value = `Pr(>F)`)
}

run_beta_set <- function(interaction) {
  bind_rows(
    lapply(steroids, function(st) run_steroid_beta(
      dist_bray, md, st, 'Bray-Curtis', interaction = interaction)),
    lapply(steroids, function(st) run_steroid_beta(
      dist_uwuf, md, st, 'Unweighted UniFrac', interaction = interaction))
  ) |>
    group_by(distance) |>
    mutate(p_adj_fdr = p.adjust(p_value, method = 'fdr')) |>
    ungroup() |>
    mutate(steroid = recode(steroid, Allo_z = 'Allo', PA_z = 'PA',
                            Epi_z = 'Epi', IsoAllo_z = 'IsoAllo')) |>
    arrange(distance, p_value)
}

beta_main <- run_beta_set(FALSE)
beta_interaction <- run_beta_set(TRUE)

write.csv(beta_main, 'beta_steroid_main_effects.csv', row.names = FALSE)
write.csv(beta_interaction, 'beta_steroid_PMDD_interactions.csv', row.names = FALSE)

# ============================================================
# 4. TAXA
# ============================================================
# ps_all is constructed in the separate microbiome analysis script.
ps_taxa <- prune_samples(
  sample_names(ps_all) %in% md$`sample-id`,
  ps_all
)
md_taxa <- md |>
  filter(`sample-id` %in% sample_names(ps_taxa)) |>
  arrange(match(`sample-id`, sample_names(ps_taxa))) |>
  as.data.frame()
rownames(md_taxa) <- md_taxa$`sample-id`
stopifnot(identical(rownames(md_taxa), sample_names(ps_taxa)))
sample_data(ps_taxa) <- sample_data(md_taxa)

OUTDIR <- 'steroid_taxa_ancombc'
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)


run_steroid_taxa <- function(steroid, label) {
  fit <- ancombc2(
    data = ps_taxa,
    assay_name = 'counts',
    tax_level = NULL,
    fix_formula = paste0(steroid, ' * PMDD + Phase'),
    rand_formula = '(1 | enrollment_id)',
    p_adj_method = 'fdr',
    prv_cut = 0.10,
    lib_cut = 0,
    group = 'PMDD',
    struc_zero = TRUE,
    neg_lb = TRUE,
    global = FALSE
  )
  saveRDS(fit, file.path(OUTDIR, paste0('taxa_', label, '_PMDD.rds')))
  res <- fit$res
  main_suffix <- steroid
  interaction_suffix <- paste0(steroid, ':PMDDPMDD')
  if (!paste0('p_', interaction_suffix) %in% names(res)) {
    possible <- grep(paste0('^p_', steroid, ':PMDD'), names(res), value = TRUE)
    if (length(possible) != 1) stop('Unexpected interaction columns for ', label,
                                   ': ', paste(possible, collapse = ', '))
    interaction_suffix <- sub('^p_', '', possible)
  }
  for (effect in c('main', 'interaction')) {
    suffix <- if (effect == 'main') main_suffix else interaction_suffix
    p_col <- paste0('p_', suffix)
    diff_col <- paste0('diff_', suffix)
    stopifnot(p_col %in% names(res), diff_col %in% names(res))
    # FDR-significant taxa
    sig_fdr <- res |> filter(.data[[diff_col]] %in% TRUE)
    # Nominally significant taxa
    sig_p <- res |> filter(!is.na(.data[[p_col]]) & .data[[p_col]] < 0.05)
    tag <- if (effect == 'main') label else paste0(label, 'xPMDD')
    write.csv(sig_fdr, file.path(OUTDIR, paste0('sig_taxa_', tag, '.csv')),
              row.names = FALSE)
    write.csv(sig_p, file.path(OUTDIR, paste0('p05_taxa_', tag, '.csv')),
              row.names = FALSE)
    cat(label, effect, ': FDR significant =', nrow(sig_fdr),
        '; nominal p < 0.05 =', nrow(sig_p), '\n')
  }
  invisible(fit)
}

fit_allo <- run_steroid_taxa('Allo_z', 'Allo')
fit_PA <- run_steroid_taxa('PA_z', 'PA')
fit_epi <- run_steroid_taxa('Epi_z', 'Epi')
fit_isoallo <- run_steroid_taxa('IsoAllo_z', 'IsoAllo')