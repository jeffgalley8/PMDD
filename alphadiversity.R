# ============================================================
# PMDD MICROBIOME - ALPHA DIVERSITY ANALYSIS
# Longitudinal mixed-effects models: all metrics and models
# Input: metadata.tsv
# ============================================================

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

# 1. Read metadata (alpha metrics already included)
md <- read.delim(INPUT, header = TRUE, sep = "\t",
                 check.names = FALSE, stringsAsFactors = FALSE)

metrics <- c("shannon", "observed", "faith", "evenness")

# The supplied metadata contains numeric alpha metrics and covariates.
required_vars <- c("enrollment_id", "PMDD", "Phase", metrics,
                   "demos_age", "bmi", "Kcal", "Fibe", "omega")

md <- md |>
  mutate(
    PMDD = factor(as.character(PMDD), levels = c("0", "2"),
                  labels = c("Control", "PMDD")),
    Phase = factor(Phase, levels = c("F", "L")),
    enrollment_id = factor(enrollment_id)
  )

# 2. Define the original M0/M1/M2 models
make_formulas <- function(outcome) {
  list(
    M0 = as.formula(paste0(outcome, " ~ PMDD * Phase + (1 | enrollment_id)")),
    M1 = as.formula(paste0(outcome, " ~ PMDD * Phase + demos_age + bmi + (1 | enrollment_id)")),
    M2 = as.formula(paste0(outcome, " ~ PMDD * Phase + demos_age + bmi + Kcal + Fibe + omega + (1 | enrollment_id)"))
  )
}

# 3. Fit all 12 models
alpha_models <- setNames(lapply(metrics, function(metric) {
  lapply(make_formulas(metric), function(f) lmer(f, data = md, REML = FALSE))
}), metrics)

# 4. Extract every fixed-effect coefficient from every model
extract_coefficients <- function(fit, metric, model_name) {
  s <- summary(fit)$coefficients
  data.frame(
    metric = metric,
    model = model_name,
    n_obs = nobs(fit),
    n_subjects = nlevels(getME(fit, "flist")$enrollment_id),
    term = rownames(s),
    estimate = s[, "Estimate"],
    std_error = s[, "Std. Error"],
    df = s[, "df"],
    t_value = s[, "t value"],
    p_value = s[, "Pr(>|t|)"],
    row.names = NULL
  )
}
alpha_results <- bind_rows(lapply(names(alpha_models), function(metric) {
  bind_rows(lapply(names(alpha_models[[metric]]), function(model_name) {
    extract_coefficients(alpha_models[[metric]][[model_name]], metric, model_name)
  }))
}))
write.csv(alpha_results, file.path(OUTDIR, "alpha_ALL_coefficients.csv"),
          row.names = FALSE)

# 5. Focused PMDD, Phase, and PMDD-by-Phase coefficient table
key_terms <- c("PMDDPMDD", "PhaseL", "PMDDPMDD:PhaseL")
alpha_key_results <- alpha_results |>
  filter(term %in% key_terms) |>
  arrange(metric, model, term)
write.csv(alpha_key_results,
          file.path(OUTDIR, "alpha_PMDD_Phase_key_effects.csv"),
          row.names = FALSE)

# 6. Within-group F vs L contrasts for all 12 models
extract_phase_contrasts <- function(fit, metric, model_name) {
  em <- emmeans(fit, ~ Phase | PMDD, lmer.df = "satterthwaite")
  out <- as.data.frame(pairs(em, adjust = "none"))
  out$metric <- metric
  out$model <- model_name
  out
}
alpha_phase_contrasts <- bind_rows(lapply(names(alpha_models), function(metric) {
  bind_rows(lapply(names(alpha_models[[metric]]), function(model_name) {
    extract_phase_contrasts(alpha_models[[metric]][[model_name]], metric, model_name)
  }))
})) |>
  select(metric, model, PMDD, contrast, estimate, SE, df, t.ratio, p.value) |>
  arrange(metric, model, PMDD)
write.csv(alpha_phase_contrasts,
          file.path(OUTDIR, "alpha_ALL_phase_contrasts_unadjusted.csv"),
          row.names = FALSE)

# 7. Save fitted models and session information
saveRDS(alpha_models, file.path(OUTDIR, "alpha_models.rds"))
writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, "session_info.txt"))