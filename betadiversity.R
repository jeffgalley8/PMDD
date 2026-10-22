
# ============================================================
# PMDD MICROBIOME - BETA DIVERSITY (PMDD x PHASE)
# Bray-Curtis and unweighted UniFrac PERMANOVA
# ============================================================

library(vegan)
library(dplyr)

# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

# Distance matrices are exported from the original QIIME2

md <- read.delim(
  "metadata.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

read_distance <- function(path) {
  mat <- as.matrix(
    read.delim(
      path,
      row.names = 1,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )

  storage.mode(mat) <- "numeric"

  stopifnot(
    nrow(mat) == ncol(mat),
    identical(rownames(mat), colnames(mat))
  )

  as.dist(mat)
}

dist_bray <- read_distance("bray_curtis.tsv")
dist_uwuf <- read_distance("unweighted_unifrac.tsv")

OUTDIR <- "results/beta"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Metadata preparation
# ------------------------------------------------------------

md <- md |>
  mutate(
    PMDD = factor(
      PMDD,
      levels = c(0, 2),
      labels = c("Control", "PMDD")
    ),
    Phase = factor(Phase, levels = c("F", "L")),
    enrollment_id = factor(enrollment_id)
  )

# ------------------------------------------------------------
# 3. Align distance matrices with the analytic cohort
# ------------------------------------------------------------

study_ids <- md$`sample-id`

stopifnot(
  !anyDuplicated(study_ids),
  !anyNA(study_ids)
)

subset_distance <- function(d, ids) {
  stopifnot(all(ids %in% attr(d, "Labels")))

  result <- as.dist(
    as.matrix(d)[ids, ids, drop = FALSE]
  )

  stopifnot(
    identical(attr(result, "Labels"), ids)
  )

  result
}

dist_bray_72 <- subset_distance(dist_bray, study_ids)
dist_uwuf_72 <- subset_distance(dist_uwuf, study_ids)

# ------------------------------------------------------------
# 4. Define the original models
# ------------------------------------------------------------

model_vars <- list(
  M0 = c(
    "PMDD", "Phase", "enrollment_id"
  ),
  M1 = c(
    "PMDD", "Phase", "demos_age",
    "bmi", "enrollment_id"
  ),
  M2 = c(
    "PMDD", "Phase", "demos_age",
    "bmi", "Kcal", "Fibe", "omega",
    "enrollment_id"
  )
)

model_rhs <- list(
  M0 = "PMDD * Phase",
  M1 = "PMDD * Phase + demos_age + bmi",
  M2 = paste(
    "PMDD * Phase + demos_age + bmi +",
    "Kcal + Fibe + omega"
  )
)

# ------------------------------------------------------------
# 5. Run PERMANOVA
# ------------------------------------------------------------

run_beta <- function(distance, metadata, model_name) {

  keep <- complete.cases(
    metadata[, model_vars[[model_name]], drop = FALSE]
  )

  dat <- droplevels(
    metadata[keep, , drop = FALSE]
  )

  dsub <- subset_distance(
    distance,
    dat$`sample-id`
  )

  form <- as.formula(
    paste("dsub ~", model_rhs[[model_name]])
  )

  fit <- adonis2(
    form,
    data = dat,
    permutations = 9999,
    strata = dat$enrollment_id,
    by = "terms"
  )

  list(
    fit = fit,
    n_samples = nrow(dat),
    n_subjects = n_distinct(dat$enrollment_id),
    sample_ids = dat$`sample-id`
  )
}

# Bray-Curtis

beta_bray <- lapply(
  names(model_rhs),
  function(m) run_beta(dist_bray_72, md, m)
)

names(beta_bray) <- names(model_rhs)

# Unweighted UniFrac

beta_uwuf <- lapply(
  names(model_rhs),
  function(m) run_beta(dist_uwuf_72, md, m)
)

names(beta_uwuf) <- names(model_rhs)

beta_models <- list(
  Bray = beta_bray,
  Unweighted_UniFrac = beta_uwuf
)

# ------------------------------------------------------------
# 6. Combine and export results
# ------------------------------------------------------------

beta_results <- bind_rows(
  lapply(names(beta_models), function(metric) {

    bind_rows(
      lapply(names(beta_models[[metric]]), function(model_name) {

        obj <- beta_models[[metric]][[model_name]]

        tab <- as.data.frame(obj$fit)

        tab$term <- rownames(tab)
        tab$distance <- metric
        tab$model <- model_name
        tab$n_samples <- obj$n_samples
        tab$n_subjects <- obj$n_subjects

        rownames(tab) <- NULL

        tab
      })
    )
  })
) |>
  relocate(
    distance, model, n_samples,
    n_subjects, term
  )

write.csv(
  beta_results,
  file.path(OUTDIR, "beta_PERMANOVA_ALL_models.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Print all six models
# ------------------------------------------------------------

for (metric in names(beta_models)) {

  for (model_name in names(beta_models[[metric]])) {

    cat(
      "\n====",
      metric,
      model_name,
      "====\n"
    )

    print(
      beta_models[[metric]][[model_name]]$fit
    )
  }
}

# ------------------------------------------------------------
# 8. Save fitted models
# ------------------------------------------------------------

saveRDS(
  beta_models,
  file.path(OUTDIR, "beta_models.rds")
)
