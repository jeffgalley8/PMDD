# ============================================================
# CYTOKINE - ALPHA DIVERSITY
# ============================================================

library(dplyr)
library(lmerTest)

# 1. Import data

md <- read.delim(
  "metadata.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metrics <- c("shannon", "observed", "faith", "evenness")

cytokines <- c(
  "il1_z", "il6_z", "il8_z", "tnfa_z"
)

numeric_vars <- c(
  metrics,
  "blood_il1", "blood_il6",
  "blood_il8", "blood_tnfa"
)

md[numeric_vars] <- lapply(md[numeric_vars], as.numeric)

md <- md |>
  mutate(
    PMDD = factor(
      as.character(PMDD),
      levels = c("0", "2"),
      labels = c("Control", "PMDD")
    ),
    Phase = factor(Phase, levels = c("F", "L")),
    enrollment_id = factor(enrollment_id)
  )

# 2. Standardize cytokines
# Luteal models: within luteal cohort
# Longitudinal models: across both phases

md_L <- md |>
  filter(Phase == "L") |>
  mutate(
    il1_z = as.numeric(scale(blood_il1)),
    il6_z = as.numeric(scale(blood_il6)),
    il8_z = as.numeric(scale(blood_il8)),
    tnfa_z = as.numeric(scale(blood_tnfa))
  )

md_long <- md |>
  mutate(
    timepoint_f = Phase,
    il1_z = as.numeric(scale(blood_il1)),
    il6_z = as.numeric(scale(blood_il6)),
    il8_z = as.numeric(scale(blood_il8)),
    tnfa_z = as.numeric(scale(blood_tnfa))
  )

# 3. Luteal main effects, adjusted for PMDD
#    alpha ~ cytokine + PMDD (one luteal sample per participant)

run_alpha_main <- function(metric, cytokine) {

  dat <- md_L |>
    filter(
      !is.na(.data[[metric]]),
      !is.na(.data[[cytokine]]),
      !is.na(PMDD)
    ) |>
    droplevels()

  fit <- lm(
    reformulate(c(cytokine, "PMDD"), response = metric),
    data = dat
  )

  s <- summary(fit)$coefficients

  data.frame(
    metric = metric,
    cytokine = cytokine,
    n = nobs(fit),
    estimate = s[cytokine, "Estimate"],
    std_error = s[cytokine, "Std. Error"],
    t_value = s[cytokine, "t value"],
    p_value = s[cytokine, "Pr(>|t|)"]
  )
}

alpha_main <- bind_rows(
  lapply(metrics, function(m) {
    bind_rows(lapply(cytokines, function(cy) {
      run_alpha_main(m, cy)
    }))
  })
) |>
  group_by(metric) |>
  mutate(
    p_adj_fdr_metric = p.adjust(p_value, method = "fdr")
  ) |>
  ungroup() |>
  mutate(
    p_adj_fdr_global = p.adjust(p_value, method = "fdr")
  )

write.csv(
  alpha_main,
  "cytokine_alpha_luteal.csv",
  row.names = FALSE
)

# 4. Luteal cytokine x PMDD interactions

run_alpha_int <- function(metric, cytokine) {

  dat <- md_L |>
    filter(
      !is.na(.data[[metric]]),
      !is.na(.data[[cytokine]]),
      !is.na(PMDD)
    ) |>
    droplevels()

  fit <- lm(
    as.formula(paste(metric, "~", cytokine, "* PMDD")),
    data = dat
  )

  term <- paste0(cytokine, ":PMDDPMDD")
  s <- summary(fit)$coefficients

  data.frame(
    metric = metric,
    cytokine = cytokine,
    n = nobs(fit),
    estimate = s[term, "Estimate"],
    std_error = s[term, "Std. Error"],
    t_value = s[term, "t value"],
    p_value = s[term, "Pr(>|t|)"]
  )
}

alpha_int <- bind_rows(
  lapply(metrics, function(m) {
    bind_rows(lapply(cytokines, function(cy) {
      run_alpha_int(m, cy)
    }))
  })
) |>
  group_by(metric) |>
  mutate(
    p_adj_fdr_metric = p.adjust(p_value, method = "fdr")
  ) |>
  ungroup() |>
  mutate(
    p_adj_fdr_global = p.adjust(p_value, method = "fdr")
  )

write.csv(
  alpha_int,
  "cytokine_x_PMDD_alpha_luteal.csv",
  row.names = FALSE
)

# 5. Longitudinal cytokine x phase interactions

run_alpha_phase <- function(metric, cytokine) {

  dat <- md_long |>
    filter(
      !is.na(.data[[metric]]),
      !is.na(.data[[cytokine]]),
      !is.na(timepoint_f),
      !is.na(enrollment_id)
    )

  fit <- lmer(
    as.formula(paste0(
      metric, " ~ timepoint_f * ", cytokine,
      " + (1 | enrollment_id)"
    )),
    data = dat,
    REML = FALSE
  )

  term <- paste0("timepoint_fL:", cytokine)
  s <- summary(fit)$coefficients

  data.frame(
    metric = metric,
    cytokine = cytokine,
    n = nobs(fit),
    n_subjects = n_distinct(dat$enrollment_id),
    estimate = s[term, "Estimate"],
    std_error = s[term, "Std. Error"],
    df = s[term, "df"],
    t_value = s[term, "t value"],
    p_value = s[term, "Pr(>|t|)"]
  )
}

alpha_phase <- bind_rows(
  lapply(metrics, function(m) {
    bind_rows(lapply(cytokines, function(cy) {
      run_alpha_phase(m, cy)
    }))
  })
) |>
  group_by(metric) |>
  mutate(
    p_adj_fdr_metric = p.adjust(p_value, method = "fdr")
  ) |>
  ungroup() |>
  mutate(
    p_adj_fdr_global = p.adjust(p_value, method = "fdr")
  )

write.csv(
  alpha_phase,
  "cytokine_x_phase_alpha.csv",
  row.names = FALSE
)






# ============================================================
# CYTOKINE - BETA DIVERSITY
# Luteal phase | Bray-Curtis and Unweighted UniFrac
# ============================================================

library(dplyr)
library(vegan)

# 1. Import metadata

md <- read.delim(
  "metadata.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

numeric_vars <- c(
  "blood_il1", "blood_il6",
  "blood_il8", "blood_tnfa"
)

md[numeric_vars] <- lapply(md[numeric_vars], as.numeric)

md <- md |>
  mutate(
    PMDD = factor(
      as.character(PMDD),
      levels = c("0", "2"),
      labels = c("Control", "PMDD")
    ),
    Phase = factor(Phase, levels = c("F", "L")),
    enrollment_id = factor(enrollment_id)
  )

md_L <- md |>
  filter(Phase == "L") |>
  mutate(
    il1_z = as.numeric(scale(blood_il1)),
    il6_z = as.numeric(scale(blood_il6)),
    il8_z = as.numeric(scale(blood_il8)),
    tnfa_z = as.numeric(scale(blood_tnfa))
  )

cytokines <- c(
  "il1_z", "il6_z", "il8_z", "tnfa_z"
)

# 2. Import original QIIME2 distance matrices

read_distance <- function(path) {

  mat <- as.matrix(read.delim(
    path,
    row.names = 1,
    check.names = FALSE
  ))

  storage.mode(mat) <- "numeric"

  as.dist(mat)
}

dist_bray <- read_distance("bray_curtis.tsv")
dist_uwuf <- read_distance("unweighted_unifrac.tsv")

distances <- list(
  Bray = dist_bray,
  Unweighted_UniFrac = dist_uwuf
)

# 3. Luteal main cytokine effects

run_beta_main <- function(distance, cytokine) {

  dat <- md_L |>
    filter(!is.na(.data[[cytokine]]))

  ids <- dat$`sample-id`

  dsub <- as.dist(
    as.matrix(distance)[ids, ids, drop = FALSE]
  )

  set.seed(12345)

  fit <- adonis2(
    reformulate(cytokine, response = "dsub"),
    data = dat,
    permutations = 9999,
    by = "terms"
  )

  term <- as.data.frame(fit)[cytokine, , drop = FALSE]

  data.frame(
    cytokine = cytokine,
    n_samples = nrow(dat),
    n_subjects = n_distinct(dat$enrollment_id),
    R2 = term$R2,
    F_value = term$F,
    p_value = term[["Pr(>F)"]]
  )
}

beta_main <- bind_rows(
  lapply(names(distances), function(metric) {

    bind_rows(lapply(cytokines, function(cy) {
      run_beta_main(distances[[metric]], cy)
    })) |>
      mutate(distance = metric)

  })
) |>
  group_by(distance) |>
  mutate(
    p_adj_fdr = p.adjust(p_value, method = "fdr")
  ) |>
  ungroup() |>
  select(
    distance, cytokine, everything()
  )

write.csv(
  beta_main,
  "cytokine_beta_luteal.csv",
  row.names = FALSE
)

# 4. Cytokine x PMDD interactions

run_beta_int <- function(distance, cytokine) {

  dat <- md_L |>
    filter(
      !is.na(.data[[cytokine]]),
      !is.na(PMDD)
    ) |>
    droplevels()

  ids <- dat$`sample-id`

  dsub <- as.dist(
    as.matrix(distance)[ids, ids, drop = FALSE]
  )

  set.seed(12345)

  fit <- adonis2(
    as.formula(paste0(
      "dsub ~ ", cytokine, " * PMDD"
    )),
    data = dat,
    permutations = 9999,
    by = "terms"
  )

  term_name <- paste0(cytokine, ":PMDD")

  term <- as.data.frame(fit)[
    term_name, ,
    drop = FALSE
  ]

  data.frame(
    cytokine = cytokine,
    n_samples = nrow(dat),
    n_subjects = n_distinct(dat$enrollment_id),
    R2 = term$R2,
    F_value = term$F,
    p_value = term[["Pr(>F)"]]
  )
}

beta_int <- bind_rows(
  lapply(names(distances), function(metric) {

    bind_rows(lapply(cytokines, function(cy) {
      run_beta_int(distances[[metric]], cy)
    })) |>
      mutate(distance = metric)

  })
) |>
  group_by(distance) |>
  mutate(
    p_adj_fdr = p.adjust(p_value, method = "fdr")
  ) |>
  ungroup() |>
  select(
    distance, cytokine, everything()
  )

write.csv(
  beta_int,
  "cytokine_x_PMDD_beta_luteal.csv",
  row.names = FALSE
)

# 5. Export IL-8 and IL-8 x PMDD from the same Bray-Curtis model
#    Luteal phase only; one sample per participant.

md_il8 <- md_L |>
  filter(!is.na(il8_z), !is.na(PMDD)) |>
  droplevels()

il8_ids <- md_il8$`sample-id`
stopifnot(
  !anyDuplicated(il8_ids),
  all(il8_ids %in% attr(dist_bray, "Labels")),
  n_distinct(md_il8$enrollment_id) == nrow(md_il8)
)

d_il8 <- as.dist(
  as.matrix(dist_bray)[il8_ids, il8_ids, drop = FALSE]
)

# Same specification as run_beta_int(dist_bray, "il8_z").
fit_bray_il8 <- adonis2(
  d_il8 ~ il8_z * PMDD,
  data = md_il8,
  permutations = 9999,
  by = "terms"
)

il8_beta_table <- as.data.frame(fit_bray_il8) |>
  tibble::rownames_to_column("term") |>
  filter(term %in% c("il8_z", "il8_z:PMDD")) |>
  transmute(
    distance = "Bray-Curtis",
    term = recode(
      term,
      il8_z = "IL-8",
      `il8_z:PMDD` = "IL-8 x PMDD"
    ),
    n_samples = nrow(md_il8),
    n_subjects = n_distinct(md_il8$enrollment_id),
    Df,
    SumOfSqs,
    R2,
    F,
    p_value = `Pr(>F)`
  )

stopifnot(nrow(il8_beta_table) == 2)

write.csv(
  il8_beta_table,
  "Table_BrayCurtis_IL8_PMDD_interaction.csv",
  row.names = FALSE
)




# ============================================================
# CYTOKINE–TAXA ANALYSIS
# Luteal phase | ANCOM-BC2
# ============================================================

library(phyloseq)
library(ANCOMBC)
library(dplyr)

# ------------------------------------------------------------
# 1. Import genus counts and corrected metadata
# ------------------------------------------------------------

GENUS_TSV <- "genus_table_export/genus_table.tsv"
META_TSV  <- "metadata.tsv"

otu_genus <- read.delim(
  GENUS_TSV,
  skip = 1,
  row.names = 1,
  check.names = FALSE,
  comment.char = ""
)

if ("taxonomy" %in% colnames(otu_genus)) {
  otu_genus$taxonomy <- NULL
}

otu_genus <- as.matrix(otu_genus)
storage.mode(otu_genus) <- "numeric"

md <- read.delim(
  META_TSV,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

md <- md |>
  mutate(
    PMDD = factor(
      as.character(PMDD),
      levels = c("0", "2"),
      labels = c("Control", "PMDD")
    ),
    Phase = factor(
      as.character(Phase),
      levels = c("F", "L")
    ),
    enrollment_id = factor(enrollment_id)
  )

# Match count-table columns to metadata sample order
otu_genus <- otu_genus[
  ,
  md$`sample-id`,
  drop = FALSE
]

md_ps <- as.data.frame(md)
rownames(md_ps) <- md_ps$`sample-id`

tax_genus <- matrix(
  rownames(otu_genus),
  ncol = 1,
  dimnames = list(
    rownames(otu_genus),
    "Genus"
  )
)

ps_all <- phyloseq(
  otu_table(otu_genus, taxa_are_rows = TRUE),
  sample_data(md_ps),
  tax_table(tax_genus)
)

# ------------------------------------------------------------
# 2. Restrict to luteal phase and standardize cytokines
# ------------------------------------------------------------

ps_cyto_L <- subset_samples(ps_all, Phase == "L")

meta_L <- as.data.frame(sample_data(ps_cyto_L))

meta_L <- meta_L |>
  mutate(
    blood_il1 = as.numeric(blood_il1),
    blood_il6 = as.numeric(blood_il6),
    blood_il8 = as.numeric(blood_il8),
    blood_tnfa = as.numeric(blood_tnfa),

    il1_z = as.numeric(scale(blood_il1)),
    il6_z = as.numeric(scale(blood_il6)),
    il8_z = as.numeric(scale(blood_il8)),
    tnfa_z = as.numeric(scale(blood_tnfa))
  )

sample_data(ps_cyto_L) <- sample_data(meta_L)

cytokines_z <- c(
  "il1_z",
  "il6_z",
  "il8_z",
  "tnfa_z"
)

# ------------------------------------------------------------
# 3. Main cytokine associations
# ------------------------------------------------------------

run_cyto_taxa_main <- function(ps, cytokine) {

  meta <- as.data.frame(sample_data(ps))

  keep_ids <- rownames(meta)[
    !is.na(meta[[cytokine]]) &
    !is.na(meta$PMDD)
  ]

  ps_sub <- prune_samples(keep_ids, ps)

  fit <- ancombc2(
    data = ps_sub,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = cytokine,
    rand_formula = NULL,
    p_adj_method = "fdr",
    prv_cut = 0.10,
    lib_cut = 1000,
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

cyto_taxa_main <- lapply(
  cytokines_z,
  function(cyt) {
    run_cyto_taxa_main(ps_cyto_L, cyt)
  }
)

names(cyto_taxa_main) <- cytokines_z

# ------------------------------------------------------------
# 4. Extract main-effect results
# ------------------------------------------------------------

extract_cyto_main <- function(obj, cytokine) {

  res <- obj$fit$res

  data.frame(
    cytokine = cytokine,
    n = obj$n,
    taxon = res$taxon,
    lfc = res[[paste0("lfc_", cytokine)]],
    se = res[[paste0("se_", cytokine)]],
    W = res[[paste0("W_", cytokine)]],
    p_value = res[[paste0("p_", cytokine)]],
    q_value = res[[paste0("q_", cytokine)]],
    significant = res[[paste0("diff_", cytokine)]],
    robust = res[[paste0("diff_robust_", cytokine)]]
  )
}

cyto_taxa_main_results <- bind_rows(
  lapply(cytokines_z, function(cyt) {
    extract_cyto_main(cyto_taxa_main[[cyt]], cyt)
  })
) |>
  mutate(
    cytokine = recode(
      cytokine,
      il1_z = "IL-1β",
      il6_z = "IL-6",
      il8_z = "IL-8",
      tnfa_z = "TNFα"
    )
  )

cyto_taxa_main_sig <- cyto_taxa_main_results |>
  filter(significant == TRUE)

write.csv(
  cyto_taxa_main_results,
  "cytokine_taxa_main_ALL.csv",
  row.names = FALSE
)

write.csv(
  cyto_taxa_main_sig,
  "cytokine_taxa_main_SIG.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Cytokine × PMDD interaction models
# ------------------------------------------------------------

run_cyto_taxa_interaction <- function(ps, cytokine) {

  meta <- as.data.frame(sample_data(ps))

  keep_ids <- rownames(meta)[
    !is.na(meta[[cytokine]]) &
    !is.na(meta$PMDD)
  ]

  ps_sub <- prune_samples(keep_ids, ps)

  fit <- ancombc2(
    data = ps_sub,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = paste0(cytokine, " * PMDD"),
    rand_formula = NULL,
    p_adj_method = "fdr",
    prv_cut = 0.10,
    lib_cut = 1000,
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

cyto_taxa_int <- lapply(
  cytokines_z,
  function(cyt) {
    run_cyto_taxa_interaction(ps_cyto_L, cyt)
  }
)

names(cyto_taxa_int) <- cytokines_z

# ------------------------------------------------------------
# 6. Extract interaction results
# ------------------------------------------------------------

extract_cyto_int <- function(obj, cytokine) {

  res <- obj$fit$res

  int_term <- paste0(cytokine, ":PMDDPMDD")

  data.frame(
    cytokine = cytokine,
    n = obj$n,
    taxon = res$taxon,
    lfc = res[[paste0("lfc_", int_term)]],
    se = res[[paste0("se_", int_term)]],
    W = res[[paste0("W_", int_term)]],
    p_value = res[[paste0("p_", int_term)]],
    q_value = res[[paste0("q_", int_term)]],
    significant = res[[paste0("diff_", int_term)]],
    robust = res[[paste0("diff_robust_", int_term)]]
  )
}

cyto_taxa_int_results <- bind_rows(
  lapply(cytokines_z, function(cyt) {
    extract_cyto_int(cyto_taxa_int[[cyt]], cyt)
  })
) |>
  mutate(
    cytokine = recode(
      cytokine,
      il1_z = "IL-1β",
      il6_z = "IL-6",
      il8_z = "IL-8",
      tnfa_z = "TNFα"
    )
  )

cyto_taxa_int_sig <- cyto_taxa_int_results |>
  filter(significant == TRUE)

write.csv(
  cyto_taxa_int_results,
  "cytokine_x_PMDD_taxa_ALL.csv",
  row.names = FALSE
)

write.csv(
  cyto_taxa_int_sig,
  "cytokine_x_PMDD_taxa_SIG.csv",
  row.names = FALSE
)