
# ============================================================
# PMDD MICROBIOME - GENUS DIFFERENTIAL ABUNDANCE
# ANCOM-BC2: PMDD x Phase
# ============================================================

library(phyloseq)
library(ANCOMBC)
library(dplyr)

# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

GENUS_TSV <- "genus_table_export/genus_table.tsv"
META_TSV <- "metadata.tsv"

OUTDIR <- "results/taxa_pmdd_phase"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Import genus-level count table
# ------------------------------------------------------------

otu_genus <- read.delim(
  GENUS_TSV,
  skip = 1,
  row.names = 1,
  check.names = FALSE,
  comment.char = ""
)

# Remove taxonomy column if present
if ("taxonomy" %in% colnames(otu_genus)) {
  otu_genus$taxonomy <- NULL
}

otu_genus <- as.matrix(otu_genus)
storage.mode(otu_genus) <- "numeric"


# ------------------------------------------------------------
# 3. Import and prepare corrected metadata
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 4. Match samples and construct phyloseq
# ------------------------------------------------------------

study_ids <- md$`sample-id`

# Restrict and order genus counts to the study cohort
otu_genus <- otu_genus[
  ,
  study_ids,
  drop = FALSE
]

# Assign sample IDs as metadata row names
md_ps <- as.data.frame(md)
rownames(md_ps) <- md_ps$`sample-id`

# Construct genus-level taxonomy table
tax_genus <- matrix(
  rownames(otu_genus),
  ncol = 1,
  dimnames = list(
    rownames(otu_genus),
    "Genus"
  )
)

# Construct phyloseq object
ps_all <- phyloseq(
  otu_table(otu_genus, taxa_are_rows = TRUE),
  sample_data(md_ps),
  tax_table(tax_genus)
)


# ------------------------------------------------------------
# 5. Original ANCOM-BC2 model settings
# ------------------------------------------------------------

run_taxa <- function(ps, fixed_formula) {

  ancombc2(
    data = ps,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = fixed_formula,
    rand_formula = "(1 | enrollment_id)",
    p_adj_method = "fdr",
    prv_cut = 0.10,
    lib_cut = 1000,
    group = "PMDD",
    struc_zero = TRUE,
    neg_lb = TRUE,
    global = FALSE
  )
}

# ------------------------------------------------------------
# 6. M0: PMDD x Phase
# ------------------------------------------------------------

fit_taxa_m0 <- run_taxa(
  ps_all,
  "PMDD * Phase"
)

saveRDS(
  fit_taxa_m0,
  file.path(OUTDIR, "taxa_M0.rds")
)

# ------------------------------------------------------------
# 7. M1: PMDD x Phase + age + BMI
# ------------------------------------------------------------

fit_taxa_m1 <- run_taxa(
  ps_all,
  "PMDD * Phase + demos_age + bmi"
)

saveRDS(
  fit_taxa_m1,
  file.path(OUTDIR, "taxa_M1.rds")
)

# ------------------------------------------------------------
# 8. M2: Additional dietary covariates
# ------------------------------------------------------------

m2_vars <- c(
  "PMDD",
  "Phase",
  "demos_age",
  "bmi",
  "Kcal",
  "Fibe",
  "omega",
  "enrollment_id"
)

m2_ids <- md$`sample-id`[
  complete.cases(
    md[, m2_vars, drop = FALSE]
  )
]

ps_m2 <- prune_samples(
  sample_names(ps_all) %in% m2_ids,
  ps_all
)

fit_taxa_m2 <- run_taxa(
  ps_m2,
  paste(
    "PMDD * Phase + demos_age + bmi +",
    "Kcal + Fibe + omega"
  )
)

saveRDS(
  fit_taxa_m2,
  file.path(OUTDIR, "taxa_M2.rds")
)

# ------------------------------------------------------------
# 9. Export full ANCOM-BC2 results
# ------------------------------------------------------------

taxa_models <- list(
  M0 = fit_taxa_m0,
  M1 = fit_taxa_m1,
  M2 = fit_taxa_m2
)

for (model_name in names(taxa_models)) {

  write.csv(
    taxa_models[[model_name]]$res,
    file.path(
      OUTDIR,
      paste0("taxa_", model_name, "_ALL.csv")
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 10. Summarize FDR-significant taxa
# ------------------------------------------------------------

significance_summary <- bind_rows(
  lapply(names(taxa_models), function(model_name) {

    res <- taxa_models[[model_name]]$res

    diff_cols <- grep(
      "^diff_",
      names(res),
      value = TRUE
    )

    data.frame(
      model = model_name,
      term = diff_cols,
      n_significant = vapply(
        diff_cols,
        function(x) {
          sum(res[[x]], na.rm = TRUE)
        },
        numeric(1)
      )
    )
  })
)

write.csv(
  significance_summary,
  file.path(OUTDIR, "taxa_significance_counts.csv"),
  row.names = FALSE
)

print(significance_summary)
