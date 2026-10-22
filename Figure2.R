# ============================================================
# PMDD MICROBIOME - BETA DIVERSITY PCoA FIGURE
# ============================================================

library(dplyr)
library(ggplot2)
library(ape)
library(patchwork)

# ------------------------------------------------------------
# 1. Read the metadata and distance matrices
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 2. Metadata coding: identical to finalized beta analysis
# ------------------------------------------------------------
md_plot <- md |>
  mutate(
    PMDD = factor(PMDD, levels = c(0, 2), labels = c("Control", "PMDD")),
    Phase = factor(Phase, levels = c("F", "L")),
    enrollment_id = factor(enrollment_id)
  )

study_ids <- md_plot$`sample-id`
stopifnot(
  !anyDuplicated(study_ids),
  !anyNA(study_ids),
  !anyNA(md_plot$PMDD),
  !anyNA(md_plot$Phase)
)

# ------------------------------------------------------------
# 3. Align both distances to exact metadata sample order
# ------------------------------------------------------------
subset_distance <- function(d, ids) {
  stopifnot(all(ids %in% attr(d, "Labels")))
  result <- as.dist(as.matrix(d)[ids, ids, drop = FALSE])
  stopifnot(identical(attr(result, "Labels"), ids))
  result
}

dist_bray_72 <- subset_distance(dist_bray, study_ids)
dist_uwuf_72 <- subset_distance(dist_uwuf, study_ids)

# ------------------------------------------------------------
# 4. Compute PCoA scores and percent variance for plotting
# ------------------------------------------------------------
make_pcoa_df <- function(dist_obj, md, label) {
  samps <- attr(dist_obj, "Labels")
  stopifnot(identical(samps, md$`sample-id`))

  pcoa <- ape::pcoa(dist_obj)
  scores <- as.data.frame(pcoa$vectors[, 1:2, drop = FALSE])
  scores$`sample-id` <- rownames(scores)

  var1 <- round(pcoa$values$Relative_eig[1] * 100, 1)
  var2 <- round(pcoa$values$Relative_eig[2] * 100, 1)

  scores |>
    left_join(md, by = "sample-id") |>
    mutate(
      distance = label,
      axis1_lab = paste0("PCoA1 (", var1, "%)"),
      axis2_lab = paste0("PCoA2 (", var2, "%)")
    )
}

# ------------------------------------------------------------
# 5. Original figure styling, unchanged
# ------------------------------------------------------------
plot_pcoa <- function(pcoa_df, title) {
  xlab <- unique(pcoa_df$axis1_lab)
  ylab <- unique(pcoa_df$axis2_lab)

  ggplot(
    pcoa_df,
    aes(x = Axis.1, y = Axis.2, color = PMDD, shape = Phase)
  ) +
    geom_point(size = 3, alpha = 0.9) +
    theme_classic(base_size = 14) +
    labs(
      title = title,
      x = xlab,
      y = ylab,
      color = "Group",
      shape = "Phase"
    ) +
    scale_color_manual(
      values = c("Control" = "#999999", "PMDD" = "#D55E00"),
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c("F" = 16, "L" = 17),
      labels = c("F" = "Follicular", "L" = "Luteal"),
      drop = FALSE
    )
}

bray_df <- make_pcoa_df(dist_bray_72, md_plot, "Bray-Curtis")
uwuf_df <- make_pcoa_df(dist_uwuf_72, md_plot, "Unweighted UniFrac")

p_bray <- plot_pcoa(bray_df, "Bray-Curtis")
p_uwuf <- plot_pcoa(uwuf_df, "Unweighted UniFrac")
p_beta <- p_bray + p_uwuf
print(p_beta)

# ------------------------------------------------------------
# 6. Export figure, same dimensions and filenames
# ------------------------------------------------------------
ggsave(
  "Figure_beta_PMDD_Phase_PCoA.png",
  p_beta,
  width = 12,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Figure_beta_PMDD_Phase_PCoA.pdf",
  p_beta,
  width = 12,
  height = 5
)