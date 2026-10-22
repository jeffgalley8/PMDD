library(vegan)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(patchwork)

# -------------------------
# 1) Load metadata
# -------------------------
md <- read_tsv(
  "FinalMap.tsv",
  show_col_types = FALSE
)

if ("#SampleID" %in% names(md)) {
  md <- md |>
    rename(`sample-id` = `#SampleID`)
}

md <- md |>
  mutate(
    PMDD = case_when(
      PMDD == "C" ~ "Control",
      PMDD == "P" ~ "PMDD",
      TRUE ~ as.character(PMDD)
    ),
    PMDD = factor(PMDD, levels = c("Control", "PMDD")),
    DRSP_tot = as.numeric(drsp_tot),
    DRSP_mood = as.numeric(drsp_mood)
  )

# -------------------------
# 2) Load Bray-Curtis distance matrix
# -------------------------
read_qiime_dist <- function(path) {
  m <- read.delim(
    path,
    check.names = FALSE,
    row.names = 1
  )
  m <- as.matrix(m)
  storage.mode(m) <- "numeric"
  as.dist(m)
}

dist_bray <- read_qiime_dist(
  "braycurtis/distance-matrix.tsv"
)

# -------------------------
# 3) Align samples
# -------------------------
samps <- intersect(
  labels(dist_bray),
  md$`sample-id`
)

dist_sub <- as.dist(
  as.matrix(dist_bray)[samps, samps]
)

md_sub <- md |>
  filter(`sample-id` %in% samps) |>
  arrange(match(`sample-id`, samps))

stopifnot(all(md_sub$`sample-id` == samps))

# -------------------------
# 4) Calculate 2D centroid-distance proxy
# -------------------------
coords <- cmdscale(dist_sub, k = 2)
centroid <- colMeans(coords)

dist_to_centroid <- sqrt(
  rowSums(
    (
      coords -
        matrix(
          centroid,
          nrow(coords),
          2,
          byrow = TRUE
        )
    )^2
  )
)

# Preserve the literal sample-id column name.
cd_df <- data.frame(
  `sample-id` = rownames(coords),
  dist_to_centroid = dist_to_centroid,
  check.names = FALSE
)

plot_df <- cd_df |>
  left_join(md_sub, by = "sample-id")

# -------------------------
# 5) Make DRSP total and mood plots
# -------------------------
p_drsp_tot <- plot_df |>
  filter(!is.na(DRSP_tot)) |>
  ggplot(aes(x = DRSP_tot, y = dist_to_centroid)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#0072B2") +
  theme_classic(base_size = 14) +
  labs(
    x = "DRSP Total",
    y = "Distance to centroid (Bray-Curtis)",
    title = "DRSP total"
  )

p_drsp_mood <- plot_df |>
  filter(!is.na(DRSP_mood)) |>
  ggplot(aes(x = DRSP_mood, y = dist_to_centroid)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#0072B2") +
  theme_classic(base_size = 14) +
  labs(
    x = "DRSP Mood",
    y = "Distance to centroid (Bray-Curtis)",
    title = "DRSP mood"
  )

p_drsp_beta <- p_drsp_tot + p_drsp_mood
print(p_drsp_beta)

# -------------------------
# 6) Save figure
# -------------------------
ggsave(
  "DRSP_BrayCentroid.png",
  p_drsp_beta,
  width = 11,
  height = 5,
  dpi = 600,
  bg = "white"
)