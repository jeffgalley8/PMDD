# IL-8 versus Bray-Curtis PCoA1

library(ape)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)

# 1. Load metadata
md <- read_tsv("metadata.tsv", show_col_types = FALSE)
if ("#SampleID" %in% names(md)) {
  md <- md |> rename(`sample-id` = `#SampleID`)
}
md <- md |>
  mutate(
    PMDD = case_when(
      PMDD == "C" ~ "Control",
      PMDD == "P" ~ "PMDD",
      TRUE ~ as.character(PMDD)
    ),
    PMDD = factor(PMDD, levels = c("Control", "PMDD")),
    IL8 = as.numeric(blood_il8)
  )

# 2. Load Bray-Curtis distance matrix
read_qiime_dist <- function(path) {
  m <- read.delim(path, check.names = FALSE, row.names = 1)
  m <- as.matrix(m)
  storage.mode(m) <- "numeric"
  as.dist(m)
}
dist_bray <- read_qiime_dist("braycurtis/distance-matrix.tsv")

# 3. Align samples
samps <- intersect(labels(dist_bray), md$`sample-id`)
cat("Overlapping samples:", length(samps), "\n")
dist_sub <- as.dist(as.matrix(dist_bray)[samps, samps])
md_sub <- md |>
  filter(`sample-id` %in% samps) |>
  arrange(match(`sample-id`, samps))
stopifnot(all(md_sub$`sample-id` == samps))

# 4. PCoA
pcoa <- ape::pcoa(dist_sub)
scores <- as.data.frame(pcoa$vectors[, 1:2])
scores$`sample-id` <- rownames(scores)
plot_df <- scores |>
  left_join(md_sub, by = "sample-id") |>
  filter(!is.na(PMDD), !is.na(IL8))
var1 <- round(pcoa$values$Relative_eig[1] * 100, 1)

# 5. Plot
p_il8 <- ggplot(plot_df, aes(x = IL8, y = Axis.1, color = PMDD)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = c("Control" = "#999999", "PMDD" = "#D55E00")) +
  labs(
    x = "IL-8 (pg/mL)",
    y = paste0("PCoA1 (Bray-Curtis, ", var1, "%)"),
    color = "Group",
    title = str_wrap("Association between IL-8 and microbiome composition", 45)
  )
print(p_il8)

# 6. Export
ggsave(
  "IL8_PCoA_scatter.png",
  plot = p_il8,
  width = 6,
  height = 5,
  dpi = 600,
  bg = "white"
)