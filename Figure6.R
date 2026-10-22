############################################################
# KEGG IL-8 x PMDD interaction figure
# Run after the statistical analysis has generated its CSV.
############################################################

library(dplyr)
library(ggplot2)
library(tidyr)

outdir <- "KEGG_cytokine_interactions"
input_csv <- file.path(outdir, "KEGG_cytokine_PMDD_interactions_SIGNIFICANT.csv")
stopifnot(file.exists(input_csv))

kegg_int_sig <- read.csv(input_csv, check.names = FALSE)
required <- c("cytokine", "pathway", "pathway_name", "pathway_class", "lfc")
stopifnot(all(required %in% names(kegg_int_sig)))

plot_il8_int <- kegg_int_sig |>
  filter(cytokine == "IL-8") |>
  mutate(
    pathway_label = ifelse(is.na(pathway_name) | pathway_name == "", pathway, pathway_name),
    pathway_class = replace_na(pathway_class, "Unclassified")
  )

if (nrow(plot_il8_int) == 0) {
  stop("No significant IL-8 x PMDD interactions in the input CSV.")
}

# Keep pathway IDs unique even when multiple pathways share a name.
plot_il8_int <- plot_il8_int |>
  arrange(lfc) |>
  mutate(pathway_ordered = factor(pathway, levels = unique(pathway)))

p_il8_int <- ggplot(
  plot_il8_int,
  aes(x = lfc, y = pathway_ordered, fill = pathway_class)
) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  scale_y_discrete(labels = setNames(plot_il8_int$pathway_label, plot_il8_int$pathway)) +
  labs(
    x = "IL-8 × PMDD interaction log-fold change",
    y = NULL,
    fill = "KEGG pathway class",
    title = "KEGG pathways exhibiting IL-8 × PMDD interactions"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_il8_int)

ggsave(file.path(outdir, "KEGG_IL8_PMDD_interaction_barplot.pdf"),
       plot = p_il8_int, width = 10, height = 8, units = "in")
ggsave(file.path(outdir, "KEGG_IL8_PMDD_interaction_barplot.png"),
       plot = p_il8_int, width = 10, height = 8, units = "in", dpi = 600)