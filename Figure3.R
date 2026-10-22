# Cytokine–taxa recurrent heatmap: FIGURE 3
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)

main_results <- read.csv("cytokine_taxa_main_ALL.csv", check.names = FALSE)
int_results <- read.csv("cytokine_x_PMDD_taxa_ALL.csv", check.names = FALSE)

main_plot <- main_results |> mutate(effect = cytokine, analysis = "Main")
int_plot <- int_results |> mutate(effect = paste0(cytokine, " × PMDD"), analysis = "Interaction")
all_plot <- bind_rows(main_plot, int_plot)

effect_order <- c(
  "IL-1β", "IL-6", "IL-8", "TNFα",
  "IL-1β × PMDD", "IL-6 × PMDD", "IL-8 × PMDD", "TNFα × PMDD"
)
stopifnot(all(as.character(all_plot$effect) %in% effect_order))

# Significant in at least two of the eight analyses (49 taxa in original run).
recurrent_taxa <- all_plot |>
  filter(significant == TRUE) |>
  count(taxon, name = "n_sig") |>
  filter(n_sig >= 2)

heat_dat <- all_plot |>
  semi_join(recurrent_taxa, by = "taxon") |>
  left_join(recurrent_taxa, by = "taxon")

# Use the lowest resolved taxonomic rank as a display label.
get_taxon_label <- function(x) {
  vapply(x, function(tx) {
    ranks <- str_split(tx, ";", simplify = FALSE)[[1]]
    cleaned <- str_remove(ranks, "^[dkpcofgs]__")
    cleaned <- cleaned[!is.na(cleaned) & cleaned != "" & cleaned != "__"]
    if (length(cleaned) == 0) return("Unclassified")
    tail(cleaned, 1)
  }, character(1), USE.NAMES = FALSE)
}

heat_dat <- heat_dat |> mutate(taxon_label = get_taxon_label(taxon))

# Disambiguate duplicate display labels without combining distinct taxa.
label_map <- heat_dat |>
  distinct(taxon, taxon_label) |>
  group_by(taxon_label) |>
  mutate(
    dup_id = row_number(),
    dup_n = n(),
    taxon_label_unique = ifelse(
      dup_n > 1, paste0(taxon_label, "_", dup_id), taxon_label
    )
  ) |>
  ungroup() |>
  select(taxon, taxon_label_unique)

heat_dat <- heat_dat |>
  select(-taxon_label) |>
  left_join(label_map, by = "taxon") |>
  rename(taxon_label = taxon_label_unique) |>
  mutate(
    effect = factor(effect, levels = effect_order),
    plot_lfc = ifelse(significant == TRUE, lfc, NA_real_)
  )

# Extract phylum/class/order/family, tolerating incomplete taxonomy strings.
extract_rank <- function(x, prefix) {
  vapply(x, function(tx) {
    ranks <- str_split(tx, ";", simplify = FALSE)[[1]]
    hit <- ranks[startsWith(ranks, paste0(prefix, "__"))]
    if (length(hit) == 0) return(NA_character_)
    value <- sub(paste0("^", prefix, "__"), "", hit[1])
    if (is.na(value) || value %in% c("", "__")) NA_character_ else value
  }, character(1), USE.NAMES = FALSE)
}

tax_order_df <- heat_dat |>
  distinct(taxon, taxon_label) |>
  mutate(
    Phylum = extract_rank(taxon, "p"),
    Class = extract_rank(taxon, "c"),
    Order = extract_rank(taxon, "o"),
    Family = extract_rank(taxon, "f")
  ) |>
  arrange(Phylum, Class, Order, Family, taxon_label)

heat_mat <- heat_dat |>
  select(taxon_label, effect, plot_lfc) |>
  pivot_wider(names_from = effect, values_from = plot_lfc) |>
  as.data.frame()
rownames(heat_mat) <- heat_mat$taxon_label
heat_mat$taxon_label <- NULL
heat_mat <- as.matrix(heat_mat)
heat_mat <- heat_mat[, effect_order, drop = FALSE]

tax_meta <- tax_order_df |>
  select(taxon_label, Phylum, Class, Order, Family) |>
  distinct()
tax_meta <- tax_meta[match(rownames(heat_mat), tax_meta$taxon_label), ]

tax_ranks <- tax_meta |>
  select(Phylum, Class, Order, Family) |>
  mutate(across(everything(), ~ ifelse(is.na(.x), "Unclassified", .x)))

n <- nrow(tax_ranks)
tax_dist_mat <- matrix(0, nrow = n, ncol = n)
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    same <- as.character(tax_ranks[i, ]) == as.character(tax_ranks[j, ])
    tax_dist_mat[i, j] <- sum(!same)
  }
}
rownames(tax_dist_mat) <- tax_meta$taxon_label
colnames(tax_dist_mat) <- tax_meta$taxon_label
tax_dist <- as.dist(tax_dist_mat)
tax_hclust <- hclust(tax_dist, method = "complete")
tax_dend <- as.dendrogram(tax_hclust)

col_fun <- colorRamp2(
  c(min(heat_mat, na.rm = TRUE), 0, max(heat_mat, na.rm = TRUE)),
  c("#2166AC", "white", "#B2182B")
)

ht <- Heatmap(
  heat_mat,
  name = "LFC",
  col = col_fun,
  na_col = "white",
  cluster_rows = tax_dend,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_dend_side = "left",
  show_row_dend = TRUE,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 9),
  column_names_gp = grid::gpar(fontsize = 11),
  column_names_rot = 45,
  border = TRUE
)

draw(ht)

png(
  "Cytokine_Taxa_Recurrent_Heatmap_Dendrogram.png",
  width = 10, height = 15, units = "in", res = 600
)
draw(ht)
dev.off()

pdf(
  "Cytokine_Taxa_Recurrent_Heatmap_Dendrogram.pdf",
  width = 10, height = 15
)
draw(ht)
dev.off()
