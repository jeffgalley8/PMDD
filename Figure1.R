# Alpha diversity figure: Control vs PMDD by follicular/luteal phase

library(dplyr)
library(tidyr)
library(ggplot2)

# Rebuild plotting data from corrected metadata
# Assumes PMDD is coded C/P and Phase is coded F/L.
df_long <- md |>
  select(
    `sample-id`,
    enrollment_id,
    PMDD,
    Phase,
    shannon,
    observed,
    faith,
    evenness
  ) |>
  pivot_longer(
    cols = c(shannon, observed, faith, evenness),
    names_to = "metric",
    values_to = "value"
  )

df_long <- df_long |>
  mutate(
    metric = recode(
      metric,
      shannon = "Shannon",
      observed = "Observed features",
      faith = "Faith's PD",
      evenness = "Pielou's evenness"
    ),
    PMDD = factor(
      PMDD,
      levels = c("C", "P"),
      labels = c("Control", "PMDD")
    ),
    Phase = factor(Phase, levels = c("F", "L"))
  )

p1 <- ggplot(
  df_long,
  aes(x = PMDD, y = value, fill = Phase)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.6,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color = Phase),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 1.5,
    alpha = 0.7
  ) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_classic(base_size = 14) +
  labs(
    x = "Group",
    y = "Alpha diversity",
    fill = "Phase",
    color = "Phase"
  ) +
  scale_fill_manual(
    values = c("F" = "#0072B2", "L" = "#D55E00"),
    labels = c("F" = "Follicular", "L" = "Luteal")
  ) +
  scale_color_manual(
    values = c("F" = "#0072B2", "L" = "#D55E00"),
    labels = c("F" = "Follicular", "L" = "Luteal")
  )

print(p1)

ggsave(
  "Figure1_alpha_phase_PMDD.png",
  p1,
  width = 10,
  height = 6,
  dpi = 300
)