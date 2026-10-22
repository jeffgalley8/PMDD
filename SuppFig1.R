# Supplemental Figure 1: Alpha diversity model coefficient heatmap
# Uses metadata/alpha object `md` already loaded in R.
# Fits standardized-outcome models for the figure's comparable beta colors.

library(dplyr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggplot2)
library(purrr)

# Ensure reference levels and coding match the analysis.
# Expected PMDD levels C (Control), P (PMDD); Phase levels F, L.
stopifnot(all(c("PMDD", "Phase", "enrollment_id", "shannon", "observed", "faith", "evenness", "demos_age", "bmi", "Kcal", "Fibe", "omega") %in% names(md)))
df_alpha_z <- md |>
  mutate(
    PMDD = factor(as.character(PMDD), levels = c("C", "P")),
    Phase = factor(as.character(Phase), levels = c("F", "L")),
    shannon_z = as.numeric(scale(shannon)),
    observed_z = as.numeric(scale(observed)),
    faith_z = as.numeric(scale(faith)),
    evenness_z = as.numeric(scale(evenness))
  )

fit_alpha_models <- function(outcome) {
  m0 <- lmer(
    as.formula(paste0(outcome, " ~ PMDD * Phase + (1 | enrollment_id)")),
    data = df_alpha_z, REML = FALSE
  )
  m1 <- lmer(
    as.formula(paste0(outcome, " ~ PMDD * Phase + demos_age + bmi + (1 | enrollment_id)")),
    data = df_alpha_z, REML = FALSE
  )
  m2 <- lmer(
    as.formula(paste0(outcome, " ~ PMDD * Phase + demos_age + bmi + Kcal + Fibe + omega + (1 | enrollment_id)")),
    data = df_alpha_z, REML = FALSE
  )
  list(M0 = m0, M1 = m1, M2 = m2)
}

alpha_models_z <- list(
  Shannon = fit_alpha_models("shannon_z"),
  Observed = fit_alpha_models("observed_z"),
  Faith = fit_alpha_models("faith_z"),
  Evenness = fit_alpha_models("evenness_z")
)

coef_df <- imap_dfr(alpha_models_z, function(metric_models, metric_name) {
  imap_dfr(metric_models, function(mod, model_name) {
    broom.mixed::tidy(mod, effects = "fixed") |>
      mutate(metric = metric_name, model = model_name)
  })
})

coef_df <- coef_df |>
  filter(term %in% c("PhaseL", "PMDDP", "PMDDP:PhaseL", "demos_age", "bmi", "Kcal", "Fibe", "omega")) |>
  mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    term_clean = recode(
      term,
      PhaseL = "Phase",
      PMDDP = "PMDD",
      `PMDDP:PhaseL` = "PMDD × Phase",
      demos_age = "Age",
      bmi = "BMI",
      Kcal = "Energy intake",
      Fibe = "Fiber",
      omega = "Omega-3"
    ),
    metric = factor(metric, levels = c("Shannon", "Observed", "Faith", "Evenness")),
    model = factor(model, levels = c("M0", "M1", "M2"))
  )

p_heat <- ggplot(coef_df, aes(x = model, y = metric, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), size = 5) +
  facet_wrap(~ term_clean) +
  scale_fill_gradient2(
    low = "#4575B4", mid = "white", high = "#D73027",
    midpoint = 0, name = "Standardized β"
  ) +
  theme_classic(base_size = 14) +
  labs(x = "Model", y = "Alpha diversity metric") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_heat)

ggsave(
  "Figure2_alpha_model_heatmap.png",
  p_heat,
  width = 11,
  height = 6,
  dpi = 600,
  bg = "white"
)