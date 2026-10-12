library(tidyverse)
library(survey)
library(gt)
library(gtsummary)
library(broom)

## Setup

df <- readRDS("data/nhanes_L_analytic.rds")
df$WTDR2D[is.na(df$WTDR2D)] <- 0

## Build variables

df <- df |>
  mutate(
    choline_dens = 1000 * choline_mean / energy_mean,
    log_alt = log(LBXSATSI),
    log_ast = log(LBXSASSI),
    log_ggt = log(LBXSGTSI),
    race3 = fct_collapse(
      factor(RIDRETH3),
      white = "Non-Hispanic White",
      black = "Non-Hispanic Black",
      other_level = "other"
    ),
    chol_q = cut(
      choline_dens,
      breaks = quantile(
        choline_dens[in_aim2 & bmi_cat != "Underweight"],
        probs = seq(0, 1, 0.25),
        na.rm = TRUE
      ),
      include.lowest = TRUE,
      labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)")
    ),
    low_chol_obese = as.numeric(chol_q == "Q1 (Lowest)" & bmi_cat == "Obese"),
    bmi_chol_group = factor(
      paste0(bmi_cat, ", ", 
             if_else(chol_q == "Q1 (Lowest)", "low choline", "higher choline")),
      levels = c("Normal, higher choline", "Normal, low choline",
                 "Overweight, higher choline", "Overweight, low choline",
                 "Obese, higher choline", "Obese, low choline")
    )
  ) 
## Quartile cutoffs
quantile(df$choline_dens[df$in_aim2 & df$bmi_cat != "Underweight"],
         probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

## Survey design

des_full <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDR2D,
  data = df,
  nest = TRUE
)

des2 <- subset(
  des_full,
  in_aim2 %in% TRUE & bmi_cat != "Underweight" & !is.na(choline_dens)
)
des2 <- update(des2, bmi_cat = fct_drop(bmi_cat))

n_total <- nrow(des2$variables)
cat("Design df:", degf(des2), "\n")

## Sample sizes

des2$variables |>
  count(bmi_cat, chol_q) |>
  pivot_wider(names_from = chol_q, values_from = n) |>
  print()

## Table 1: liver enzymes by BMI category and choline quartile

geo_means <- function(outcome) {
  svyby(reformulate(outcome), ~bmi_cat + chol_q, des2, svymean,
        na.rm = TRUE, vartype = "ci") |>
    as.data.frame() |>
    mutate(value = sprintf("%.1f (%.1f, %.1f)",
                           exp(.data[[outcome]]), exp(ci_l), exp(ci_u))) |>
    select(bmi_cat, chol_q, value)
}

tbl_1 <- bind_rows(
  ALT = geo_means("log_alt"),
  AST = geo_means("log_ast"),
  GGT = geo_means("log_ggt"),
  .id = "enzyme"
) |>
  pivot_wider(names_from = chol_q, values_from = value) |>
  gt(groupname_col = "enzyme", rowname_col = "bmi_cat") |>
  tab_spanner("Choline quartile (mg/1,000 kcal)", columns = starts_with("Q")) |>
  tab_source_note("Survey-weighted geometric mean (95% CI), U/L.")

## Main test: low choline and obese vs. all
## Running for each liver enzyme

fit_enzyme <- function(outcome, exposure = "low_chol_obese") {
  f <- as.formula(
    paste(outcome, "~", exposure, "+ RIDAGEYR + sex + race3 + INDFMPIR")
  )
  svyglm(f, design = des2)
}

m_alt <- fit_enzyme("log_alt")
m_ast <- fit_enzyme("log_ast")
m_ggt <- fit_enzyme("log_ggt")

summary(m_alt)
summary(m_ast)
summary(m_ggt)

## sample sizes for each model
nrow(model.frame(m_alt))
nrow(model.frame(m_ast))
nrow(model.frame(m_ggt))

## BMI category by low vs. higher choline

m_alt_groups <- fit_enzyme("log_alt", "bmi_chol_group")
m_ast_groups <- fit_enzyme("log_ast", "bmi_chol_group")
m_ggt_groups <- fit_enzyme("log_ggt", "bmi_chol_group")

## Enzymes differ across the six groups?

regTermTest(m_alt_groups, ~bmi_chol_group)
regTermTest(m_ast_groups, ~bmi_chol_group)
regTermTest(m_ggt_groups, ~bmi_chol_group)

## Within each BMI category: low vs. higher choline
des2 <- update(des2, chol_low = as.numeric(chol_q == "Q1 (Lowest)"))

fit_within <- function(outcome) {
  f <- as.formula(paste(outcome, "~ bmi_cat * chol_low + RIDAGEYR + sex + race3 + INDFMPIR"))
  svyglm(f, design = des2)
}

m_alt_within <- fit_within("log_alt")
m_ast_within <- fit_within("log_ast")
m_ggt_within <- fit_within("log_ggt")

## cholinen between each BMI category?
regTermTest(m_alt_within, ~chol_low + bmi_cat:chol_low)
regTermTest(m_ast_within, ~chol_low + bmi_cat:chol_low)
regTermTest(m_ggt_within, ~chol_low + bmi_cat:chol_low)

## Tables 2 and 3: adjusted ratios for each enzyme

enzyme_table <- function(models, var, label) {
  models |>
    map(\(m) tbl_regression(m, include = all_of(var), exponentiate = TRUE,
                            label = setNames(list(label), var),
                            pvalue_fun = function(x) style_pvalue(x, digits = 2)) |>
          modify_header(estimate ~ "**Ratio**")) |>
    tbl_merge(tab_spanner = c("**ALT**", "**AST**", "**GGT**"))
}

tbl_2 <- enzyme_table(list(m_alt, m_ast, m_ggt),
                      "low_chol_obese", "Low choline and obesity")

tbl_3 <- enzyme_table(list(m_alt_groups, m_ast_groups, m_ggt_groups),
                      "bmi_chol_group", "BMI category and choline intake")

## Save

saveRDS(tbl_1, "output/aim2_table_enzymes.rds")
saveRDS(tbl_2, "output/aim2_table_contrast.rds")
saveRDS(tbl_3, "output/aim2_table_groups.rds")

gtsave(tbl_1, "output/aim2_table_enzymes.png", vwidth = 1000, zoom = 2)
gtsave(as_gt(tbl_2), "output/aim2_table_contrast.png", vwidth = 1000, zoom = 2)
gtsave(as_gt(tbl_3), "output/aim2_table_groups.png", vwidth = 1000, zoom = 2)

#forest plot figure liver enzymes across six groups

get_ratios <- function(model, enzyme) {
  est <- coef(model)
  ci <- confint(model)
  keep <- str_starts(names(est), "bmi_chol_group")
  tibble(
    enzyme = enzyme,
    group = str_remove(names(est)[keep], "bmi_chol_group"),
    ratio = exp(est[keep]),
    low = exp(ci[keep, 1]),
    high = exp(ci[keep, 2])
  )
}

## reference group
reference <- tibble(
  enzyme = c("ALT", "AST", "GGT"),
  group = "Normal, higher choline",
  ratio = 1, low = 1, high = 1
)

## labels for the legend
low_label <- "Low choline: lowest quarter of intake"
high_label <- "Higher choline: top three quarters of intake"

plot_data <- bind_rows(
  get_ratios(m_alt_groups, "ALT"),
  get_ratios(m_ast_groups, "AST"),
  get_ratios(m_ggt_groups, "GGT"),
  reference
) |>
  mutate(
    bmi = factor(str_extract(group, "^[A-Za-z]+"),
                 levels = c("Obese", "Overweight", "Normal")),
    choline = factor(if_else(str_detect(group, "low"), low_label, high_label),
                     levels = c(low_label, high_label))
  )

p_groups <- ggplot(plot_data,
                   aes(x = bmi, y = ratio, ymin = low, ymax = high,
                       color = choline, shape = choline)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_pointrange(position = position_dodge(width = 0.5),
                  size = 0.5, linewidth = 0.7) +
  coord_flip() +
  facet_wrap(~enzyme, nrow = 1) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4)) +
  scale_color_manual(values = setNames(c("#eb6834", "#2a78d6"),
                                       c(low_label, high_label))) +
  labs(
    x = "BMI category",
    y = "Liver enzyme levels compared with normal weight, higher choline (ratio)",
    color = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    panel.spacing = unit(1.5, "lines"),
    panel.grid.minor = element_blank()
  )

ggsave("output/aim2_plot_groups.png", p_groups,
       width = 9, height = 4, dpi = 300, bg = "white")