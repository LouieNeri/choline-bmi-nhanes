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

## choline between each BMI category?
regTermTest(m_alt_within, ~chol_low + bmi_cat:chol_low)
regTermTest(m_ast_within, ~chol_low + bmi_cat:chol_low)
regTermTest(m_ggt_within, ~chol_low + bmi_cat:chol_low)

## Save objects for 05_tables_figures.R
save(des2, m_alt, 
     m_ast, m_ggt, 
     m_alt_groups, 
     m_ast_groups, 
     m_ggt_groups,
     file = "output/aim2_objects.RData")
