library(tidyverse)
library(survey)

# Setup and read data

df <- readRDS("data/nhanes_L_analytic.rds")

# Survey design on the full sample then subset for aim 1
df$WTDR2D[is.na(df$WTDR2D)] <- 0

des_full <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDR2D,
  data = df,
  nest = TRUE
)

des <- subset(des_full,
              in_aim1 %in% TRUE & bmi_cat 
              != "Underweight")

des <- update(des,
              bmi_cat = fct_drop(bmi_cat),
              DMDEDUC2 = fct_drop(DMDEDUC2),
              RIDRETH3 = fct_drop(RIDRETH3))

ana <- des$variables

# Overall prevalence below Ai

overall_prev <- svymean(~below_ai, 
                        des, 
                        na.rm = TRUE)

overall_prev
confint(overall_prev)

# Weighted population size
sum(weights(des))

# proportions AI differ by BMI category
svychisq(~below_ai + bmi_cat, des)
svychisq(~below_ai + bmi_cat, subset(des, sex == "Female"))
svychisq(~below_ai + bmi_cat, subset(des, sex == "Male"))


# Sex and BMI category

# Counts below the AI by sex and BMI category
sex_counts <- ana |>
  filter(!is.na(below_ai), 
         !is.na(sex), 
         !is.na(bmi_cat)) |>
  group_by(sex, bmi_cat) |>
  summarise(n_below_ai = sum(below_ai == 1), .groups = "drop")

# Total sample size by sex
sex_totals <- ana |>
  filter(!is.na(below_ai), 
         !is.na(sex), 
         !is.na(bmi_cat)) |>
  count(sex, name = "N")

Female_N <- sex_totals |> 
  filter(sex == "Female") |> 
  pull(N)

Male_N <- sex_totals |> 
  filter(sex == "Male") |> 
  pull(N)

# Weighted prevalence and confidence intervals
prev_sex_bmi <- svyby(
  ~below_ai,
  ~sex + bmi_cat,
  des,
  svymean,
  na.rm = TRUE,
  vartype = "ci"
) |>
  as.data.frame() |>
  transmute(
    sex,
    bmi_cat,
    prevalence = 100 * as.numeric(below_ai),
    ci_lower = 100 * as.numeric(ci_l),
    ci_upper = 100 * as.numeric(ci_u)
  ) |>
  left_join(sex_counts, by = c("sex", "bmi_cat")) |>
  mutate(
    cell = sprintf(
      "%d, %.1f%% (%.1f–%.1f)",
      n_below_ai,
      prevalence,
      ci_lower,
      ci_upper
    )
  )



# By BMI category

# Counts below the AI by BMI category
bmi_counts <- ana |>
  filter(!is.na(below_ai), 
         !is.na(bmi_cat)) |>
  group_by(bmi_cat) |>
  summarise(n_below_ai = sum(below_ai == 1), 
            .groups = "drop")

# Total analytic sample size
bmi_N <- ana |>
  filter(!is.na(below_ai), 
         !is.na(bmi_cat)) |>
  nrow()

prev_bmi <- svyby(
  ~below_ai,
  ~bmi_cat,
  des,
  svymean,
  na.rm = TRUE,
  vartype = "ci"
) |>
  as.data.frame() |>
  transmute(
    bmi_cat,
    prevalence = 100 * as.numeric(below_ai),
    ci_lower = 100 * as.numeric(ci_l),
    ci_upper = 100 * as.numeric(ci_u)
  ) |>
  left_join(bmi_counts, by = "bmi_cat") |>
  mutate(
    overall = sprintf(
      "%d, %.1f%% (%.1f–%.1f)",
      n_below_ai,
      prevalence,
      ci_lower,
      ci_upper
    )
  )

# Add an overall row for total BMI
overall_row <- tibble(
  bmi_cat = "All",
  prevalence = 100 * as.numeric(coef(overall_prev)),
  ci_lower = 100 * confint(overall_prev)[1],
  ci_upper = 100 * confint(overall_prev)[2],
  n_below_ai = sum(bmi_counts$n_below_ai)
) |>
  mutate(
    overall = sprintf(
      "%d, %.1f%% (%.1f–%.1f)",
      n_below_ai,
      prevalence,
      ci_lower,
      ci_upper
    )
  )

prev_bmi_all <- prev_bmi |>
  mutate(bmi_cat = as.character(bmi_cat)) |>
  bind_rows(overall_row)

# Aim 1: Women by menopause and BMI

des_f <- subset(des, sex == "Female" & menopause %in% c("pre", "post"))

ana_f <- des_f$variables

# BMI differ by menopausal status?
m_menopause <- svyglm(below_ai ~ bmi_cat * menopause,
                      design = des_f,
                      family = quasibinomial())
regTermTest(m_menopause, ~bmi_cat:menopause)

# Counts below the AI by menopause status and BMI category
menopause_counts <- ana_f |>
  filter(!is.na(below_ai), 
         !is.na(menopause), 
         !is.na(bmi_cat)) |>
  group_by(menopause, 
           bmi_cat) |>
  summarise(n_below_ai = sum(below_ai == 1), 
            .groups = "drop")

# Total sample size by menopause status
menopause_totals <- ana_f |>
  filter(!is.na(below_ai), 
         !is.na(menopause), 
         !is.na(bmi_cat)) |>
  count(menopause, name = "N")

pre_N <- menopause_totals |> 
  filter(menopause == "pre") |> 
  pull(N)

post_N <- menopause_totals |> 
  filter(menopause == "post") |> 
  pull(N)

prev_menopause <- svyby(
  ~below_ai,
  ~menopause + bmi_cat,
  des_f,
  svymean,
  na.rm = TRUE,
  vartype = "ci"
) |>
  as.data.frame() |>
  transmute(
    menopause,
    bmi_cat,
    prevalence = 100 * as.numeric(below_ai),
    ci_lower = 100 * as.numeric(ci_l),
    ci_upper = 100 * as.numeric(ci_u)
  ) |>
  left_join(menopause_counts, by = c("menopause", "bmi_cat")) |>
  mutate(
    cell = sprintf(
      "%d, %.1f%% (%.1f–%.1f)",
      n_below_ai,
      prevalence,
      ci_lower,
      ci_upper
    )
  )

# overall for menopause
menopause_all <- svyby(
  ~below_ai,
  ~menopause,
  des_f,
  svymean,
  na.rm = TRUE,
  vartype = "ci"
) |>
  as.data.frame() |>
  transmute(
    menopause,
    bmi_cat = "All",
    prevalence = 100 * as.numeric(below_ai),
    ci_lower = 100 * as.numeric(ci_l),
    ci_upper = 100 * as.numeric(ci_u)
  ) |>
  left_join(
    menopause_counts |>
      group_by(menopause) |>
      summarise(n_below_ai = sum(n_below_ai)),
    by = "menopause"
  ) |>
  mutate(
    cell = sprintf(
      "%d, %.1f%% (%.1f–%.1f)",
      n_below_ai,
      prevalence,
      ci_lower,
      ci_upper
    )
  )

prev_menopause_all <- prev_menopause |>
  mutate(bmi_cat = as.character(bmi_cat)) |>
  bind_rows(menopause_all)

# Save objects for 05_tables_figures.R
save(des, ana, prev_sex_bmi, 
     Female_N, Male_N, prev_bmi, prev_bmi_all, bmi_N,
     prev_menopause, prev_menopause_all, pre_N, post_N,
     file = "output/aim1_objects.RData")
