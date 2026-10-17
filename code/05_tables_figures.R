## Tables and figures for Aims 1 and 2
## Run 03 and 04 first

library(tidyverse)
library(survey)
library(gt)
library(gtsummary)
library(ggthemes)

load("output/aim1_objects.RData")
load("output/aim2_objects.RData")

# BMI category labels
bmi_labels <- c(
  "Normal" = "Normal weight (18.5–24.9)",
  "Overweight" = "Overweight (25.0–29.9)",
  "Obese" = "Obesity (≥30.0)",
  "All" = "All BMI categories"
)

# fixing numbers for tables
comma <- scales::comma
fmt_ci <- function(est,
                   low,
                   high) sprintf("%.1f (%.1f, %.1f)",
                                 est,
                                 low,
                                 high)
fmt_p <- function(p) style_pvalue(p, digits = 2)

## Table 1: baseline characteristics by sex

tbl_1 <- des |>
  tbl_svysummary(
    by = sex,
    include = c(RIDAGEYR,
                RIDRETH3,
                DMDEDUC2,
                INDFMPIR,
                bmi_cat,
                energy_mean,
                choline_mean, below_ai),
    label = list(
      RIDAGEYR ~ "Age, years",
      RIDRETH3 ~ "Race and ethnicity",
      DMDEDUC2 ~ "Education",
      INDFMPIR ~ "Income-to-poverty ratio",
      bmi_cat ~ "BMI category",
      energy_mean ~ "Energy intake, kcal/day",
      choline_mean ~ "Choline intake, mg/day",
      below_ai ~ "Below choline AI"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n_unweighted} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    missing = "no"
  ) |>
  add_overall(last = TRUE) |>
  add_p(test = all_continuous() ~ "svy.t.test") |> # t-test to match mean (SD)
  modify_header(
    label ~ "Characteristic",
    stat_1 ~ paste0("Female (n = ", comma(sum(ana$sex == "Female")), ")"),
    stat_2 ~ paste0("Male (n = ", comma(sum(ana$sex == "Male")), ")"),
    stat_0 ~ paste0("Total (n = ", comma(nrow(ana)), ")"),
    p.value ~ "P value"
  ) |>
  modify_spanning_header(c(stat_1, stat_2) ~ "Sex") |>
  modify_footnote(
    all_stat_cols() ~ "Unweighted n (survey-weighted %) or survey-weighted mean (SD)."
  ) |>
  as_gt() |>
  tab_source_note("AI, Adequate Intake (550 mg/day for men, 425 mg/day for women); BMI, body mass index.")

## Table 2: % below the AI by sex and BMI

# % below AI for each sex, all BMI categories combined
sex_all <- svyby(~below_ai,
                 ~sex,
                 des,
                 svymean,
                 na.rm = TRUE,
                 vartype = "ci") |>
  as.data.frame() |>
  transmute(sex,
            bmi_cat = "All",
            prevalence = 100 * below_ai,
            ci_lower = 100 * ci_l,
            ci_upper = 100 * ci_u) |>
  left_join(count(ana,
                  sex,
                  wt = below_ai,
                  name = "n_below_ai"),
            by = "sex")

# participants in each cell
n_cells <- bind_rows(
  count(ana,
        sex,
        bmi_cat),
  count(ana,
        sex) |>
    mutate(bmi_cat = "All"),
  count(ana,
        bmi_cat) |>
    mutate(sex = "Total"),
  tibble(sex = "Total",
         bmi_cat = "All",
         n = nrow(ana))
) |>
  mutate(bmi_cat = as.character(bmi_cat))

# Rao-Scott chi-square across BMI categories (same tests as 03)
p_female <- svychisq(~below_ai +
                       bmi_cat,
                     subset(des,
                            sex == "Female"))$p.value
p_male <- svychisq(~below_ai +
                     bmi_cat,
                   subset(des, sex == "Male"))$p.value
p_total <- svychisq(~below_ai +
                      bmi_cat, des)$p.value

tbl_2 <- bind_rows(
  prev_sex_bmi |>
    mutate(bmi_cat =
             as.character(bmi_cat)), sex_all,
  prev_bmi_all |>
    mutate(sex = "Total")
) |>
  left_join(n_cells, by = c("sex", "bmi_cat")) |>
  mutate(
    bmi_cat = factor(bmi_cat,
                     levels = names(bmi_labels),
                     labels = bmi_labels),
    n = paste0(comma(n_below_ai), "/", comma(n)),
    pct = fmt_ci(prevalence,
                 ci_lower,
                 ci_upper)
  ) |>
  arrange(bmi_cat) |>
  select(sex,
         bmi_cat,
         n,
         pct) |>
  pivot_wider(names_from = sex,
              values_from = c(n, pct),
              names_vary = "slowest") |>
  mutate(bmi_cat = as.character(bmi_cat)) |>
  bind_rows(tibble(bmi_cat = "P value",
                   pct_Female = fmt_p(p_female),
                   pct_Male = fmt_p(p_male),
                   pct_Total = fmt_p(p_total))) |>
  gt() |>
  cols_label(
    bmi_cat = "BMI category (kg/m²)",
    n_Female = "n/N",
    pct_Female = "% (95% CI)",
    n_Male = "n/N",
    pct_Male = "% (95% CI)",
    n_Total = "n/N",
    pct_Total = "% (95% CI)"
  ) |>
  tab_spanner(paste0("Female (n = ", comma(Female_N), ")"),
              columns = c(n_Female, pct_Female)) |>
  tab_spanner(paste0("Male (n = ",
                     comma(Male_N), ")"), columns = c(n_Male, pct_Male)) |>
  tab_spanner(paste0("Total (n = ",
                     comma(bmi_N), ")"), columns = c(n_Total, pct_Total)) |>
  sub_missing(missing_text = "") |>
  tab_footnote("Unweighted number below the AI / number of participants.",
               locations = cells_column_labels(columns = starts_with("n_"))) |>
  tab_footnote("Survey-weighted percentage.",
               locations = cells_column_labels(columns = starts_with("pct_"))) |>
  tab_footnote("Rao-Scott chi-square test across BMI categories.",
               locations = cells_body(columns = bmi_cat, rows = bmi_cat == "P value")) |>
  tab_source_note("AI, Adequate intake (550 mg/day for men, 425 mg/day for women); BMI, body mass index; CI, confidence interval.")

## Table 3: % below the AI by menopausal status and BMI

# participants in each cell
ana_f <- filter(ana,
                sex == "Female",
                menopause %in% c("pre", "post"))

n_cells_f <- bind_rows(
  count(ana_f,
        menopause,
        bmi_cat),
  count(ana_f,
        menopause) |>
    mutate(bmi_cat = "All")) |>
  mutate(bmi_cat = as.character(bmi_cat))

# BMI × menopause interaction (same test as 03)
des_f <- subset(des,
                sex == "Female" & menopause %in% c("pre", "post"))
m_menopause <- svyglm(below_ai ~ bmi_cat *
                        menopause,
                      design = des_f,
                      family = quasibinomial())
p_interaction <- as.numeric(regTermTest(m_menopause, ~bmi_cat:menopause)$p)

tbl_3 <- prev_menopause_all |>
  left_join(n_cells_f, by = c("menopause", "bmi_cat")) |>
  mutate(
    bmi_cat = factor(bmi_cat, levels = names(bmi_labels), labels = bmi_labels),
    n = paste0(comma(n_below_ai), "/", comma(n)),
    pct = fmt_ci(prevalence, ci_lower, ci_upper)
  ) |>
  arrange(bmi_cat) |>
  select(menopause, bmi_cat, n, pct) |>
  pivot_wider(names_from = menopause, values_from = c(n, pct), names_vary = "slowest") |>
  mutate(bmi_cat = as.character(bmi_cat)) |>
  bind_rows(tibble(bmi_cat = "P for interaction",
                   pct_post = fmt_p(p_interaction))) |>
  gt() |>
  cols_label(
    bmi_cat = "BMI category (kg/m²)",
    n_pre = "n/N", pct_pre = "% (95% CI)",
    n_post = "n/N", pct_post = "% (95% CI)"
  ) |>
  tab_spanner(paste0("Premenopausal (n = ", comma(pre_N), ")"), columns = c(n_pre, pct_pre)) |>
  tab_spanner(paste0("Postmenopausal (n = ", comma(post_N), ")"), columns = c(n_post, pct_post)) |>
  sub_missing(missing_text = "") |>
  tab_footnote("Unweighted number below the AI / number of participants.",
               locations = cells_column_labels(columns = starts_with("n_"))) |>
  tab_footnote("Survey-weighted percentage.",
               locations = cells_column_labels(columns = starts_with("pct_"))) |>
  tab_footnote("BMI × menopausal status interaction from survey-weighted logistic regression. Women with a hysterectomy or unknown status were excluded.",
               locations = cells_body(columns = bmi_cat, rows = bmi_cat == "P for interaction")) |>
  tab_source_note("AI, Adequate intake (425 mg/day for women); BMI, body mass index; CI, confidence interval.")

## Table 4: liver enzymes by BMI category and choline quartile

# quartile cutoffs for the column headers
cuts <- quantile(des2$variables$choline_dens, probs = c(0.25, 0.5, 0.75))

geo_means <- function(outcome) {
  svyby(reformulate(outcome), ~bmi_cat + chol_q, des2, svymean,
        na.rm = TRUE, vartype = "ci") |>
    as.data.frame() |>
    mutate(value = fmt_ci(exp(.data[[outcome]]), exp(ci_l), exp(ci_u))) |>
    select(bmi_cat, chol_q, value)
}

# participants in each cell
cell_n <- des2$variables |>
  count(bmi_cat, chol_q) |>
  mutate(value = comma(n)) |>
  select(bmi_cat, chol_q, value)

tbl_4 <- bind_rows(
  "ALT, U/L" = geo_means("log_alt"),
  "AST, U/L" = geo_means("log_ast"),
  "GGT, U/L" = geo_means("log_ggt"),
  "Participants, n" = cell_n,
  .id = "enzyme"
) |>
  mutate(bmi_cat = unname(bmi_labels[as.character(bmi_cat)])) |>
  pivot_wider(names_from = chol_q, values_from = value) |>
  gt(groupname_col = "enzyme") |>
  cols_label(
    bmi_cat = "BMI category (kg/m²)",
    "Q1 (Lowest)" = sprintf("Q1 (<%.1f)", cuts[1]),
    "Q2" = sprintf("Q2 (%.1f–%.1f)", cuts[1], cuts[2]),
    "Q3" = sprintf("Q3 (%.1f–%.1f)", cuts[2], cuts[3]),
    "Q4 (Highest)" = sprintf("Q4 (>%.1f)", cuts[3])
  ) |>
  tab_spanner("Choline quartile (mg/1,000 kcal)", columns = starts_with("Q")) |>
  tab_footnote("Survey-weighted geometric mean (95% CI).",
               locations = cells_column_spanners()) |>
  tab_source_note("ALT, alanine aminotransferase; AST, aspartate aminotransferase; BMI, body mass index; CI, confidence interval; GGT, gamma-glutamyl transferase.")

## Tables 5 and 6: adjusted ratios for each enzyme

enzyme_table <- function(models, var, label) {
  # n used in each model
  n <- map_int(models, \(m) nrow(model.frame(m)))
  
  models |>
    map(\(m) tbl_regression(m, include = all_of(var), exponentiate = TRUE,
                            label = setNames(list(label), var),
                            pvalue_fun = fmt_p) |>
          modify_header(estimate ~ "Ratio", p.value ~ "P value")) |>
    tbl_merge(tab_spanner = paste0(c("ALT", "AST", "GGT"), " (n = ", comma(n), ")")) |>
    modify_footnote(starts_with("ci") ~ NA)
}

tbl_5 <- enzyme_table(list(m_alt, m_ast, m_ggt),
                      "low_chol_obese", "Low choline and obesity (vs. all others)") |>
  modify_footnote(
    starts_with("estimate") ~ "Ratio of geometric means, adjusted for age, sex, race and ethnicity, and income-to-poverty ratio."
  ) |>
  as_gt() |>
  tab_source_note("Low choline = lowest quartile of choline density.") |>
  tab_source_note("ALT, alanine aminotransferase; AST, aspartate aminotransferase; CI, confidence interval; GGT, gamma-glutamyl transferase.")

# joint test across the six groups (same as 04)
p_joint <- c(
  ALT = as.numeric(regTermTest(m_alt_groups, ~bmi_chol_group)$p),
  AST = as.numeric(regTermTest(m_ast_groups, ~bmi_chol_group)$p),
  GGT = as.numeric(regTermTest(m_ggt_groups, ~bmi_chol_group)$p)
)

tbl_6 <- enzyme_table(list(m_alt_groups, m_ast_groups, m_ggt_groups),
                      "bmi_chol_group", "BMI category and choline intake") |>
  modify_footnote(
    starts_with("estimate") ~ "Ratio of geometric means, adjusted for age, sex, race and ethnicity, and income-to-poverty ratio.",
    starts_with("p.value") ~ paste0(
      "Joint test across the six groups: ",
      paste0(names(p_joint), " P = ", fmt_p(p_joint), collapse = "; "), "."
    )
  ) |>
  as_gt() |>
  tab_source_note("Low choline = lowest quartile of choline density; higher choline = quartiles 2 to 4.") |>
  tab_source_note("ALT, alanine aminotransferase; AST, aspartate aminotransferase; BMI, body mass index; CI, confidence interval; GGT, gamma-glutamyl transferase.")

## Table 7: low vs. higher choline within each BMI category
## Estimates come from 04 (within_bmi); this section only formats

tbl_7 <- within_bmi |>
  mutate(bmi_cat = unname(bmi_labels[bmi_cat]),
         ci = sprintf("%.2f, %.2f", ci_low, ci_high),
         ratio = sprintf("%.2f", ratio),
         p = fmt_p(p)) |>
  select(enzyme, 
         bmi_cat, 
         ratio, 
         ci, 
         p) |>
  pivot_wider(names_from = enzyme, 
              values_from = c(ratio, ci, p), 
              names_vary = "slowest") |>
  gt() |>
  cols_label(
    bmi_cat = "BMI category (kg/m²)",
    ratio_ALT = "Ratio", 
    ci_ALT = "95% CI", 
    p_ALT = "P value",
    ratio_AST = "Ratio", 
    ci_AST = "95% CI", 
    p_AST = "P value",
    ratio_GGT = "Ratio", 
    ci_GGT = "95% CI", 
    p_GGT = "P value"
  ) |>
  tab_spanner("ALT", 
              columns = ends_with("_ALT")) |>
  tab_spanner("AST", 
              columns = ends_with("_AST")) |>
  tab_spanner("GGT", 
              columns = ends_with("_GGT")) |>
  tab_footnote("Ratio of geometric means for low vs. higher choline within the same BMI category, adjusted for age, sex, race and ethnicity, and income-to-poverty ratio.",
               locations = cells_column_labels(columns = starts_with("ratio_"))) |>
  tab_source_note("Low choline = lowest quartile of choline density; 
                  higher choline = quartiles 2 to 4.") |>
  tab_source_note("ALT, alanine aminotransferase; 
                  AST, aspartate aminotransferase; 
                  BMI, body mass index; 
                  CI, confidence interval; 
                  GGT, gamma-glutamyl transferase.")

## Figure 1: % below AI by sex and BMI

fig_1 <- bind_rows(
  prev_sex_bmi |> mutate(bmi_cat = as.character(bmi_cat)),
  prev_bmi_all |> mutate(sex = "Total")
) |>
  filter(bmi_cat != "All") |>
  mutate(bmi_cat = factor(bmi_cat,
                          levels = names(bmi_labels),
                          labels = bmi_labels)) |>
  ggplot(aes(x = bmi_cat,
             y = prevalence,
             ymin = ci_lower,
             ymax = ci_upper,
             color = sex,
             shape = sex)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  labs(title = "Percentage Below Choline Adequate Intake (AI) by Sex and Body Mass Index (BMI)",
       subtitle = "U.S. Adults aged 20 years and older, NHANES August 2021-2023",
       x = "BMI category (kg/m²)",
       y = "Below the choline AI, % (95% CI)",
       caption = "Survey-weighted percentages with 95% Cls. 
       Al, Adequate intake; CI, confidence interval
       Al = 550 mg/d for men and 425 mg/d for women.",
       color = "Sex",
       shape = "Sex") +
  scale_color_few()

## Figure 2: % below AI by menopausal status and BMI

fig_2 <- prev_menopause_all |>
  filter(bmi_cat != "All") |>
  mutate(bmi_cat = factor(bmi_cat, levels = names(bmi_labels), labels = bmi_labels),
         menopause = factor(menopause, levels = c("pre", "post"),
                            labels = c("Premenopausal", "Postmenopausal"))) |>
  ggplot(aes(x = bmi_cat,
             y = prevalence,
             ymin = ci_lower,
             ymax = ci_upper,
             color = menopause,
             shape = menopause)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  labs(title = "Percentage below the choline Adequate Intake by Menopausal Status and BMI",
       subtitle = "U.S. Women aged 20 years and older, NHANES August 2021-August 2023",
       x = "BMI category (kg/m²)",
       y = "Below the choline AI, % (95% CI)",
       caption = "Survey-weighted percentages with 95% Cls. 
       Al = 425 mg/d. 
       Women who had undergone a hysterectomy and those of unknown status were excluded. 
       Al, Adequate intake; CI, confidence interval",
       color = "Menopausal Status", 
       shape = "Menopausal Status") +
  scale_color_few()

## Figure 3: liver enzymes across the six groups

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

# reference group
reference <- tibble(
  enzyme = c("ALT", "AST", "GGT"),
  group = "Normal, higher choline",
  ratio = 1, low = 1, high = 1
)

fig_3 <- bind_rows(
  get_ratios(m_alt_groups, "ALT"),
  get_ratios(m_ast_groups, "AST"),
  get_ratios(m_ggt_groups, "GGT"),
  reference
) |>
  mutate(
    bmi_cat = factor(str_extract(group, "^[A-Za-z]+"),
                     levels = rev(names(bmi_labels)[1:3]),
                     labels = rev(bmi_labels[1:3])),
    choline = factor(if_else(str_detect(group, "low"), "Low choline", "Higher choline"),
                     levels = c("Low choline", "Higher choline"))
  ) |>
  ggplot(aes(x = bmi_cat, 
             y = ratio, 
             ymin = low, 
             ymax = high,
             color = choline, 
             shape = choline)) +
  geom_hline(yintercept = 1, 
             linetype = 1, 
             color = "grey") +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  coord_flip() +
  facet_wrap(~enzyme, 
             nrow = 1) +
  scale_y_log10(breaks = c(0.8, 1, 1.25, 1.5, 2)) + # log scale for ratios
  labs(title = "Adjusted liver enzyme levels by BMI category and choline intake",
       subtitle = "US adults aged 20 years and older, NHANES August 2021-August 2023",
       x = "BMI category (kg/m²)",
       y = "Ratio vs. normal weight, higher choline (log scale)",
       color = "Choline intake",
       shape = "Choline intake",
       caption = "Ratios of geometric means with 95% CIs, adjusted for age, gender, race/ethnicity, and income-to-poverty ratio.
       Low Choline: lowest quartile of choline density (<121.4 mg/1,000 kcal); 
       High Choline: quartiles 2 through 4.
       ALT: alanine aminotransferase; 
       AST: aspartate aminotransferase; 
       BMI: body mass index; 
       CI: confidence interval; 
       GGT: gamma-glutamyl transferase") +
  scale_colour_few()

## Save tables

# tables as rds
saveRDS(tbl_1, "output/table1_baseline.rds")
saveRDS(tbl_2, "output/table2_sex_bmi.rds")
saveRDS(tbl_3, "output/table3_menopause.rds")
saveRDS(tbl_4, "output/table4_enzymes.rds")
saveRDS(tbl_5, "output/table5_contrast.rds")
saveRDS(tbl_6, "output/table6_groups.rds")
saveRDS(tbl_7, "output/table7_within.rds")

# tables as png
gtsave(tbl_1, "output/table1_baseline.png", vwidth = 1100, zoom = 2)
gtsave(tbl_2, "output/table2_sex_bmi.png", vwidth = 1100, zoom = 2)
gtsave(tbl_3, "output/table3_menopause.png", vwidth = 1100, zoom = 2)
gtsave(tbl_4, "output/table4_enzymes.png", vwidth = 1100, zoom = 2)
gtsave(tbl_5, "output/table5_contrast.png", vwidth = 1100, zoom = 2)
gtsave(tbl_6, "output/table6_groups.png", vwidth = 1100, zoom = 2)
gtsave(tbl_7, "output/table7_within.png", vwidth = 1100, zoom = 2)

# figures
ggsave("output/figure1_sex_bmi.png", fig_1, width = 8, height = 4.5, dpi = 300, bg = "white")
ggsave("output/figure2_menopause.png", fig_2, width = 8, height = 4.5, dpi = 300, bg = "white")
ggsave("output/figure3_groups.png", fig_3, width = 8, height = 4.5, dpi = 300, bg = "white")