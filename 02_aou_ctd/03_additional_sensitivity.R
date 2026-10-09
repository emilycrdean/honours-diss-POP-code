# =====================================================================
# All of Us (Aim 1): additional sensitivity analyses (Supplementary Table Y)
#
# Run in the EDS notebook after 01_build_analysis_dataset.R (needs eds_combined,
# match_data and surveyOccurrence_df). Covers:
#   - full unmatched cohort: Fisher's exact test, and logistic regression adjusted for
#     age, BMI and deprivation, then + menopause code, then + smoking
#   - matched sample: Fisher's exact test; conditional logistic regression + race/ethnicity
#   - matching on age only; exact matching on race/ethnicity; + matching on deprivation
# =====================================================================
library(tidyverse); library(survival); library(MatchIt); library(broom)

# ---- smoking status (cases from this notebook, controls from the saved file) ----
system2("gsutil", c("cp", "gs://rerun-controls-YOUR-WORKSPACE-ID/control_smoking.rds", "/tmp/"),
        stdout = TRUE, stderr = TRUE)
control_smoking <- readRDS("/tmp/control_smoking.rds")
stopifnot(nrow(control_smoking) > 200000)          # guards against an empty file

case_smoking <- surveyOccurrence_df %>%
  filter(question == 1585857 | question == "1585857") %>%
  transmute(person_id = as.character(person_id), ans = as.character(answer)) %>%
  group_by(person_id) %>%
  summarise(smoking_status = case_when(any(ans == "1585858") ~ "Ever",
                                       any(ans == "1585859") ~ "Never",
                                       TRUE ~ NA_character_), .groups = "drop")

eds_combined <- eds_combined %>%
  select(-any_of("smoking_status")) %>%
  left_join(bind_rows(case_smoking, control_smoking), by = "person_id")

# ---- helpers ----
row_tidy <- function(label, fit) {
  t <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "is_case")
  tibble(analysis = label, OR = t$estimate, lower = t$conf.low, upper = t$conf.high, p = t$p.value)
}
row_fisher <- function(label, tab) {
  f <- fisher.test(tab)
  tibble(analysis = label, OR = unname(f$estimate), lower = f$conf.int[1], upper = f$conf.int[2], p = f$p.value)
}

# ---- matched samples (1:4, caliper 0.2) ----
m  <- matchit(is_case ~ age + bmi_value, data = match_data, method = "nearest", ratio = 4, caliper = 0.2)
md <- match.data(m)                                              # primary: age + BMI
m2  <- matchit(is_case ~ age, data = eds_combined, method = "nearest", ratio = 4, caliper = 0.2)
md2 <- match.data(m2)                                            # age only, missing BMI retained
m3  <- matchit(is_case ~ age + bmi_value, data = match_data, method = "nearest", ratio = 4, caliper = 0.2,
               exact = ~ race_eth)
md3 <- match.data(m3)                                            # exact on race/ethnicity
md_dep_data <- match_data %>% filter(!is.na(deprivation_index))
m4  <- matchit(is_case ~ age + bmi_value + deprivation_index, data = md_dep_data, method = "nearest",
               ratio = 4, caliper = 0.2, exact = ~ race_eth)
md4 <- match.data(m4)                                            # + deprivation, exact on race/ethnicity

# ---- full-cohort logistic regression ----
glm_a <- glm(prolapse_case ~ is_case + age + bmi_value + deprivation_index, data = eds_combined, family = binomial)
glm_b <- glm(prolapse_case ~ is_case + age + bmi_value + deprivation_index + menopause_code, data = eds_combined, family = binomial)
glm_c <- glm(prolapse_case ~ is_case + age + bmi_value + deprivation_index + smoking_status, data = eds_combined, family = binomial)

# ---- results table ----
eds_results <- bind_rows(
  row_fisher("Full cohort, unadjusted (Fisher's exact)", table(eds_combined$is_case, eds_combined$prolapse_case)),
  row_tidy("Full cohort, adjusted: age, BMI, deprivation", glm_a),
  row_tidy("Full cohort, adjusted: + menopause code", glm_b),
  row_tidy("Full cohort, adjusted: + smoking", glm_c),
  row_fisher("Matched (age + BMI), unadjusted (Fisher's exact)", table(md$is_case, md$prolapse_case)),
  row_tidy("Matched (age + BMI), conditional logistic [primary]",
           clogit(prolapse_case ~ is_case + strata(subclass), data = md, method = "exact")),
  row_tidy("Matched (age + BMI), conditional logistic + race/ethnicity",
           clogit(prolapse_case ~ is_case + race_eth + strata(subclass), data = md, method = "exact")),
  row_tidy("Matched (age only), conditional logistic",
           clogit(prolapse_case ~ is_case + strata(subclass), data = md2, method = "exact")),
  row_tidy("Matched (age + BMI, exact on race/ethnicity), conditional logistic",
           clogit(prolapse_case ~ is_case + strata(subclass), data = md3, method = "exact")),
  row_tidy("Matched (age + BMI + deprivation, exact on race/ethnicity), conditional logistic",
           clogit(prolapse_case ~ is_case + strata(subclass), data = md4, method = "exact"))
) %>% mutate(across(c(OR, lower, upper), ~ round(.x, 2)),
             p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# ---- sample sizes and balance (inside the Workbench only) ----
rbind(primary = summary(m, un = FALSE)$nn["Matched", ], age_only = summary(m2, un = FALSE)$nn["Matched", ],
      exact_eth = summary(m3, un = FALSE)$nn["Matched", ], dep_eth = summary(m4, un = FALSE)$nn["Matched", ])
summary(m, un = FALSE)$sum.matched[, c("Means Treated", "Means Control", "Std. Mean Diff.")]
c(glm_a = nobs(glm_a), glm_b = nobs(glm_b), glm_c = nobs(glm_c))

eds_results                                         # no counts: safe to export
write_csv(eds_results, "/tmp/eds_results_table_v4.csv")
system2("gsutil", c("cp", "/tmp/eds_results_table_v4.csv", "gs://rerun-eds-YOUR-WORKSPACE-ID/"),
        stdout = TRUE, stderr = TRUE)
