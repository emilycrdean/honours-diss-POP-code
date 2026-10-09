# =====================================================================
# All of Us: descriptive tables + covariate analysis
# Paste into each condition notebook BELOW the cell that defines run_all
# (it reuses the same data objects). Counts of 1-19 are replaced with "<20"
# automatically, so the printed output is safe to copy out.
# =====================================================================
library(tidyverse); library(MatchIt)

describe_all <- function(combined, exposure, covariate_model = FALSE) {

  sup  <- function(x) ifelse(x > 0 & x < 20, "<20", format(x, big.mark = ",", trim = TRUE))
  pct  <- function(k, n) ifelse(k > 0 & k < 20, "-", sprintf("%.1f", 100 * k / n))

  # ---- same exclusions, smoking and race/ethnicity coding as run_all ----
  # there are lots of "your-workspace-id"s here -- they were used to call up different ids for different cohorts
  # i.e., one was the control cohort, OI cohort etc - they have the data project (CDR) in AoU (something silly usually)
  # then the specification - so like gs://controls-wb-silky-bagel-3498 or something equally strange
  get_ids <- function(bucket) {
    files <- system2("gsutil", c("ls", paste0(bucket, "/**")), stdout = TRUE, stderr = TRUE)
    f <- grep("person[^/]*\\.csv(\\.gz)?$", files, value = TRUE)
    unique(as.character(unlist(lapply(f, function(x) {
      tmp <- file.path("/tmp", basename(x))
      system2("gsutil", c("cp", x, tmp), stdout = TRUE, stderr = TRUE)
      read_csv(tmp, col_select = "person_id", show_col_types = FALSE)$person_id
    }))))
  }
  exclude <- union(get_ids("gs://YOUR-WORKSPACE-ID"),
                   get_ids("gs://YOUR-WORKSPACE-ID"))
  system2("gsutil", c("cp", "gs://YOUR-WORKSPACE-ID/control_smoking.rds", "/tmp/"),
          stdout = TRUE, stderr = TRUE)
  case_smoking <- surveyOccurrence_df %>%
    filter(question == 1585857 | question == "1585857") %>%
    transmute(person_id = as.character(person_id), ans = as.character(answer)) %>%
    group_by(person_id) %>%
    summarise(smoking_status = case_when(any(ans == "1585858") ~ "Ever",
                                         any(ans == "1585859") ~ "Never",
                                         TRUE ~ NA_character_), .groups = "drop")
  smoking <- bind_rows(case_smoking, readRDS("/tmp/control_smoking.rds")) %>%
    distinct(person_id, .keep_all = TRUE)

  d <- combined %>%
    select(-any_of(c("smoking_status", "race_eth", "is_case"))) %>%
    filter(!(cohort == "control" & person_id %in% exclude)) %>%
    mutate(is_case = as.integer(cohort != "control"),
           race_eth = case_when(
             grepl("^hispanic", ethnicity, ignore.case = TRUE) ~ "Hispanic/Latino",
             race == "White" ~ "Non-Hispanic White",
             grepl("black", race, ignore.case = TRUE) ~ "Non-Hispanic Black",
             race == "Asian" ~ "Non-Hispanic Asian",
             grepl("more than one|another single|american indian|native hawaiian|middle eastern", race, ignore.case = TRUE) ~ "Other/multiple",
             TRUE ~ "Unknown/not reported"),
           race_eth = relevel(factor(race_eth), ref = "Non-Hispanic White")) %>%
    left_join(smoking, by = "person_id") %>%
    mutate(smoking_status = factor(coalesce(smoking_status, "Unknown"), levels = c("Never", "Ever", "Unknown")))
  stopifnot(!any(duplicated(d$person_id)))

  # ---- A. POP in the full, unmatched cohort ----
  cat("\n=== A. POP in the full cohort:", exposure, "===\n")
  d %>% mutate(group = ifelse(is_case == 1, exposure, "Controls")) %>%
    group_by(group) %>% summarise(n_raw = n(), pop_raw = sum(prolapse_case), .groups = "drop") %>%
    transmute(group, n = sup(n_raw), POP = sup(pop_raw), `POP %` = pct(pop_raw, n_raw)) %>%
    print()

  # ---- B. Characteristics of the primary matched sample (Table 1) ----
  md_data <- d %>% filter(!is.na(age), !is.na(bmi_value))
  m  <- matchit(is_case ~ age + bmi_value, data = md_data, method = "nearest", ratio = 4, caliper = 0.2)
  md <- match.data(m)
  n1 <- sum(md$is_case == 1); n0 <- sum(md$is_case == 0)

  cont_row <- function(v, label) {
    x1 <- md[[v]][md$is_case == 1]; x0 <- md[[v]][md$is_case == 0]
    tibble(variable = label,
           cases    = sprintf("%.1f (%.1f)", mean(x1, na.rm = TRUE), sd(x1, na.rm = TRUE)),
           controls = sprintf("%.1f (%.1f)", mean(x0, na.rm = TRUE), sd(x0, na.rm = TRUE)),
           SMD      = sprintf("%.3f", (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) /
                                       sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)),
           missing  = paste(sup(sum(is.na(x1))), "/", sup(sum(is.na(x0)))))
  }
  cat_rows <- function(v, label) {
    lv <- levels(factor(md[[v]]))
    bind_rows(tibble(variable = label, cases = "", controls = "", SMD = "", missing = ""),
      map_dfr(lv, function(l) {
        k1 <- sum(md[[v]][md$is_case == 1] == l, na.rm = TRUE)
        k0 <- sum(md[[v]][md$is_case == 0] == l, na.rm = TRUE)
        p1 <- k1 / n1; p0 <- k0 / n0
        tibble(variable = paste0("  ", l),
               cases    = paste0(sup(k1), " (", pct(k1, n1), ")"),
               controls = paste0(sup(k0), " (", pct(k0, n0), ")"),
               SMD      = sprintf("%.3f", (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)),
               missing  = "")
      }))
  }
  tab1 <- bind_rows(
    tibble(variable = "n", cases = sup(n1), controls = sup(n0), SMD = "", missing = ""),
    cont_row("age", "Age, mean (SD)"),
    cont_row("bmi_value", "BMI, mean (SD)"),
    cont_row("deprivation_index", "Deprivation index, mean (SD)"),
    cat_rows("race_eth", "Race and ethnicity, n (%)"),
    cat_rows("smoking_status", "Smoking, n (%)"),
    if ("menopause_code" %in% names(md)) cat_rows("menopause_code", "Menopause code, n (%)")
  )
  cat("\n=== B. Matched sample characteristics:", exposure, "===\n")
  print(as.data.frame(tab1), right = FALSE, row.names = FALSE)

  # ---- C. Covariate associations with POP in the control population ----
  if (covariate_model) {
    ctrl <- d %>% filter(is_case == 0, !is.na(age), !is.na(bmi_value), !is.na(deprivation_index)) %>%
      mutate(age10  = age / 10,
             bmi5   = bmi_value / 5,
             dep_sd = (deprivation_index - mean(deprivation_index)) / sd(deprivation_index))
    fit_tab <- function(fit) {
      co <- summary(fit)$coefficients
      tibble(term = rownames(co), OR = exp(co[, 1]),
             lower = exp(co[, 1] - 1.96 * co[, 2]), upper = exp(co[, 1] + 1.96 * co[, 2]),
             p = co[, 4]) %>%
        filter(term != "(Intercept)") %>%
        mutate(across(c(OR, lower, upper), ~ round(.x, 2)),
               p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
    }
    f1 <- glm(prolapse_case ~ age10 + bmi5 + dep_sd + race_eth + smoking_status,
              data = ctrl, family = binomial)
    cat("\n=== C1. POP risk factors in controls (n =", sup(nrow(ctrl)), "; POP =", sup(sum(ctrl$prolapse_case)), ") ===\n")
    cat("Age per 10 years, BMI per 5 units, deprivation per SD\n")
    print(fit_tab(f1), n = Inf)
    cat("POP events by race/ethnicity level:\n")
    print(ctrl %>% group_by(race_eth) %>% summarise(n = sup(n()), POP = sup(sum(prolapse_case))))
    if ("menopause_code" %in% names(ctrl)) {
      f2 <- update(f1, . ~ . + menopause_code)
      cat("\n=== C2. As C1, plus menopause code ===\n")
      print(fit_tab(f2), n = Inf)
    }
  }
  invisible(NULL)
}

# ---- run: covariate_model = TRUE only in the EDS notebook (same control pool everywhere) ----
describe_all(eds_combined, "EDS", covariate_model = TRUE)
# describe_all(hyp_combined, "Hypermobility")
# describe_all(marfan_combined, "Marfan")
# describe_all(bind_rows(control_df %>% filter(!person_id %in% cases_df$person_id) %>% mutate(cohort = "control"), cases_df), "OI")
