library(tidyverse); library(survival); library(MatchIt); library(broom)

run_all <- function(combined, exposure) {

  # ---- 1. remove hypermobility + OI or whatever extra people from the controls ----
  get_ids <- function(bucket) {
    files <- system2("gsutil", c("ls", paste0(bucket, "/**")), stdout = TRUE, stderr = TRUE)
    f <- grep("person[^/]*\\.csv(\\.gz)?$", files, value = TRUE)
    unique(as.character(unlist(lapply(f, function(x) {
      tmp <- file.path("/tmp", basename(x))
      system2("gsutil", c("cp", x, tmp), stdout = TRUE, stderr = TRUE)
      read_csv(tmp, col_select = "person_id", show_col_types = FALSE)$person_id
    }))))
  }
  exclude <- union(get_ids("gs://hypermobility-syndrome-YOUR-WORKSPACE-ID"),
                   get_ids("gs://oi-YOUR-WORKSPACE-ID"))

  # ---- 2. smoking (cases from this notebook, controls from saved file) ----
  system2("gsutil", c("cp", "gs://YOUR-WORKSPACE-ID/controls_smoking.rds", "/tmp/"),
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

  # ---- 3. analysis dataset ----
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
    mutate(smoking_status = factor(coalesce(smoking_status, "Unknown"),     # missing kept as its own level
                                   levels = c("Never", "Ever", "Unknown")))
  stopifnot(!any(duplicated(d$person_id)))

  # ---- 4. primary match (age + BMI) and age-only match ----
  md_data <- d %>% filter(!is.na(age), !is.na(bmi_value))
  m   <- matchit(is_case ~ age + bmi_value, data = md_data, method = "nearest", ratio = 4, caliper = 0.2)
  md  <- match.data(m)
  m2  <- matchit(is_case ~ age, data = d, method = "nearest", ratio = 4, caliper = 0.2)
  md2 <- match.data(m2)

  ct <- function(fit) tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "is_case")
  p1 <- ct(clogit(prolapse_case ~ is_case + strata(subclass), data = md,  method = "exact"))
  p2 <- ct(clogit(prolapse_case ~ is_case + strata(subclass), data = md2, method = "exact"))

  # ---- 5. matched + adjusted: SAME matched sets, extra covariates ----
  md_adj <- md %>% filter(!is.na(deprivation_index))
  # keep only matched sets that still contain the case and >= 1 control
  md_adj <- md_adj %>% group_by(subclass) %>%
    filter(sum(is_case) == 1, sum(is_case == 0) >= 1) %>% ungroup() %>% droplevels()
  n_params <- 1 + (nlevels(md_adj$race_eth) - 1) + (nlevels(md_adj$smoking_status) - 1)
  n_events <- sum(md_adj$prolapse_case)
  enough   <- n_events >= 10 * n_params              # ~10 events per parameter

  if (enough) {
    p3 <- ct(clogit(prolapse_case ~ is_case + deprivation_index + race_eth + smoking_status + strata(subclass),
                    data = md_adj, method = "exact"))
  } else {
    p3 <- tibble(estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_)
  }

  # ---- 6. crude ----
  f <- fisher.test(table(d$is_case, d$prolapse_case))

  out <- tibble(
    exposure = exposure,
    analysis = c("Matched, age + BMI (primary)", "Matched, age only",
                 "Matched + adjusted (deprivation, race/ethnicity, smoking)", "Crude (unmatched)"),
    n_used   = c(paste(sum(md$is_case), "/", sum(md$is_case == 0)),
                 paste(sum(md2$is_case), "/", sum(md2$is_case == 0)),
                 paste(sum(md_adj$is_case), "/", sum(md_adj$is_case == 0)),
                 as.character(nrow(d))),
    OR    = round(c(p1$estimate, p2$estimate, p3$estimate, unname(f$estimate)), 2),
    lower = round(c(p1$conf.low, p2$conf.low, p3$conf.low, f$conf.int[1]), 2),
    upper = round(c(p1$conf.high, p2$conf.high, p3$conf.high, f$conf.int[2]), 2),
    p     = c(p1$p.value, p2$p.value, p3$p.value, f$p.value)
  ) %>% mutate(p = case_when(is.na(p) ~ "not estimable (too few events)",
                             p < 0.001 ~ "<0.001",
                             TRUE ~ sprintf("%.3f", p)))

  print(table(d$cohort))
  print(summary(m, un = FALSE)$sum.matched[, "Std. Mean Diff."])
  cat("Adjusted model: events =", n_events, "| parameters =", n_params,
      "| needed ~", 10 * n_params, "| run =", enough, "\n")
  out
}

# ---- change this last line per notebook ----
run_all(cases_combined, "cases")
# run_all(marfan_combined, "Marfan")
# run_all(hyp_combined, "Hypermobility")
# run_all(bind_rows(control_df %>% filter(!person_id %in% cases_df$person_id) %>% mutate(cohort = "control"), cases_df), "OI")
