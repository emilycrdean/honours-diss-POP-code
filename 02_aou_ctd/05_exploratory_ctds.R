# =====================================================================
# Extra CTD analyses (osteogenesis imperfecta, cutis laxa, hypermobility syndrome)
# One notebook per condition. Run the Dataset Builder loading cells first
# (person, conditionOccurrence, bmi, zipcodeSocioeconomic), then these two cells.
# =====================================================================


# ---- CELL A: build the case data frame ----------------------------------
# >>> change this one line per notebook: "oi_case", "eds_case" or "hypermobility_case"
case_label <- "condition_case"

library(tidyverse)

cond_flags <- conditionOccurrence_df %>%
  mutate(person_id = as.character(person_id),
         code  = T_DISP_source_concept_code,
         vocab = T_DISP_source_vocabulary,
         name  = coalesce(T_DISP_source_concept_name, T_DISP_standard_concept_name),
         is_pop = case_when(
           vocab == "ICD10CM" ~ grepl("^N81", code) | code == "N99.3",
           vocab == "ICD9CM"  ~ grepl("^618", code),
           TRUE ~ grepl("prolapse|cystocele|rectocele|enterocele|urethrocele|herniation of rectum into vagina",
                        name, ignore.case = TRUE) &
                  !grepl("rectal prolapse|anal prolapse|pregnan|gravid|maternal|ovary|fallopian",
                         name, ignore.case = TRUE)),
         is_meno = !is_pop & (grepl("^N95|^627", code) | grepl("menopaus", name, ignore.case = TRUE)))

cases_df <- person_df %>%
  transmute(person_id = as.character(person_id),
            age = as.numeric(Sys.Date() - as.Date(birth_datetime)) / 365.25,
            race = T_DISP_race, ethnicity = T_DISP_ethnicity) %>%
  distinct(person_id, .keep_all = TRUE) %>%
  left_join(cond_flags %>% filter(is_pop) %>%
              group_by(person_id) %>%
              summarise(first_prolapse_date = suppressWarnings(min(as.Date(condition_start_datetime), na.rm = TRUE)),
                        pop_any = 1L),
            by = "person_id") %>%
  left_join(cond_flags %>% filter(is_meno) %>% distinct(person_id) %>% mutate(menopause_code = 1L),
            by = "person_id") %>%
  left_join(bmi_df %>%
              group_by(person_id = as.character(person_id)) %>%
              summarise(bmi_value = median(value_numeric, na.rm = TRUE), n_bmi = n()),
            by = "person_id") %>%
  left_join(zipcodeSocioeconomic_df %>%
              distinct(person_id = as.character(person_id), .keep_all = TRUE) %>%
              select(person_id, deprivation_index),
            by = "person_id") %>%
  mutate(cohort = case_label,
         prolapse_case  = coalesce(pop_any, 0L),
         menopause_code = coalesce(menopause_code, 0L)) %>%
  select(-pop_any)

# checks (inside the Workbench only)
cond_flags %>% filter(is_pop) %>% distinct(person_id, code, name) %>% count(code, name, sort = TRUE) %>% print(n = Inf)
cases_df %>% summarise(n = n(), pop = sum(prolapse_case), missing_bmi = sum(is.na(bmi_value)))
# NO saveRDS / gsutil lines here, so the controls file can't be overwritten


# ---- CELL B: full analysis ---------------------------------------------------
library(survival)
library(MatchIt)
library(broom)

system("gsutil cp gs://YOUR-WORKSPACE-ID/control_analysis_df_v3.rds /tmp/")
control_df <- readRDS("/tmp/control_analysis_df_v3.rds")
stopifnot(nrow(control_df) > 250000, all(control_df$cohort == "control"))   # guards against a wrong file

run_ctd <- function(cases, controls, label) {
  combined <- bind_rows(
      controls %>% filter(!person_id %in% cases$person_id) %>% mutate(cohort = "control"),
      cases    %>% mutate(cohort = label)) %>%
    mutate(is_case = as.integer(cohort != "control"))

  safe <- function(expr) tryCatch(expr, error = function(e) NULL)
  fmt  <- function(lbl, est, lo, hi, p) tibble(analysis = lbl, OR = est, lower = lo, upper = hi, p = p)
  from_fisher <- function(lbl, tab) {
    f <- safe(fisher.test(tab))
    if (is.null(f)) return(fmt(lbl, NA, NA, NA, NA))
    fmt(lbl, unname(f$estimate), f$conf.int[1], f$conf.int[2], f$p.value)
  }
  from_fit <- function(lbl, fit) {
    if (is.null(fit)) return(fmt(lbl, NA, NA, NA, NA))
    t <- safe(tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "is_case"))
    if (is.null(t) || nrow(t) == 0) return(fmt(lbl, NA, NA, NA, NA))
    fmt(lbl, t$estimate, t$conf.low, t$conf.high, t$p.value)
  }

  # primary: 1:4 on age + BMI
  md_data <- combined %>% filter(!is.na(age), !is.na(bmi_value))
  m  <- matchit(is_case ~ age + bmi_value, data = md_data, method = "nearest", ratio = 4, caliper = 0.2)
  md <- match.data(m)
  clr <- safe(clogit(prolapse_case ~ is_case + strata(subclass), data = md, method = "exact"))

  # sensitivity: age only (keeps those missing BMI)
  m2  <- matchit(is_case ~ age, data = combined, method = "nearest", ratio = 4, caliper = 0.2)
  md2 <- match.data(m2)
  clr2 <- safe(clogit(prolapse_case ~ is_case + strata(subclass), data = md2, method = "exact"))

  # full cohort, adjusted
  glm_full <- safe(glm(prolapse_case ~ is_case + age + bmi_value + deprivation_index,
                       data = combined, family = binomial))

  results <- bind_rows(
    from_fisher("Full cohort, unadjusted (Fisher's exact)", table(combined$is_case, combined$prolapse_case)),
    from_fit("Full cohort, adjusted: age, BMI, deprivation", glm_full),
    from_fisher("Matched (age + BMI), unadjusted (Fisher's exact)", table(md$is_case, md$prolapse_case)),
    from_fit("Matched (age + BMI), conditional logistic [primary]", clr),
    from_fit("Matched (age only), conditional logistic", clr2)
  ) %>% mutate(exposure = label,
               across(c(OR, lower, upper), ~ round(.x, 2)),
               p = ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))) %>%
    relocate(exposure)

  # counts for YOU only (inside the Workbench) - suppress anything 1-20 before exporting
  counts <- list(
    cases_total      = sum(combined$is_case),
    cases_matched    = sum(md$is_case),
    controls_matched = sum(md$is_case == 0),
    pop_cases_full   = sum(combined$prolapse_case[combined$is_case == 1]),
    pop_cases_matched = sum(md$prolapse_case[md$is_case == 1]),
    pop_controls_matched = sum(md$prolapse_case[md$is_case == 0])
  )
  balance <- summary(m, un = FALSE)$sum.matched[, c("Means Treated", "Means Control", "Std. Mean Diff.")]

  list(results = results, counts = counts, balance = balance)
}

out <- run_ctd(cases_df, control_df, case_label)
out$counts      # Workbench only
out$balance
out$results     # no counts in it: safe to export if the case total is >20
write_csv(out$results, paste0(case_label, "_results_table.csv"))
