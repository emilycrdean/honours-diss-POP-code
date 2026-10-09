# =====================================================================
# All of Us (Aim 1): build the analysis dataset for one cohort
#
# Run in each cohort's notebook (controls, EDS, Marfan) BELOW the cells the
# Dataset Builder generates. Those cells load person_df, conditionOccurrence_df,
# bmi_df, zipcodeSocioeconomic_df and surveyOccurrence_df; they are generated per
# workspace and are not included here (concept sets: Supplementary Table S2).
#
# Set this_cohort below, run the notebook, and the dataset is saved to that
# cohort's workspace bucket. Intermediate files go to /tmp, never the notebook folder.
# =====================================================================
library(tidyverse)

this_cohort <- "eds_case"      # "control" in the controls notebook; "marfan_case" for Marfan
bucket      <- "gs://rerun-eds-YOUR-WORKSPACE-ID"    # this cohort's bucket (controls: gs://rerun-controls-YOUR-WORKSPACE-ID)

# ---- 1. classify each condition record as POP, menopause or neither ----
# POP: ICD-10-CM N81.x or N99.3, ICD-9-CM 618.x, or SNOMED genital prolapse concepts;
# rectal/anal prolapse, prolapse in pregnancy and ovary/fallopian tube prolapse removed.
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

# ---- 2. one row per participant ----
# Age at data extraction; BMI = median of all recorded values (limits erroneous entries).
analysis_df <- person_df %>%
  transmute(person_id = as.character(person_id),
            age = as.numeric(Sys.Date() - as.Date(birth_datetime)) / 365.25,
            race = T_DISP_race, ethnicity = T_DISP_ethnicity) %>%
  left_join(cond_flags %>% filter(is_pop) %>%
              group_by(person_id) %>%
              summarise(first_prolapse_date = min(as.Date(condition_start_datetime), na.rm = TRUE),
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
  mutate(cohort = this_cohort,
         prolapse_case  = coalesce(pop_any, 0L),
         menopause_code = coalesce(menopause_code, 0L)) %>%
  select(-pop_any)

# ---- 3. checks (inside the Workbench only: these print small counts) ----
cond_flags %>% filter(is_pop) %>% distinct(person_id, code, name) %>% count(code, name, sort = TRUE) %>% print(n = Inf)
analysis_df %>% summarise(n = n(), pop = sum(prolapse_case), pct_pop = mean(prolapse_case) * 100,
                          meno = sum(menopause_code))

# ---- 4. save to the cohort's bucket ----
out <- file.path("/tmp", paste0(sub("_case$", "", this_cohort), "_analysis_df_v3.rds"))
saveRDS(analysis_df, out)
system2("gsutil", c("cp", out, paste0(bucket, "/")), stdout = TRUE, stderr = TRUE)

# Controls notebook only: smoking status (Lifestyle survey, "smoked at least 100 cigarettes",
# question 1585857), derived the same way as for cases in 02_run_all.R
if (this_cohort == "control") {
  control_smoking <- surveyOccurrence_df %>%
    filter(question == 1585857 | question == "1585857") %>%
    transmute(person_id = as.character(person_id), ans = as.character(answer)) %>%
    group_by(person_id) %>%
    summarise(smoking_status = case_when(any(ans == "1585858") ~ "Ever",
                                         any(ans == "1585859") ~ "Never",
                                         TRUE ~ NA_character_), .groups = "drop")
  saveRDS(control_smoking, "/tmp/control_smoking.rds")
  system2("gsutil", c("cp", "/tmp/control_smoking.rds", paste0(bucket, "/")), stdout = TRUE, stderr = TRUE)
}

# ---- 5. case notebooks: combine with the control pool ----
# (in the Marfan notebook, use marfan_df / "marfan_case" and call the result marfan_combined)
if (this_cohort != "control") {
  system2("gsutil", c("cp", "gs://rerun-controls-YOUR-WORKSPACE-ID/control_analysis_df_v3.rds", "/tmp/"),
          stdout = TRUE, stderr = TRUE)
  control_df <- readRDS("/tmp/control_analysis_df_v3.rds")
  sum(control_df$person_id %in% analysis_df$person_id)     # cases sitting in the control pool

  eds_combined <- bind_rows(
      control_df %>% filter(!person_id %in% analysis_df$person_id) %>% mutate(cohort = "control"),
      analysis_df
    ) %>%
    mutate(is_case = if_else(cohort == "control", 0, 1),
           race_eth = case_when(
             grepl("^hispanic", ethnicity, ignore.case = TRUE)                 ~ "Hispanic/Latino",
             race == "White"                                                   ~ "Non-Hispanic White",
             grepl("black", race, ignore.case = TRUE)                          ~ "Non-Hispanic Black",
             race == "Asian"                                                   ~ "Non-Hispanic Asian",
             grepl("more than one|another single|american indian|native hawaiian|middle eastern",
                   race, ignore.case = TRUE)                                   ~ "Other/multiple",
             TRUE                                                              ~ "Unknown/not reported"))

  match_data <- eds_combined %>% filter(!is.na(age), !is.na(bmi_value))   # primary match excludes missing BMI
}
