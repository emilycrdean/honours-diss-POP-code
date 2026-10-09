# =====================================================================
# Evidence for results section 6.1 (UK Biobank Marfan).
# Add as a new chunk at the end of marfan_analysis.Rmd, after matching
# (needs ukb_clean, matched_data, ps_formula from earlier chunks).
# =====================================================================
library(dplyr)
library(MatchIt)
library(survival)

# ---- 1. Cases vs matched controls vs whole control pool --------------
summ <- function(d, label) data.frame(
  group            = label,
  n                = nrow(d),
  pop_n            = sum(d$POP == 1),
  pop_pct          = round(100 * mean(d$POP == 1), 1),
  birth_year       = round(mean(d$DOB), 1),
  live_births_mean = round(mean(d$X._live_birth), 2),
  pct_2plus_births = round(100 * mean(d$X._live_birth >= 2), 1),
  bmi              = round(mean(d$BMI), 1),
  townsend         = round(mean(d$townsend_deprivation), 2),
  pct_ever_smoked  = round(100 * mean(d$ever_smoked == "Yes"), 1),
  pct_menopause    = round(100 * mean(d$had_menopause == "Yes"), 1)
)
evidence_groups <- bind_rows(
  summ(filter(matched_data, cohort == "cases"),    "Marfan cases"),
  summ(filter(matched_data, cohort == "controls"), "Matched controls"),
  summ(filter(ukb_clean,    cohort == "controls"), "Control pool")
)
print(evidence_groups)

# ---- 2. POP in the control pool by number of live births and birth cohort ----
pool <- filter(ukb_clean, cohort == "controls")
pop_by_births <- pool %>%
  mutate(live_births = ifelse(X._live_birth >= 3, "3+", as.character(X._live_birth))) %>%
  group_by(live_births) %>%
  summarise(n = n(), pop_n = sum(POP == 1), pop_pct = round(100 * mean(POP == 1), 1))
print(pop_by_births)

pop_by_cohort <- pool %>%
  group_by(DOB_bin) %>%
  summarise(n = n(), pop_n = sum(POP == 1), pop_pct = round(100 * mean(POP == 1), 1))
print(pop_by_cohort)

# ---- 3. Repeated matching: 500 random draws of controls within the caliper ----
# Same approach as the All of Us Marfan analysis. Nearest-neighbour matching always picks the
# same closest controls, so re-running it in a different order barely changes anything. Instead,
# each case gets 4 controls drawn at random from every control that would have been an
# acceptable match (same 5-year birth cohort, propensity score within 0.2 SD).
set.seed(2026)
ps_vals <- glm(ps_formula, data = ukb_clean, family = binomial)$fitted.values
cal     <- 0.2 * sd(ps_vals)
case_i  <- which(ukb_clean$marfan == 1)
ctrl_i  <- which(ukb_clean$marfan == 0)
eligible <- lapply(case_i, function(i)
  ctrl_i[ukb_clean$DOB_bin[ctrl_i] == ukb_clean$DOB_bin[i] & abs(ps_vals[ctrl_i] - ps_vals[i]) <= cal])
cat("Eligible controls per case: min", min(lengths(eligible)), "median", median(lengths(eligible)), "\n")

rep_df <- bind_rows(lapply(1:500, function(r) {
  used <- integer(0); rows <- list()
  for (k in sample(seq_along(case_i))) {                 # random case order, no control reused
    avail <- setdiff(eligible[[k]], used)
    pick  <- if (length(avail) > 4) sample(avail, 4) else avail
    used  <- c(used, pick)
    rows[[k]] <- data.frame(POP = ukb_clean$POP[c(case_i[k], pick)],
                            case = c(1, rep(0, length(pick))), set = k)
  }
  d   <- bind_rows(rows)
  fit <- clogit(POP ~ case + strata(set), data = d, method = "exact")
  data.frame(or = unname(exp(coef(fit))),
             p  = unname(summary(fit)$coefficients[, "Pr(>|z|)"]),
             ctrl_pop = sum(d$POP[d$case == 0]))
}))
repeat_summary <- c(
  median_or = round(median(rep_df$or), 2),
  or_2.5    = round(unname(quantile(rep_df$or, 0.025)), 2),
  or_97.5   = round(unname(quantile(rep_df$or, 0.975)), 2),
  pct_or_above_1 = round(100 * mean(rep_df$or > 1)),
  pct_p_below_05 = round(100 * mean(rep_df$p < 0.05))
)
print(repeat_summary)

# ---- 4. Save everything for the thesis ------------------------------------
write.csv(evidence_groups, "ukb_evidence_groups.csv",   row.names = FALSE)
write.csv(pop_by_births,   "ukb_pop_by_births.csv",     row.names = FALSE)
write.csv(pop_by_cohort,   "ukb_pop_by_birth_cohort.csv", row.names = FALSE)
write.csv(rep_df,          "ukb_repeated_matching.csv", row.names = FALSE)
