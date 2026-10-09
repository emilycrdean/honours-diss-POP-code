# =====================================================================
# All of Us (Aim 1): Marfan syndrome resampling sensitivity analysis
#
# The Marfan estimate rests on few POP events, so this checks how much it depends
# on which controls were selected. For each case, 4 controls are drawn at random from
# every control within the caliper; logistic model is then refitted, and
# this is repeated 500 times. Run in the Marfan notebook after
# 01_build_analysis_dataset.R (needs match_data built from marfan_combined).
# Output is a summary of the OR distribution only (no counts) due to AoU data privacy
# =====================================================================
library(tidyverse); library(survival)

set.seed(2026)

# propensity score and caliper, same as the primary match
ps_fit <- glm(is_case ~ age + bmi_value, data = match_data, family = binomial)
match_data$ps <- predict(ps_fit, type = "response")
cal <- 0.2 * sd(match_data$ps)

cases    <- match_data %>% filter(is_case == 1)
controls <- match_data %>% filter(is_case == 0)

# every control within the caliper of each case
pools <- lapply(cases$ps, function(p) which(abs(controls$ps - p) <= cal))
pool_size <- lengths(pools)
summary(pool_size)                          # eligible controls per case
stopifnot(all(pool_size >= 4))

one_draw <- function() {
  sets <- lapply(seq_len(nrow(cases)), function(i) {
    picked <- controls[sample(pools[[i]], 4), ]
    bind_rows(cases[i, ], picked) %>% mutate(set_id = i)
  })
  d <- bind_rows(sets)
  fit <- tryCatch(clogit(prolapse_case ~ is_case + strata(set_id), data = d, method = "exact"),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(OR = NA, p = NA))
  s <- summary(fit)$coefficients
  c(OR = exp(s["is_case", "coef"]), p = s["is_case", "Pr(>|z|)"])
}

draws <- t(replicate(500, one_draw()))
draws <- as.data.frame(draws) %>% filter(!is.na(OR), is.finite(OR))

c(draws_used     = nrow(draws),
  median_OR      = median(draws$OR),
  OR_2.5pct      = unname(quantile(draws$OR, 0.025)),
  OR_97.5pct     = unname(quantile(draws$OR, 0.975)),
  pct_draws_p05  = mean(draws$p < 0.05) * 100)

hist(draws$OR, breaks = 40, main = "Marfan: OR across 500 random control sets",
     xlab = "Odds ratio", col = "grey80", border = "white")
abline(v = c(1, median(draws$OR)), lty = c(2, 1), col = c("black", "firebrick"))
