# =====================================================================
# Aim 1 figures (UK Biobank + All of Us). Run thesis_theme.R first.
# All numbers are final results; nothing here needs Workbench data.
# UK Biobank Marfan numbers: corrected-data run (kms.html, 2,491 control records fixed).
# =====================================================================
out <- "thesis_figures"
dir.create(out, showWarnings = FALSE)

# ---- Figures: one condition each, case vs matched control, p from primary analysis ----
save_thesis(plot_case_control("Ehlers-Danlos syndrome", 210, 1816, 250, 7264, "p < 0.001"),
            "fig_eds_bar", "half", height_cm = 8, dir = out)
save_thesis(plot_case_control("Hypermobility syndrome", 88, 1001, 158, 4004, "p < 0.001"),
            "fig_hyp_bar", "half", height_cm = 8, dir = out)
save_thesis(plot_case_control("Marfan syndrome", 6, 37, 14, 148, "p = 0.24"),
            "fig_ukb_marfan_bar", "half", height_cm = 8, dir = out)

# ---- Figure: All of Us conditions together (Marfan and OI cannot be shown: counts <20) ----
prev <- tribble(
  ~condition,                ~group,             ~events, ~n,   ~p_lab,
  "Ehlers-Danlos syndrome",  "Cases",            210,     1816, "p < 0.001",
  "Ehlers-Danlos syndrome",  "Matched controls", 250,     7264, "p < 0.001",
  "Hypermobility syndrome",  "Cases",            88,      1001, "p < 0.001",
  "Hypermobility syndrome",  "Matched controls", 158,     4004, "p < 0.001"
)
save_thesis(plot_prevalence(prev), "fig_aou_prevalence", "full", height_cm = 8, dir = out)

# ---- Figure: primary matched analysis, both cohorts ----
main <- tribble(
  ~panel,        ~label,                                                            ~OR,  ~lo,  ~hi,
  "All of Us",   "Ehlers-Danlos syndrome\n1,816 cases / 7,264 controls",            3.78, 3.11, 4.61,
  "All of Us",   "Marfan syndrome\n124 cases / 496 controls",                       1.96, 0.90, 4.23,
  "All of Us",   "Hypermobility syndrome*\n1,001 cases / 4,004 controls",           2.38, 1.81, 3.13,
  "All of Us",   "Osteogenesis imperfecta*\n73 cases / 292 controls",               1.07, 0.34, 3.39,
  "UK Biobank",  "Marfan syndrome\n37 cases / 148 controls",                        1.89, 0.65, 5.47
)
fig_main <- plot_forest(main, breaks = c(0.25, 0.5, 1, 2, 4, 8))
save_thesis(fig_main, "fig_aim1_forest", "full", height_cm = 9, dir = out)

# ---- Supplementary figure: sensitivity analyses, both cohorts ----
sens <- tribble(
  ~panel,                                  ~label,                          ~OR,  ~lo,  ~hi,   ~bold,
  "Ehlers-Danlos syndrome (All of Us)",    "Matched, age + BMI (primary)",  3.78, 3.11, 4.61,  TRUE,
  "Ehlers-Danlos syndrome (All of Us)",    "Matched, age only",             3.49, 2.89, 4.21,  FALSE,
  "Ehlers-Danlos syndrome (All of Us)",    "Matched + adjusted",            3.65, 2.95, 4.50,  FALSE,
  "Ehlers-Danlos syndrome (All of Us)",    "Unmatched",                     2.48, 2.15, 2.86,  FALSE,
  "Marfan syndrome (All of Us)",           "Matched, age + BMI (primary)",  1.96, 0.90, 4.23,  TRUE,
  "Marfan syndrome (All of Us)",           "Matched, age only",             2.10, 0.95, 4.62,  FALSE,
  "Marfan syndrome (All of Us)",           "Matched + adjusted",            NA,   NA,   NA,    FALSE,
  "Marfan syndrome (All of Us)",           "Unmatched",                     1.65, 0.77, 3.14,  FALSE,
  "Marfan syndrome (UK Biobank)",          "Matched 1:4 (primary)",         1.89, 0.65, 5.47,  TRUE,
  "Marfan syndrome (UK Biobank)",          "Matched 1:4, Fisher's exact",   1.85, 0.54, 5.63,  FALSE,
  "Marfan syndrome (UK Biobank)",          "Matched 1:4 + first birth <25", 2.27, 0.74, 6.96,  FALSE,
  "Marfan syndrome (UK Biobank)",          "Matched 1:10",                  2.70, 1.02, 7.16,  FALSE,
  "Marfan syndrome (UK Biobank)",          "Full cohort, adjusted",         2.68, 1.10, 6.57,  FALSE,
  "Hypermobility syndrome (All of Us)",    "Matched, age + BMI (primary)",  2.38, 1.81, 3.13,  TRUE,
  "Hypermobility syndrome (All of Us)",    "Matched, age only",             2.62, 2.00, 3.44,  FALSE,
  "Hypermobility syndrome (All of Us)",    "Matched + adjusted",            2.55, 1.88, 3.45,  FALSE,
  "Hypermobility syndrome (All of Us)",    "Unmatched",                     1.89, 1.51, 2.34,  FALSE,
  "Osteogenesis imperfecta (All of Us)",   "Matched, age + BMI (primary)",  1.07, 0.34, 3.39,  TRUE,
  "Osteogenesis imperfecta (All of Us)",   "Matched, age only",             1.37, 0.47, 3.96,  FALSE,
  "Osteogenesis imperfecta (All of Us)",   "Matched + adjusted",            NA,   NA,   NA,    FALSE,
  "Osteogenesis imperfecta (All of Us)",   "Unmatched",                     1.38, 0.44, 3.37,  FALSE
)
fig_sens <- plot_forest(sens, breaks = c(0.25, 0.5, 1, 2, 4, 8, 16))
save_thesis(fig_sens, "fig_aim1_sensitivity", "full", height_cm = 15, dir = out)

# ---- Supplementary figure: UK Biobank covariate balance ----
bal <- tribble(
  ~group,               ~variable,                     ~before,  ~after,
  NA,                   "Townsend deprivation",        0.4322,   0.0759,
  NA,                   "BMI",                        -0.2747,   0.1320,
  NA,                   "Number of live births",      -0.4834,   0.0059,
  NA,                   "Birth cohort (largest band)", 0.1567,   0.0000,
  "Smoking status",     "Never",                      -0.0728,   0.0550,
  "Smoking status",     "Ever",                        0.0244,  -0.0818,
  "Smoking status",     "Unknown",                     0.1460,   0.0833,
  "Menopausal status",  "No",                          0.0082,  -0.0164,
  "Menopausal status",  "Yes",                        -0.1153,   0.0136,
  "Menopausal status",  "Hysterectomy",                0.1151,   0.0183,
  "Menopausal status",  "Unsure/unknown",              0.0501,  -0.0299
)
fig_bal <- plot_balance(bal)
save_thesis(fig_bal, "fig_ukb_balance", "full", height_cm = 10, dir = out)

# ---- Figure: Marfan syndrome, both cohorts and pooled (fixed-effect, inverse variance) ----
marfan <- tribble(
  ~label,                                         ~OR,  ~lo,  ~hi,   ~pooled, ~bold,
  "All of Us\n124 cases / 496 controls",          1.96, 0.90, 4.23,  FALSE,   FALSE,
  "UK Biobank\n37 cases / 148 controls",          1.89, 0.65, 5.47,  FALSE,   FALSE,
  "Pooled (fixed-effect)\nI\u00b2 = 0%",        1.94, 1.03, 3.62,  TRUE,    TRUE
)
fig_marfan <- plot_forest(marfan, breaks = c(0.5, 1, 2, 4, 8), label_frac = 0.75)
save_thesis(fig_marfan, "fig_marfan_pooled", "full", height_cm = 4.8, dir = out)

# ---- Figure: all conditions by OR (the only way Marfan/OI can be shown in All of Us) ----
allc <- tribble(
  ~panel,                     ~label,                                              ~OR,  ~lo,  ~hi,   ~pooled,
  "Primary conditions",       "Ehlers-Danlos syndrome (All of Us)",               3.78, 3.11, 4.61,  FALSE,
  "Primary conditions",       "Marfan syndrome (All of Us)",                      1.96, 0.90, 4.23,  FALSE,
  "Primary conditions",       "Marfan syndrome (UK Biobank)",                     1.89, 0.65, 5.47,  FALSE,
  "Primary conditions",       "Marfan syndrome (pooled)",                         1.94, 1.03, 3.62,  TRUE,
  "Exploratory conditions",   "Hypermobility syndrome (All of Us)",               2.38, 1.81, 3.13,  FALSE,
  "Exploratory conditions",   "Osteogenesis imperfecta (All of Us)",              1.07, 0.34, 3.39,  FALSE
)
fig_all <- plot_forest(allc, breaks = c(0.25, 0.5, 1, 2, 4, 8, 16), label_frac = 1.0)
save_thesis(fig_all, "fig_aim1_all_conditions", "full", height_cm = 8, dir = out)

# ---- Figure: All of Us only, all four conditions by OR ----
aou_all <- tribble(
  ~panel,                   ~label,                                                        ~OR,  ~lo,  ~hi,
  "Primary conditions",     "Ehlers-Danlos syndrome\n1,816 cases / 7,264 controls",         3.78, 3.11, 4.61,
  "Primary conditions",     "Marfan syndrome\n124 cases / 496 controls",                    1.96, 0.90, 4.23,
  "Exploratory conditions", "Hypermobility syndrome\n1,001 cases / 4,004 controls",         2.38, 1.81, 3.13,
  "Exploratory conditions", "Osteogenesis imperfecta\n73 cases / 292 controls",             1.07, 0.34, 3.39
)
fig_aou_all <- plot_forest(aou_all, breaks = c(0.25, 0.5, 1, 2, 4, 8), label_frac = 0.95)
save_thesis(fig_aou_all, "fig_aou_all_conditions", "full", height_cm = 8.5, dir = out)

# ---- Figure: osteogenesis imperfecta, all analyses (All of Us) ----
oi <- tribble(
  ~label,                                         ~OR,  ~lo,  ~hi,  ~bold,
  "Matched, age + BMI (primary)\n73 cases / 292 controls", 1.07, 0.34, 3.39, TRUE,
  "Matched, age only\n78 cases / 312 controls",            1.37, 0.47, 3.96, FALSE,
  "Matched + adjusted",                                     NA,   NA,   NA,   FALSE,
  "Unmatched",                                              1.38, 0.44, 3.37, FALSE
)
save_thesis(plot_forest(oi, breaks = c(0.25, 0.5, 1, 2, 4, 8), label_frac = 0.85),
            "fig_oi_analyses", "full", height_cm = 6.5, dir = out)

# ---- Figure: Marfan syndrome, all analyses (All of Us) ----
mf <- tribble(
  ~label,                                                   ~OR,  ~lo,  ~hi,  ~bold,
  "Matched, age + BMI (primary)\n124 cases / 496 controls", 1.96, 0.90, 4.23, TRUE,
  "Matched, age only",                                      2.10, 0.95, 4.62, FALSE,
  "Matched + adjusted",                                     NA,   NA,   NA,   FALSE,
  "Unmatched",                                              1.65, 0.77, 3.14, FALSE
)
save_thesis(plot_forest(mf, breaks = c(0.25, 0.5, 1, 2, 4, 8), label_frac = 0.85),
            "fig_marfan_aou_analyses", "full", height_cm = 6.5, dir = out)
