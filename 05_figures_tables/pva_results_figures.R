# =====================================================================
# Aim 2 (rare-variant analysis) figures. Run from the repository root.
# Numbers typed from the SHAREABLE Workbench output only (CELL 6 / 6b / 7 / 8 printouts,
# every count shown is 0 or >= 20). Nothing from the "FULL results" view is used.
# =====================================================================
theme_src <- readLines("thesis_theme.R"); theme_src <- sub("library(tidyverse)", "suppressMessages({library(ggplot2); library(dplyr); library(tibble); library(tidyr); library(purrr)})", theme_src, fixed = TRUE); eval(parse(text = theme_src))
out <- "thesis_figures"

grp_levels <- c("Early-onset POP\n(cases)", "Complete\ncomparison group", "Aged 18-35,\nno menopause",
                "Aged >40,\nno menopause", "Aged >40,\nmenopause")

# ---- Figure: carrier frequency, 167-gene panel vs negative-control genes (primary tier) ----
carr <- tribble(
  ~panel,                                   ~group,        ~events, ~n,
  "Connective tissue panel (167 genes)",    grp_levels[1], 135,     1330,
  "Connective tissue panel (167 genes)",    grp_levels[2], 30532,   320028,
  "Connective tissue panel (167 genes)",    grp_levels[3], 4803,    50191,
  "Connective tissue panel (167 genes)",    grp_levels[4], 17651,   185861,
  "Connective tissue panel (167 genes)",    grp_levels[5], 4858,    50301,
  "Negative-control genes (10 genes)",      grp_levels[1], 319,     1330,
  "Negative-control genes (10 genes)",      grp_levels[2], 68472,   320028,
  "Negative-control genes (10 genes)",      grp_levels[3], 10627,   50191,
  "Negative-control genes (10 genes)",      grp_levels[4], 39190,   185861,
  "Negative-control genes (10 genes)",      grp_levels[5], 11305,   50301
) %>% rowwise() %>%
  mutate(pct = 100 * events / n,
         lo  = 100 * prop.test(events, n)$conf.int[1],
         hi  = 100 * prop.test(events, n)$conf.int[2]) %>% ungroup() %>%
  mutate(group = factor(group, levels = grp_levels),
         role  = ifelse(group == grp_levels[1], "Early-onset POP (cases)", "Comparison groups"),
         role  = factor(role, levels = c("Early-onset POP (cases)", "Comparison groups")),
         panel = factor(panel, levels = unique(panel)))

top <- max(carr$hi)
fig_carr <- ggplot(carr, aes(group, pct, fill = role)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, colour = thesis_col[["ns"]], linewidth = 0.4) +
  geom_text(aes(y = hi + top * 0.03, label = sprintf("%.1f%%", pct)), size = 2.8, vjust = 0,
            colour = thesis_col[["ink"]]) +
  facet_wrap(~panel, ncol = 1) +
  scale_fill_manual(values = c(`Early-onset POP (cases)` = thesis_col[["case"]],
                               `Comparison groups` = thesis_col[["comparison"]])) +
  scale_y_continuous(limits = c(0, top * 1.15), expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Carriers of a qualifying variant") +
  theme_thesis() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(colour = thesis_col[["ink"]], size = rel(0.85), lineheight = 0.95),
        legend.position = "none")
save_thesis(fig_carr, "fig_pva_carrier_frequency", "full", height_cm = 13, dir = out)

# ---- Supplementary figure: ORs across tiers, panels and sensitivity analyses ----
cmp <- c("Complete comparison group", "Aged 18-35, no menopause", "Aged >40, no menopause", "Aged >40, menopause")
mk <- function(panel, vals) tibble(panel = panel, label = cmp,
                                   OR = vals[c(1, 4, 7, 10)], lo = vals[c(2, 5, 8, 11)], hi = vals[c(3, 6, 9, 12)])
sens <- bind_rows(
  mk("167-gene panel, primary tier (main analysis)",
     c(1.07, 0.90, 1.28,  1.07, 0.89, 1.28,  1.08, 0.90, 1.29,  1.06, 0.88, 1.27)),
  mk("167-gene panel, strict tier",
     c(1.01, 0.72, 1.42,  1.09, 0.77, 1.54,  1.01, 0.72, 1.42,  0.94, 0.67, 1.33)),
  mk("167-gene panel, extended tier",
     c(0.93, 0.83, 1.03,  0.91, 0.82, 1.02,  0.92, 0.82, 1.02,  0.95, 0.85, 1.06)),
  mk("167-gene panel, primary tier, White participants only",
     c(0.99, 0.79, 1.25,  0.95, 0.75, 1.20,  1.00, 0.80, 1.25,  1.00, 0.80, 1.26)),
  mk("Original 146-gene panel, primary tier",
     c(1.07, 0.89, 1.28,  1.07, 0.89, 1.29,  1.08, 0.90, 1.29,  1.05, 0.88, 1.26)),
  mk("Negative-control genes, primary tier",
     c(1.16, 1.02, 1.32,  1.17, 1.03, 1.33,  1.18, 1.04, 1.34,  1.09, 0.96, 1.24))
)
save_thesis(plot_forest(sens, breaks = c(0.5, 0.75, 1, 1.5, 2), label_frac = 1.25, value_frac = 0.7),
            "fig_pva_sensitivity_forest", "full", height_cm = 17, dir = out)

# ---- Supplementary figure: mechanisms with >= 20 case carriers (primary tier) ----
mech <- bind_rows(
  mk("Structural collagens and collagen-associated matrix",
     c(1.05, 0.74, 1.49,  1.01, 0.71, 1.45,  1.04, 0.73, 1.48,  1.10, 0.77, 1.57)),
  mk("Elastic fibres and microfibrils",
     c(1.06, 0.80, 1.42,  1.12, 0.84, 1.50,  1.08, 0.81, 1.43,  1.00, 0.75, 1.34)),
  mk("Syndromic or other, no direct connective tissue mechanism",
     c(1.19, 0.79, 1.78,  1.25, 0.83, 1.88,  1.18, 0.79, 1.77,  1.13, 0.75, 1.70))
)
save_thesis(plot_forest(mech, breaks = c(0.5, 0.75, 1, 1.5, 2), label_frac = 1.25, value_frac = 0.7),
            "fig_pva_mechanisms_forest", "full", height_cm = 9.5, dir = out)
