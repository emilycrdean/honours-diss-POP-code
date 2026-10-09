# All of Us: odds ratio bars for all four conditions (compliant: no counts or prevalence shown)
d <- tribble(
  ~condition,                    ~OR,  ~lo,  ~hi,  ~p_lab,
  "Ehlers-Danlos\nsyndrome",     3.78, 3.11, 4.61, "p < 0.001",
  "Hypermobility\nsyndrome*",    2.38, 1.81, 3.13, "p < 0.001",
  "Marfan\nsyndrome",            1.96, 0.90, 4.23, "p = 0.088",
  "Osteogenesis\nimperfecta*",   1.07, 0.34, 3.39, "p = 0.91"
) %>% mutate(condition = factor(condition, levels = condition),
             x = as.numeric(condition), sig = lo > 1)
fig <- ggplot(d) +
  geom_hline(yintercept = 1, colour = thesis_col[["muted"]], linetype = "dashed", linewidth = 0.4) +
  geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = 1, ymax = OR, fill = sig)) +
  geom_errorbar(aes(x = x, ymin = lo, ymax = hi), width = 0.12, colour = thesis_col[["ns"]], linewidth = 0.4) +
  geom_text(aes(x, hi * 1.22, label = sprintf("%.2f", OR)), size = 3.2, fontface = "bold", colour = thesis_col[["ink"]]) +
  geom_text(aes(x, hi * 1.6, label = p_lab), size = 3, fontface = "italic", colour = thesis_col[["ink"]]) +
  scale_fill_manual(values = c(`TRUE` = thesis_col[["case"]], `FALSE` = "#c9c9c4"), guide = "none") +
  scale_x_continuous(breaks = d$x, labels = levels(d$condition), expand = expansion(add = 0.5)) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = c("0.25", "0.5", "1", "2", "4", "8"),
                limits = c(0.25, 9.5), expand = c(0, 0)) +
  labs(x = NULL, y = "Odds ratio for POP vs matched controls (log scale)") +
  theme_thesis() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(colour = thesis_col[["ink"]], size = rel(1), lineheight = 1.05))
save_thesis(fig, "fig_aou_or_bars", "full", height_cm = 9, dir = "thesis_figures")
