# =====================================================================
# Participant flow diagrams (STROBE). Run thesis_theme.R first.
# Fill in the two starting counts below, then re-run.
# =====================================================================
aou_total <- "293,202"   # All of Us participants with female sex at birth and EHR data
ukb_total <- "271,791"   # women in the UK Biobank extract (thesis presentation figure)
ukb_excl  <- "54,924"    # 271,754 - 216,830

out <- "thesis_figures"
ink <- thesis_col[["ink"]]; line <- thesis_col[["ns"]]

boxes <- function(b, size = 2.5) list(
  geom_rect(data = b, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
                          fill = kind, linetype = kind), colour = line, linewidth = 0.35),
  geom_text(data = b, aes(x, y, label = label, fontface = ifelse(kind == "main", "bold", "plain")),
            size = size, lineheight = 1.0, colour = ink),
  scale_fill_manual(values = c(main = "#e6f2f3", step = "white", excl = "#f7f7f5"), guide = "none"),
  scale_linetype_manual(values = c(main = "solid", step = "solid", excl = "dashed"), guide = "none"))
arrows <- function(a) geom_segment(data = a, aes(x = x, y = y, xend = xend, yend = yend),
                                   colour = line, linewidth = 0.35,
                                   arrow = if (TRUE) NULL)
heads  <- function(a) geom_segment(data = a, aes(x = x, y = y, xend = xend, yend = yend),
                                   colour = line, linewidth = 0.35,
                                   arrow = arrow(length = unit(1.6, "mm"), type = "closed"))
canvas <- function(xl, yl) list(coord_cartesian(xlim = xl, ylim = yl, expand = FALSE),
                                theme_void(), theme(plot.margin = margin(4, 4, 4, 4)))

# ---------------- All of Us ----------------
cx <- c(1, 2.05, 3.1, 4.15); px <- 5.3          # condition columns, control pool column
cond <- c("Ehlers-Danlos\nsyndrome\nn = 2,015", "Marfan\nsyndrome",
          "Hypermobility\nsyndrome*\nn = 1,086", "Osteogenesis\nimperfecta*")
excl <- c("No BMI recorded\nn = 199", "No BMI recorded\nn < 20",
          "No BMI recorded\nn = 85", "No BMI recorded\nn < 20")
mtch <- c("1,816 cases\n7,264 controls", "124 cases\n496 controls",
          "1,001 cases\n4,004 controls", "73 cases\n292 controls")
b <- bind_rows(
  tibble(x = 3.15, y = 5.0, w = 4.2, h = 0.55, kind = "main",
         label = paste0("All of Us participants with female sex at birth and EHR data\nn = ", aou_total)),
  tibble(x = cx, y = 3.75, w = 0.95, h = 0.75, kind = "step", label = cond),
  tibble(x = px, y = 3.75, w = 1.2, h = 0.75, kind = "step",
         label = "Control pool\nno EDS, Marfan, HS or OI\nn = 289,864"),
  tibble(x = cx, y = 2.65, w = 0.95, h = 0.55, kind = "excl", label = excl),
  tibble(x = cx, y = 1.55, w = 0.95, h = 0.6, kind = "main", label = mtch)
)
a_lines <- bind_rows(
  tibble(x = 3.15, y = 4.725, xend = 3.15, yend = 4.4),              # source down
  tibble(x = cx[1], y = 4.4, xend = px, yend = 4.4),                  # split bar
  tibble(x = px, y = 3.375, xend = px, yend = 0.95),                  # control pool down
  tibble(x = px, y = 0.95, xend = cx[1], yend = 0.95)                 # bus under matched row
)
a_heads <- bind_rows(
  tibble(x = c(cx, px), y = 4.4, xend = c(cx, px), yend = 4.125),     # into condition / pool boxes
  tibble(x = cx, y = 3.375, xend = cx, yend = 2.925),                 # condition -> exclusion
  tibble(x = cx, y = 2.375, xend = cx, yend = 1.85),                  # exclusion -> matched
  tibble(x = cx, y = 0.95, xend = cx, yend = 1.25)                    # control pool -> matched (from below)
)
fig_aou_flow <- ggplot() + geom_segment(data = a_lines, aes(x = x, y = y, xend = xend, yend = yend),
                                        colour = line, linewidth = 0.35) +
  heads(a_heads) + boxes(b) +
  annotate("text", x = 3.15, y = 0.72, label = "Controls drawn from the control pool for each condition, matched 1:4 on age and BMI",
           size = 2.5, fontface = "italic", colour = line) +
  canvas(c(0.4, 6.0), c(0.5, 5.35))
save_thesis(fig_aou_flow, "fig_flow_aou", "full", height_cm = 10, dir = out)

# ---------------- UK Biobank ----------------
ukb_nonmarfan <- "271,754"   # 271,791 - 37
b2 <- bind_rows(
  tibble(x = 2, y = 5.0, w = 3.2, h = 0.55, kind = "main",
         label = paste0("Women in the UK Biobank extract\nn = ", ukb_total)),
  tibble(x = 1, y = 3.85, w = 1.5, h = 0.6, kind = "step", label = "Marfan syndrome (ICD-10 Q87.4)\nn = 37"),
  tibble(x = 3, y = 3.85, w = 1.5, h = 0.6, kind = "step",
         label = paste0("Women without Marfan syndrome\nn = ", ukb_nonmarfan)),
  tibble(x = 3, y = 2.8, w = 1.5, h = 0.55, kind = "excl",
         label = paste0("Excluded: other CTD diagnosis or\nincomplete matching data, n = ", ukb_excl)),
  tibble(x = 3, y = 1.8, w = 1.5, h = 0.55, kind = "step", label = "Control pool\nn = 216,830"),
  tibble(x = 2, y = 0.75, w = 2.2, h = 0.6, kind = "main", label = "Matched 1:4\n37 cases, 148 controls")
)
b2_lines <- tibble(x = c(2, 1), y = c(4.725, 4.4), xend = c(2, 3), yend = c(4.4, 4.4))
b2_heads <- bind_rows(
  tibble(x = c(1, 3), y = 4.4, xend = c(1, 3), yend = 4.15),
  tibble(x = 1, y = 3.55, xend = 1, yend = 1.05),
  tibble(x = 3, y = 3.55, xend = 3, yend = 3.075),
  tibble(x = 3, y = 2.525, xend = 3, yend = 2.075),
  tibble(x = 3, y = 1.525, xend = 3, yend = 1.05)
)
fig_ukb_flow <- ggplot() + geom_segment(data = b2_lines, aes(x = x, y = y, xend = xend, yend = yend),
                                        colour = line, linewidth = 0.35) +
  heads(b2_heads) + boxes(b2, size = 2.7) +
  canvas(c(0.1, 3.9), c(0.35, 5.35))
save_thesis(fig_ukb_flow, "fig_flow_ukb", "full", height_cm = 10, dir = out)
