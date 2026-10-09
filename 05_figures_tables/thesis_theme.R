# =====================================================================
# thesis_theme.R  -  one look for every results figure
# Paste this as the FIRST cell of any notebook that makes a figure
# (All of Us, UK Biobank, assay data), then use the helpers below.
# =====================================================================
library(tidyverse)

# ---- Colours: one meaning each, used the same way in every figure ----
thesis_col <- c(
  case       = "#089099",  # teal: cases / exposure group / estimate whose CI excludes 1
  comparison = "#EEB479",  # peach: controls / comparison group (bar fills)
  ns         = "#55575a",  # estimate whose CI includes 1; secondary text
  extra1     = "#80407A",  # soft plum: 3rd group, ONLY when a figure needs one
  extra2     = "#CF597E",  # dusty rose: 4th group, ONLY with labels on the plot (close to teal for some colourblind readers)
  ink        = "#0b0b0b",  # main text
  muted      = "#898781",  # axis ticks, reference line
  grid       = "#e6e6e3",  # gridlines
  axis       = "#c3c3be"   # axis line
)

# ---- Sizes: save at the size it will be printed (A4, 2.5 cm margins) ----
thesis_width_cm <- c(full = 16, half = 8)

# ---- Theme ----
theme_thesis <- function(base_size = 10, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text               = element_text(colour = thesis_col[["ink"]]),
      axis.text          = element_text(colour = thesis_col[["ns"]], size = rel(0.9)),
      axis.title         = element_text(colour = thesis_col[["ns"]], size = rel(0.9)),
      axis.line.x        = element_line(colour = thesis_col[["axis"]], linewidth = 0.4),
      axis.ticks.x       = element_line(colour = thesis_col[["axis"]], linewidth = 0.4),
      panel.grid.major   = element_line(colour = thesis_col[["grid"]], linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      legend.position    = "top",
      legend.justification = "left",
      legend.title       = element_blank(),
      legend.text        = element_text(size = rel(0.9)),
      strip.text         = element_text(face = "bold", hjust = 0, size = rel(0.95)),
      strip.placement    = "outside",
      panel.spacing      = unit(8, "pt"),
      plot.title         = element_blank(),   # titles go in the Word figure legend
      plot.caption       = element_blank(),   # so do notes and abbreviations
      plot.margin        = margin(6, 6, 6, 6)
    )
}

# ---- Save: PNG (600 dpi, for Word) + PDF, at print size ----
save_thesis <- function(p, name, width = c("full", "half"), height_cm = 8,
                        dir = "/tmp", bucket = NULL) {
  width <- match.arg(width)
  w <- thesis_width_cm[[width]]
  paths <- file.path(dir, paste0(name, c(".png", ".pdf")))
  ggsave(paths[1], p, width = w, height = height_cm, units = "cm", dpi = 600, bg = "white")
  ggsave(paths[2], p, width = w, height = height_cm, units = "cm", bg = "white")
  if (!is.null(bucket)) for (f in paths)
    system2("gsutil", c("cp", f, bucket), stdout = TRUE, stderr = TRUE)
  invisible(paths)
}

# ---- Shared layout: label column | plot | value column, all on one canvas ----
# Rows are laid out top to bottom; an optional grouping column becomes bold header rows,
# with the rows beneath indented. Everything is drawn in data units, so labels can never
# be squeezed or clipped by facet strips.
.layout_rows <- function(df, group_col = NULL, gap = 0.5) {
  if (is.null(group_col) || !group_col %in% names(df)) {
    return(df %>% mutate(.y = -(seq_len(n()) - 1), .type = "row", .text = label))
  }
  out <- list(); y <- 0
  for (g in unique(df[[group_col]])) {
    if (y < 0) y <- y - gap
    out[[length(out) + 1]] <- tibble(.y = y, .type = "header", .text = g)
    sub <- df[df[[group_col]] == g, ]
    out[[length(out) + 1]] <- sub %>% mutate(.y = y - seq_len(n()), .type = "row", .text = label)
    y <- y - nrow(sub) - 1
  }
  bind_rows(out)
}

.axis_layers <- function(breaks_x, break_labels, ymin, ymax, title, ref = NULL, pad = 0) {
  ax_y <- ymin - 0.7
  list(
    geom_segment(data = tibble(x = breaks_x), aes(x = x, xend = x, y = ax_y, yend = ymax + 0.5),
                 colour = thesis_col[["grid"]], linewidth = 0.3, inherit.aes = FALSE),
    if (!is.null(ref)) geom_segment(aes(x = ref, xend = ref, y = ax_y, yend = ymax + 0.5),
                                    linetype = "dashed", colour = thesis_col[["muted"]], linewidth = 0.4,
                                    inherit.aes = FALSE),
    annotate("segment", x = min(breaks_x) - pad, xend = max(breaks_x) + pad, y = ax_y, yend = ax_y,
             colour = thesis_col[["axis"]], linewidth = 0.4),
    annotate("text", x = breaks_x, y = ax_y - 0.55, label = break_labels,
             size = 2.9, colour = thesis_col[["ns"]]),
    annotate("text", x = mean(range(breaks_x)), y = ax_y - 1.45, label = title,
             size = 3, colour = thesis_col[["ns"]])
  )
}

# ---- Forest plot (odds ratios) ----
# df columns: label, OR, lo, hi
# optional:   panel (grouping -> bold header rows), bold (TRUE = bold estimate text)
# Filled teal = CI excludes 1; hollow grey = CI includes 1; NA rows = "not estimable".
# label_frac / value_frac: width of the label and value columns, as multiples of the plot width.
plot_forest <- function(df, xlab = "Odds ratio (95% CI, log scale)",
                        breaks = c(0.25, 0.5, 1, 2, 4, 8),
                        label_frac = 1.15, value_frac = 0.62, text_size = 3) {
  if (!"bold" %in% names(df)) df$bold <- FALSE
  grouped <- "panel" %in% names(df)
  lo_b <- log10(min(breaks)); hi_b <- log10(max(breaks)); span <- hi_b - lo_b
  pad   <- 0.04 * span                    # breathing room inside the axis at both ends
  x_lab <- lo_b - pad - label_frac * span
  x_val <- hi_b + pad + 0.07 * span
  x_end <- hi_b + value_frac * span
  d <- .layout_rows(df, if (grouped) "panel") %>%
    mutate(sig    = .type == "row" & !is.na(OR) & (lo > 1 | hi < 1),
           lOR = log10(OR), llo = pmax(log10(lo), lo_b), lhi = pmin(log10(hi), hi_b),
           ortext = ifelse(is.na(OR), "not estimable", sprintf("%.2f (%.2f-%.2f)", OR, lo, hi)),
           x_text = ifelse(.type == "header", x_lab, x_lab + ifelse(grouped, 0.05 * span, 0)))
  if (!"pooled" %in% names(d)) d$pooled <- FALSE
  d$pooled[is.na(d$pooled)] <- FALSE
  rows <- d %>% filter(.type == "row")
  est  <- rows %>% filter(!is.na(OR), !pooled)
  dia  <- rows %>% filter(!is.na(OR), pooled) %>%
    mutate(id = row_number()) %>%
    reframe(id = rep(id, each = 4), sig = rep(sig, each = 4),
            x = as.vector(t(cbind(llo, lOR, lhi, lOR))),
            y = as.vector(t(cbind(.y, .y + 0.32, .y, .y - 0.32))))
  ymin <- min(d$.y); ymax <- max(d$.y)
  ggplot() +
    .axis_layers(log10(breaks), format(breaks, drop0trailing = TRUE, trim = TRUE),
                 ymin, ymax, xlab, ref = 0, pad = pad) +
    geom_text(data = d, aes(x_text, .y, label = .text,
                            fontface = ifelse(.type == "header", "bold", "plain")),
              hjust = 0, size = text_size, lineheight = 0.95, colour = thesis_col[["ink"]]) +
    geom_segment(data = est, aes(x = llo, xend = lhi, y = .y, yend = .y, colour = sig), linewidth = 0.7) +
    geom_point(data = est, aes(lOR, .y, fill = sig, colour = sig), shape = 21, size = 2.4, stroke = 0.8) +
    geom_polygon(data = dia, aes(x, y, group = id, fill = sig, colour = sig), linewidth = 0.5) +
    geom_text(data = rows, aes(x_val, .y, label = ortext, fontface = ifelse(bold, "bold", "plain")),
              hjust = 0, size = text_size, colour = thesis_col[["ink"]]) +
    scale_colour_manual(values = c(`TRUE` = thesis_col[["case"]], `FALSE` = thesis_col[["ns"]]), guide = "none") +
    scale_fill_manual(values = c(`TRUE` = thesis_col[["case"]], `FALSE` = "white"), guide = "none") +
    coord_cartesian(xlim = c(x_lab, x_end), ylim = c(ymin - 2.6, ymax + 0.6), expand = FALSE, clip = "off") +
    theme_void() + theme(plot.margin = margin(6, 6, 6, 6))
}

# ---- Proportions with 95% CI: cases vs comparison group, several conditions ----
# df columns: condition, group ("Cases" / comparison label), events, n
# optional:   p_lab (one per condition, e.g. "p < 0.001") -> bracket over the pair
plot_prevalence <- function(df, ylab = "Women with a POP diagnosis") {
  df <- df %>% rowwise() %>%
    mutate(pct = 100 * events / n,
           lo  = 100 * prop.test(events, n)$conf.int[1],
           hi  = 100 * prop.test(events, n)$conf.int[2]) %>% ungroup()
  groups <- unique(df$group)
  if (length(groups) > 4) stop("Max 4 groups per figure - split into facets instead")
  cols <- setNames(c(thesis_col[["case"]], thesis_col[["comparison"]], thesis_col[["extra1"]],
                     thesis_col[["extra2"]])[seq_along(groups)], groups)
  df <- df %>% mutate(group = factor(group, levels = groups),
                      condition = factor(condition, levels = unique(condition)))
  top <- max(df$hi)
  p <- ggplot(df, aes(condition, pct, fill = group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.8),
                  width = 0.12, colour = thesis_col[["ns"]], linewidth = 0.4) +
    geom_text(aes(y = hi + top * 0.03, label = sprintf("%.1f%%", pct)),
              position = position_dodge(width = 0.8), size = 3, vjust = 0, colour = thesis_col[["ink"]])
  if ("p_lab" %in% names(df)) {
    br <- df %>% group_by(condition, p_lab) %>% summarise(y = max(hi) + top * 0.2, .groups = "drop") %>%
      mutate(x = as.numeric(condition))
    p <- p +
      geom_segment(data = br, aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y), inherit.aes = FALSE,
                   colour = thesis_col[["ink"]], linewidth = 0.4) +
      geom_segment(data = br, aes(x = x - 0.2, xend = x - 0.2, y = y, yend = y - top * 0.03), inherit.aes = FALSE,
                   colour = thesis_col[["ink"]], linewidth = 0.4) +
      geom_segment(data = br, aes(x = x + 0.2, xend = x + 0.2, y = y, yend = y - top * 0.03), inherit.aes = FALSE,
                   colour = thesis_col[["ink"]], linewidth = 0.4) +
      geom_text(data = br, aes(x, y + top * 0.025, label = p_lab), inherit.aes = FALSE,
                size = 3, vjust = 0, fontface = "italic", colour = thesis_col[["ink"]])
  }
  p + scale_fill_manual(values = cols) +
    scale_y_continuous(limits = c(0, top * 1.38), expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = ylab) +
    theme_thesis() +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(colour = thesis_col[["ink"]], size = rel(1), lineheight = 1.05))
}

# ---- One condition: case bar vs control bar, p-value bracket ----
plot_case_control <- function(case_name, case_events, case_n, ctrl_events, ctrl_n, p_lab,
                              ctrl_name = "Matched controls", ylab = "Women with a POP diagnosis") {
  d <- tibble(group = c(ctrl_name, case_name), events = c(ctrl_events, case_events), n = c(ctrl_n, case_n)) %>%
    rowwise() %>%
    mutate(pct = 100 * events / n,
           lo  = 100 * prop.test(events, n)$conf.int[1],
           hi  = 100 * prop.test(events, n)$conf.int[2]) %>% ungroup() %>%
    mutate(xlab = sprintf("%s\nn = %s", sapply(group, function(g) paste(strwrap(g, 16), collapse = "\n")),
                         format(n, big.mark = ",", trim = TRUE)),
           xlab = factor(xlab, levels = xlab))
  top <- max(d$hi); yb <- top * 1.2
  ggplot(d, aes(xlab, pct, fill = group)) +
    geom_col(width = 0.55) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.1, colour = thesis_col[["ns"]], linewidth = 0.4) +
    geom_text(aes(y = hi + top * 0.03, label = sprintf("%.1f%%", pct)), size = 3.2, vjust = 0,
              fontface = "bold", colour = thesis_col[["ink"]]) +
    annotate("segment", x = 1, xend = 2, y = yb, yend = yb, colour = thesis_col[["ink"]], linewidth = 0.4) +
    annotate("segment", x = c(1, 2), xend = c(1, 2), y = yb, yend = yb - top * 0.035,
             colour = thesis_col[["ink"]], linewidth = 0.4) +
    annotate("text", x = 1.5, y = yb + top * 0.03, label = p_lab, size = 3.1, vjust = 0,
             fontface = "italic", colour = thesis_col[["ink"]]) +
    scale_fill_manual(values = setNames(c(thesis_col[["comparison"]], thesis_col[["case"]]),
                                        c(ctrl_name, case_name)), guide = "none") +
    scale_y_continuous(limits = c(0, top * 1.38), expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = ylab) +
    theme_thesis() +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(colour = thesis_col[["ink"]], size = rel(1), lineheight = 1.05))
}

# ---- Covariate balance (matching): absolute SMD before vs after ----
# df columns: variable, before, after;  optional: group (bold header rows, e.g. "Smoking status")
plot_balance <- function(df, threshold = 0.1, label_frac = 0.62, text_size = 3) {
  names(df)[names(df) == "variable"] <- "label"
  grouped <- "group" %in% names(df)
  if (grouped) df$group[is.na(df$group)] <- ""
  top <- ceiling(max(abs(c(df$before, df$after))) * 10) / 10
  if (top - max(abs(c(df$before, df$after))) < 0.05) top <- top + 0.1   # never let a point sit on the axis end
  breaks <- seq(0, top, by = 0.1)
  pad  <- 0.03 * top
  xmax <- top + pad
  x_lab <- -pad - label_frac * top
  d <- .layout_rows(df, if (grouped) "group") %>%
    filter(!(.type == "header" & .text == "")) %>%
    mutate(x_text = ifelse(.type == "header", x_lab,
                    ifelse(grouped & !is.na(group) & group != "", x_lab + 0.04 * xmax, x_lab)))
  rows <- d %>% filter(.type == "row") %>% mutate(b = abs(before), a = abs(after))
  ymin <- min(d$.y); ymax <- max(d$.y)
  key_y <- ymax + 1.3
  ggplot() +
    .axis_layers(breaks, format(breaks, nsmall = 1), ymin, ymax,
                 "Absolute standardised mean difference", ref = threshold, pad = pad) +
    geom_text(data = d, aes(x_text, .y, label = .text,
                            fontface = ifelse(.type == "header", "bold", "plain")),
              hjust = 0, size = text_size, colour = thesis_col[["ink"]]) +
    geom_segment(data = rows, aes(x = pmin(a, b), xend = pmax(a, b), y = .y, yend = .y),
                 colour = thesis_col[["grid"]], linewidth = 0.7) +
    geom_point(data = rows, aes(b, .y), shape = 21, size = 2.4, stroke = 0.8,
               fill = "white", colour = thesis_col[["ns"]]) +
    geom_point(data = rows, aes(a, .y), shape = 21, size = 2.4, stroke = 0.8,
               fill = thesis_col[["case"]], colour = thesis_col[["case"]]) +
    # key, above the plot
    annotate("point", x = 0.02 * xmax, y = key_y, shape = 21, size = 2.4, stroke = 0.8,
             fill = "white", colour = thesis_col[["ns"]]) +
    annotate("text", x = 0.06 * xmax, y = key_y, label = "Before matching", hjust = 0,
             size = text_size, colour = thesis_col[["ink"]]) +
    annotate("point", x = 0.42 * xmax, y = key_y, shape = 21, size = 2.4, stroke = 0.8,
             fill = thesis_col[["case"]], colour = thesis_col[["case"]]) +
    annotate("text", x = 0.46 * xmax, y = key_y, label = "After matching", hjust = 0,
             size = text_size, colour = thesis_col[["ink"]]) +
    coord_cartesian(xlim = c(x_lab, xmax), ylim = c(ymin - 2.6, key_y + 0.6), expand = FALSE, clip = "off") +
    theme_void() + theme(plot.margin = margin(6, 10, 6, 6))
}
