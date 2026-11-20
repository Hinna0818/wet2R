#' Plot qPCR fold changes and optional statistical comparisons
#'
#' Draw geometric-mean bars, back-transformed dCt SD ranges, and one point
#' per biological replicate. The returned ggplot can be customized or saved
#' with [ggplot2::ggsave()].
#'
#' @param result A result from [prepare_qpcr()].
#' @param gene Target gene name, used as the centered plot title.
#' @param group_labels Group labels in `group_order`, or a vector named by group.
#' @param colors Fill colors in `group_order`, or a vector named by group.
#'   Supply one per group. Defaults to distinct gray shades, from dark to light.
#' @param font Font family. Defaults to `"sans"`. Unavailable named fonts cause
#'   an error; generic `"sans"`, `"serif"`, and `"mono"` are accepted.
#' @param tests Optional table from [test_qpcr()] for these biological data.
#'   Subsetting rows does not recompute the original correction family.
#' @param significance_label `"stars"` (default) or `"p_adj"` for adjusted P values.
#'   Stars use adjusted P <= 0.05, 0.01, 0.001, and 0.0001; otherwise `"ns"`.
#' @param base_size Base font size in points.
#' @param bar_width Bar width, greater than zero and at most one.
#' @param point_size Biological replicate point size.
#' @param label_angle X-axis label angle in degrees, between -90 and 90. Default 45.
#' @param show_legend Whether to show a legend. X-axis labels are always shown.
#'
#' @return A ggplot object with a gene title and no caption. Test settings remain
#'   available in the table returned by [test_qpcr()].
#' @examples
#' path <- system.file("extdata", "qpcr_hsp12.csv", package = "wet2R")
#' result <- prepare_qpcr(path, control = "C-")
#' tests <- test_qpcr(result, comparisons = list(c("P-", "C-"), c("T-", "C-")))
#' p <- plot_qpcr(result, gene = "HSP12", tests = tests)
#' p + ggplot2::theme(legend.position = "bottom")
#' @export
plot_qpcr <- function(result, gene, group_labels = NULL, colors = NULL, font = "sans",
                      tests = NULL, significance_label = c("stars", "p_adj"),
                      base_size = 12, bar_width = 0.6, point_size = 2.5,
                      label_angle = 45, show_legend = FALSE) {
  if (!inherits(result, "wet2r_qpcr")) stop("result must come from prepare_qpcr().")
  qpcr_string(gene, "gene")
  qpcr_string(font, "font")
  if (!font %in% c("sans", "serif", "mono") && !font %in% systemfonts::system_fonts()$family) {
    stop("Font is not installed: ", font)
  }
  significance_label <- match.arg(significance_label)
  qpcr_positive(base_size, "base_size")
  qpcr_positive(bar_width, "bar_width")
  qpcr_positive(point_size, "point_size")
  if (bar_width > 1) stop("bar_width must be at most 1.")
  if (!is.numeric(label_angle) || length(label_angle) != 1L ||
      !is.finite(label_angle) || abs(label_angle) > 90) stop("label_angle must be between -90 and 90.")
  if (!is.logical(show_legend) || length(show_legend) != 1L || is.na(show_legend)) {
    stop("show_legend must be TRUE or FALSE.")
  }
  groups <- result$group_order
  if (is.null(group_labels)) group_labels <- groups
  group_labels <- qpcr_group_values(group_labels, groups, "group_labels")
  if (is.null(colors)) {
    colors <- grDevices::gray.colors(length(groups), start = 0.45, end = 0.85)
  }
  colors <- qpcr_group_values(colors, groups, "colors")
  grDevices::col2rgb(colors)
  d <- result$biological_data
  s <- result$summary
  s$x <- match(s$group, groups)
  d$x <- match(d$group, groups)
  p <- ggplot2::ggplot(s, ggplot2::aes(x = x, y = geometric_mean, fill = group)) +
    ggplot2::geom_col(width = bar_width, colour = "black", linewidth = 0.8) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = error_lower, ymax = error_upper),
                           width = bar_width / 2, linewidth = 0.45) +
    ggplot2::geom_point(data = d, ggplot2::aes(x = x, y = relative_expression),
                        inherit.aes = FALSE, shape = 16, size = point_size,
                        colour = "black",
                        position = ggplot2::position_jitter(width = bar_width / 5, height = 0, seed = 1)) +
    ggplot2::scale_fill_manual(values = colors, limits = groups, breaks = groups,
                               labels = unname(group_labels)) +
    ggplot2::scale_x_continuous(breaks = seq_along(groups), labels = unname(group_labels),
                                limits = c(0.4, length(groups) + 0.6), expand = ggplot2::expansion()) +
    ggplot2::scale_y_continuous(limits = c(0, NA), breaks = scales::breaks_pretty(n = 5),
                                expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::theme_classic(base_size = base_size, base_family = font) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                     margin = ggplot2::margin(b = 12)),
                    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
                    axis.text = ggplot2::element_text(colour = "black", face = "bold"),
                    axis.text.x = ggplot2::element_text(angle = label_angle,
                                                       vjust = if (label_angle == 0) 0.5 else 1,
                                                       hjust = if (label_angle == 0) 0.5 else if (label_angle > 0) 1 else 0),
                    axis.line = ggplot2::element_line(linewidth = 0.65),
                    axis.ticks = ggplot2::element_line(linewidth = 0.65),
                    legend.position = if (show_legend) "right" else "none",
                    plot.margin = ggplot2::margin(8, 10, 8, 8))
  if (!is.null(tests)) {
    required <- c("group1", "group2", "p_adj", "method", "p_adjust_method", "n_comparisons")
    if (!is.data.frame(tests) || !all(required %in% names(tests)) ||
        !identical(attr(tests, "biological_data"), result$biological_data)) {
      stop("tests must come from test_qpcr() for these biological data.")
    }
    if (!nrow(tests) || anyNA(tests[required]) || !is.numeric(tests$p_adj) ||
        any(!is.finite(tests$p_adj) | tests$p_adj < 0 | tests$p_adj > 1) ||
        !all(c(tests$group1, tests$group2) %in% groups)) stop("Invalid comparison table.")
    top <- max(s$error_upper, d$relative_expression)
    step <- top * 0.14
    annotations <- data.frame(x = match(tests$group1, groups),
                              xend = match(tests$group2, groups),
                              y = top + step * seq_len(nrow(tests)))
    annotations$label <- if (significance_label == "p_adj") {
      paste0(ifelse(tests$p_adjust_method == "none", "P ", "adj. P "),
              ifelse(tests$p_adj < 0.001, "< 0.001",
                      paste0("= ", format(tests$p_adj, digits = 3, trim = TRUE))))
    } else {
      as.character(cut(tests$p_adj, c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                        labels = c("****", "***", "**", "*", "ns")))
    }
    segments <- rbind(
      transform(annotations, yend = y),
      transform(annotations, xend = x, yend = y - step / 4),
      transform(annotations, x = xend, yend = y - step / 4)
    )
    p <- p + ggplot2::geom_segment(data = segments,
                                    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                                    inherit.aes = FALSE, linewidth = 0.65) +
      ggplot2::geom_text(data = annotations,
                         ggplot2::aes(x = (x + xend) / 2, y = y + step / 8, label = label),
                         inherit.aes = FALSE, family = font,
                         fontface = if (significance_label == "stars") "bold" else "plain",
                         size = base_size / 3, vjust = 0)
  }
  p + ggplot2::labs(title = gene, x = NULL, y = "Relative expression", fill = NULL)
}

qpcr_group_values <- function(values, groups, argument) {
  if (!is.character(values) || length(values) != length(groups) ||
      anyNA(values) || any(!nzchar(trimws(values)))) {
    stop(argument, " must contain one nonempty value per group.")
  }
  if (is.null(names(values))) {
    names(values) <- groups
  } else {
    if (anyDuplicated(names(values)) || !setequal(names(values), groups)) {
      stop(argument, " names must match all groups exactly once.")
    }
    values <- values[groups]
  }
  values
}
