#' wet2R: qPCR analysis and visualization utilities
#'
#' Tools for qPCR data analysis. Relative expression is calculated with the
#' 2^-ddCt method, and results are plotted as bar charts with error bars and
#' biological replicate points.
#'
#' The typical workflow has three steps: \code{\link{prepare_qpcr}} reads and
#' calculates the data, \code{\link{plot_qpcr}} creates the plot, and
#' \code{\link{save_qpcr}} exports the figures and tables.
#'
#' @keywords internal
#' @importFrom utils read.csv write.csv globalVariables
#' @importFrom stats aggregate sd
#' @importFrom grid unit
#' @importFrom scales breaks_pretty
#' @importFrom systemfonts system_fonts
#' @importFrom ragg agg_png
#' @importFrom svglite svglite
#' @importFrom ggplot2 ggplot aes geom_col geom_errorbar geom_point
#'   scale_fill_manual scale_x_continuous scale_y_continuous labs
#'   theme_classic theme element_text element_blank element_line
#'   expansion margin ggsave
"_PACKAGE"

# Variables used in ggplot2::aes() must be declared to avoid R CMD check notes.
utils::globalVariables(c("x", "geometric_mean", "group", "error_lower",
                         "error_upper", "relative_expression"))

# Default categorical palette (light mode, from the dataviz reference palette).
# The first three colors pass the colorblind and contrast checks.
qpcr_palette <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                  "#e87ba4", "#008300", "#6250d6", "#e34948")

#' Read qPCR Ct values and calculate relative expression
#'
#' Reads a CSV file of Ct values, validates the format, and calculates relative
#' expression with the 2^-ddCt method. Each biological replicate is averaged over
#' its technical replicates and contributes one dCt value; technical replicates
#' are not weighted by their number.
#'
#' The CSV file must contain the following columns:
#' \describe{
#'   \item{group}{Group name, e.g. WT or KO.}
#'   \item{bio_rep}{Biological replicate ID.}
#'   \item{tech_rep}{Technical replicate ID.}
#'   \item{target_ct}{Ct value of the target gene.}
#'   \item{reference_ct}{Ct value of the reference (housekeeping) gene.}
#' }
#'
#' @param csv_file Path to the input CSV file.
#' @param control Name of the control group. Must be one of the groups in the
#'   CSV file. Default is "WT".
#' @param group_order Character vector giving the plotting order of the groups.
#'   It must contain every group exactly once. Default puts the control group
#'   first, followed by the other groups in order of appearance.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{raw}{Raw data with the required columns only.}
#'   \item{biological_data}{Per-biological-replicate mean Ct, dCt, ddCt, and
#'     relative expression.}
#'   \item{summary}{Per-group mean, standard deviation, geometric mean, and
#'     error range.}
#'   \item{technical_qc}{Range of technical replicate Ct values within each
#'     biological replicate, for quality control.}
#'   \item{control}{Name of the control group.}
#'   \item{group_order}{Plotting order of the groups.}
#' }
#'
#' @details
#' Each group needs at least 2 biological replicates to calculate the standard
#' deviation. Ct values must be finite positive numbers. Undetected reactions
#' must be checked before the data are used.
#'
#' @examples
#' \dontrun{
#' result <- prepare_qpcr("data/example_qpcr.csv", control = "WT")
#' result$summary
#' }
#' @export
prepare_qpcr <- function(csv_file, control = "WT", group_order = NULL) {
  raw <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = c(group = "character", bio_rep = "character",
                                 tech_rep = "character", target_ct = "numeric",
                                 reference_ct = "numeric"))
  required <- c("group", "bio_rep", "tech_rep", "target_ct", "reference_ct")
  if (!all(required %in% names(raw))) {
    stop("CSV must contain the columns: ", paste(required, collapse = ", "))
  }
  raw <- raw[required]
  if (!nrow(raw) || anyNA(raw)) stop("Input must not be empty or contain missing values.")
  for (column in c("group", "bio_rep", "tech_rep")) {
    if (any(!nzchar(trimws(raw[[column]])))) stop("Empty IDs found in column: ", column)
  }
  if (anyDuplicated(raw[c("group", "bio_rep", "tech_rep")])) {
    stop("group + bio_rep + tech_rep must be unique; check for duplicate entries.")
  }
  if (any(!is.finite(raw$target_ct) | raw$target_ct <= 0 |
          !is.finite(raw$reference_ct) | raw$reference_ct <= 0)) {
    stop("Ct values must be finite positive numbers; verify undetected reactions first.")
  }

  groups <- unique(raw$group)
  if (length(control) != 1L || !control %in% groups) {
    stop("control must match one control group in the CSV file.")
  }
  if (is.null(group_order)) group_order <- c(control, groups[groups != control])
  if (anyDuplicated(group_order) || !setequal(group_order, groups)) {
    stop("group_order must contain every group name exactly once.")
  }

  # Average technical replicates within each biological replicate
  bio <- aggregate(cbind(target_ct, reference_ct) ~ group + bio_rep, raw, mean)
  if (any(table(bio$group) < 2L)) {
    stop("Each group needs at least 2 biological replicates to calculate the SD.")
  }
  bio <- bio[order(match(bio$group, group_order), bio$bio_rep), ]
  rownames(bio) <- NULL
  bio$delta_ct <- bio$target_ct - bio$reference_ct
  control_mean <- mean(bio$delta_ct[bio$group == control])
  bio$delta_delta_ct <- bio$delta_ct - control_mean
  bio$relative_expression <- 2^(-bio$delta_delta_ct)

  group_summary <- do.call(rbind, lapply(group_order, function(group) {
    d <- bio[bio$group == group, ]
    center <- mean(d$delta_delta_ct)
    spread <- sd(d$delta_ct)
    data.frame(group = group, n_bio = nrow(d),
               mean_target_ct = mean(d$target_ct),
               mean_reference_ct = mean(d$reference_ct),
               mean_delta_ct = mean(d$delta_ct),
               mean_delta_delta_ct = center, sd_delta_ct = spread,
               geometric_mean = 2^(-center),
               error_lower = 2^(-center - spread),
               error_upper = 2^(-center + spread))
  }))

  technical_qc <- do.call(rbind, lapply(seq_len(nrow(bio)), function(i) {
    d <- raw[raw$group == bio$group[i] & raw$bio_rep == bio$bio_rep[i], ]
    data.frame(group = bio$group[i], bio_rep = bio$bio_rep[i], n_technical = nrow(d),
               target_ct_range = diff(range(d$target_ct)),
               reference_ct_range = diff(range(d$reference_ct)))
  }))

  if (any(!is.finite(bio$relative_expression)) ||
      any(!is.finite(group_summary$geometric_mean)) ||
      any(!is.finite(group_summary$error_upper))) {
    stop("Exponential transformation produced non-finite values; check the input Ct values.")
  }

  list(raw = raw, biological_data = bio, summary = group_summary,
       technical_qc = technical_qc, control = control, group_order = group_order)
}

#' Plot relative expression as a bar chart
#'
#' Bars show the geometric mean, error bars show the range derived from the
#' standard deviation of dCt, and points show individual biological replicates.
#' With 4 or fewer groups, group names are shown directly on the x axis. With
#' more than 4 groups, a legend is shown instead.
#'
#' @param result A list returned by \code{\link{prepare_qpcr}}.
#' @param gene Name of the target gene, used in the y-axis title.
#' @param group_labels Character vector of x-axis labels, one per group. Must
#'   have the same length as the number of groups. Default uses the group names.
#' @param colors Character vector of fill colors, one per group. Default uses
#'   the built-in palette; it must be specified when there are more than 8 groups.
#' @param font Font family name. If the font is not installed, a warning is
#'   given and the default sans-serif font is used.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' result <- prepare_qpcr("data/example_qpcr.csv", control = "WT")
#' p <- plot_qpcr(result, gene = "Il6", group_labels = c("Control", "Treatment"))
#' print(p)
#' }
#' @export
plot_qpcr <- function(result, gene, group_labels = NULL, colors = NULL, font = "Arial") {
  groups <- result$group_order
  if (length(gene) != 1L || !nzchar(gene)) stop("gene must be a single gene name.")
  if (!font %in% system_fonts()$family) {
    warning("Font not installed: ", font, ". Using the default font instead.", call. = FALSE)
    font <- "sans"
  }
  if (is.null(group_labels)) group_labels <- groups
  if (length(group_labels) != length(groups)) {
    stop("group_labels must have the same length as the number of groups.")
  }
  if (is.null(colors)) {
    if (length(groups) > length(qpcr_palette)) {
      stop("For more than 8 groups, specify the colors argument.")
    }
    colors <- qpcr_palette[seq_along(groups)]
  }
  if (length(colors) != length(groups)) {
    stop("colors must have the same length as the number of groups.")
  }
  names(colors) <- groups

  d <- result$biological_data
  s <- result$summary
  s$x <- match(s$group, groups)
  d$x <- match(d$group, groups)
  # Spread replicate points horizontally within each group to avoid overlap
  for (group in groups) {
    index <- which(d$group == group)
    offset <- if (length(index) > 1L) seq(-0.12, 0.12, length.out = length(index)) else 0
    d$x[index] <- d$x[index] + offset
  }
  upper <- max(s$error_upper, d$relative_expression) * 1.13
  show_legend <- length(groups) > 4L

  ggplot(s, aes(x = x, y = geometric_mean, fill = group)) +
    geom_col(width = 0.60, colour = "#333333", linewidth = 0.35) +
    geom_errorbar(aes(ymin = error_lower, ymax = error_upper),
                  width = 0.18, linewidth = 0.40) +
    geom_point(data = d, aes(x = x, y = relative_expression),
               inherit.aes = FALSE, shape = 21, size = 1.9, stroke = 0.45,
               colour = "#222222", fill = "white") +
    scale_fill_manual(values = colors, labels = group_labels,
                      guide = if (show_legend) "legend" else "none") +
    scale_x_continuous(breaks = seq_along(groups), labels = group_labels,
                       limits = c(0.45, length(groups) + 0.55),
                       expand = expansion(mult = 0)) +
    scale_y_continuous(limits = c(0, upper), expand = expansion(mult = 0),
                       breaks = breaks_pretty(n = 4)) +
    labs(x = NULL, y = paste("Relative", gene, "expression"), fill = NULL) +
    theme_classic(base_size = 10, base_family = font) +
    theme(axis.title.y = element_text(size = 10, margin = margin(r = 7)),
          axis.text = element_text(colour = "black", size = 9),
          axis.text.x = element_text(margin = margin(t = 5), lineheight = 0.95),
          axis.ticks.x = element_blank(),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks.y = element_line(linewidth = 0.35),
          axis.ticks.length = unit(2.5, "pt"),
          legend.position = if (show_legend) "right" else "none",
          plot.margin = margin(t = 5, r = 7, b = 5, l = 5))
}

#' Export a qPCR figure and summary tables
#'
#' Saves the plot as a 600 dpi PNG and as an SVG, and saves the biological
#' replicate data, group summary, and technical QC table as CSV files. The output
#' directory is created if it does not exist.
#'
#' @param result A list returned by \code{\link{prepare_qpcr}}.
#' @param plot A \code{ggplot} object returned by \code{\link{plot_qpcr}}.
#' @param output_prefix Output path prefix, e.g. "results/Il6".
#' @param width Figure width in inches. Default is 3.4.
#' @param height Figure height in inches. Default is 3.5.
#' @param dpi Resolution of the PNG file. Default is 600.
#'
#' @return The input \code{plot}, returned invisibly.
#'
#' @examples
#' \dontrun{
#' save_qpcr(result, p, output_prefix = "results/Il6")
#' }
#' @export
save_qpcr <- function(result, plot, output_prefix, width = 3.4, height = 3.5, dpi = 600) {
  if (any(!is.finite(c(width, height, dpi))) || any(c(width, height, dpi) <= 0)) {
    stop("width, height, and dpi must be positive numbers.")
  }
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(output_prefix, ".png"), plot = plot, device = agg_png,
         width = width, height = height, units = "in", dpi = dpi, bg = "white")
  ggsave(paste0(output_prefix, ".svg"), plot = plot, device = svglite,
         width = width, height = height, units = "in", bg = "white")
  write.csv(result$biological_data, paste0(output_prefix, "_biological.csv"),
            row.names = FALSE)
  write.csv(result$summary, paste0(output_prefix, "_summary.csv"), row.names = FALSE)
  write.csv(result$technical_qc, paste0(output_prefix, "_technical_qc.csv"),
            row.names = FALSE)
  invisible(plot)
}
