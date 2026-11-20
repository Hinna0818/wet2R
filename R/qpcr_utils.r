#' qPCR analysis and visualization utilities
#'
#' Use [prepare_qpcr()] to calculate expression and technical QC,
#' [test_qpcr()] to compare biological replicates, and [plot_qpcr()]
#' to return an editable ggplot object.
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c("x", "xend", "y", "yend", "label", "group",
                         "geometric_mean", "error_lower", "error_upper",
                         "relative_expression"))

#' Prepare single-target qPCR data
#'
#' Average technical replicates within each biological sample and calculate
#' relative expression using the 2^-ddCt method.
#'
#' @param csv_file CSV path or a data frame containing `group`, `bio_rep`,
#'   `tech_rep`, `target_ct`, and `reference_ct`. Identifiers must be nonempty.
#'   Each row contains one measurement of each assay from the same sample.
#' @param control Control group name. Defaults to `"WT"`.
#' @param group_order All group names, exactly once, in plotting order.
#' @param qc_range_max Optional positive Ct range threshold. Values above the
#'   threshold are flagged, never removed. `NULL` disables this threshold.
#'
#' @details
#' Input must contain one target and one reference assay, with equal numbers
#' of measured technical replicates per sample. Do not duplicate reference
#' values to fill rows. Undetected reactions must be resolved before analysis.
#' Extra columns are sample metadata: they must be complete and constant within
#' each `group`/`bio_rep`. They are retained in `raw` and `biological_data`.
#' Columns named `gene`, `target_gene`, or `reference_gene` must each contain
#' a single value across the input. Sample IDs must distinguish independent
#' samples across plates and batches.
#'
#' Each group needs at least two biological replicates. The control mean dCt
#' is subtracted from every sample dCt. Bars use the geometric mean of fold
#' changes; error bounds back-transform mean ddCt +/- SD(dCt). These describe
#' sample dispersion conditional on the observed control mean, not a confidence
#' interval for a treatment effect. Validate reference stability and comparable,
#' near-doubling amplification efficiencies before using this method.
#'
#' @return A `wet2r_qpcr` list with `raw`, `biological_data`, `summary`,
#'   `technical_qc`, `control`, `group_order`, and `qc_range_max`.
#'   QC includes counts, SDs, ranges, and `qc_flag`/`qc_reason`. A single
#'   technical replicate has undefined SD and is flagged as insufficient QC.
#' @examples
#' path <- system.file("extdata", "qpcr_hsp12.csv", package = "wet2R")
#' result <- prepare_qpcr(path, control = "C-")
#' result$summary
#' result$technical_qc
#' @export
prepare_qpcr <- function(csv_file, control = "WT", group_order = NULL,
                         qc_range_max = NULL) {
  required <- c("group", "bio_rep", "tech_rep", "target_ct", "reference_ct")
  if (is.data.frame(csv_file)) {
    raw <- as.data.frame(csv_file)
  } else {
    qpcr_string(csv_file, "csv_file")
    raw <- utils::read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE,
                           colClasses = "character")
    for (column in intersect(c("target_ct", "reference_ct"), names(raw))) {
      raw[[column]] <- utils::type.convert(raw[[column]], as.is = TRUE)
    }
  }
  if (anyDuplicated(names(raw)) || any(!nzchar(names(raw)))) {
    stop("Input column names must be nonempty and unique.")
  }
  if (!all(required %in% names(raw))) {
    stop("Input must contain: ", paste(required, collapse = ", "))
  }
  if (!nrow(raw) || anyNA(raw)) stop("Input must not be empty or contain missing values.")
  for (column in c("group", "bio_rep", "tech_rep")) {
    if (!is.atomic(raw[[column]])) stop("Invalid identifier column: ", column)
    raw[[column]] <- as.character(raw[[column]])
    if (any(!nzchar(trimws(raw[[column]]))) ||
        any(raw[[column]] != trimws(raw[[column]]))) {
      stop("Identifiers must be nonempty with no surrounding whitespace: ", column)
    }
  }
  if (anyDuplicated(raw[c("group", "bio_rep", "tech_rep")])) {
    stop("group + bio_rep + tech_rep must be unique.")
  }
  for (column in c("target_ct", "reference_ct")) {
    if (!is.numeric(raw[[column]]) ||
        any(!is.finite(raw[[column]]) | raw[[column]] <= 0)) {
      stop(column, " must contain finite positive numeric Ct values.")
    }
  }
  for (column in intersect(c("gene", "target_gene", "reference_gene"), names(raw))) {
    if (length(unique(raw[[column]])) != 1L) stop("Use one assay per input: ", column)
  }
  metadata <- setdiff(names(raw), required)
  if (any(metadata %in% c("delta_ct", "delta_delta_ct", "relative_expression"))) {
    stop("Input contains reserved calculated column names.")
  }
  if (!is.null(qc_range_max)) qpcr_positive(qc_range_max, "qc_range_max")
  qpcr_string(control, "control")
  groups <- unique(raw$group)
  if (!control %in% groups) stop("control must match an input group.")
  if (is.null(group_order)) group_order <- c(control, groups[groups != control])
  if (!is.character(group_order) || anyNA(group_order) ||
      anyDuplicated(group_order) || !setequal(group_order, groups)) {
    stop("group_order must contain every group name exactly once.")
  }

  samples <- split(raw, interaction(match(raw$group, groups),
                                    match(raw$bio_rep, unique(raw$bio_rep)), drop = TRUE))
  bio <- do.call(rbind, lapply(samples, function(d) {
    for (column in metadata) {
      if (!is.atomic(d[[column]]) || length(unique(d[[column]])) != 1L) {
        stop("Metadata must be constant within each biological sample: ", column)
      }
    }
    out <- d[1L, c("group", "bio_rep", metadata), drop = FALSE]
    out$target_ct <- mean(d$target_ct)
    out$reference_ct <- mean(d$reference_ct)
    out
  }))
  if (any(table(bio$group) < 2L)) stop("Each group needs at least 2 biological replicates.")
  bio <- bio[order(match(bio$group, group_order), bio$bio_rep), , drop = FALSE]
  rownames(bio) <- NULL
  bio$delta_ct <- bio$target_ct - bio$reference_ct
  bio$delta_delta_ct <- bio$delta_ct - mean(bio$delta_ct[bio$group == control])
  bio$relative_expression <- 2^(-bio$delta_delta_ct)

  summary <- do.call(rbind, lapply(group_order, function(group) {
    d <- bio[bio$group == group, , drop = FALSE]
    center <- mean(d$delta_delta_ct)
    spread <- stats::sd(d$delta_ct)
    data.frame(group = group, n_bio = nrow(d),
               mean_target_ct = mean(d$target_ct),
               mean_reference_ct = mean(d$reference_ct),
               mean_delta_ct = mean(d$delta_ct), mean_delta_delta_ct = center,
               sd_delta_ct = spread, geometric_mean = 2^(-center),
               error_lower = 2^(-center - spread), error_upper = 2^(-center + spread))
  }))
  qc <- do.call(rbind, lapply(samples, function(d) {
    target_range <- diff(range(d$target_ct))
    reference_range <- diff(range(d$reference_ct))
    reasons <- c(if (nrow(d) < 2L) "fewer than two technical replicates",
                 if (!is.null(qc_range_max) && target_range > qc_range_max) "target Ct range",
                 if (!is.null(qc_range_max) && reference_range > qc_range_max) "reference Ct range")
    data.frame(group = d$group[1L], bio_rep = d$bio_rep[1L],
               n_technical = nrow(d),
               target_ct_sd = stats::sd(d$target_ct), reference_ct_sd = stats::sd(d$reference_ct),
               target_ct_range = target_range, reference_ct_range = reference_range,
               qc_flag = length(reasons) > 0L, qc_reason = paste(reasons, collapse = "; "))
  }))
  qc <- qc[order(match(qc$group, group_order), qc$bio_rep), , drop = FALSE]
  rownames(qc) <- NULL
  transformed <- c(bio$relative_expression, summary$geometric_mean,
                   summary$error_lower, summary$error_upper)
  if (any(!is.finite(transformed) | transformed <= 0)) {
    stop("Exponential transformation overflowed or underflowed; check Ct values.")
  }
  structure(list(raw = raw, biological_data = bio, summary = summary,
                 technical_qc = qc, control = control, group_order = group_order,
                 qc_range_max = qc_range_max), class = "wet2r_qpcr")
}

qpcr_string <- function(x, argument) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))) {
    stop(argument, " must be one nonempty character string.")
  }
}

qpcr_positive <- function(x, argument) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(argument, " must be one finite positive number.")
  }
}
