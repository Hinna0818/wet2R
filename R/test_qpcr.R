#' Compare biological-replicate delta Ct values
#'
#' Perform prespecified pairwise tests and adjust all P values in one family.
#'
#' @param result A result from [prepare_qpcr()].
#' @param method One of `"welch"` (default), `"student"`, `"wilcox"`,
#'   `"paired_t"`, or `"paired_wilcox"`. Student's test assumes equal variances.
#' @param comparisons List of two-element character vectors, each ordered
#'   `c(group1, group2)`. `NULL` compares every treatment against the control.
#' @param p_adjust_method A method in [stats::p.adjust.methods], default `"holm"`.
#' @param pair_by Column in `biological_data` identifying matched samples.
#'   Required for paired tests; must be `NULL` for independent tests.
#' @param alternative `"two.sided"`, `"less"`, or `"greater"`, on the dCt scale:
#'   `"greater"` means group1 has higher dCt (lower relative expression).
#' @param conf_level Confidence level for t-test intervals, strictly between 0 and 1.
#'
#' @details
#' The unit of inference is the biological replicate. Pairing is explicit and
#' matched by identifier, never by row order. Incomplete or duplicated pairs
#' cause an error. Confirm pairing from the experimental design.
#'
#' Correction applies to all requested comparisons before plotting. Holm
#' controls family-wise error; BH controls false discovery rate under its
#' dependence assumptions. Separate calls are separate correction families;
#' coordinate cross-gene correction outside this single-target function.
#' No normality-based automatic test selection is performed. Batch effects,
#' factorial interactions and repeated measures beyond pairs require a model
#' appropriate to the experiment.
#'
#' Wilcoxon tests use [stats::wilcox.test()] defaults for exact/approximate
#' P values and continuity correction; warnings about ties are retained.
#' They test rank/location hypotheses. Descriptive mean dCt differences and
#' fold changes remain available, but their confidence intervals and degrees
#' of freedom are `NA` for Wilcoxon tests.
#'
#' @return A data frame with groups, sample sizes, method, statistic, `df`,
#'   `p_value`, `p_adj`, `p_adjust_method`, `n_comparisons`, `alternative`,
#'   `pair_by`, `mean_delta_ct_difference` (group1 minus group2), `log2_fold_change`,
#'   `fold_change` (group1/group2), `conf_level`, `delta_ct_low`, `delta_ct_high`,
#'   `fold_change_low`, and `fold_change_high`. T-test intervals are pointwise,
#'   unadjusted intervals, even when P values are adjusted. Source biological
#'   data are attached for consistency checks by [plot_qpcr()].
#' @examples
#' path <- system.file("extdata", "qpcr_hsp12.csv", package = "wet2R")
#' result <- prepare_qpcr(path, control = "C-")
#' test_qpcr(result, comparisons = list(c("P-", "C-"), c("T-", "C-")))
#' @export
test_qpcr <- function(result,
                      method = c("welch", "student", "wilcox", "paired_t", "paired_wilcox"),
                      comparisons = NULL, p_adjust_method = "holm", pair_by = NULL,
                      alternative = c("two.sided", "less", "greater"), conf_level = 0.95) {
  if (!inherits(result, "wet2r_qpcr")) stop("result must come from prepare_qpcr().")
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)
  qpcr_positive(conf_level, "conf_level")
  if (conf_level >= 1) stop("conf_level must be less than 1.")
  paired <- method %in% c("paired_t", "paired_wilcox")
  rank_test <- method %in% c("wilcox", "paired_wilcox")
  d <- result$biological_data
  if (paired) {
    qpcr_string(pair_by, "pair_by")
    if (!pair_by %in% names(d)) stop("pair_by must name a biological_data column.")
    if (anyNA(d[[pair_by]]) || any(!nzchar(trimws(as.character(d[[pair_by]]))))) {
      stop("Pair identifiers must be complete and nonempty.")
    }
  } else if (!is.null(pair_by)) {
    stop("pair_by is only valid for paired tests.")
  }
  if (is.null(comparisons)) {
    comparisons <- lapply(setdiff(result$group_order, result$control),
                          function(group) c(group, result$control))
  }
  if (!is.list(comparisons) || !length(comparisons) ||
      !all(vapply(comparisons, function(pair) {
        is.character(pair) && length(pair) == 2L && !anyNA(pair) &&
          !anyDuplicated(pair) && all(pair %in% result$group_order)
      }, logical(1)))) {
    stop("comparisons must be a nonempty list of two distinct existing group names.")
  }
  if (anyDuplicated(lapply(comparisons, sort))) stop("Duplicate or reversed comparisons found.")

  tests <- do.call(rbind, lapply(comparisons, function(pair) {
    a <- d[d$group == pair[1L], , drop = FALSE]
    b <- d[d$group == pair[2L], , drop = FALSE]
    if (paired) {
      if (anyDuplicated(a[[pair_by]]) || anyDuplicated(b[[pair_by]]) ||
          !setequal(a[[pair_by]], b[[pair_by]])) {
        stop("Pairs must be unique and complete for ", paste(pair, collapse = " vs "), ".")
      }
      b <- b[match(a[[pair_by]], b[[pair_by]]), , drop = FALSE]
    }
    x <- a$delta_ct
    y <- b$delta_ct
    if (length(x) < 2L || length(y) < 2L || any(!is.finite(c(x, y)))) {
      stop("Each comparison needs at least two finite biological observations per group.")
    }
    if (rank_test) {
      fit <- stats::wilcox.test(x, y, paired = paired, alternative = alternative)
      ci <- c(NA_real_, NA_real_)
      df <- NA_real_
    } else {
      fit <- stats::t.test(x, y, paired = paired, var.equal = method == "student",
                           alternative = alternative, conf.level = conf_level)
      ci <- unname(fit$conf.int)
      df <- unname(fit$parameter)
    }
    if (!is.finite(fit$p.value)) stop("Test produced an undefined P value for ", paste(pair, collapse = " vs "))
    difference <- mean(x) - mean(y)
    data.frame(group1 = pair[1L], group2 = pair[2L], n1 = length(x), n2 = length(y),
               method = method, statistic = unname(fit$statistic), df = df,
               p_value = fit$p.value, alternative = alternative,
               pair_by = if (paired) pair_by else NA_character_,
               mean_delta_ct_difference = difference, log2_fold_change = -difference,
               fold_change = 2^(-difference), conf_level = if (rank_test) NA_real_ else conf_level,
               delta_ct_low = ci[1L], delta_ct_high = ci[2L],
               fold_change_low = 2^(-ci[2L]), fold_change_high = 2^(-ci[1L]))
  }))
  tests$p_adj <- stats::p.adjust(tests$p_value, method = p_adjust_method)
  tests$p_adjust_method <- p_adjust_method
  tests$n_comparisons <- nrow(tests)
  attr(tests, "biological_data") <- d
  tests
}
