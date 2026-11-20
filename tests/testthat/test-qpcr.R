test_that("measured Ct values retain biological sample size and metadata", {
  raw <- qpcr_measured_data()
  result <- prepare_qpcr(raw, control = "C-", qc_range_max = 0.5)
  expect_s3_class(result, "wet2r_qpcr")
  expect_equal(result$raw, raw)
  csv_result <- prepare_qpcr(system.file("extdata", "qpcr_hsp12.csv", package = "wet2R"), "C-")
  expect_type(csv_result$raw$plate, "character")
  expect_type(csv_result$raw$bio_rep, "character")
  expect_equal(csv_result$summary, result$summary)
  expect_equal(nrow(result$biological_data), 6L)
  expect_equal(result$summary$n_bio, c(2L, 2L, 2L))
  expect_equal(result$summary$mean_delta_ct, c(0.533333333333333, 3.786666666666667, 2.871666666666667))
  expect_equal(result$summary$geometric_mean[1], 1)
  expect_equal(result$biological_data$plate, rep(1:2, 3))
  expect_equal(result$technical_qc$n_technical, rep(3L, 6))
  expect_true(any(result$technical_qc$qc_flag))
  expect_true(all(nzchar(result$technical_qc$qc_reason[result$technical_qc$qc_flag])))
  expect_true(all(result$summary$error_lower < result$summary$geometric_mean))
  expect_true(all(result$summary$error_upper > result$summary$geometric_mean))

  reordered <- prepare_qpcr(raw[nrow(raw):1, ], "C-")
  expect_equal(reordered$biological_data[order(reordered$biological_data$group, reordered$biological_data$bio_rep), ],
               result$biological_data[order(result$biological_data$group, result$biological_data$bio_rep), ],
               ignore_attr = TRUE)
})

test_that("incomplete and ambiguous input stops before aggregation", {
  raw <- qpcr_measured_data()
  expect_error(prepare_qpcr(raw[, -5], "C-"), "Input must contain")
  expect_error(prepare_qpcr(raw[FALSE, ], "C-"), "empty")
  expect_error(prepare_qpcr(rbind(raw, raw[1, ]), "C-"), "unique")
  expect_error(prepare_qpcr(raw[raw$bio_rep == "1", ], "C-"), "2 biological")
  expect_error(prepare_qpcr(raw, "unknown"), "match an input group")
  expect_error(prepare_qpcr(raw, "C-", group_order = c("C-", "P-")), "every group")
  expect_error(prepare_qpcr(raw, "C-", qc_range_max = 0), "positive")
  duplicate_names <- raw
  names(duplicate_names)[ncol(duplicate_names)] <- "group"
  expect_error(prepare_qpcr(duplicate_names, "C-"), "unique")
  one_well <- prepare_qpcr(raw[raw$tech_rep == "1", ], "C-")
  expect_true(all(one_well$technical_qc$qc_flag))
  expect_true(all(is.na(one_well$technical_qc$target_ct_sd)))
})

test_that("independent methods agree with R stats on biological dCt", {
  result <- prepare_qpcr(qpcr_measured_data(), "C-")
  d <- result$biological_data
  x <- d$delta_ct[d$group == "P-"]
  y <- d$delta_ct[d$group == "C-"]
  for (method in c("welch", "student", "wilcox")) {
    tests <- test_qpcr(result, method = method)
    expected <- if (method == "wilcox") stats::wilcox.test(x, y) else
      stats::t.test(x, y, var.equal = method == "student")
    expect_equal(tests$p_value[1], expected$p.value)
    expect_equal(tests$n1, c(2L, 2L))
    expect_equal(tests$p_adj, stats::p.adjust(tests$p_value, "holm"))
    expect_equal(tests$fold_change[1], 2^(-3.253333333333333))
    if (method == "wilcox") {
      expect_true(all(is.na(tests$fold_change_low)))
    } else {
      expect_equal(tests$fold_change_low[1], unname(2^(-expected$conf.int[2])))
      expect_equal(tests$fold_change_high[1], unname(2^(-expected$conf.int[1])))
    }
  }
  for (adjustment in stats::p.adjust.methods) {
    tests <- test_qpcr(result, p_adjust_method = adjustment)
    expect_equal(tests$p_adj, stats::p.adjust(tests$p_value, adjustment))
  }
  greater <- test_qpcr(result, alternative = "greater")
  expect_equal(greater$p_value[1], stats::t.test(x, y, alternative = "greater")$p.value)
  expect_true(is.infinite(greater$delta_ct_high[1]))
  expect_equal(greater$fold_change_low[1], 0)
  expect_error(test_qpcr(result, comparisons = list(c("P-", "C-"), c("C-", "P-"))), "Duplicate")
  expect_error(test_qpcr(result, comparisons = list(c("C-", "C-"))), "distinct")
  expect_error(test_qpcr(result, comparisons = list()), "nonempty")
  expect_error(test_qpcr(result, method = "anova"), "arg")
  expect_error(test_qpcr(result, p_adjust_method = "invalid"), "arg")
  expect_error(test_qpcr(result, conf_level = 1), "less than 1")
})

test_that("paired methods match explicit identifiers and reject ambiguous pairs", {
  result <- prepare_qpcr(qpcr_measured_data(), "C-")
  # plate 仅用于验证匹配算法；实验配对依据由分析者确认。
  result$biological_data <- result$biological_data[c(1, 2, 4, 3, 5, 6), ]
  d <- result$biological_data
  a <- d[d$group == "P-", ]
  b <- d[d$group == "C-", ]
  b <- b[match(a$plate, b$plate), ]
  for (method in c("paired_t", "paired_wilcox")) {
    tests <- test_qpcr(result, method = method, pair_by = "plate")
    expected <- if (method == "paired_t") stats::t.test(a$delta_ct, b$delta_ct, paired = TRUE) else
      stats::wilcox.test(a$delta_ct, b$delta_ct, paired = TRUE)
    expect_equal(tests$p_value[1], expected$p.value)
  }
  expect_error(test_qpcr(result, method = "paired_t"), "pair_by")
  expect_error(test_qpcr(result, pair_by = "plate"), "only valid")
  expect_error(test_qpcr(result, method = "paired_t", pair_by = "group"), "unique and complete")
  expect_error(test_qpcr(result, method = "paired_t", pair_by = "target_ct"), "unique and complete")
})

test_that("plots preserve named colors, comparisons and all biological points", {
  raw <- qpcr_measured_data()
  result <- prepare_qpcr(raw, "C-", group_order = c("C-", "T-", "P-"))
  colors <- c("P-" = "red", "C-" = "grey", "T-" = "blue")
  labels <- c("T-" = "Thiolutin", "P-" = "Phenanthroline", "C-" = "Control")
  tests <- test_qpcr(result)
  p <- plot_qpcr(result, "HSP12", colors = colors, group_labels = labels,
                  tests = tests[1, , drop = FALSE], show_legend = TRUE)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  fill <- built$plot$scales$get_scales("fill")
  expect_equal(fill$get_breaks(), result$group_order, ignore_attr = TRUE)
  expect_equal(fill$get_labels(), unname(labels[result$group_order]))
  expect_equal(fill$map(result$group_order), unname(colors[result$group_order]))
  expect_equal(nrow(built$data[[3]]), 6L)
  expect_equal(built$data[[3]]$y, result$biological_data$relative_expression)
  expect_equal(built$data[[3]]$x, ggplot2::ggplot_build(p)$data[[3]]$x)
  expect_equal(nrow(built$data[[4]]), 3L)
  expect_true(all(is.finite(built$data[[5]]$y)))
  expect_null(p$labels$caption)
  expect_equal(p$labels$title, "HSP12")
  expect_equal(p$labels$y, "Relative expression")
  expect_equal(tests$n_comparisons, c(2L, 2L))
  expect_equal(tests$p_adj, stats::p.adjust(tests$p_value, "holm"))
  expect_no_warning(ggplot2::ggplotGrob(p))
  stars <- plot_qpcr(result, "HSP12", tests = tests, significance_label = "stars")
  expect_equal(ggplot2::ggplot_build(stars)$data[[5]]$label, c("ns", "ns"))
  unadjusted <- plot_qpcr(result, "HSP12", tests = test_qpcr(result, p_adjust_method = "none"),
                          significance_label = "p_adj")
  expect_true(all(startsWith(ggplot2::ggplot_build(unadjusted)$data[[5]]$label, "P = ")))
  expect_error(plot_qpcr(result, "HSP12", colors = c("red", "blue")), "per group")
  expect_error(plot_qpcr(result, "HSP12", font = "wet2R_missing_font"), "not installed")
  expect_error(plot_qpcr(result, "HSP12", tests = test_qpcr(prepare_qpcr(raw[raw$tech_rep == "1", ], "C-"))),
                "these biological data")
})
