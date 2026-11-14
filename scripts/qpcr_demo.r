# qPCR analysis demo
#
# Workflow: simulate data -> read and calculate -> plot -> export
# Run from the package root:  Rscript examples/qpcr_demo.r
# For real data, set example_path to your own CSV file.

pkgload::load_all(".")

# ---- 1. Prepare example data (skip for real analysis) -----------------------
example_path <- "data/example_qpcr.csv"
dir.create(dirname(example_path), recursive = TRUE, showWarnings = FALSE)

set.seed(42)
n_bio  <- 3   # biological replicates per group
n_tech <- 3   # technical replicates per biological replicate
sim <- expand.grid(group = c("WT", "KO"), bio_rep = 1:n_bio,
                   tech_rep = 1:n_tech, stringsAsFactors = FALSE)
# Each biological replicate shares a base Ct; technical replicates get small noise
bio_key <- paste(sim$group, sim$bio_rep)
bio_target <- setNames(rnorm(length(unique(bio_key)), 0, 0.3), unique(bio_key))
bio_ref    <- setNames(rnorm(length(unique(bio_key)), 0, 0.2), unique(bio_key))
sim$target_ct <- round(24 + ifelse(sim$group == "KO", 2.5, 0) +
                         bio_target[bio_key] + rnorm(nrow(sim), 0, 0.1), 2)
sim$reference_ct <- round(18 + bio_ref[bio_key] + rnorm(nrow(sim), 0, 0.1), 2)
write.csv(sim, example_path, row.names = FALSE)

# ---- 2. Read data and calculate relative expression -------------------------
result <- prepare_qpcr(example_path, control = "WT", group_order = c("WT", "KO"))

cat("== Group summary ==\n")
print(result$summary[, c("group", "n_bio", "geometric_mean",
                         "error_lower", "error_upper")], row.names = FALSE)

cat("\n== Technical QC (Ct range) ==\n")
print(result$technical_qc, row.names = FALSE)

# ---- 3. Plot ------------------------------------------------------------------
p <- plot_qpcr(result, gene = "Il6", group_labels = c("WT", "KO"))
dir.create("results", showWarnings = FALSE)
# Preview with the ragg device: R's default PDF device does not support
# the Arial font selected through systemfonts
if (interactive()) {
  print(p)
} else {
  ragg::agg_png("results/Il6_preview.png", width = 3.4, height = 3.5,
                units = "in", res = 150)
  print(p)
  invisible(grDevices::dev.off())
}

# ---- 4. Export figures and tables -------------------------------------------
save_qpcr(result, p, output_prefix = "results/Il6")
cat("\nResults saved to results/\n")
