# wet2R 0.2.0

- `test_qpcr()` supports Welch, Student, Wilcoxon, and explicit paired tests
  on biological-replicate dCt, with configurable P-value adjustment.
- `prepare_qpcr()` accepts CSV paths or data frames, preserves sample metadata,
  validates single-assay inputs, and reports technical SD and optional QC flags.
- `plot_qpcr()` returns an editable ggplot with reproducible replicate positions,
  named group colors, optional adjusted-P annotations, and generic fonts.
- Use `ggplot2::ggsave()` and `utils::write.csv()` for exports. The `save_qpcr()`
  export is removed; update existing callers accordingly.
- Prepared results now have class `wet2r_qpcr`. Recreate previously saved lists
  with `prepare_qpcr()` before passing them to analysis or plotting functions.
- Installed examples use a documented subset of publicly available measured
  HSP12/PGK1 data via `system.file("extdata", ..., package = "wet2R")`.
