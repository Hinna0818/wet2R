# HSP12 / PGK1 measured qPCR example

Measurements: Edward Wallace, 15 June 2018, Saccharomyces cerevisiae.
Source: rOpenSci tidyqpcr, `multifactor_vignette` and its accompanying
LightCycler exports. Retrieved 20 November 2026.

Experiment description:
https://docs.ropensci.org/tidyqpcr/articles/multifactor_vignette.html

Original files (preserved in this repository's `data-raw/`):

- https://raw.githubusercontent.com/ropensci/tidyqpcr/main/inst/extdata/Edward_qPCR_TxnInhibitors_HS_2018-06-15_plate1_Cq.txt.gz
- https://raw.githubusercontent.com/ropensci/tidyqpcr/main/inst/extdata/Edward_qPCR_TxnInhibitors_HS_2018-06-15_plate2_Cq.txt.gz

The source describes two biological replicates, one per plate, and three
technical replicates. `qpcr_hsp12.csv` contains the complete no-heat subset
for HSP12 (row B) and PGK1 (row H): columns 1–3, 7–9, 13–15 on each plate.
C- is no drug; P- is phenanthroline; T- is thiolutin. Columns 1–3, 7–9,
and 13–15 are technical repeats 1, 2, and 3 respectively. Target/reference
wells are put on the same row by plate column for storage, then averaged
separately. Ct values are unchanged; no selected wells are imputed or removed.

Run `Rscript data-raw/prepare-example.R` from the repository root to reproduce
the CSV using readr. Raw files and conversion code are excluded from the
installed package; the CSV and this attribution are included.

This is a pilot dataset for software examples and regression tests. Each
group has only two biological replicates. PGK1 technical variation should be
reviewed in `technical_qc`; inclusion here does not validate reference stability,
amplification efficiency, independence, or a biological treatment conclusion.
The plate identifier describes the measurement batch; it is not evidence of
matched individuals. Paired-method tests verify identifier matching numerically.

Upstream license: Apache License 2.0; see `tidyqpcr-LICENSE.md` alongside this
file. This subset is distributed with attribution and the upstream license.
