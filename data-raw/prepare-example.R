# 原始实验与板位说明见 inst/extdata/SOURCE.md；从仓库根目录执行。
library(readr)

example_data <- do.call(rbind, lapply(1:2, function(plate) {
  raw <- read_tsv(paste0("data-raw/plate", plate, "_Cq.txt.gz"), skip = 1,
                  col_types = cols(.default = col_character(), Cp = col_double()))
  stop_for_problems(raw)
  columns <- c(1:3, 7:9, 13:15)
  target <- raw$Cp[match(paste0("B", columns), raw$Pos)]
  reference <- raw$Cp[match(paste0("H", columns), raw$Pos)]
  stopifnot(length(target) == 9L, !anyNA(target), !anyNA(reference))
  data.frame(group = rep(c("C-", "P-", "T-"), 3),
             bio_rep = as.character(plate), tech_rep = rep(1:3, each = 3),
             target_ct = target, reference_ct = reference,
             target_gene = "HSP12", reference_gene = "PGK1", plate = as.character(plate))
}))
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
write.csv(example_data, "inst/extdata/qpcr_hsp12.csv", row.names = FALSE)
