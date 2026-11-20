# 安装 wet2R 后，在任意工作目录运行本示例。
library(wet2R)

path <- system.file("extdata", "qpcr_hsp12.csv", package = "wet2R", mustWork = TRUE)
result <- prepare_qpcr(path, control = "C-")
tests <- test_qpcr(result, comparisons = list(c("P-", "C-"), c("T-", "C-")),
                   method = "welch", p_adjust_method = "holm")
print(result$technical_qc)
print(tests)
p <- plot_qpcr(result, gene = "HSP12", tests = tests)
print(p)

# 按需保存：ggplot2::ggsave("HSP12.pdf", plot = p, width = 6, height = 5)。
