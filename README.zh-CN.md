# wet2R

[English](README.md) | 中文版

**wet2R** 整理了常用的湿实验（wet lab）操作 protocol，以及用于处理和可视化实验数据的脚本与 R 包。目的是让实验室成员在日常分子生物学实验中，数据分析步骤可以被重复和追溯。

## 仓库内容

| 路径 | 说明 |
|---|---|
| [`protocols/`](protocols/) | 逐步操作的实验 protocol（见下文） |
| [`R/`](R/) | `wet2R` R 包：qPCR 分析与绘图函数 |
| [`scripts/`](scripts/) | 可直接运行的分析脚本，包括完整的 qPCR 演示流程 |
| [`data/`](data/) | 小型示例数据（大体积原始数据不要提交） |
| [`man/`](man/) | 由 roxygen2 生成的函数文档 |

### 实验 protocol

- **[RNA 提取（TRIzol 法）](protocols/RNA_extraction.md)**：TRIzol 裂解组织、氯仿分相、异丙醇沉淀 RNA、75% 乙醇洗涤。附录说明 TRIzol 各成分的作用、氯仿的功能，以及 rpm 与 ×g 的区别。
- **[离心机使用](protocols/centrifugation.md)**：转速单位的选择（rpm 与 ×g）、离心管放置、转子配平、运行前检查、运行安全，以及结束后的清理。
- **[RT-qPCR](protocols/RT-qPCR.md)**：RNA 定量与取样、去除基因组 DNA、反转录合成 cDNA、SYBR Green qPCR 体系配制，并包含对照和质控要求。
- **[小鼠基因型鉴定](protocols/genotyping.md)**：用碱裂解法从脚趾或尾部组织粗提基因组 DNA，进行基因型 PCR 扩增，并用琼脂糖凝胶电泳检测。

## 安装

从 GitHub 直接安装 R 包：

```r
# install.packages("remotes")
remotes::install_github("Hinna0818/wet2R")
```

或者克隆仓库后本地安装：

```bash
git clone https://github.com/Hinna0818/wet2R.git
cd wet2R
```

```r
install.packages(c("ggplot2", "ragg", "scales", "svglite", "systemfonts"))
devtools::install(".")
```

## 快速上手：qPCR 分析

包读取的 CSV 文件每行对应一个技术重复，需要包含以下列：

| 列名 | 含义 |
|---|---|
| `group` | 分组名称，如 `WT`、`KO` |
| `bio_rep` | 生物学重复编号 |
| `tech_rep` | 技术重复编号 |
| `target_ct` | 目标基因的 Ct 值 |
| `reference_ct` | 内参（管家）基因的 Ct 值 |

相对表达量使用 2^-ddCt 方法计算。每个生物学重复先对其技术重复取平均，再以对照组作为基线。

```r
library(wet2R)

result <- prepare_qpcr("data/example_qpcr.csv", control = "WT")
result$summary          # 各组均值、标准差及误差范围
result$technical_qc     # 每个生物学重复内技术重复的 Ct 极差

p <- plot_qpcr(result, gene = "Il6", group_labels = c("WT", "KO"))
print(p)

save_qpcr(result, p, output_prefix = "results/Il6")
# 输出 Il6.png（600 dpi）、Il6.svg 以及三个 CSV 表格
```

完整流程（包括模拟数据）见 [`scripts/qpcr_demo.r`](scripts/qpcr_demo.r)。在仓库根目录运行：

```bash
Rscript scripts/qpcr_demo.r
```

### 分析说明

- 每组至少需要 2 个生物学重复才能计算标准差。
- Ct 值必须为有限的正数。未检出的反应需要先核实再分析。
- 柱高为 2^-ddCt 的几何均值；误差线由 dCt 的标准差换算得到。
- 分组数不超过 4 时，组名直接标在 x 轴上；超过 4 组时显示图例。
- 默认配色已做色觉缺陷检查。如需更换，可通过 `colors` 参数指定。

## 引用

如果你在工作中使用了 wet2R，请引用本项目。安装包后可以运行以下命令获取引用格式：

```r
citation("wet2R")
```

BibTeX：

```bibtex
@Misc{wet2R,
  title  = {wet2R: Wet-lab protocols and qPCR analysis utilities},
  author = {Nan He},
  year   = {2026},
  note   = {R package version 0.1.0},
  url    = {https://github.com/Hinna0818/wet2R},
}
```

## 作者

**Nan He**
邮箱：[hinna01@163.com](mailto:hinna01@163.com)
GitHub：[@Hinna0818](https://github.com/Hinna0818)

## 贡献与数据规范

- 不要提交敏感数据、大体积原始文件或个人信息。
- 修改 protocol 或脚本时，请在提交信息中简要说明改动内容。
- protocol 反映的是编写时的实验室操作方式。使用前请核对试剂盒说明书和所在机构的安全规定。

## 许可证

R 包以 MIT 许可证发布，详见 [LICENSE](LICENSE)。
