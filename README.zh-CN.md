# wet2R

[English](README.md) | 中文版

**wet2R** 整理了常用的湿实验（wet lab）操作 protocol，以及用于处理和可视化实验数据的脚本与 R 包。目的是让实验室成员在日常分子生物学实验中，数据分析步骤可以被重复和追溯。

## 仓库内容

| 路径 | 说明 |
|---|---|
| [`protocols/`](protocols/) | 逐步操作的实验 protocol（见下文） |
| [`R/`](R/) | `wet2R` R 包：qPCR 分析与绘图函数 |
| [`scripts/`](scripts/) | 可直接运行的分析脚本，包括完整的 qPCR 演示流程 |
| [`inst/extdata/`](inst/extdata/) | 实测示例数据、来源说明与上游许可证 |
| [`data-raw/`](data-raw/) | 公开原始数据及可复现的整理脚本 |
| [`man/`](man/) | 由 roxygen2 生成的函数文档 |
| [`tests/testthat/`](tests/testthat/) | 使用实测数据的回归检查 |

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
install.packages("remotes")
remotes::install_local(".")
```

## qPCR 快速上手

安装后，直接运行以下示例，读取随包提供的实测数据：

```r
library(wet2R)

path <- system.file("extdata", "qpcr_hsp12.csv", package = "wet2R", mustWork = TRUE)
result <- prepare_qpcr(path, control = "C-")
tests <- test_qpcr(result)
p <- plot_qpcr(result, gene = "HSP12", tests = tests)
p
```

默认计算 2^-ddCt 相对表达量，在生物学重复的 dCt 上使用 Welch t 检验，将各处理组与 `C-` 比较，并采用 Holm 校正。图中使用不同深浅的灰色柱形，星号表示校正后的显著性。`p` 为 ggplot 对象；汇总、技术质控和检验结果分别保存在 `result$summary`、`result$technical_qc` 和 `tests` 中。

HSP12/PGK1 示例包含三个未热休克条件：对照（`C-`）、phenanthroline（`P-`）和 thiolutin（`T-`）。每组两个生物学重复，每个样本三个技术重复。来源见[数据说明](inst/extdata/SOURCE.md)，完整调用见[示例脚本](scripts/qpcr_demo.r)。解释这个预实验的结果前，需要检查技术变异和样本量。

### 使用自己的数据

将 `path` 替换为自己的 CSV 路径，`control` 替换为对照组名称，`gene` 替换为目标基因名称。输入也支持 data.frame，每行对应一个技术重复，包含以下必需列：

| 列名 | 含义 |
|---|---|
| `group` | 分组名称，如 `WT`、`KO` |
| `bio_rep` | 生物学重复编号 |
| `tech_rep` | 技术重复编号 |
| `target_ct` | 目标基因的 Ct 值 |
| `reference_ct` | 内参基因的 Ct 值 |

每次输入包含一个目标基因和一个内参。额外的样本信息列会保留，在同一个 `group`/`bio_rep` 内必须一致。目标基因和内参需要具有相同数量的有效技术测量；未检出反应需要在分析前核实处理。

### 分析说明

- 每组至少需要 2 个生物学重复才能计算标准差。
- Ct 值必须为有限的正数。未检出的反应需要先核实再分析。
- 使用前确认内参表达稳定、目标基因与内参扩增效率相近且接近每轮倍增。
- 柱高为 2^-ddCt 的几何均值；误差线将平均 ddCt ± SD(dCt) 转换到表达量尺度，对照均值视为已观测基线。误差线表示样本离散程度，处理效应置信区间另见检验结果。
- `qc_range_max` 可指定技术重复 Ct 极差阈值。质控只标记、保留原因，不自动删除数据。只有一个技术重复时，SD 为 NA，并标记质控信息不足。
- 检验使用生物学重复的 dCt。`method` 支持 `welch`、`student`、`wilcox`、`paired_t`、`paired_wilcox`；根据实验设计选择，Student t 检验要求等方差。
- 配对检验必须提供 `pair_by`，并具有完整且唯一的匹配关系。相同的重复编号或实验板编号本身不能证明生物学配对。
- `alternative` 指 dCt 的方向：`greater` 表示 group1 的 dCt 更高、表达量更低。
- `p_adjust_method` 接受 `stats::p.adjust.methods` 中的方法，默认 `holm`。BH 在相应依赖假设下控制 FDR。同一次调用的全部比较构成校正范围；跨基因校正需统一汇总处理。
- 返回原始 P 值、校正 P 值及描述性表达倍数。t 检验的置信区间为逐个比较的未校正区间；Wilcoxon 不提供均值差的置信区间及自由度。
- 显著性标注使用校正 P 值，默认 `significance_label = "stars"`，也支持 `"p_adj"`。筛选检验表中的部分行用于绘图时，保留原校正范围。
- 默认布局使用居中的基因标题、由深到浅的不同灰色柱形、黑色实心散点和倾斜 45 度的组名。统计方法和误差线含义见检验表及帮助文档。
- 命名配色和标签按组名匹配。可调整 `base_size`、`bar_width`、`point_size`、`label_angle`、`show_legend`，也可继续添加 ggplot2 图层。默认字体为 `sans`，显式指定未安装字体会报错。
- 因素交互作用、批次调整和超过两个条件的重复测量需要另行建立适合实验设计的统计模型。

### R 包开发与检查

```r
install.packages(c("roxygen2", "pkgload", "testthat", "readr"))
pkgload::load_all(".")
testthat::test_local(".")
roxygen2::roxygenise(".")
```

```bash
R CMD build .
R CMD check --no-manual wet2R_0.2.0.tar.gz
```

帮助示例通过 `system.file()` 读取安装包中的数据。测试核对计算结果、样本分组、配对匹配、P 值校正、颜色对应及图形构建。原始数据保存在 `data-raw/`，执行 `Rscript data-raw/prepare-example.R` 可重新生成示例 CSV。0.2.0 的接口变化见 [NEWS.md](NEWS.md)。

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
  note   = {R package version 0.2.0},
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
