# 项目说明

## 定位

这是一个**正常房腔 bulk 基线分析**，目标是建立健康 `LA vs RA` 的 donor 配对参考。

这一步只回答两件事：

- 哪些基因在正常左房和右房之间存在稳定差异
- 哪些程序在正常左房和右房之间存在稳定偏向

这一步**不是**：

- 机制证明
- 细胞来源判断
- AF rewiring 结论

## 分析原则

- 统计单位是 `donor`，不是单个样本点
- 比较方式是 donor 内 `LA vs RA` 配对比较
- 主分析依赖 `~ donor + chamber`
- PCA/QC 只作为背景检查，不作为主结论依据

## 输入文件

- `GSE112339/GSE112339_HTSeq_hg19_Ensembl.txt.gz`
- `GSE112339/metadata.xlsx`
- `gene_sets/h.all.v2026.1.Hs.symbols.gmt`
- `gene_sets/c2.cp.reactome.v2026.1.Hs.symbols.gmt`

## 运行方式

如果终端里没有 `Rscript`，直接用：

```powershell
.\run_normal_chamber_analysis.ps1 `
  --count=GSE112339/GSE112339_HTSeq_hg19_Ensembl.txt.gz `
  --metadata=GSE112339/metadata.xlsx `
  --hallmark_gmt=gene_sets/h.all.v2026.1.Hs.symbols.gmt `
  --reactome_gmt=gene_sets/c2.cp.reactome.v2026.1.Hs.symbols.gmt `
  --out_dir=.
```

## 脚本结构

- `normal_chamber_analysis.R`
  主分析：metadata、配对差异、FGSEA、leave-one-donor-out
- `normal_chamber_qc.R`
  QC 和作图：PCA、MDS、mean-variance、Figure 1
- `normal_chamber_export.R`
  结果整理：full/core/core_strict/program/prior

## 主分析主线

1. 用 `~ donor + chamber` 做 gene-level 差异分析
2. 用完整排序结果做 pre-ranked FGSEA
3. 做 leave-one-donor-out，检查 gene/program 是否稳定

## 关键输出

### 1. 全量基因表

- `normal_chamber_genes.tsv`
- `normal_chamber_genes_full.tsv`

作用：

- 保留完整 gene-level 结果
- 供追溯和后续筛选使用

### 2. 核心基因 prior

- `normal_chamber_genes_core.tsv`

筛选规则：

- `stable_gene == TRUE`
- `FDR < 0.05`
- `abs(logFC) >= 0.5`

作用：

- 作为 baseline chamber-associated gene prior

### 3. 严格核心基因 prior

- `normal_chamber_genes_core_strict.tsv`

筛选规则：

- 满足 `normal_chamber_genes_core.tsv` 条件
- 额外要求 `is_core_gene == TRUE`

作用：

- 作为更严格的 program-anchored 基因集合
- 适合做 highlight 或后续重点展示

### 4. 程序表

- `normal_chamber_programs.tsv`

包含：

- `program_id`
- `program`
- `program_class`
- `direction`
- `NES`
- `FDR`
- `top_core_genes`
- `leading_edge_genes`
- `stable_program`

作用：

- 作为 baseline chamber-associated program prior

### 5. 整理版派生表

- `bulk_chamber_DE_full.tsv`
- `bulk_chamber_stability.tsv`
- `bulk_chamber_core_genes.tsv`
- `bulk_chamber_fgsea.tsv`
- `NormalChamberPrior.tsv`

这些文件主要用于整理和对外输出，不是主分析逻辑本身。

## QC 输出

- `qc/qc_pca.pdf`
- `qc/qc_mds.pdf`
- `qc/qc_logcpm_mean_variance.pdf`
- `qc/qc_pca_chamber_lm.txt`
- `figures/figure1A_gse112339_pca_la_vs_ra.pdf`
- `figures/figure1B_gse112339_volcano_heatmap.pdf`

## 解释边界

这批结果只能解释为：

- 正常房腔 bulk 层面的基线差异
- chamber-associated genes / programs

不能直接解释为：

- 因果机制
- 细胞类型来源
- AF 中已经发生的 rewiring

后续如果进入 AF 分析，应把这些结果当作**基线参照系**，而不是最终结论。
