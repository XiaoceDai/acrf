规则

==============================
全局叙事依赖硬规则（Figure stage lock）
==============================

目的：
保证主文图逻辑严格循序渐进，避免“前图调用后图结果”造成叙事穿越。

硬规则 1（前向依赖）：
Figure n 只能使用以下来源：
1) 原始数据/元数据
2) 同图内当步计算结果
3) Figure 1..(n-1) 已冻结输出

硬规则 2（禁止后向依赖）：
Figure n 严禁调用 Figure (n+1..end) 的任何对象，包括但不限于：
- 后续图的结果表（results/*.tsv）
- 后续图脚本产物（中间文件/缓存）
- 后续图结论性阈值或筛选名单
- 从后续图“反提取”出来的 genes/modules/anchors

阶段锁定（主文）：
- Figure 1：只讲 healthy LA/RA baseline construction（不允许 AF/SR 结果注入）
- Figure 2：可使用 Figure 1 冻结对象；禁止使用 Figure 3+
- Figure 3：可使用 Figure 1-2；禁止使用 Figure 4+
- Figure 4：可使用 Figure 1-3；禁止使用 Figure 5+
- Figure 5：可使用 Figure 1-4；禁止使用 Figure 6+
- Figure 6：可使用 Figure 1-5

违规处理：
- 若某元素依赖未来图信息：必须后移到对应后续图，或下沉至 Supplementary；
- 不允许为满足当前图叙事而提前“借用”后续图对象。

执行检查（每次改图前自检）：
- 当前图使用的每个输入文件，是否来自 <= 当前图阶段？
- 当前图标题/图注是否出现未来图结论性措辞？
- 当前图脚本是否读取了未来图专用输出路径？
不做重聚类
不做重新整合
不做 composition 推断
不做 cell–cell communication
不做 trajectory
不做新的 signature
不做 fancy 亚群命名\

只提取必要列

打开 adata.obs.columns 后，只找这几类列，并统一命名：

donor
region
cell_type
modality 或 assay（如果有）
n_counts / n_genes（如果已有，也可保留）

最终你需要的 obs 最少只保留：

cell_id
donor
region
cell_type
modality

这一步的目标不是整理得全面，而是整理得足够做 pseudobulk。

3.2 先限制分析范围

只保留：

区域
LA
RA
细胞类型

第一轮只做 4 类主细胞：

atrial_cardiomyocyte
fibroblast
endothelial_or_vascular
immune

这里的关键不是注释名长什么样，而是你要把原始注释映射成统一的上层类。
比如：

aCM1 / aCM2 / cycling CM → atrial_cardiomyocyte
FB1 / FB2 / activated FB → fibroblast
EC / lymphatic / vascular EC → 先可并成 endothelial_or_vascular
macrophage / monocyte / T / B / mast → 先并成 immune

理由很简单：
你当前目的是做 bulk 归因，不是单细胞精细图谱文章。
上层类足够稳，也更适合 donor 配对。

三、最关键的原则：不看逐细胞，只看 pseudobulk

这一步的灵魂就是：

推断单位是 donor，不是 cell。

3.3 pseudobulk 的最小粒度

按：

donor × chamber × cell_type

聚合。

例如：

D1_LA_atrialCM
D1_RA_atrialCM
D1_LA_fibro
D1_RA_fibro

每个 pseudobulk 只存 3 类信息：

聚合后的 counts
细胞数 n_cells
metadata：donor / chamber / cell_type
3.4 过滤规则

为了避免把稀疏 pseudobulk 也拿去建模，规则尽量简单：

组合层面过滤

每个 donor × chamber × cell_type 组合：

n_cells >= 30

低于 30，先丢掉。

细胞类型层面过滤

某个 cell type 进入主分析必须满足：

至少有 4 个 donor
且这些 donor 同时有 LA 和 RA

如果配对 donor < 4：

不做主分析
只做探索性补充

这样可以避免某个 cell type 看起来有细胞，其实没有足够配对信息。

四、每个 cell type 只做两件事

这是奥卡姆剃刀最重要的一刀：

每个 cell type 只产出 2 张表。

4.1 gene-level chamber DE

对每个 cell type 单独建模：

design_ct <- model.matrix(~ donor + chamber, data = meta_ct)

其中：

donor 是固定效应
chamber 只比较 LA vs RA

建议和 bulk 一样，用：

edgeR + voom + limma
或
直接 limma-trend

保持和前面方法一致，不要换工具链。

输出：

logFC
moderated_t
P.Value
FDR
direction

文件名：

atrialCM_chamber_DE.tsv
fibro_chamber_DE.tsv
endo_chamber_DE.tsv
immune_chamber_DE.tsv
4.2 program-level chamber FGSEA

和 bulk 保持同一套 gene sets：

Hallmark
Reactome

仍然用 gene-level 排名结果，比如 moderated_t 做 pre-ranked FGSEA。

输出：

pathway
NES
padj
direction
leading_edge_genes

文件名：

atrialCM_chamber_fgsea.tsv
fibro_chamber_fgsea.tsv
endo_chamber_fgsea.tsv
immune_chamber_fgsea.tsv

不单独解释 cell type，而是给 bulk 结果做归因

这里是这一步最核心的部分。

不是要写：

“atrial CM 有哪些 LA-high genes，fibro 有哪些 RA-high genes。”

而是要回答：

bulk 的 LA/RA 差异，主要由哪类细胞承载。

所以归因应该围绕 bulk 结果展开。