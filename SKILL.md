---
name: "scrna-geo-pipeline"
description: "GEO单细胞数据集→发表级图表全流程：自动格式识别、分步审查点、CNS可视化、HTML报告 (R/Seurat)"
---

# scrna-geo-pipeline

GEO下载的单细胞RNA-seq数据 → 发表级图表 + Seurat对象 + HTML综合报告。
R/Seurat 驱动，中文交互，每步出图可审查。

---

## 触发条件

用户给出单细胞数据集路径，关键词含：`单细胞`、`scRNA`、`GSE`、`10X`、`Seurat对象`、`降维聚类`、`细胞注释`

典型触发语：
- "这是单细胞数据，路径在 ~/Desktop/GSE12345/"
- "帮我分析这个 scRNA-seq 数据集"
- "对这个单细胞数据做全套分析"

---

## 工作流总览

```
Phase 0: 数据侦探 → 格式识别 + ENSEMBL转换 → 初检报告
Phase 1: 质量控制 → MAD阈值 + per-sample统计 + 双细胞检测 → 【审查点①】
Phase 2: 标准化   → SCT + 细胞周期 + 批次评估 → 【审查点②】
Phase 3: 降维     → Elbow + Harmony + LISI → UMAP → 【审查点③】
Phase 4: 聚类     → 多分辨率画廊 + Silhouette → 【审查点④】
Phase 5: 注释     → SingleR + FeaturePlot + 比例图 → 【审查点⑤】
Phase 6: 差异分析 → 增强火山 + 热图 + GO网络图 → 【审查点⑥】
Phase 7: 报告打包 → HTML报告 + 对象保存
```

**核心原则**：每完成一个 Phase，必须停下来展示结果 + 给出推荐，等用户确认后进入下一步。
**绝不一次性跑完全流程**，审查点是这个 skill 的灵魂。

---

## Phase 0: 数据侦探

### 目标
自动识别数据格式，读入或转换为标准 Seurat 对象。自动处理 ENSEMBL ID 转换。

### 执行流程

1. 用 `list.files(path, recursive=TRUE)` 扫描用户给的目录
2. 按优先级判定格式：

| 优先级 | 检测条件 | 格式 | 读取方法 |
|--------|---------|------|---------|
| 1 | 目录下存在 `matrix.mtx.gz` | 10X标准三件套 | `Read10X(dir)` |
| 2 | 存在 `.h5` 且非 `.h5ad` | 10X HDF5 | `Read10X_h5(file)` |
| 3 | 存在 `.rds` | Seurat对象 | `readRDS(file)` |
| 4 | 存在 `.h5ad` | Scanpy/AnnData | `SeuratDisk::Convert()` → `LoadH5Seurat()` |
| 5 | 存在 `.csv.gz` / `.tsv.gz` / `.txt.gz` | 纯计数矩阵 | `read.csv()` + 侦探自检 |

3. **矩阵文件侦探自检**（仅格式5）：

```r
# 报告以下信息：
# - 维度：N基因 × M细胞
# - 是否整数？→ 原始counts OR 已标准化
# - 最大值 → 辅助判断（>100可能是counts，<20可能是log）
# - 是否有负值？→ 已标准化（z-score等）
# - 基因名格式：Symbol / ENSEMBL ID / 混合
# - 计数稀疏度：zero_proportion
```

根据侦探结果自动决定：
- 整数且无负值 → 当作原始counts，直建Seurat对象
- 小数且有负值 → 警告用户数据已预处理过，询问是继续（跳过SCT）还是放弃
- 小数但全是正值 → 可能是TPM/FPKM，警告并用 `LogNormalize` 代替SCT

4. **RDS 状态侦查**（仅格式3）：
```r
# 检查 Seurat 对象已处理到哪步
# - 有无 counts assay？default assay 是什么？
# - 是否已 NormalizeData？已 FindVariableFeatures？已 ScaleData？
# - 是否已 RunPCA？已 RunUMAP？已 FindClusters？
# 报告给用户：这个对象已经被处理到了XX步骤
```

5. **ENSEMBL ID → Symbol 自动转换**

如果侦探发现基因名是 ENSEMBL 格式（`^ENSG` 开头占比 > 80%）：
```r
# 自动转换
library(org.Hs.eg.db)
gene_ids <- gsub("\\..*", "", rownames(counts))  # 去掉版本号后缀
symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, 
                  column = "SYMBOL", keytype = "ENSEMBL", 
                  multiVals = "first")
# 去重：ENSEMBL→相同Symbol的取平均值
# 保留有Symbol的行，报告转换率和丢失数
```

在初检报告中单独说明：
```
🔤 ENSEMBL转换: 25,032 → 23,891 基因 (95.4%)
   丢失: 1,141 (无对应Symbol，可能是非编码或假基因)
   去重合并: 342 个多对一映射已处理
```

### 输出
- 格式检测报告 + 基本统计（基因数、细胞数、sparsity）
- ENSEMBL转换报告（如适用）
- 已构建的 Seurat 对象

---

## Phase 1: 质量控制 + 双细胞检测

### 目标
计算QC指标，MAD自动推荐过滤阈值，检测双细胞，per-sample QC统计，等用户拍板。

### Step A: QC 指标计算

```r
# 计算QC指标
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[ABDEGQZ]")

# 基础统计
obj[["log10GenesPerUMI"]] <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
```

### Step B: MAD 自动阈值推荐

**核心方法：MAD（Median Absolute Deviation）**

```r
# nFeature_RNA: 下界 = max(200, median - 3*MAD), 上界 = median + 3*MAD
# nCount_RNA: 同上
# percent.mt: 上界 = min(median + 3*MAD, 20%), 以20%为天花板
# 如果 MAD 推导出的上界 < 5%，则用 5%

# 自检1：推荐阈值过滤>50%细胞 → 警告
# 自检2：上界推荐<500 → 可能有问题，提示
```

### Step C: Per-Sample QC 统计表

**重要：多数 GEO 数据集是多样本的，这个表是 paper Table S1 标配**

```r
# 生成 per-sample 表格
sample_qc <- obj@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    n_cells_raw = n(),
    mean_nCount = mean(nCount_RNA),
    median_nCount = median(nCount_RNA),
    mean_nFeature = mean(nFeature_RNA),
    median_nFeature = median(nFeature_RNA),
    mean_pct_mt = mean(percent.mt),
    median_pct_mt = median(percent.mt),
    .groups = "drop"
  )
```

生成 per-sample 小提琴图（每个样本一条，fill by sample），用 cns_palette 着色。

### Step D: 双细胞检测 (scDblFinder)

```r
library(scDblFinder)

# 按样本独立检测双细胞（避免样本间差异被误判为doublet）
set.seed(42)
sce <- as.SingleCellExperiment(obj)
sce <- scDblFinder(sce, samples = "orig.ident")
obj$doublet_score <- sce$scDblFinder.score
obj$doublet_class <- sce$scDblFinder.class  # "singlet" / "doublet"
```

生成双细胞可视化：
1. **UMAP按doublet_score着色**（渐变色，高亮潜在doublet）
2. **双细胞统计表**：per-sample doublet rate

```r
doublet_summary <- obj@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    n_cells = n(),
    n_doublets = sum(doublet_class == "doublet"),
    doublet_rate = round(100 * n_doublets / n_cells, 2)
  )
```

**自动判定**：
- doublet_rate < 1%：报告"双细胞率极低，数据质量良好"
- doublet_rate 1-5%：正常范围
- doublet_rate > 5%：警告"双细胞率偏高，检查数据质量"

### 审查点① —— 展示内容

**必须展示 4 组内容，一次呈现后一起提问：**

1. **QC 小提琴图**（5面板：nFeature, nCount, %MT, %RB, log10GenesPerUMI，水平虚线标推荐阈值）
2. **QC 散点图**（nCount vs nFeature + nCount vs %MT，标过滤边界）
3. **双细胞 UMAP**（doublet_score 着色）+ **per-sample doublet 统计表**
4. **Per-sample QC 统计表**（markdown 表格）
5. **推荐阈值文本摘要**

```
📊 质控报告

推荐过滤：
- nFeature_RNA: [下界] ~ [上界]
- nCount_RNA:   [下界] ~ [上界]  
- percent.mt:   < [值]%
- percent.rb:   [参考]
- 综合保留: [N]/[M] 细胞 ([X]%)

双细胞检测 (scDblFinder)：
- 检出双细胞: [N] / [M] ([X]%)
- Singlets保留: [Y]

QC图 + per-sample统计表见附件。

用推荐阈值继续？还是需要调整？
- 回复「继续」→ 用推荐阈值 + 移除双细胞
- 回复「nFeature 300-5000, mt<10」→ 自定义阈值
- 回复「保留双细胞」→ 不过滤双细胞（不推荐）
```

---

## Phase 2: 标准化 + 细胞周期 + 批次评估

### 目标
SCTransform 标准化，评估细胞周期效应，评估批次效应。

### Step A: 细胞周期评分

```r
# 使用 Seurat 内置的细胞周期基因集
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

obj <- CellCycleScoring(obj, s.features = s.genes, 
                         g2m.features = g2m.genes, set.ident = FALSE)
```

**细胞周期评估**：
```r
# 检查 S.Score 和 G2M.Score 的分布
# 如果细胞周期相关基因在 HVG 中占比 >15%，提示 regress
hvgs <- VariableFeatures(obj)
cc_genes_all <- c(s.genes, g2m.genes)
cc_in_hvg <- intersect(cc_genes_all, hvgs)
cc_hvg_pct <- length(cc_in_hvg) / length(hvgs) * 100
```

生成：
1. **细胞周期 UMAP**（按 Phase 着色，G1/S/G2M 三色）
2. **细胞周期分布饼图或bar图**（per-sample）

### Step B: SCTransform

```r
# 如果细胞周期效应明显 (cc_hvg_pct > 15%)，回归掉 S.Score 和 G2M.Score 的差异
# 注意：不要回归 Phase（分类变量会让SCT崩溃），只回归连续分数差
if (cc_hvg_pct > 15) {
  obj <- SCTransform(obj, vst.flavor = "v2", 
                     vars.to.regress = c("S.Score", "G2M.Score"),
                     verbose = FALSE)
} else {
  obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE)
}

# 多样本/多数据集时按 orig.ident 分别SCT
```

### Step C: 批次评估

- `length(unique(obj$orig.ident))` = 1 → 无批次，跳过
- ≥2 → 快速PCA看样本分离情况
- PCA按样本着色 + 每个样本的细胞分布密度

### 审查点② —— 展示内容

1. **HVG散点图**（前3000高变基因标记）
2. **细胞周期 UMAP**（Phase着色 + per-sample比例图）
3. **细胞周期评估文本**（cc_in_hvg占比，是否回归）
4. **批次PCA图**（如有多个样本，按orig.ident着色）

```
📊 标准化完成 (SCTransform v2)

细胞周期评估：
- 细胞周期基因在HVG中占比: [X]%（[是否超过15%]）
- [若超过]已在SCT中回归S/G2M分数差异

批次效应评估：
- [N]个样本/批次
- [有/无]明显批次分离

- 回复「继续」→ 进入降维
- 回复「不回归细胞周期」→ 重新SCT但不回归
```

---

## Phase 3: 降维 + 批次校正

### 目标
PCA降维 → 自动选PC数 → Harmony校正（如有批次） → UMAP可视化 → LISI量化评估

### Step A: PCA + 自动选PC

```r
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

# Elbow图 + 自动选PC（两条规则取max）：
# 规则A：连续两点方差下降<0.1%的第一个点
# 规则B：累积方差>80%所需最少PC数
```

### Step B: Harmony 批次校正（如有批次效应）

```r
obj <- RunHarmony(obj, group.by.vars = "orig.ident", 
                  reduction = "pca", assay.use = "SCT",
                  dims.use = 1:npc_auto, max.iter.harmony = 20,
                  theta = 2)  # theta=2适合同种组织不同样本
```

### Step C: LISI 量化评估 ⭐ 新增

**比肉眼对比更客观，paper审稿人认可**

```r
library(lisi)

# Harmony 前
lisi_before <- compute_lisi(
  Embeddings(obj, "pca")[, 1:npc_auto],
  obj@meta.data,
  "orig.ident"
)

# Harmony 后
lisi_after <- compute_lisi(
  Embeddings(obj, "harmony")[, 1:npc_auto],
  obj@meta.data,
  "orig.ident"
)

# LISI = 1 表示完美混合，LISI = N（样本数）表示完全分离
# 报告 mean(LISI_before) vs mean(LISI_after)
```

生成：
1. **LISI 分布对比小提琴图**（before vs after，一条水平虚线=1）
2. **Harmony 前后 UMAP 对比**（按样本着色，并排）

### Step D: UMAP

```r
obj <- RunUMAP(obj, reduction = if(has_batch) "harmony" else "pca", 
               dims = 1:npc_auto, verbose = FALSE)
```

### 审查点③ —— 展示内容

1. **Elbow图**（标注推荐PC数 + 两条规则说明）
2. **UMAP**（按样本着色）
3. 如有Harmony：**Harmony前后 UMAP 对比 + LISI 小提琴图**

```
📊 降维完成

推荐PC数: [N]（Elbow拐点=[X], 80%方差需要=[Y]）
UMAP使用: [pca / harmony] reduction

[如有Harmony]
批次校正效果 (LISI)：
  Harmony前: [mean_before] → Harmony后: [mean_after]
  (LISI≈1表示完美混合)

- 回复「继续」→ 用[N]个PC进入聚类
- 回复「PC=30」→ 自定义PC数
```

---

## Phase 4: 聚类

### 目标
多分辨率聚类 → Silhouette自动推荐 → 等用户选择。

### 执行流程

```r
# 构建KNN图
obj <- FindNeighbors(obj, reduction = reduction_name, dims = 1:npc)

# 多分辨率聚类
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (res in resolutions) {
  obj <- FindClusters(obj, resolution = res, algorithm = 1, verbose = FALSE)
}

# Silhouette score 对每个分辨率
# 使用前npc个harmony/pca维度计算距离
```

### 自动推荐策略

三条规则综合：
1. **Silhouette Score**：最高优先
2. **合理簇数范围**：√(细胞数)/2 ~ √(细胞数)×2
3. **稳定性**：分辨率跃迁时cluster数量变化不过于剧烈

如果最高silhouette在合理范围内 → 推荐。
否则推荐合理范围内silhouette最高的。

### 审查点④ —— 展示内容

1. **多分辨率 UMAP 画廊**（2×3面板）
2. **Silhouette趋势图**（X=分辨率，Y=silhouette，标注推荐）
3. **Clustree 树图**（如有clustree包）

```
📊 聚类完成

| 分辨率 | 簇数 | Silhouette |
|--------|------|------------|
| 0.2    | 8    | 0.23       |
| 0.4    | 14   | 0.27 ★     |
| ...    | ...  | ...        |

推荐: res = 0.4 (14簇)

- 回复「继续」→ 用推荐分辨率
- 回复「res=0.8」→ 自选
```

---

## Phase 5: 细胞注释 + Feature Plot + 比例图

### 目标
SingleR自动注释 → FeaturePlot marker验证 → 细胞比例图 → 用户审核修正。

### Step A: 参考数据库选择

```
请选择参考数据库：
1. HumanPrimaryCellAtlasData (HPCA) — 713种细胞，推荐首选
2. BlueprintEncodeData — 造血/基质
3. MonacoImmuneData — 免疫精细
4. DatabaseImmuneCellExpressionData — 免疫
5. NovershternHematopoieticData — 造血

你的数据类型？[回复数字]
```

如果用户不确定：
- NPC/HNSCC/实体瘤 → HPCA
- 血液/免疫 → Blueprint + Monaco

### Step B: SingleR 自动注释

```r
pred <- SingleR(test = GetAssayData(obj, assay = "SCT", layer = "data"),
                ref = ref,
                labels = ref$label.main,
                clusters = obj$seurat_clusters)
obj$singleR_label <- pred$labels[match(obj$seurat_clusters, rownames(pred))]
```

### Step C: Feature Plot 验证 ⭐

```r
# 每个推测的细胞类型，选择2-3个经典marker做FeaturePlot
# 例：预测是T cells → FeaturePlot(c("CD3D", "CD3E", "CD2"))
# 生成 Feature Plot 画廊（多行多列，每行一种细胞类型）

# 经典marker映射表（内置知识）：
feature_map <- list(
  "T cells" = c("CD3D", "CD3E", "CD2"),
  "CD4+ T cells" = c("CD3D", "CD4", "IL7R"),
  "CD8+ T cells" = c("CD3D", "CD8A", "CD8B"),
  "NK cells" = c("NKG7", "GNLY", "KLRD1"),
  "B cells" = c("CD79A", "MS4A1", "CD19"),
  "Plasma cells" = c("MZB1", "SDC1", "JCHAIN"),
  "Monocytes" = c("CD14", "LYZ", "S100A9"),
  "Macrophages" = c("CD68", "CD163", "CSF1R"),
  "Dendritic cells" = c("FCER1A", "CST3", "CLEC10A"),
  "Endothelial cells" = c("PECAM1", "CDH5", "VWF"),
  "Fibroblasts" = c("COL1A1", "COL1A2", "DCN"),
  "Epithelial cells" = c("EPCAM", "KRT19", "CDH1"),
  "Myeloid cells" = c("ITGAM", "CD14", "LYZ")
)
```

### Step D: 细胞比例图 ⭐

```r
# 1. Per-sample cell type proportion bar plot (堆叠柱状图)
# 2. Per-sample box plot (每种细胞类型一个box，X轴细胞类型，Y轴比例)
#    这对比较组间细胞组成差异非常关键

# 配色使用 cns_palette(n_cell_types_数量)
```

### Step E: Marker 验证表

```r
all_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25)
top3 <- all_markers %>% group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC)

# 对照表：Cluster | SingleR | Top Markers | Known Markers | Match
```

### 审查点⑤ —— 展示内容（最多的一组图）

1. **UMAP标注图**（SingleR 预测的细胞类型，cns_palette着色）
2. **Feature Plot画廊**（每种细胞类型2-3个经典marker）
3. **细胞比例堆叠柱状图**（per-sample）
4. **细胞比例箱线图**（per-celltype across samples）
5. **Marker验证对照表**

```
📊 细胞注释完成

SingleR 自动注释: [N]种细胞类型
Feature Plot验证 + 比例图见附件。
Marker验证对照表见附件。

- 回复「继续」→ 注释OK，进入差异分析
- 回复「Cluster X 改成 Y细胞」→ 手动修正
- 回复「重新注释，用 Blueprint」→ 换参考数据库
```
手动修正后重新出UMAP标注图确认。

---

## Phase 6: 差异分析 + 增强可视化

### 目标
按用户指定分组做差异分析，出增强火山图、热图、GO网络图。

### Step A: 确认比较组

```
要对哪些组做差异分析？

当前可用分组变量：
- seurat_clusters / cell_type
- orig.ident（样本来源）
- [其他meta列]

请指定：分组变量 + 比较哪两个组？
也支持"每种细胞类型 vs 其余所有"
```

### Step B: 差异分析

```r
deg <- FindMarkers(obj, ident.1 = group_a, ident.2 = group_b,
                   min.pct = 0.25, logfc.threshold = 0.25)
```

### Step C: 增强火山图 ⭐

比基础版多：
1. **基因名标注**：top10上调和top10下调基因标名称（ggrepel避重）
2. **双色着色**：上调=红色(#FE0500)，下调=蓝色(#306AF0)，不显著=灰色
3. **统计标注**：显著性阈值虚线 + 图内标注"N上调=[X], N下调=[Y]"
4. **FC阈值线**：|logFC|=0.5 垂直虚线
5. **自定义标题**："[A组] vs [B组]"

### Step D: 热图

- Top 50 DEGs（上调和下调各取top25），z-score标准化
- 按细胞类型/组别分组annotation bar
- 蓝白红色阶

### Step E: GO/KEGG cnetplot 网络图 ⭐

```r
library(clusterProfiler)

ego <- enrichGO(gene = deg_genes, OrgDb = org.Hs.eg.db, 
                ont = "BP", pAdjustMethod = "BH")

# cnetplot: 基因-通路网络图
# 展示哪些基因富集到哪些通路，比气泡图更直观
cnetplot(ego, showCategory = 8, 
         node_label = "category",  # 或 "all" 
         circular = FALSE,
         colorEdge = TRUE,
         cex_label_category = 0.8,
         cex_label_gene = 0.6)
```

同时出：
1. **cnetplot**（基因-通路网络，主图）
2. **气泡图**（富集通路概览，补充）
3. **上调+下调分开做富集**

### 审查点⑥ —— 展示内容

1. **增强火山图**（标注top基因，双色，统计信息）
2. **热图**（top50 DEGs）
3. **GO/KEGG cnetplot**（基因-通路网络）
4. **GO/KEGG气泡图**（富集概览）

```
📊 差异分析完成

[A组] vs [B组]
- 上调: [N] genes  |  下调: [M] genes (|logFC|>0.5, p_adj<0.05)

增强火山图 + 热图 + GO网络图见附件。

- 回复「继续」→ 生成最终报告
- 回复「再加一组比较」→ 再跑一次Phase 6
- 回复「下调FC阈值到0.25」→ 调整参数重跑
```

---

## Phase 7: 报告打包

### 目标
生成HTML综合报告 + 保存最终Seurat对象 + 打包所有图表。

### 执行流程

1. **保存最终 Seurat 对象**
```r
saveRDS(obj, file = file.path(output_dir, "seurat_final.rds"))
```

2. **图表命名规范**
```
figures/
├── 00_format_report.png
├── 01_qc_violin.png
├── 01_qc_scatter.png
├── 01_doublet_umap.png
├── 01_per_sample_qc.png
├── 02_hvg.png
├── 02_cellcycle_umap.png
├── 02_cellcycle_bar.png
├── 02_batch_pca.png
├── 03_elbow.png
├── 03_lisi_violin.png
├── 03_harmony_before_umap.png
├── 03_harmony_after_umap.png
├── 04_umap_multires.png
├── 04_silhouette.png
├── 05_umap_annotation.png
├── 05_feature_plots.png
├── 05_cell_proportion_bar.png
├── 05_cell_proportion_box.png
├── 06_volcano.png
├── 06_heatmap.png
├── 06_go_cnet.png
└── 06_go_dotplot.png
```

3. **HTML报告生成**（RMarkdown）包含：
- 数据集概览
- 各Phase参数选择
- 全部图表（内嵌）
- Per-sample QC统计表
- Marker验证对照表
- 差异基因结果表
- 分析日志

4. **分析摘要**
```
📋 分析摘要
━━━━━━━━━━━━━━━━━━━━━
数据集: GSE12345
细胞数: 14,230 → 13,850 (QC) → 13,520 (去doublet)
使用PC: 15 | 聚类分辨率: 0.4
细胞类型数: 14 | 差异比较组: [N]
输出: ~/Desktop/GSE12345_analysis/
  ├── seurat_final.rds
  ├── analysis_report.html
  └── figures/
```

---

## 🎨 CNS 可视化标准

### 配色方案

```r
cns_colors <- c(
  "#670073", "#306AF0", "#54F90B", "#FFBE03", "#FE0500", "#AA0C00"
)

# 多细胞类型扩展
extended_palette <- c(
  "#670073", "#306AF0", "#54F90B", "#FFBE03", "#FE0500", "#AA0C00",
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
  "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
  "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"
)

theme_cns <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )
```

### 图表规范

| 图类型 | 基本设置 |
|--------|---------|
| UMAP/tSNE | pt.size=0.3, raster=TRUE（>1万细胞）, order=shuffle |
| Feature Plot | order=TRUE（高表达在上面）, pt.size=0.3 |
| 小提琴图 | pt.size=0, 填充用cns_palette |
| 增强火山图 | |logFC|>0.5虚线, p_adj=0.05虚线, ggrepel标注top10± |
| 热图 | z-score, 蓝白红色阶, complexheatmap或pheatmap |
| 气泡图 | 点大小=表达比例, 颜色=平均表达 |
| 比例图 | 按cns_palette着色, 排序按丰度降序 |
| **所有图表** | **300 DPI**, `ggsave(width=8, height=6)` 默认 |
| KM曲线 | 中蓝#5B7BA7 / 中红#D05555, surv.median.line, Number at risk |

---

## ⚠️ 内存管理

- MacBook Air M5 16GB RAM
- **细胞数阈值**：
  - ≤50,000：正常运行
  - 50,000-100,000：提示用户建议降采样到50k
  - >100,000：强制降采样到50k

```r
if (ncol(obj) > 100000) {
  set.seed(42)
  obj <- subset(obj, cells = sample(Cells(obj), 50000))
} else if (ncol(obj) > 50000) {
  # 询问用户是否降采样
}
```

- **并行计算**：`future::plan("multisession", workers = min(4, parallel::detectCores()))`

---

## 📦 依赖管理

### 必需
`Seurat` (≥5.0), `tidyverse`, `patchwork`, `SingleR`, `celldex`, `scDblFinder`

### 按需加载
`harmony`, `lisi`, `clustree`, `cluster`, `clusterProfiler`, `org.Hs.eg.db`, `SeuratDisk`, `ggrepel`, `pheatmap`/`ComplexHeatmap`

### 缺少包时
1. 告知用户缺少哪个包
2. 执行 `install.packages()` / `BiocManager::install()`
3. 失败则给出手动命令

---

## 🔧 常见问题排查

| 症状 | 可能原因 | 处理 |
|------|---------|------|
| QC后细胞几乎全丢 | 阈值太严 | 放宽MAD倍数3→2.5，或手动设 |
| 双细胞率>10% | 10X过载或组织解离差 | 警告用户检查原始数据 |
| 细胞周期基因在HVG>25% | 增殖信号太强 | 确认SCT回归S/G2M |
| LISI Harmony后仍>2 | 批次效应太强 | 换CCA-LIGER或报告给用户 |
| UMAP全混在一起 | 数据可能已标准化 | 检查Phase 0侦探结果 |
| FindMarkers全不显著 | 组间差异小 | 降低logfc.threshold |
| 内存不足 | 数据太大 | 降采样 |
| SCT报错 | 基因太少(已过滤) | 用LogNormalize+ScaleData替代 |

---

## 📂 输出目录结构

```
{输入目录}_analysis/
├── seurat_final.rds
├── analysis_report.html
├── figures/                  # 所有图表（见Phase 7命名规范）
└── tables/
    ├── per_sample_qc.csv
    ├── doublet_summary.csv
    ├── marker_validation.csv
    └── deg_results.csv
```
