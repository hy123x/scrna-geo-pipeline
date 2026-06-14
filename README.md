# 🧬 scrna-geo-pipeline

> GEO 单细胞数据集 → 发表级图表 + HTML 报告，零手动参数调优

R/Seurat 驱动的单细胞 RNA-seq 全流程分析管线，专为从 GEO 下载的数据设计。

## 特性

- 🔍 **5种 GEO 格式自动识别** — 10X 三件套 / HDF5 / Seurat RDS / h5ad / 纯矩阵
- 🎯 **分步审查点** — 每步出图可审查，不盲目跑完全流程
- 🧹 **自动参数推荐** — MAD 阈值 / Elbow 选 PC / Silhouette 选分辨率
- 🎨 **CNS 发表级可视化** — 统一配色主题，300 DPI，火山图/热图/FeaturePlot/cnetplot
- 📊 **双细胞 + 细胞周期 + 批次校正** — scDblFinder + Harmony + LISI 量化
- 🏷️ **SingleR 自动注释 + 经典 marker 验证** — 内置 50+ 细胞类型 marker 映射
- 📝 **HTML 综合报告** — 打包全部图表和统计表

## 工作流

```
Phase 0  数据侦探     5种格式识别 + ENSEMBL→Symbol
Phase 1  质控         MAD阈值 + per-sample表 + 双细胞    [审查]
Phase 2  标准化       SCT + 细胞周期 + 批次评估          [审查]
Phase 3  降维         Elbow + Harmony + LISI + UMAP      [审查]
Phase 4  聚类         多分辨率画廊 + Silhouette           [审查]
Phase 5  注释         SingleR + FeaturePlot + 比例图      [审查]
Phase 6  差异分析     增强火山 + 热图 + GO cnetplot       [审查]
Phase 7  报告打包     HTML + Seurat对象 + 全部图表
```

## 使用方法

把单细胞数据集路径告诉 AI Agent：

```
这是单细胞数据，路径在 ~/Desktop/GSE12345/
```

管线会自动识别格式，逐步分析，每完成一步停下来让你审查。

## 安装

### OpenClaw 用户

```bash
git clone https://github.com/hy123x/scrna-geo-pipeline.git \
  ~/.openclaw/workspace/skills/scrna-geo-pipeline
```

安装后在对话中说「单细胞分析 + 数据路径」即可自动触发。

### 手动使用

```bash
git clone https://github.com/hy123x/scrna-geo-pipeline.git
```

参考 `SKILL.md` 中的工作流说明，在 R 中按 Phase 0-7 逐步执行。

## 依赖

```r
# 必需
Seurat (≥5.0), tidyverse, patchwork, SingleR, celldex, scDblFinder

# 按需加载
harmony, lisi, clustree, cluster, clusterProfiler, org.Hs.eg.db, 
SeuratDisk, ggrepel, pheatmap
```

## 文件结构

```
scrna-geo-pipeline/
├── SKILL.md                  # 完整工作流说明书
├── scripts/
│   ├── cns_theme.R           # CNS配色 + ggplot主题
│   ├── format_detector.R     # 格式识别 + ENSEMBL转换
│   ├── mad_threshold.R       # MAD阈值 + per-sample QC
│   └── feature_map.R         # 50+细胞类型marker映射
└── README.md
```

## 适用场景

- 从 GEO 下载的单细胞 RNA-seq 数据快速分析
- 肿瘤微环境 / 免疫细胞图谱构建
- 多样本整合 + 批次校正
- 需要发表级图表的论文准备

## 注意

- 为 MacBook 16GB RAM 优化，>10万细胞自动降采样
- 默认人类数据（MT 基因用 `^MT-` 匹配，小鼠需改）
