# MAD-based threshold recommendation for scRNA-seq QC
# Returns recommended thresholds and expected retention

recommend_thresholds <- function(obj, 
                                  nFeature_mad = 3,
                                  nCount_mad = 3,
                                  mt_ceiling = 20,
                                  mt_min = 5) {
  
  metrics <- obj@meta.data
  
  # --- nFeature_RNA ---
  nf_median <- median(metrics$nFeature_RNA)
  nf_mad <- mad(metrics$nFeature_RNA)
  nf_lower <- max(200, nf_median - nFeature_mad * nf_mad)
  nf_upper <- nf_median + nFeature_mad * nf_mad
  
  # --- nCount_RNA ---
  nc_median <- median(metrics$nCount_RNA)
  nc_mad <- mad(metrics$nCount_RNA)
  nc_lower <- max(100, nc_median - nCount_mad * nc_mad)
  nc_upper <- nc_median + nCount_mad * nc_mad
  
  # --- percent.mt ---
  if ("percent.mt" %in% colnames(metrics)) {
    mt_median <- median(metrics$percent.mt)
    mt_mad <- mad(metrics$percent.mt)
    mt_upper_raw <- mt_median + 3 * mt_mad
    mt_upper <- min(mt_ceiling, mt_upper_raw)
    mt_upper <- max(mt_min, mt_upper)
  } else {
    mt_upper <- mt_ceiling
    warning("percent.mt not found, using ceiling value")
  }
  
  # --- Retention estimate ---
  keep <- metrics$nFeature_RNA > nf_lower &
    metrics$nFeature_RNA < nf_upper &
    metrics$nCount_RNA > nc_lower &
    metrics$nCount_RNA < nc_upper
  
  if ("percent.mt" %in% colnames(metrics)) {
    keep <- keep & metrics$percent.mt < mt_upper
  }
  
  n_total <- ncol(obj)
  n_keep <- sum(keep)
  pct_keep <- round(100 * n_keep / n_total, 1)
  
  # Return
  list(
    nFeature_lower = round(nf_lower),
    nFeature_upper = round(nf_upper),
    nCount_lower = round(nc_lower),
    nCount_upper = round(nc_upper),
    mt_upper = round(mt_upper, 1),
    nFeature_median = round(nf_median),
    nCount_median = round(nc_median),
    n_total = n_total,
    n_keep = n_keep,
    pct_keep = pct_keep,
    warning = if (pct_keep < 50) {
      sprintf("警告：推荐阈值将过滤掉 %.0f%% 细胞(>50%%), 建议手动设定",
              100 - pct_keep)
    } else {
      NULL
    }
  )
}

# Apply thresholds to Seurat object
apply_thresholds <- function(obj, thresholds) {
  cells_keep <- colnames(obj)[
    obj$nFeature_RNA > thresholds$nFeature_lower &
    obj$nFeature_RNA < thresholds$nFeature_upper &
    obj$nCount_RNA > thresholds$nCount_lower &
    obj$nCount_RNA < thresholds$nCount_upper
  ]
  
  if ("percent.mt" %in% colnames(obj@meta.data)) {
    cells_keep <- intersect(cells_keep,
      colnames(obj)[obj$percent.mt < thresholds$mt_upper])
  }
  
  subset(obj, cells = cells_keep)
}

# Per-sample QC summary
per_sample_qc <- function(obj) {
  obj@meta.data %>%
    group_by(orig.ident) %>%
    summarise(
      n_cells = n(),
      mean_nCount = round(mean(nCount_RNA), 1),
      median_nCount = round(median(nCount_RNA), 1),
      mean_nFeature = round(mean(nFeature_RNA), 1),
      median_nFeature = round(median(nFeature_RNA), 1),
      mean_pct_mt = round(mean(percent.mt), 2),
      median_pct_mt = round(median(percent.mt), 2),
      .groups = "drop"
    )
}
