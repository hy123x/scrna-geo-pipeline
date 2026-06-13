# Phase 0 Helper: Detect scRNA-seq data format from GEO download
# Returns list with format type and loading instructions

detect_sc_format <- function(data_path) {
  all_files <- list.files(data_path, recursive = TRUE, full.names = TRUE)
  all_names <- basename(all_files)
  
  result <- list(
    format = NA_character_,
    method = NA_character_,
    files = character(0),
    warnings = character(0),
    n_cells = NA_integer_,
    n_genes = NA_integer_,
    sparsity = NA_real_
  )
  
  # --- Priority 1: 10X directory (barcodes + features + matrix) ---
  has_mtx <- any(grepl("matrix\\.mtx(\\.gz)?$", all_names))
  has_barcodes <- any(grepl("barcodes\\.tsv(\\.gz)?$", all_names))
  has_features <- any(grepl("features\\.tsv(\\.gz)?$", all_names))
  
  if (has_mtx && has_barcodes && has_features) {
    mtx_dir <- unique(dirname(all_files[grepl("matrix\\.mtx", all_names)]))
    result$format <- "10X_directory"
    result$method <- "Read10X"
    result$files <- mtx_dir
    return(result)
  }
  
  # --- Priority 2: 10X HDF5 ---
  h5_files <- all_files[grepl("\\.h5$", all_names)]
  if (length(h5_files) > 0) {
    result$format <- "10X_HDF5"
    result$method <- "Read10X_h5"
    result$files <- h5_files[1]
    return(result)
  }
  
  # --- Priority 3: Seurat RDS ---
  rds_files <- all_files[grepl("\\.rds$", all_names)]
  if (length(rds_files) > 0) {
    result$format <- "Seurat_RDS"
    result$method <- "readRDS"
    result$files <- rds_files[1]
    return(result)
  }
  
  # --- Priority 4: h5ad (Scanpy) ---
  h5ad_files <- all_files[grepl("\\.h5ad$", all_names)]
  if (length(h5ad_files) > 0) {
    result$format <- "h5ad"
    result$method <- "SeuratDisk"
    result$files <- h5ad_files[1]
    result$warnings <- c(result$warnings, 
      "h5ad格式需要SeuratDisk包转换")
    return(result)
  }
  
  # --- Priority 5: Flat count matrix ---
  csv_files <- all_files[grepl("\\.(csv|tsv|txt)(\\.gz)?$", all_names)]
  if (length(csv_files) > 0) {
    sizes <- file.info(csv_files)$size
    target <- csv_files[which.max(sizes)]
    result$format <- "count_matrix"
    result$method <- "read_matrix"
    result$files <- target
    result$warnings <- c(result$warnings,
      "纯计数矩阵格式，将进行侦探自检以确认数据类型")
    return(result)
  }
  
  result$format <- "unknown"
  result$warnings <- c(result$warnings,
    paste("无法识别数据格式。找到的文件：", 
          paste(all_names, collapse = ", ")))
  return(result)
}

# Matrix self-inspection for flat count matrix files
inspect_matrix <- function(mat, max_preview = 5) {
  info <- list()
  
  info$n_genes <- nrow(mat)
  info$n_cells <- ncol(mat)
  info$is_integer <- all(mat == round(mat), na.rm = TRUE)
  info$has_negative <- any(mat < 0, na.rm = TRUE)
  info$max_value <- max(mat, na.rm = TRUE)
  info$min_value <- min(mat, na.rm = TRUE)
  info$mean_value <- mean(mat, na.rm = TRUE)
  info$sparsity <- sum(mat == 0) / (nrow(mat) * ncol(mat))
  
  # Gene name format detection
  gene_names <- rownames(mat)
  if (length(gene_names) > 0) {
    ensembl_ratio <- mean(grepl("^ENSG", gene_names))
    info$gene_format <- if (ensembl_ratio > 0.8) "ENSEMBL" else "Symbol"
    info$gene_preview <- head(gene_names, max_preview)
    info$gene_suffix <- any(grepl("\\.\\d+$", gene_names))
  }
  
  # Verdict
  if (info$is_integer && !info$has_negative && info$max_value > 100) {
    info$verdict <- "raw_counts"
    info$action <- "直接创建Seurat对象，后续使用SCT"
  } else if (info$has_negative) {
    info$verdict <- "normalized"
    info$action <- "数据已标准化（含负值），无法用SCT，需用LogNormalize"
  } else if (!info$is_integer && info$max_value < 30) {
    info$verdict <- "log_transformed"
    info$action <- "数据已log变换，建议用LogNormalize"
  } else if (!info$is_integer && info$max_value > 100) {
    info$verdict <- "probable_TPM"
    info$action <- "可能是TPM/FPKM，用LogNormalize+ScaleData替代SCT"
  } else {
    info$verdict <- "uncertain"
    info$action <- "不确定，手动确认"
  }
  
  return(info)
}

# Inspect existing Seurat object state
inspect_seurat <- function(obj) {
  state <- list()
  
  state$n_cells <- ncol(obj)
  state$n_genes <- nrow(obj)
  state$assays <- names(obj@assays)
  state$default_assay <- DefaultAssay(obj)
  state$has_counts <- "RNA" %in% names(obj@assays) && 
    !is.null(obj@assays$RNA$counts)
  state$has_data <- "RNA" %in% names(obj@assays) && 
    !is.null(obj@assays$RNA$data)
  state$has_sct <- "SCT" %in% names(obj@assays)
  state$has_pca <- "pca" %in% names(obj@reductions)
  state$has_umap <- "umap" %in% names(obj@reductions)
  state$has_clusters <- "seurat_clusters" %in% colnames(obj@meta.data)
  state$n_clusters <- if (state$has_clusters) 
    length(unique(obj$seurat_clusters)) else NA_integer_
  
  if (!state$has_counts && !state$has_sct) {
    state$stage <- "empty"
  } else if (state$has_clusters && state$has_umap) {
    state$stage <- "fully_processed"
  } else if (state$has_umap) {
    state$stage <- "reduced"
  } else if (state$has_sct) {
    state$stage <- "normalized"
  } else if (state$has_counts) {
    state$stage <- "raw"
  } else {
    state$stage <- "unknown"
  }
  
  return(state)
}

# ENSEMBL ID to Symbol conversion
convert_ensembl_to_symbol <- function(counts) {
  # counts: matrix with ENSEMBL gene IDs as rownames
  library(org.Hs.eg.db)
  
  n_before <- nrow(counts)
  
  # Remove version suffix (ENSG00000123456.12 → ENSG00000123456)
  gene_ids <- gsub("\\..*", "", rownames(counts))
  
  # Map to Symbol
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db, 
    keys = gene_ids, 
    column = "SYMBOL", 
    keytype = "ENSEMBL", 
    multiVals = "first"
  )
  
  # Remove genes with no Symbol
  keep <- !is.na(symbols)
  counts <- counts[keep, ]
  symbols <- symbols[keep]
  
  # Handle many-to-one (multiple ENSEMBL → same Symbol)
  # Take the one with highest mean expression
  dup_sym <- symbols[duplicated(symbols)]
  if (length(dup_sym) > 0) {
    for (sym in unique(dup_sym)) {
      idx <- which(symbols == sym)
      if (length(idx) > 1) {
        row_means <- rowMeans(counts[idx, , drop = FALSE])
        best <- idx[which.max(row_means)]
        idx_remove <- setdiff(idx, best)
        counts <- counts[-idx_remove, ]
        symbols <- symbols[-idx_remove]
      }
    }
  }
  
  rownames(counts) <- symbols
  
  result <- list(
    counts = counts,
    n_before = n_before,
    n_after = nrow(counts),
    n_lost = n_before - nrow(counts),
    pct_kept = round(100 * nrow(counts) / n_before, 1),
    n_dups_merged = length(unique(dup_sym))
  )
  
  return(result)
}
