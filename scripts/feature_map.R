# Known cell type marker gene mapping
# Used in Phase 5 for FeaturePlot validation

known_markers <- list(
  # T cells
  "T cells" = c("CD3D", "CD3E", "CD2", "TRAC"),
  "CD4+ T cells" = c("CD3D", "CD4", "IL7R", "CCR7"),
  "CD8+ T cells" = c("CD3D", "CD8A", "CD8B", "GZMK"),
  "Treg" = c("CD3D", "CD4", "FOXP3", "IL2RA", "CTLA4"),
  "Naive T cells" = c("CD3D", "CCR7", "SELL", "LEF1"),
  "Memory T cells" = c("CD3D", "CD44", "CD69", "S100A4"),
  "Cytotoxic T cells" = c("CD3D", "CD8A", "GZMB", "PRF1", "NKG7"),
  
  # NK cells
  "NK cells" = c("NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1"),
  
  # B cells
  "B cells" = c("CD79A", "MS4A1", "CD19", "BANK1", "PAX5"),
  "Plasma cells" = c("MZB1", "SDC1", "JCHAIN", "XBP1", "IGHG1"),
  "Germinal center B" = c("MS4A1", "BCL6", "AICDA", "MKI67"),
  
  # Myeloid
  "Monocytes" = c("CD14", "LYZ", "S100A9", "S100A8", "VCAN"),
  "CD14+ Monocytes" = c("CD14", "LYZ", "S100A9", "S100A8"),
  "CD16+ Monocytes" = c("FCGR3A", "MS4A7", "LST1", "AIF1"),
  "Macrophages" = c("CD68", "CD163", "CSF1R", "MRC1", "APOE"),
  "M1 Macrophages" = c("CD68", "IL1B", "TNF", "CXCL8"),
  "M2 Macrophages" = c("CD68", "CD163", "MRC1", "MSR1"),
  "Dendritic cells" = c("FCER1A", "CST3", "CLEC10A", "CLEC9A"),
  "pDC" = c("LILRA4", "IRF7", "TCF4", "IL3RA"),
  "cDC1" = c("CLEC9A", "XCR1", "CADM1", "BATF3"),
  "cDC2" = c("FCER1A", "CLEC10A", "CD1C", "SIRPA"),
  "Myeloid DC" = c("ITGAX", "ITGAM", "HLA-DRA"),
  
  # Granulocytes
  "Neutrophils" = c("CSF3R", "FCGR3B", "CXCR2", "S100A8", "S100A9"),
  "Mast cells" = c("KIT", "TPSAB1", "CPA3", "HDC"),
  "Eosinophils" = c("CCR3", "IL5RA", "PRG2", "EPX"),
  "Basophils" = c("FCER1A", "CCR3", "IL3RA", "GATA2"),
  
  # Stromal
  "Endothelial cells" = c("PECAM1", "CDH5", "VWF", "CLDN5", "ENG"),
  "Lymphatic endothelial" = c("PROX1", "LYVE1", "PDPN", "CCL21"),
  "Fibroblasts" = c("COL1A1", "COL1A2", "DCN", "LUM", "FAP"),
  "CAF" = c("COL1A1", "ACTA2", "FAP", "PDPN", "MMP11"),
  "Myofibroblasts" = c("COL1A1", "ACTA2", "TAGLN", "MYH11"),
  "Pericytes" = c("RGS5", "PDGFRB", "CSPG4", "ACTA2"),
  "Smooth muscle" = c("ACTA2", "MYH11", "TAGLN", "CNN1"),
  
  # Epithelial
  "Epithelial cells" = c("EPCAM", "KRT19", "CDH1", "KRT18", "CLDN4"),
  "Basal epithelial" = c("KRT5", "KRT14", "TP63", "ITGA6"),
  "Luminal epithelial" = c("KRT8", "KRT18", "EPCAM", "CDH1"),
  "Glandular epithelial" = c("EPCAM", "KRT7", "KRT18", "MUC1"),
  "Squamous epithelial" = c("KRT5", "KRT14", "DSP", "PKP1"),
  
  # Tumor
  "Malignant cells" = c("EPCAM", "KRT19", "MKI67", "PCNA"),
  
  # Other
  "Erythrocytes" = c("HBA1", "HBA2", "HBB", "HBD", "GYPA"),
  "Platelets/Megakaryocytes" = c("PPBP", "PF4", "ITGA2B", "GP9"),
  "Melanocytes" = c("PMEL", "TYRP1", "MITF", "MLANA")
)

# Get feature set for a predicted cell type
# Falls back to partial match or returns NULL
get_markers_for_type <- function(cell_type, fallback = TRUE) {
  # Exact match
  if (cell_type %in% names(known_markers)) {
    return(known_markers[[cell_type]])
  }
  
  # Partial match (e.g., "T cells (CD4+)" matches "T cells")
  if (fallback) {
    for (nm in names(known_markers)) {
      if (grepl(nm, cell_type, ignore.case = TRUE) ||
          grepl(cell_type, nm, ignore.case = TRUE)) {
        return(known_markers[[nm]])
      }
    }
  }
  
  return(NULL)
}
