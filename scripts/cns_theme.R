# CNS Visual Theme for scRNA-seq pipeline
# Follows 老大's color standards from MEMORY.md

library(ggplot2)

# Core palette
cns_colors <- c(
  "#670073",  # deep purple
  "#306AF0",  # bright blue
  "#54F90B",  # neon green
  "#FFBE03",  # gold
  "#FE0500",  # bright red
  "#AA0C00"   # dark red
)

# Extended palette for many clusters/cell types
extended_palette <- c(
  "#670073", "#306AF0", "#54F90B", "#FFBE03", "#FE0500", "#AA0C00",
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
  "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
  "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"
)

# CNS ggplot theme
theme_cns <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    plot.margin = margin(10, 10, 10, 10),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA)
  )

# Save with CNS defaults
save_cns <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(filename, plot, width = width, height = height, dpi = dpi, 
         bg = "white", create.dir = TRUE)
}

# Color mapper: generate n colors from extended palette
cns_palette <- function(n) {
  if (n <= length(extended_palette)) {
    extended_palette[1:n]
  } else {
    colorRampPalette(extended_palette)(n)
  }
}
