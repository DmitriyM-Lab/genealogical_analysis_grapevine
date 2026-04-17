library(SNPRelate)
library(ggplot2)
library(ggrepel)

# --- Workspace Setup ---
cat("A window will open: please select the input clean GDS file...\n")
gds_clean.fn <- file.choose() 
work_dir <- dirname(gds_clean.fn)
setwd(work_dir)
genofile <- snpgdsOpen(gds_clean.fn)

# --- Run PCA ---
cat("\nRunning PCA...\n")
pca <- snpgdsPCA(genofile, autosome.only = FALSE, num.thread = 2, verbose = TRUE)

# Get proportion of genetic variance
pc_percent <- round(pca$varprop * 100, 2)

# Create a base dataframe for plotting
pca_tab <- data.frame(
  Sample_ID = pca$sample.id,
  PC1 = pca$eigenvect[, 1], 
  PC2 = pca$eigenvect[, 2],
  stringsAsFactors = FALSE
)

# --- Visualization & Export ---
cat("\nGenerating plot and saving as TIFF (600 DPI)...\n")

output_tiff <- "PCA_Plot.tiff"

# Create the plot using ggplot2
pca_plot <- ggplot(pca_tab, aes(x = PC1, y = PC2, label = Sample_ID)) +
  # Add dashed axes passing through zero
  geom_hline(yintercept = 0, color = "gray60", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray60", linetype = "dashed") +
  # Draw points
  geom_point(color = "forestgreen", size = 2.5) +
  # Smart labels using ggrepel (prevents overlapping)
  geom_text_repel(
    size = 3, 
    color = "gray20",
    max.overlaps = Inf,       # Allow drawing ALL labels
    box.padding = 0.5,        # Padding between text boxes
    point.padding = 0.3       # Padding between text and point
  ) +
  # Axis labels with variance percentage
  labs(
    x = paste0("PC 1 (", pc_percent[1], "%)"), 
    y = paste0("PC 2 (", pc_percent[2], "%)")
  ) +
  theme_bw()

# Export plot to high-resolution TIFF
ggsave(
  filename = output_tiff, 
  plot = pca_plot, 
  device = "tiff", 
  width = 8, 
  height = 6, 
  units = "in", 
  dpi = 600, 
  compression = "lzw"
)

snpgdsClose(genofile)

cat("\nSUCCESS! Files saved to:", work_dir, "\n")
cat("- Plot:", output_tiff, "\n")
