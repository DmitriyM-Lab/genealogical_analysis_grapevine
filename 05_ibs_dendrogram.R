library(SNPRelate)
library(openxlsx)

# --- 1. Workspace Setup & File Selection ---
cat("Step 1: Please select the Clean GDS file...\n")
gds_clean.fn <- file.choose() 

cat("Step 2: Please select the IBD Results file (Excel)...\n")
ibd_results.fn <- file.choose()

work_dir <- dirname(gds_clean.fn)
setwd(work_dir)
cat("Working directory set to:", work_dir, "\n")

# --- 2. Calculate IBS Distance Matrix ---
cat("Calculating Identity-by-State (IBS) Distance...\n")
genofile <- snpgdsOpen(gds_clean.fn)
ibs <- snpgdsIBS(genofile, autosome.only = FALSE, num.thread = 2, verbose = FALSE)
sample_ids <- ibs$sample.id

# Convert IBS to genetic distance (1 - IBS)
ibs_dist_mat <- 1 - ibs$ibs
colnames(ibs_dist_mat) <- rownames(ibs_dist_mat) <- sample_ids

snpgdsClose(genofile)

# --- 3. UPGMA Hierarchical Clustering ---
cat("Running UPGMA clustering...\n")
GeneticDistance <- as.dist(ibs_dist_mat)
hc_upgma <- hclust(GeneticDistance, method = "average")
hc_upgma$labels <- sample_ids

# --- 4. Load & Filter IBD Kinship Data ---
ibd_data <- read.xlsx(ibd_results.fn)

rel_colors <- c(
  "Monozygotic twin / Clone" = "purple",
  "Parent-offspring" = "red",
  "Full sib" = "blue",
  "2nd Degree" = "green",
  "3rd Degree" = "orange"
)

final_links <- ibd_data[ibd_data$Relationship %in% names(rel_colors), ]
final_links <- unique(final_links)

# --- 5. Coordinate Preparation for Kinship Arrows ---
dend <- as.dendrogram(hc_upgma)
label_order <- labels(dend)
y_coords <- setNames(seq_along(label_order), label_order)

if (nrow(final_links) > 0) {
  final_links$y_node1 <- y_coords[final_links$ID1]
  final_links$y_node2 <- y_coords[final_links$ID2]
  
  valid_links <- final_links[!is.na(final_links$y_node1) & !is.na(final_links$y_node2), ]
  
  valid_links$y_dist <- abs(valid_links$y_node1 - valid_links$y_node2)
  valid_links$y_min <- pmin(valid_links$y_node1, valid_links$y_node2)
  valid_links$Rel_Factor <- factor(valid_links$Relationship, levels = names(rel_colors))
  
  # Sort to prevent overlapping (short arcs on bottom)
  valid_links <- valid_links[order(valid_links$Rel_Factor, valid_links$y_dist, valid_links$y_min), ]
} else {
  valid_links <- data.frame()
}

max_d <- max(ibs_dist_mat)

# --- 6. Dynamic Plot Width Calculation ---
max_char_len <- max(nchar(sample_ids))
est_text_width <- max_char_len * (max_d * 0.015) 
base_offset <- -est_text_width - (max_d * 0.08)  

cascade_step <- -max_d * 0.02     
cat_gap <- -max_d * 0.08          

start_x <- list()
curr_x <- base_offset

for (rel in names(rel_colors)) {
  start_x[[rel]] <- curr_x
  n_arrows <- sum(valid_links$Relationship == rel)
  if (n_arrows > 0) {
    curr_x <- curr_x + (n_arrows * cascade_step) + cat_gap
  }
}

x_limit_min <- max_d * 1.05  
x_limit_max <- curr_x + (-max_d * 0.05) 

# --- 7. Plotting Function ---
draw_dendrogram <- function() {
  par(mar = c(8, 4, 2, 6)) 
  
  plot(dend, horiz = TRUE, xlim = c(x_limit_min, x_limit_max), 
       main = "", axes = FALSE, leaflab = "perpendicular")
  
  axis_ticks <- seq(0, max_d, by = 0.05)
  axis(1, at = axis_ticks, labels = axis_ticks, cex.axis = 1.1)
  mtext("Genetic Dissimilarity (1 - IBS)", side = 1, line = 3, at = max_d / 2, cex = 1.3)
  
  if (nrow(valid_links) > 0) {
    shift_counters <- c("Monozygotic twin / Clone" = 0, "Parent-offspring" = 0, "Full sib" = 0, "2nd Degree" = 0, "3rd Degree" = 0)
    
    for(i in 1:nrow(valid_links)) {
      rel <- valid_links$Relationship[i]
      color <- rel_colors[rel]
      
      x_pos <- start_x[[rel]] + (cascade_step * shift_counters[[rel]])
      shift_counters[[rel]] <- shift_counters[[rel]] + 1
      
      text_padding <- max_d * 0.015
      w1 <- abs(strwidth(as.character(valid_links$ID1[i]), cex = 1))
      w2 <- abs(strwidth(as.character(valid_links$ID2[i]), cex = 1))
      
      s1_width <- -w1 - text_padding
      s2_width <- -w2 - text_padding
      
      segments(x0 = s1_width, y0 = valid_links$y_node1[i], x1 = x_pos, y1 = valid_links$y_node1[i], col = "gray65", lty = 3, lwd = 1.2)
      segments(x0 = s2_width, y0 = valid_links$y_node2[i], x1 = x_pos, y1 = valid_links$y_node2[i], col = "gray65", lty = 3, lwd = 1.2)
      
      arrows(x0 = x_pos, y0 = valid_links$y_node1[i], 
             x1 = x_pos, y1 = valid_links$y_node2[i], 
             col = color, lwd = 2.5, length = 0.07, angle = 22, code = 3)
    }
    
    tick_pos <- c()
    tick_labs <- c()
    
    kinship_text <- c(
      "Monozygotic twin / Clone" = "Mt/C\n(\u03D5\u22480.50)",
      "Parent-offspring" = "PO\n(\u03D5\u22480.25)",
      "Full sib" = "FS\n(\u03D5\u22480.25)",
      "2nd Degree" = "2nd Deg\n(\u03D5\u22480.12)",
      "3rd Degree" = "3rd Deg\n(\u03D5\u22480.06)"
    )
    
    for (rel in names(rel_colors)) {
      n_arrows <- sum(valid_links$Relationship == rel)
      if (n_arrows > 0) {
        center_x <- start_x[[rel]] + ((n_arrows - 1) * cascade_step / 2)
        tick_pos <- c(tick_pos, center_x)
        tick_labs <- c(tick_labs, kinship_text[rel])
      }
    }
    
    if (length(tick_pos) > 0) {
      axis(1, at = tick_pos, labels = tick_labs, padj = 0.5, cex.axis = 0.95, gap.axis = -1)
      axis_center <- (min(tick_pos) + max(tick_pos)) / 2
      mtext("Kinship Coefficients (\u03D5)", side = 1, line = 4.0, at = axis_center, cex = 1.3)
    }
  }
}

# --- 8. Export Plot ---
output_png <- "Dendrogram_IBS_Kinship.png"
png(output_png, width = 20, height = 14, units = "in", res = 600, bg = "white")
draw_dendrogram()
dev.off()

cat("SUCCESS! Dendrogram saved as:", output_png, "\n")
