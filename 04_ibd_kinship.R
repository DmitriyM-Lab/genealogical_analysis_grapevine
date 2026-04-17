library(SNPRelate)
library(openxlsx)

# --- Workspace Setup ---
cat("A window will open: please select the input Clean GDS file...\n")
gds_clean.fn <- file.choose() 

work_dir <- dirname(gds_clean.fn)
setwd(work_dir)
cat("Working directory set to:", work_dir, "\n")

# --- 1. Calculate IBD (KING-robust) ---
# The KING-robust method is specifically designed for samples 
# with population stratification (different varieties/geography)
cat("\n--- Calculating IBD (KING-robust) & Kinship ---\n")
genofile <- snpgdsOpen(gds_clean.fn)

ibd <- snpgdsIBDKING(genofile, autosome.only = FALSE, type = "KING-robust", verbose = TRUE)

# --- 2. Prepare Results Table ---
# KING-robust returns IBS0 and kinship instead of standard k0 and k1
ibd_results <- snpgdsIBDSelection(ibd)

# --- 3. Relationship Classification Logic ---
# Based on theoretical thresholds from Manichaikul et al., 2010 (Table 1)
ibd_results$Relationship <- NA

# 1. Monozygotic twin / Clone / Duplicates
ibd_results$Relationship[ibd_results$kinship > 0.3535] <- "Monozygotic twin / Clone"

# 2. Parent-Offspring 
# (1st-degree kinship [0.1767 - 0.3535] + IBS0 near zero. < 0.005 allows for genotyping errors)
ibd_results$Relationship[ibd_results$kinship >= 0.1767 & ibd_results$kinship <= 0.3535 & ibd_results$IBS0 < 0.005] <- "Parent-offspring"

# 3. Full sibs (Remaining 1st-degree pairs with IBS0 > 0.005)
ibd_results$Relationship[ibd_results$kinship >= 0.1767 & ibd_results$kinship <= 0.3535 & is.na(ibd_results$Relationship)] <- "Full sib"

# 4. 2nd Degree
ibd_results$Relationship[ibd_results$kinship >= 0.0884 & ibd_results$kinship < 0.1767 & is.na(ibd_results$Relationship)] <- "2nd Degree"

# 5. 3rd Degree
ibd_results$Relationship[ibd_results$kinship >= 0.0442 & ibd_results$kinship < 0.0884 & is.na(ibd_results$Relationship)] <- "3rd Degree"

# 6. Unrelated
# Kinship below 0.0442 (including negative values).
# Strongly negative values indicate samples originate from different populations.
ibd_results$Relationship[ibd_results$kinship < 0.0442 & is.na(ibd_results$Relationship)] <- "Unrelated"

# Sort table by descending kinship
ibd_results <- ibd_results[order(ibd_results$kinship, decreasing = TRUE), ]

# Reorder columns (removed k0, k1, and LTR)
ibd_results <- ibd_results[, c("ID1", "ID2", "IBS0", "kinship", "Relationship")]

# --- 4. Export and Cleanup ---
output_excel <- "IBD_KING_Robust_Results.xlsx"
write.xlsx(ibd_results, output_excel, overwrite = TRUE)

# Safe close
if (exists("genofile") && !is.null(genofile)) {
  try({ snpgdsClose(genofile) }, silent = TRUE)
}

cat("\nSUCCESS! Results saved to:", output_excel, "\n")
