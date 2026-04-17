library(SNPRelate)
library(writexl)

# --- Workspace Setup ---
cat("A window will open: please select the input VCF file...\n")
vcf.fn <- file.choose() 
work_dir <- dirname(vcf.fn)
setwd(work_dir)
cat("Working directory set to:", work_dir, "\n")

raw_gds <- "vcf_data.gds"
output_excel <- "Individual_Heterozygosity_Results.xlsx"

if (file.exists(raw_gds)) file.remove(raw_gds)

# Convert VCF to GDS
snpgdsVCF2GDS(
  vcf.fn = vcf.fn,
  out.fn = raw_gds,
  method = "biallelic.only",
  ignore.chr.prefix = "chr",
  verbose = TRUE
)

cat("\n--- Calculating Individual Heterozygosity ---\n")

genofile <- snpgdsOpen(raw_gds)
all_samples <- read.gdsn(index.gdsn(genofile, "sample.id"))

# LD Pruning
snpset_init <- snpgdsLDpruning(
  gdsobj = genofile,
  autosome.only = FALSE,
  maf = 0.05,            
  missing.rate = 0.10,   
  method = "composite",
  ld.threshold = 0.8,
  verbose = TRUE,
  start.pos = "first"
)
pruned_snp_init <- unlist(unname(snpset_init))

# Calculate Inbreeding Coefficient (F)
inb_res_init <- snpgdsIndInb(
  genofile, 
  snp.id = pruned_snp_init, 
  method = "mom.weir", 
  verbose = FALSE
)

# Extract genotype matrix and allele frequencies
geno_mat_init <- snpgdsGetGeno(genofile, snp.id = pruned_snp_init, snpfirstdim = FALSE)
af_init <- snpgdsSNPRateFreq(genofile, snp.id = pruned_snp_init)$AlleleFreq

# Calculate He locus
He_locus <- 2 * af_init * (1 - af_init)

# Calculate individual Ho and He
Ho_init <- rowMeans(geno_mat_init == 1, na.rm = TRUE)
is_typed_init <- !is.na(geno_mat_init)
He_init <- as.numeric(is_typed_init %*% He_locus) / rowSums(is_typed_init)

# Export results to dataframe
df_ind <- data.frame(
  Sample = all_samples,
  F = inb_res_init$inbreeding,
  Ho = Ho_init,
  He = He_init
)

# Sort by inbreeding level and save
df_ind <- df_ind[order(df_ind$F), ]
write_xlsx(df_ind, path = output_excel)

snpgdsClose(genofile)

cat("\nResults successfully saved to:", output_excel, "\n")
