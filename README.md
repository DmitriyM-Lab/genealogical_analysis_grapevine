# genealogical_analysis_R_scripts
R pipeline for analysis of SNP data (heterozygosity, PCA, IBD, clustering) based on the SNPRelate package (Zheng et al., 2012)

To run this scripts, you will need R installed along with the following packages: BiocManager from 
## 01_heterozygosity_calc 
Converts VCF data into GDS format, performs LD pruning, and calculates individual observed (Ho) and expected (He) heterozygosity.
## 02_filter_gds: 
Cleans the GDS dataset based on the results from 01_heterozygosity_calc (e.g., removing samples with abnormal heterozygosity or high inbreeding). If any samples are removed, it automatically performs a new LD pruning for the remaining dataset.


