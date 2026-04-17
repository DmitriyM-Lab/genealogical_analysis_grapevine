# Genealogical Analysis of Grapevine Cultivars R scripts
R pipeline for analysis of SNP data (heterozygosity, PCA, IBD, clustering) based on the SNPRelate package (Zheng et al., 2012).
#### To run this scripts, you will need R installed along with the following packages:
#### install.packages("BiocManager")
#### BiocManager::install("SNPRelate")
#### install.packages("writexl")
## 01_heterozygosity_calc 
Converts VCF data into GDS format, performs LD pruning, and calculates individual observed (Ho) and expected (He) heterozygosity.
## 02_filter_gds
Cleans the GDS dataset based on the results from 01_heterozygosity_calc (e.g., removing samples with abnormal heterozygosity or high inbreeding). If any samples are removed, it automatically performs a new LD pruning for the remaining dataset.
## 03_pca_analysis
Performs Principal Component Analysis (PCA) on the cleaned GDS dataset.
## 04_ibd_kinship
Calculates Identity-By-Descent (IBD) and kinship coefficients using the KING-robust method. It automatically classifies the relationships between sample pairs (e.g., Clones, Parent-offspring, Full sibs, Unrelated) based on theoretical thresholds (Manichaikul et al. 2010).
## 05_ibs_dendrogram
Generates a UPGMA hierarchical clustering dendrogram based on Identity-by-State (IBS) genetic distances. It visually overlays the kinship relationships calculated in the previous step using colored arrows (e.g., clones, parent-offspring).
## 06_connection_analysis
Parses the IBD relationship dataset to count the number of valid genetic connections (edges) for each cultivar. It visually represents this data as a custom stacked bar chart using ggplot2, displaying total connections on the left and the proportion of specific relationship types (clones, full sibs, etc.) on the right.

### References

Manichaikul, A., Mychaleckyj, J. C., Rich, S. S., Daly, K., Sale, M., & Chen, W.-M. (2010). Robust relationship inference in genome-wide association studies. Bioinformatics, 26(22), 2867–2873. https://doi.org/10.1093/bioinformatics/btq559

Zheng, X., Levine, D., Shen, J., Gogarten, S. M., Laurie, C., & Weir, B. S. (2012). A high-performance computing toolset for relatedness and principal component analysis of SNP data. Bioinformatics, 28(24), 3326–3328. https://doi.org/10.1093/bioinformatics/bts606
