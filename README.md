# Genealogical Analysis of Grapevine Cultivars
## Bioinformatics Pipeline for Variant Calling
This is the bioinformatics pipeline used for mapping reads to a reference genome, post-processing, variant calling, and subsequent filtration.
## 1. Reference Genome Indexing
`bwa-mem2 index REFERENCE.fasta`
`samtools faidx REFERENCE.fasta`
## 2. Read Mapping and Post-processing
### 2.1. Read Mapping
`bwa-mem2 mem -M -t <THREADS> -R '@RG\tID:<SAMPLE_ID>\tSM:<SAMPLE_NAME>\tPL:ILLUMINA' REFERENCE.fasta SAMPLE_R1_clean.fastq.gz SAMPLE_R2_clean.fastq.gz | samtools view -@ <THREADS> -b > SAMPLE.bam`
### 2.2. BAM File Preparation and Sorting
`samtools collate -@ <THREADS> -Ou SAMPLE.bam | samtools fixmate -@ <THREADS> -m - - | samtools sort -@ <THREADS> -o SAMPLE_sorted.bam -`
### 2.3. PCR Duplicate Marking
`samtools markdup -@ <THREADS> -d 100 SAMPLE_sorted.bam SAMPLE_markdup.bam`
### 2.4. Indexing and Mapping Statistics
`samtools index -@ <THREADS> -b SAMPLE_markdup.bam`
`samtools flagstat -@ <THREADS> SAMPLE_markdup.bam > SAMPLE_mapping_stats.txt`
## 3. Variant Calling
`bcftools mpileup --threads <THREADS> -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD -f REFERENCE.fasta -R TARGET_PANEL_REGIONS.txt *_markdup.bam | bcftools call --threads <THREADS> -f GQ,GP -m -Oz -o ALL_CULTIVARS_RAW.vcf.gz`
## 4. Variant Filtering
`bcftools filter -Ou -e "INFO/DP>2.5*AVG(INFO/DP) || MQ<30 || QUAL<30 || SP>10" ALL_CULTIVARS_RAW.vcf.gz | bcftools view -Ou -M 2 -V indels | bcftools filter -Ou -S . -e "FMT/GQ<20 || FMT/DP<5" | bcftools filter -Oz -e "F_MISSING>0.2" > ALL_CULTIVARS_FILTERED.vcf.gz`
`bcftools stats -s - ALL_CULTIVARS_FILTERED.vcf.gz > ALL_CULTIVARS_FILTERED_STATS.txt`





An interactive R-based workflow for analyzing SNP profiles of grapevine cultivars. Built on top of the SNPRelate package (Zheng et al., 2012), this pipeline covers:
* VCF to GDS conversion & LD pruning
* Individual Heterozygosity calculation (Ho & He)
* Principal Component Analysis (PCA)
* IBD & Kinship estimation (KING-robust)
* UPGMA Clustering & Network visualization

To run this scripts, you will need R installed along with the following packages:

`install.packages("BiocManager")`

`BiocManager::install("SNPRelate")`

`install.packages("writexl")`
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
