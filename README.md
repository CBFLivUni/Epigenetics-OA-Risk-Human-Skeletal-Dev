# Differential methylation and mQTL analysis of cartilage tissue in the developing knee joint
<p align="center">
<img src="https://github.com/CBFLivUni/SarahRice_mQTL/assets/8311721/7200daa9-99a6-4db5-8f9a-6c48971f405b" width=70%>
</p>


## Workflow
This repository contains the code required to reproduce the bioinformatics analysis of the methylation and genotyping arrays.

**01_Methyl_Preprocess.Rmd**  
Quality control and normalisation of the EPIC methylation arrays

**02_Methyl_EDA.Rmd**  
Principal component analysis of the methylation data

**03_Methyl_Differential.Rmd**  
Differential methylated site and region analysis

**04_DMR_mFuzzClustering.Rmd**  
Mfuzz clustering of CpGs and enrichment analysis

**05_processROADMAP.Rmd**  
Overlap of DMRs with ROADMAP chondrocyte chromatin states and ATAC-seq data

**06_DMR_TFMotifs.qmd**  
Enrichment anaylsis of transcription factor motifs in DMRs

**07_genotype_arrays.Rmd**  
Genotyping array processing (using bash scripts), quality control and identification of GWAS OA SNP LD proxies

**08_matrixQTL_prep.Rmd**  
Preparation of methylation and genotyping data for mQTL analysis

**09_matrixQTL_analysis.Rmd**  
Identification and plotting of mQTLs

**10_mQTLComparisons.qmd**  
Comparison of the mQTLs against existing mQTLs from foetal brain and comparison of CpGs in the developmental DMRs vs the mQTLs.

**11_Colocalisation.qmd**  
Colocalisation of the mQTL SNPs with osteoarthritis GWAS data



All `Rmd` and `qmd` files can be rendered within RStudio.

## Installation
Install the needed R packages
```
RScript install/installRLibs.R
```
