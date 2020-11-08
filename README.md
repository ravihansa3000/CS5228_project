# Exploring Gene Expression Data Using Weighted Gene Co-expression Network Analysis Across Multiple Cancer Types

## CS5228 Knowledge Discovery and Data Mining
A.V. AKILA RAVIHANSA PERERA – A0212216X

TRAN KHANH HUNG – A0212253W


## Method

- WGCNA

Weighted gene correlation network analysis (WGCNA) is a powerful method that uses a topological overlap module approach 
for constructing co-expression networks based on gene expression data. This method involves reconstructing 
gene co-expression modules and summarizing modules using module eigengenes (ME) and intramodular hub genes.

- Gene Enrichment and Pathway Analysis

Biologically interesting modules were identified using Fisher's exact test. The overlapping and union sets of 
genes from theses interesting gene module pairs were subjected to Gene Set Enrichment Analysis using topGO package).

## Dataset

Three gene expression datasets for three cancer types (GBM, OV, BRCA) were selected from TCGA (The Cancer Genome Atlas).

 - Glioblastoma Multiforme (GBM) gene expression by RNAseq
 https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/HiSeqV2_PANCAN.gz
 
 - Ovarian Serous Cystadenocarcinoma (OV) gene expression by RNAseq
 https://tcga.xenahubs.net/download/TCGA.OV.sampleMap/HiSeqV2_PANCAN.gz
 
 - Breast Invasive Carcinoma (BRCA) gene expression by RNAseq
https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2_PANCAN.gz


## Analysis

###  Data Exploration 

Scale Free Topology Model       |  Mean Connectivity      | Selected Soft Threshold
:-------------------------:|:-------------------------:|:-------------------------:
![](results/2_SFTM_Fit_GBM.png)  |  ![](results/2_Mean_Connectivity_GBM.png) | 9
![](results/2_SFTM_Fit_OV.png)  |  ![](results/2_Mean_Connectivity_OV.png) | 6
![](results/2_SFTM_Fit_BRCA.png)  |  ![](results/2_Mean_Connectivity_BRCA.png) | 12


### Clustering Tree

Genes       |      Module Eigengenes
:-------------------------:|:-------------------------:
![GBM](results/2_Clustering_Tree_Genes_GBM.png) | ![GBM](results/2_Clustering_Tree_ME_GBM.png)
![GBM](results/2_Clustering_Tree_Genes_OV.png) | ![GBM](results/2_Clustering_Tree_ME_OV.png)
![GBM](results/2_Clustering_Tree_Genes_BRCA.png) | ![GBM](results/2_Clustering_Tree_ME_BRCA.png)

<br/><br/>

### Gene Expression Network

Network Heatmap       |      Eigengene Adjacency Heatmap      |      Eigengene Dendrogram
:-------------------------:|:-------------------------:|:-------------------------:
![](results/3_Network_heatmap_GBM.png) |  ![](results/3_Eigengene_heatmap_GBM.png)  |  ![](results/3_Eigengene_dendrogram_GBM.png)
![](results/3_Network_heatmap_OV.png) |  ![](results/3_Eigengene_heatmap_OV.png)  |  ![](results/3_Eigengene_dendrogram_OV.png)
![](results/3_Network_heatmap_BRCA.png) |  ![](results/3_Eigengene_heatmap_BRCA.png) |  ![](results/3_Eigengene_dendrogram_BRCA.png)


### Pairwise Analysis of Gene Modules

## Overview

Item                     | BRCA and GBM               |  GBM and OV               | OV and BRCA
:-----------------------:|:-------------------------:|:-------------------------:|:-------------------------:
|||
--- Intersection --- |||
Name | [black_lightcyan_genes.txt](results/4_BRCA_GBM_lowP_Intersection_P0_black_lightcyan_genes.txt) | [lightcyan_red_genes.txt](results/4_GBM_OV_lowP_Intersection_P0_lightcyan_red_genes.txt) | [grey_grey_genes.txt](results/4_OV_BRCA_lowP_Intersection_P0_grey_grey_genes.txt) 
Lowest p-value | 0 | 0 | 0
Gene count in the most significant (lowest pvalue) module pair | 557 | 864 | 1577
topGO plot | [topGO.pdf](results/4_BRCA_GBM_lowP_Intersection_black_lightcyan_topGOPlot_fullnames.pdf) | [topGO.pdf](results/4_GBM_OV_lowP_Intersection_lightcyan_red_topGOPlot_fullnames.pdf) | [topGO.pdf](results/4_OV_BRCA_lowP_Intersection_grey_grey_topGOPlot_fullnames.pdf)
topGO analysis | [topGO.csv](results/4_BRCA_GBM_lowP_Intersection_black_lightcyan_summary_topGO_analysis.csv) | [topGO.csv](results/4_GBM_OV_lowP_Intersection_lightcyan_red_summary_topGO_analysis.csv) | [topGO.csv](results/4_OV_BRCA_lowP_Intersection_grey_grey_summary_topGO_analysis.csv)
|||
|||
--- Union --- |||
Name | [black_lightcyan_genes.txt](results/4_BRCA_GBM_lowP_Union_P0_black_lightcyan_genes.txt) | [lightcyan_red_genes.txt](results/4_GBM_OV_lowP_Union_P0_lightcyan_red_genes.txt) | [grey_grey_genes.txt](results/4_OV_BRCA_lowP_Union_P0_grey_grey_genes.txt)
Lowest p-value | 0 | 0 | 0
Gene count in top module pair | 1546 | 3377 | 6766
topGO plot | [topGP.pdf](results/4_BRCA_GBM_lowP_Union_black_lightcyan_topGOPlot_fullnames.pdf) | [topGO.pdf](results/4_GBM_OV_lowP_Union_lightcyan_red_topGOPlot_fullnames.pdf) | [topGO.pdf](results/results/4_OV_BRCA_lowP_Union_grey_grey_topGOPlot_fullnames.pdf)
topGO analysis | [topGO.csv](results/4_BRCA_GBM_lowP_Union_black_lightcyan_summary_topGO_analysis.csv) | [topGO.csv](results/4_GBM_OV_lowP_Union_lightcyan_red_summary_topGO_analysis.csv) | [topGO.csv](results/4_OV_BRCA_lowP_Union_grey_grey_summary_topGO_analysis.csv)


## Heatmaps of Overlap (Across Cancer Types)

![BRCA_GBM](results/5_BRCA_GBM_heatmap_gene_module_pairs.png)
<br/><br/>

![GBM_OV](results/5_GBM_OV_heatmap_gene_module_pairs.png)
<br/><br/>

![OV_BRCA](results/5_OV_BRCA_heatmap_gene_module_pairs.png)
<br/><br/>


## Project Structure and Run Instructions

 - Download and extract datasets to `./data` directory
 - Install dependencies (R packages)
    - Run `0_install_dependencies.R`
    
 - Gene filtering
    - Run `1_wgcna_cluster.R`
    
 - Build gene expression network and identify gene modules
    - Run `2_module_detection.R`
    
 - Generate gene network plots 
    - Run `3_network_visualization.R`

 - Gene Enrichment Analysis
    - Run `4_gene_enrichment.R`
    
