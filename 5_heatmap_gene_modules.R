library(gplots)

#### load data from previous steps ####
lnames = load(file = "./Routput/4-geneOverlap.RData");

dim(brca_gbm_gene_overlap)
dim(gbm_ov_gene_overlap)
dim(ov_brca_gene_overlap)

########### Generate heatmaps for each pairwise analysis of gene modules ###########
generateGeneModuleHeatmap <- function(pval_matrix, title, prefix)
{
  myCol <- c("gray15", "gray25", "blue", "green", "yellow", "orange", "gray25", "gray15")
  # Defining breaks for the color scale
  myBreaks <- c(-1, -0.06, -0.05, -0.001, 0, 0.001, 0.05, 0.06, 1)
  png(file = paste0(prefix, "_heatmap_gene_module_pairs.png"), width = 1200, height = 1000);
  par(cex.main = 0.8)
  heatmap.2(
    pval_matrix, scale = "none", Rowv = T, Colv = T,
    col = myCol,
    breaks = myBreaks,
    symkey = FALSE,
    main = title,
    key.title = "",
    margins = c(3, 3), cexRow = 1, cexCol = 1, key = FALSE, keysize = 1,
    trace = "none")
  legend("topleft", fill = myCol, cex = 1, xpd = TRUE, lty = 1, lwd = 1,
         legend = c(">0.6", "0.6 to 0.05", "0.049 to 0.001 (OR <1)", "0.001 to 0 (OR <1)", "0 to 0.001 (OR >1)",
                    "0.001 to 0.049 (OR >1)", "0.05 to 0.6", ">0.6"))
  dev.off()
}


brca_gbm_pval_matrix = brca_gbm_gene_overlap$modulePairMatrix
generateGeneModuleHeatmap(brca_gbm_pval_matrix, "Overlap of Gene Modules - BRCA and GBM", "./results/5_BRCA_GBM")

gbm_ov_pval_matrix = gbm_ov_gene_overlap$modulePairMatrix
generateGeneModuleHeatmap(gbm_ov_pval_matrix, "Overlap of Gene Modules - GBM and OV", "./results/5_GBM_OV")

ov_brca_pval_matrix = ov_brca_gene_overlap$modulePairMatrix
generateGeneModuleHeatmap(ov_brca_pval_matrix, "Overlap of Gene Modules - OV and BRCA", "./results/5_OV_BRCA")
