#' Perform gene set enrichment analysis

#' Reference: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
#'
#' This script identifies gene or protein classes that are significantly overrepresented 
#` within a large dataset and potentially associated with specific phenotypes.
#' 
#' Outputs:
#' - gsea_{gene_id}_{type}_res:
#'   - ID: The unique identifier of the gene set
#'   - Description: A short name or definition of the gene set
#'   - setSize: The total number of genes in this gene set that are present in your input ranked gene list
#'   - enrichmentScore (ES): The running-sum statistic from the GSEA algorithm. It reflects the degree to which this gene set is overrepresented at the top or bottom of your ranked gene list. Positive ES means enrichment at the top (upregulated genes), negative means enrichment at the bottom (downregulated)
#'   - NES: Normalized enrichment score
#'   - pvalue: The nominal p-value from permutation testing
#'   - p.adjust: Adjusted p-value by Benjamini-Hochberg FDR correction to account for multiple testing
#'   - rank: The position (rank) in the input ranked gene list where the enrichment score reaches its maximum
#'   - leading_edge: A summary of the subset of genes contributing most to the ES
#'   - core_enrichment: The list of genes that form the leading edge subset
#'
#' Usage:
#'   Rscript b_GSEA.R
#'
#' Dependencies:
#'   R ≥ 4.0.x with packages: clusterProfiler, org.Hs.eg.db
#'
#' Author:
#'   Junmo Sung
#'   jsung5@bwh.harvard.edu
#'
#' Date:
#'   2025-10-02

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# Define variables
gene_id <- "ensembl_id" # ensembl_id or gene_symbol
gsea_gene_id <- "ENSEMBL" # SYMBOL or ENSEMBL or ENTREZID 
type <- "ALL" # "MF" for Molecular Function, "BP" for Biological Process, "CC" for Cellular Component , ALL(MF + BP + CC), KEGG
res_dir='./toy_example'

# Read the DESeq2 results
deseq2_results <- read.csv('deg_result.csv') # Deseq2 result as an input
deseq2_results <- na.omit(deseq2_results)
deseq2_results$ensembl_id <- sapply(strsplit(deseq2_results$X, "_"), function(x) strsplit(x[1], "\\.")[[1]][1])
deseq2_results$gene_symbol <- sub("^[^_]+_", "", deseq2_results$X)

deseq2_results$signed_pvalue<-sign(deseq2_results$log2FoldChange) * -log10(deseq2_results$pvalue) # Rank genes by p-value
gene_list <- deseq2_results$signed_pvalue
names(gene_list) <- deseq2_results[[gene_id]]
gene_list <- na.omit(gene_list)
gene_list = sort(gene_list, decreasing = TRUE) 
cat("gene_list length:", length(gene_list), "\n")
                                    
# Perform GSEA
gsea_result <- gseGO(geneList = gene_list, 
                     OrgDb = org.Hs.eg.db,
                     keyType = gsea_gene_id,
                     ont = type,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     verbose = FALSE)

# Convert GSEA results to a data frame and save as CSV
gsea_df <- as.data.frame(gsea_result)
gsea_df <- gsea_df[order(gsea_df$NES,decreasing=TRUE), ]
output_file <- paste0(res_dir,"/gsea_", gsea_gene_id, "_", type, "_res.csv")
write.csv(gsea_df, output_file, row.names = FALSE)