#' Generate differentially expressed genes 

#' Reference: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
#'
#' This script uses DESeq2 to identify differentially expressed genes 
#` between different experiment conditions using RNA-seq data. 
#' 
#' Outputs:
#' - deg_result.csv:
#'   - baseMean: The mean normalized expression of the gene across all samples (control + treatment)
#'   - log2FoldChange: The log₂ fold change in expression between the two conditions
#'   - lfcSE: Standard error of the log₂ fold change estimate
#'   - stat: Wald test statistic for the hypothesis that log₂FC ≠ 0
#'   - pvalue: Raw p-value for the Wald test statistic
#'   - padj: Adjusted p-value (FDR) using Benjamini–Hochberg
#'
#' Usage:
#'   Rscript a_DEG.R
#'
#' Dependencies:
#'   R ≥ 3.5.x with packages: DESeq2
#'
#' Author:
#'   Junmo Sung
#'   jsung5@bwh.harvard.edu
#'
#' Date:
#'   2025-10-02

# Load library
library(DESeq2)

# Load data
countData <- read.csv('./EXAMPLE_DEG/count.csv', row.names=1, check.names=FALSE) # Change the example data to your data.
colData <- read.csv('./EXAMPLE_DEG/meta.csv', row.names=1) # Change the example data to your data
countData <- round(countData)

colData$sex <- factor(colData$sex)
colData$race <- factor(colData$race)
colData$condition <- factor(colData$condition)
colData$batch <- factor(colData$batch)
colData$age <- as.numeric(colData$age)
colData$age <- scale(colData$age)  # Center and scale the age variable

# Parse data to run deseq2 with covariates added
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ age + sex + race + batch + condition)

# 10 counts in at least 3 samples as desired.
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Run differential gene expression analysis using DESeq2
dds <- DESeq(dds)
res <- results(dds)
summary(res)
write.csv(as.data.frame(res), file = "deg_result.csv")