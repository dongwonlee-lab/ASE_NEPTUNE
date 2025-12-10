#' Calculate Cell Type Fraction of Gene
#'
#' This script calculates the cell type fraction of gene using single cell RNA-seq data.
#` Basically,these fractions represent the relative expression level per single cell of each type
#' Since we used glomerular (GLOM) and tubulointerstitial (TUBE) samples, we selected
#' cell types that typically exist in either compartment. We regrouped cell types from
#' KPMP, performed pseudobulk analysis, normalized by the total expression of that cell type,
#` and calculated the proportion of normalized expression contributed by each cell type.
#' 
#' Outputs:
#' - {compartment_type}_norm_exp.tsv: Normalized single-cell RNA expression dataframe
#'   with the first column containing gene names in HUGO ID format and remaining columns
#'   representing grouped cell types (compartment_type = GLOM or TUBE)
#' - {compartment_type}_cell_fraction.tsv: Cell fraction dataframe where the cell type
#'   fraction of a given gene is calculated as the proportion of its normalized expression
#'   within a cell type divided by the total normalized expression across all grouped cell
#'   types
#`
#' Usage:
#'   Rscript calculate_cell_fraction.R
#'
#' Dependencies:
#'   R ≥ 4.0.x with packages: Seurat, SeuratDisk, dplyr
#'
#' Author:
#'   Junmo Sung
#'   jsung5@bwh.harvard.edu
#'
#' Date:
#'   2025-10-03

############################### Load Library ###############################
library(Seurat)
library(SeuratDisk)
library(dplyr)

############################### Load Single Cell RNA-seq data ###############################
# Refernece: https://atlas.kpmp.org/explorer/dataviz
seurat.object <- LoadH5Seurat('./KPMP_PREMIERE_SC_version1.5_ForExplorer_withRC.032624.h5Seurat') # KPMP single cell data with 54 cell types

######################################## Functions to Use ########################################
# Aggregate by grouped cell types
aggregate_by_group <- function(counts, groups) {
  do.call(cbind, lapply(groups, function(celltypes) {
    if (is.character(celltypes)) {
      rowSums(counts[, colnames(counts) %in% celltypes, drop = FALSE])
    } else {
      rowSums(counts[, colnames(counts) %in% unlist(celltypes), drop = FALSE])
    }
  }))
}

# Normalize (divide by total counts per group and multiply by scale factor)
normalize_counts <- function(counts, scale_factor = 10000) {
  prop_counts <- sweep(counts, 2, colSums(counts), FUN = "/")
  prop_counts * scale_factor
}

######################################## Analysis Begin! ########################################
## Returns summed counts ("pseudobulk") for each identity class (https://satijalab.org/seurat/reference/aggregateexpression) ##
DefaultAssay(seurat.object) <- 'RNA'
aggregate_expression <- AggregateExpression(seurat.object, group.by = "celltype", assays = "RNA", slot='counts', return.seurat=TRUE)

## Assigning / regrouping cell types into each tissue compartment (cycEPI excluded) ##
# Assigning cell types to each kidney compartment
celltype_unique_g <- c("POD", 
                       "PEC", 
                       "EC-GC","MC")
celltype_unique_t <- c(
  "aPT", "aTAL1", "aTAL2", "ATL", "C-TAL", "CCD-PC", "CNT", "CNT-IC-A", "CNT-PC",
  "DCT1", "DCT2", "dPC", "dPT", "DTL1", "EC-AVR", "EC-LYM", "EC-PTC", "IC-A",
  "IC-B", "IMCD", "M-TAL", "PC", "PT-S1/S2", "PT-S3", "tPC-IC"
)
celltype_both <- c(
  "aFIB", "B", "cDC", "cycMNP", "cycT", "dFIB", "dVSMC", "EC-AEA", "FIB",
  "MAC-M2", "MAST", "MDC", "MON", "MYOF", "ncMON", "NK1", "NKT", "pDC", "PL",
  "REN", "T", "T-CYT", "VSMC", "VSMC/P"
)

# Grouping cells by grouped cell types
celltype_grouped_g <- list(
    POD = c("POD"),
    PEC = c("PEC"),

    EC = c("EC-AEA", "EC-GC"),
    vSMC.MC = c("dVSMC", "MC", "REN", "VSMC", "VSMC/P"),
    FIB = c("aFIB", "dFIB", "FIB", "MYOF"),
    IMMUNE = c("B", "cDC", "cycMNP", "cycT", "MAC-M2", "MAST", "MDC", "MON", "ncMON", "NK1", "NKT", "pDC", "PL", "T", "T-CYT") 
)
celltype_grouped_t <- list(
    PT = c("aPT", "dPT", "PT-S1/S2", "PT-S3"),
    DTL.ATL = c("ATL", "DTL1"),
    TAL = c("aTAL1", "aTAL2", "C-TAL", "M-TAL"),
    DT = c("DCT1", "DCT2", "CNT"),
    PC = c("CCD-PC", "CNT-PC", "dPC", "IMCD", "PC"),
    IC = c("CNT-IC-A", "IC-A", "IC-B", "tPC-IC"),
    
    EC = c("EC-AEA", "EC-AVR", "EC-LYM", "EC-PTC"),
    vSMC.MC = c("dVSMC", "REN", "VSMC", "VSMC/P"),
    FIB = c("aFIB", "dFIB", "FIB", "MYOF"),
    IMMUNE = c("B", "cDC", "cycMNP", "cycT", "MAC-M2", "MAST", "MDC", "MON", "ncMON", "NK1", "NKT", "pDC", "PL", "T", "T-CYT")
)

## Combine the unique and shared cell types for each tissue compartment ##
selected_g <- c(celltype_unique_g, celltype_both)
selected_t <- c(celltype_unique_t, celltype_both)

## Subset the aggregate_expression matrix with cell types for each tissue compartment ##
subset_aggregate_g <- aggregate_expression@assays$RNA$counts[, colnames(aggregate_expression@assays$RNA$counts) %in% selected_g]
subset_aggregate_t <- aggregate_expression@assays$RNA$counts[, colnames(aggregate_expression@assays$RNA$counts) %in% selected_t]

## Aggregate by grouped cell types ##
subset_aggregate_g_grouped <- aggregate_by_group(subset_aggregate_g, celltype_grouped_g)
subset_aggregate_t_grouped <- aggregate_by_group(subset_aggregate_t, celltype_grouped_t)

## Normalize ##
scale_factor <- 10000
norm_subset_aggregate_g <- normalize_counts(subset_aggregate_g_grouped, scale_factor)
norm_subset_aggregate_t <- normalize_counts(subset_aggregate_t_grouped, scale_factor)

## Save normalized expression dataframe ##
norm_subset_aggregate_g <- as.data.frame(norm_subset_aggregate_g)
norm_subset_aggregate_g$gene <- rownames(norm_subset_aggregate_g) 
norm_subset_aggregate_g <- norm_subset_aggregate_g %>%
  select(gene, everything())
norm_subset_aggregate_t <- as.data.frame(norm_subset_aggregate_t)
norm_subset_aggregate_t$gene <- rownames(norm_subset_aggregate_t) 
norm_subset_aggregate_t <- norm_subset_aggregate_t %>%
  select(gene, everything())
write.table(norm_subset_aggregate_g, "./g_norm_exp.tsv", sep='\t', row.names=FALSE)
write.table(norm_subset_aggregate_t, "./t_norm_exp.tsv", sep='\t', row.names=FALSE)

## Create cell-type fraction dataframe ##
g_cell_fraction <- sweep(
  norm_subset_aggregate_g[ , -1],
  1,
  rowSums(norm_subset_aggregate_g[ , -1]),
  FUN = "/"
)

t_cell_fraction <- sweep(
  norm_subset_aggregate_t[ , -1],
  1,
  rowSums(norm_subset_aggregate_t[ , -1]),
  FUN = "/"
)

# Replace NA with 0
g_cell_fraction[is.na(g_cell_fraction)] <- 0
t_cell_fraction[is.na(t_cell_fraction)] <- 0

# Add gene column back
g_cell_fraction <- cbind(gene = norm_subset_aggregate_g$gene, g_cell_fraction)
t_cell_fraction <- cbind(gene = norm_subset_aggregate_t$gene, t_cell_fraction)

# Save
write.table(g_cell_fraction, "./g_cell_fraction.tsv", sep='\t', row.names=FALSE)
write.table(t_cell_fraction, "./t_cell_fraction.tsv", sep='\t', row.names=FALSE)
