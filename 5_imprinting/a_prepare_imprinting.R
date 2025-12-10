library(data.table)
library(tidyverse)
library(R.matlab)
library(hdf5r)
library(rhdf5)

#' Generate Baran et al. Imprinting Pipeline Input Files
#'
#' This script processes TOGA and PHASER-derived allele-specific expression (ASE)
#' data to generate input files for the Baran et al. (2015) imprinting pipeline.
#' It filters for variants within resolved COV genes, applies haplotype blacklists,
#' annotates variants with gene names, and constructs per-sample ASE matrices.
#'
#' The resulting outputs are formatted for use in the MATLAB-workflow,
#' enabling downstream imprinting analyses.
#'
#' Outputs:
#' - `../data/output/imprinting/aux_files/`:
#'   - `GTd.snps.txt`: List of all non-imputed SNP-gene pairs.
#'   - `new_ensg2name.txt`: Mapping of Ensembl IDs to gene symbols.
#'   - `known_ana_geneimprint_su.txt`: Ensembl IDs of known imprinted genes.
#'   - `ind_gene_NMD.txt`, `imputed_samples_phase_error.txt`: Empty placeholder files.
#' - `../data/output/imprinting/data/`:
#'   - `1KG.txt.sub`: List of subject IDs.
#'   - `1KG.txt.rs`: List of variant identifiers (SNP_GENE format).
#'   - `1KG.txt.chr`: Chromosome numbers.
#'   - `1KG.txt.pos`: Variant positions.
#'   - `1KG.txt.hwe.mat`: Hardy-Weinberg equilibrium p-values (MATLAB .mat format).
#'   - `1KG.txt.mat`: HDF5 file containing total and reference count matrices.
#' - `../data/output/imprinting/1KG.txt.for_query`: SNP coordinate file for bcftools query.
#' - `../data/output/imprinting/1KG.txt.queried`: Queried VCF output with HWE statistics.
#'
#' Usage:
#'   Rscript create_imprinting_files.R
#'
#' Dependencies:
#'   R ≥ 4.0.0 with packages:
#'   - data.table
#'   - tidyverse
#'   - R.matlab
#'   - hdf5r
#'   - rhdf5
#'   - system tools: bedtools, bcftools
#'
#' Author:
#'   Eric Sakkas
#'   erotokritos.sakkas@childrens.harvard.edu
#'
#' Date:
#'   2025-10-06

subject_ids <- readLines("../resource/subject_ids.tsv")

message("FILTERING TO VARIANTS IN RESOLVED COV GENES")
cat("\tDATE, TIME :", format(Sys.time(), "%Y-%m-%d, %H:%M:%S"), "\n")

message("Loading resolved COV genes...")
resolved_cov <- read.table("../data/output/toga/stats/all_samples.cov_genes.txt", sep = "\t", header = TRUE) %>%
                                filter(TAG %in% c("R", "RC", "U")) %>%
                                mutate(SUBJECT_ID = sample,
                                final_total_hap_count = ifelse(is.na(adj_total_hap_count) | adj_total_hap_count == "", total_hap_count, adj_total_hap_count)) %>%
                                dplyr::select(name, gene_symbol, variants, variants_overlapping, final_total_hap_count, SUBJECT_ID) %>%
                                separate_rows(final_total_hap_count, sep = ":") %>%
                                filter(final_total_hap_count >= 20) %>%
                                dplyr::select(name, gene_symbol, variants, variants_overlapping, SUBJECT_ID)
message(paste0("\t", nrow(resolved_cov), "resolved COV genes identified."))

message("Loading variants with ASE data...")
SAFID_allelic_list <- list()
i=1
subject_ids <- readLines("../resource/subject_ids.tsv")
for (sample in subject_ids) {    
    message(paste("Loading", sample, "allelic data..."))
    file_path <- paste0("../data/output/phaser/", sample, "/phaser_run/", sample, ".allelic_counts.txt")
    if(file.exists(file_path)) {
        SAFID_allelic_counts <- read.table(file_path, sep = "\t", header = TRUE) %>% mutate(SUBJECT_ID = sample)

        SAFID_allelic_counts <- SAFID_allelic_counts %>%
            rename(
                CHROM = contig,
                POS = position,
                SNP = variantID,
                REF = refAllele,
                ALT = altAllele,
                RNA_REF_COUNT = refCount,
                RNA_ALT_COUNT = altCount,
                RNA_TOTAL_COUNT = totalCount
                )

        SAFID_allelic_list[[i]] <- SAFID_allelic_counts
        i <- i+1
    } else {
        message(paste("File not found:", file_path))
    }
}
SAFID_allelic_counts <- do.call("rbind", SAFID_allelic_list)

message("Removing variants in haplotype blacklist...")
tempfile <- tempfile()
write.table(SAFID_allelic_counts %>% 
                ungroup() %>% 
                mutate(POS2 = POS - 1) %>% 
                dplyr::select(CHROM, POS2, POS), tempfile, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
bedtools_output <- system(paste0("bedtools intersect -a ",tempfile," -b ../resource/hg38_haplo_count_blacklist.chr.bed  -v -wa"), intern = TRUE)
intersected_df <- read.table(text = bedtools_output, header = FALSE, sep = "\t")
unlink(tempfile)

SAFID_allelic_counts <- SAFID_allelic_counts %>%
    inner_join(intersected_df %>% distinct(), by = c("CHROM" = "V1", "POS" = "V3"))

message("Adding gene names...")
tempfile <- tempfile()
write.table(SAFID_allelic_counts %>% mutate(POS2 = POS - 1) %>% dplyr::select(CHROM, POS2, POS), tempfile, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
bedtools_output <- system(paste0("bedtools intersect -a ",tempfile," -b ../resource/gencode.v.45.gene_coordinates.bed -wa -wb"), intern = TRUE)
intersected_df <- read.table(text = bedtools_output, header = FALSE, sep = "\t")
unlink(tempfile)

snp_genes <- intersected_df %>%
    dplyr::select(V1,V3,V7,V8,V5,V6)

SAFID_allelic_counts <- SAFID_allelic_counts %>%
    dplyr::select(-V2) %>%
    left_join(snp_genes, by = c("CHROM" = "V1", "POS" = "V3")) %>%
    rename("ENSEMBL" = "V7", "START" = "V5", "STOP" = "V6") %>%
    dplyr::select(-V8) %>%
    distinct()

message("Keeping only resolved COV genes...")
resolved_cov <- resolved_cov %>%
    mutate(ENSEMBL = sub("\\..*","", name))
SAFID_allelic_counts <- SAFID_allelic_counts %>%
    mutate(ENSEMBL = sub("\\..*","", ENSEMBL))

SAFID_allelic_counts_cov20 <- SAFID_allelic_counts %>%
    inner_join(resolved_cov, by = c("SUBJECT_ID", "ENSEMBL")) %>%
    dplyr::select(-variants, -variants_overlapping) %>%
    distinct()

message("GENERATING MATLAB FILES FOR BARAN ET AL PIPELINE")
cat("\tDATE, TIME :", format(Sys.time(), "%Y-%m-%d, %H:%M:%S"), "\n")

SAFID_allelic_counts_cov20$SNP_GENE <- paste(SAFID_allelic_counts_cov20$SNP, SAFID_allelic_counts_cov20$ENSEMBL, sep = "_")
SAFID_allelic_counts_cov20 <- subset(SAFID_allelic_counts_cov20, !is.na(ENSEMBL))
SAFID_allelic_counts_cov20 <- subset(SAFID_allelic_counts_cov20, RNA_TOTAL_COUNT >= 8)

outfile <- "../data/output/imprinting/all_samples.allelic_counts.cov8.tsv"
dir.create(dirname(outfile), recursive = TRUE)
write.table(SAFID_allelic_counts_cov20, outfile, sep = "\t", row.names = FALSE, col.names = TRUE)

message("Creating list of non-imputed variants:")
# In our pipeline, all imputed variants are removed by blacklist, so we include all vars:
non_imp_SNPs <- SAFID_allelic_counts_cov20 %>%
    dplyr::select(SNP_GENE) %>%
    arrange(SNP_GENE) %>%
    distinct() 
outfile <- "../data/output/imprinting/aux_files/GTd.snps.txt"
dir.create(dirname(outfile), recursive = TRUE)
write.table(non_imp_SNPs, outfile, quote = FALSE, col.names = FALSE, row.names = FALSE)
message(paste0("\t", nrow(non_imp_SNPs), " variants written to ", outfile))

message("Creating list of samples:")
subjs <- SAFID_allelic_counts_cov20 %>%
    dplyr::select(SUBJECT_ID)
message(paste("Contains", n_distinct(subjs$SUBJECT_ID), "samples."))
outfile <- "../data/output/imprinting/data/1KG.txt.sub"
dir.create(dirname(outfile), recursive = TRUE)
write.table(subjs, outfile, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
message(paste("Samples written to", outfile))

message("Creating list of variants:")
snps <- SAFID_allelic_counts_cov20 %>%
    dplyr::select(SNP_GENE)
outfile <- "../data/output/imprinting/data/1KG.txt.rs"
write.table(snps, outfile, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
message(paste("Variants written to", outfile))

message("Creating list of chr and pos:")
chr <- SAFID_allelic_counts_cov20 %>%
    dplyr::select(CHROM) %>%
    mutate(CHROM = sub("chr", "", CHROM))
pos <- SAFID_allelic_counts_cov20 %>%
    dplyr::select(POS)
outfile <- "../data/output/imprinting/data/1KG.txt.chr"
write.table(chr, outfile, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
message(paste("Chromosomes written to", outfile))
outfile <- "../data/output/imprinting/data/1KG.txt.pos"
write.table(pos, outfile, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
message(paste("Positions written to", outfile))

message("Calculating HWE probabilities:")
# Baran et al removes vars with HWE p < 0.001
SNPs_sep <- SAFID_allelic_counts_cov20 %>%
    dplyr::select(SNP, ENSEMBL) %>%
    distinct() %>%
    separate(SNP, into = c("CHROM", "POS", "REF", "ALT", "b38"), sep = "_") %>%
    arrange(CHROM, POS)
outfile <- "../data/output/imprinting/1KG.txt.for_query"
write.table(SNPs_sep, outfile, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

# Any sample's unphased VCF can be used
system("bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%HWE]' -R ../data/output/imprinting/1KG.txt.for_query ../data/HG00171/HG00171.unphased.vcf.gz -o ../data/output/imprinting/1KG.txt.queried")

hwe_df <- read.table("../data/output/imprinting/1KG.txt.queried", sep = "\t") %>% distinct()

hwep <- SNPs_sep %>%
    mutate(POS = as.integer(POS)) %>%
    left_join(hwe_df, by = c("CHROM" = "V1", "POS" = "V2", "REF" = "V3", "ALT" = "V4")) %>%
    pull(V5)

outfile <- "../data/output/imprinting/data/1KG.txt.hwe.mat"
dir.create(dirname(outfile), recursive = TRUE)
writeMat(outfile, hwep = hwep)
message(paste("HWE p-values written to", outfile))

message("Creating MATLAB matrices...")
# Matrix of total RNA count: variant x samples
for_n <- SAFID_allelic_counts_cov20 %>%
    mutate(SNP_GENE = paste(SNP, ENSEMBL, sep = "_")) %>%
    dplyr::select(SNP_GENE, SUBJECT_ID, RNA_TOTAL_COUNT) %>%
    distinct()

n <- for_n %>%
    pivot_wider(names_from = SNP_GENE, values_from = RNA_TOTAL_COUNT) %>%
    column_to_rownames(var = "SUBJECT_ID")

n_matrix <- as.matrix(n)

# Matrix of RNA ref count: variant x samples
ref <- SAFID_allelic_counts_cov20 %>%
    mutate(SNP_GENE = paste(SNP, ENSEMBL, sep = "_")) %>%
    dplyr::select(SNP_GENE, SUBJECT_ID, RNA_REF_COUNT) %>%
    distinct() %>%
    pivot_wider(names_from = SNP_GENE, values_from = RNA_REF_COUNT) %>%
    column_to_rownames(var = "SUBJECT_ID")

ref_matrix <- as.matrix(ref)

# genenames for variants (rows) in matrix
colnames_genes <- data.frame(SNP_GENE = colnames(ref)) %>% 
    left_join(SAFID_allelic_counts_cov20 %>% mutate(SNP_GENE = paste(SNP, ENSEMBL, sep = "_")) %>% dplyr::select(SNP_GENE, ENSEMBL),
             by = c("SNP_GENE")) %>%
    distinct() %>%
    pull(ENSEMBL)

h5file <- "../data/output/imprinting/data/1KG.txt.mat"
dir.create(dirname(h5file), recursive = TRUE)
h5createFile(h5file)
h5write(n_matrix, h5file, "n")
h5write(ref_matrix, h5file, "ref")
h5write(colnames_genes, h5file, "sgen")
message(paste("MATLAB matrices written to", h5file))

# Tucci+Baran et al imprinted genes corrected
imprinted_genelist <- c("MIR298", "MIR296", "MIR483", 'MIR675','H19','GNAS','DIRAS3', 'IL12RB2', 'RNU5D-1', 'AGO1', 'UTS2', 'THAP3', 'CACNA1E', 'CYP2J2', 'ACOT11', 'LINC00467', 'LRRTM1', 'TMEM247', 'THUMPD2', 'PAX8', 'PAX8-AS1', 'DNAH7', 'ICA1L', 'CMKLR2', 'CMKLR2-AS', 'ZDBF2', 'MRPL44', 'SPHKAP', 'USP4', 'SLC4A7', 'ZNF385D', 'EFCC1', 'RAB7A', 'MCCC1', 'FGF12', 'MELTF', 'GPR78', 'STX18-AS1', 'PDE6B', 'SH3BP2', 'NAP1L5', 'GRID2', 'SFRP2', 'FAM149A', 'FRG1', 'RHOBTB3', 'NUDT12', 'VTRNA2-1', 'ZNF354C', 'CUL7', 'MDGA1', 'MOCS1', 'C6orf47', 'RNF144B', 'CD83', 'FAM50B', 'CRYBG1', 'LIN28B', 'PHACTR2', 'HYMAI', 'PLAGL1', 'SLC22A2', 'SLC22A3', 'PLG', 'KIF25', 'GRB10', 'RAPGEF5', 'SCIN', 'THSD7A', 'CALCR', 'TFPI2', 'SGCE', 'PEG10', 'PDK4', 'CPA4', 'MEST', 'MESTIT1', 'COPG2IT1', 'KLF14', 'KLHDC10', 'AGBL3', 'PRKAG2', 'PTK2B', 'R3HCC1', 'CLDN23', 'DLGAP2', 'PKIA', 'ZFAT', 'ZFAT-AS1', 'PEG13', 'PSCA', 'NAPRT', 'TRAPPC9', 'KCNK9', 'DENND3', 'GLIS3', 'PGM5P3-AS1', 'EXD3', 'PTCHD3', 'ITGA8', 'PROSER2', 'PROSER2-AS1', 'JMJD1C', 'AIFM2', 'ATP5MK', 'VWA2', 'INPP5F', 'CPXM2', 'ACCS', 'ALKBH3', 'MAPK8IP1', 'WT1-AS', 'LINC00294',  'IGF2', 'IGF2-AS', 'INS', 'KCNQ1', 'KCNQ1OT1', 'KCNQ1DN', 'CDKN1C', 'PHLDA2',  'SLC22A18', 'ZNF215', 'NAV2', 'ART5', 'OVCH2', 'RNF141', 'IRF7', 'ANO1', 'PAK1', 'VSTM5', 'ZC3H12C', 'SPA17', 'NTM', 'OPCML', 'TIGAR', 'CACNA1C', 'WIF1', 'N4BP2L1', 'RB1', 'LPAR6', 'DLEU7', 'KLHL1', 'FGF14', 'PCK2', 'PAPLN-AS1', 'DLK1', 'MEG3', 'MIR337', 'RTL1', 'MEG8', 'MIR134', 'MKRN3', 'MAGEL2', 'NDN', 'NPAP1', 'SNURF', 'SNRPN', 'SNORD107', 'SNORD64', 'SNORD108', 'SNORD109A',  'SNORD109B', 'UBE3A', 'PWRN1', 'SNHG14', 'RYR3', 'DNM1P35', 'RASGRF1', 'FAM174B', 'IRAIN', 'LRRK1', 'SIAH1', 'ZNF597',  'PDPR', 'ZFP90', 'CLEC3A', 'NLGN2', 'SEPTIN4', 'ZNF714', 'AXL', 'DNMT1', 'S1PR2', 'ICAM1', 'FDX2', 'ZNF833P', 'GNG7', 'ANO8', 'CACNA1A', 'ZNF331',  'MIR512-1', 'ZIM2', 'PEG3', 'MIMT1', 'ZNF542P', 'CST1', 'MCTS2', 'ACTL10', 'NNAT', 'BLCAP', 'ZHX3', 'L3MBTL1', 'SGK2', 'CYP24A1', 'GNAS-AS1',  'GDAP1L1', 'PRMT2', 'CBR1', 'TPTEP1', 'ARVCF', 'CACNA1I', 'SNU13', 'SLC9A7')
snord116_genes <- paste0("SNORD116-", 1:30)
snord115_genes <- paste0("SNORD115-", 1:48)
baran_gene_list <- c("CPA4", "CST1", "DIRAS3", "DLK1", "FAM50B", "GRB10", "H19", "IGF2-AS", "IGF2", "INPP5F", "KCNQ1", "KIF25", "L3MBTL1", "LPAR6", "MAGEL2", "MAGI2", "MEG3", "MEG8", "MEG9", "MEST", "NAP1L5", "NDN", "NTM", "PEG10", "PEG3", "PLAGL1", "PPIEL", "PRSS50", "PWRN1", "LINC01629", "SGK2", "SNHG14", "SNRPN", "SNURF", "SYCE1", "SPMAP2L", "UBE3A", "UGT2B4", "UTS2", "ZDBF2", "ZNF331", "ZNF597")

imprinted_genelist <- c(imprinted_genelist, snord116_genes, snord115_genes, baran_gene_list)
ana_imprinted_genes <- data.frame(V5 = imprinted_genelist) %>% distinct()

ensb_to_gene_all <- read.table("../resource/gencode.v.45.gene_coordinates.bed", sep = "\t", header = FALSE)
ensb_to_gene <- ensb_to_gene_all %>% mutate(ENSEMBL = sub("\\..*", "", V4)) %>% dplyr::select(ENSEMBL, V5)

outfile <- "../data/output/imprinting/aux_files/new_ensg2name.txt"
dir.create(dirname(outfile), recursive = TRUE)
write.table(ensb_to_gene, outfile, quote = FALSE, 
            row.names = FALSE, col.names = FALSE, sep = ",")

ana_impgene_df <- ana_imprinted_genes %>%
    left_join(ensb_to_gene, by = c("V5")) %>% 
    distinct() %>%
    dplyr::select(ENSEMBL) %>%
    distinct() %>%
    filter(!is.na(ENSEMBL))

outfile <- "../data/output/imprinting/aux_files/known_ana_geneimprint_su.txt"
dir.create(dirname(outfile), recursive = TRUE)
write.table(ana_impgene_df, outfile, 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
message(paste("Known Imprinted genes written to", outfile))

# Miscellaneous empty files
vars_with_nmd_outfile <- "../data/output/imprinting/aux_files/ind_gene_NMD.txt" 
dir.create(dirname(vars_with_nmd_outfile), recursive = TRUE)
write.table(data.frame(), vars_with_nmd_outfile, row.names = FALSE, col.names = FALSE, quote = FALSE)

unphased_indiv_outfile <- "../data/output/imprinting/aux_files/imputed_samples_phase_error.txt" 
dir.create(dirname(unphased_indiv_outfile), recursive = TRUE)
write.table(data.frame(), unphased_indiv_outfile, row.names = FALSE, col.names = FALSE, quote = FALSE)

message("Done!")