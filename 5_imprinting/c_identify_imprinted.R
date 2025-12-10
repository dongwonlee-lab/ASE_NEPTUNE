library(data.table)
library(tidyverse)
library(arrow)

#' Identify Imprinted Genes and filter eQTLs
#'
#' Identifies putative imprinted genes from Baran et al pipeline
#' and tests candidate genes for whether an eQTL better explains
#' expression patterns before outputting the final list of imprinted
#' genes.
#'
#' See Supplemental Materials for information on the 3 binomial tests
#' used to determine if an eQTL better explains the imprinting-like
#' expression. A significant p (<0.05) for >=2/3 tests indicates an eGene.
#'
#' Inputs    :
#' - `../../RNA_Imprinting/results/1KG.txt.impglr.res.txt`: 
#'             Baran Imprinting likelyhood results, created by (b)
#' - `../data/output/imprinting/all_samples.allelic_counts.cov8.tsv`: 
#'              Allelic count table, created by (a)
#' - `/.../eQTLs.txt`: 
#'              per-variant eQTL data, toy file has columns (SNP, beta, gene, `t-stat`,  `p-value`, FDR)
#' - `../resource/subject_ids.tsv`: 
#'              List of sample IDs with VCF data.
#'
#' Outputs:
#' - `../data/output/imprinting/1KG.baran_imprinted.tsv`: 
#'          Initial imprinting candidate genes without removing eGenes.
#' - `../data/output/imprinting/figures/{gene}/{gene}.{eQTL}.l_plot.pdf`: 
#'          eQTL effect diagnostic plots.
#' - `../data/output/imprinting/1KG.final_imprinted.tsv`: 
#'          Final list of putatively imprinted gene candidates with eGenes removed.
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

message("\nIDENTIFYING IMPRINTED GENES")
infile <- "../../RNA_Imprinting/results/1KG.txt.impglr.res.txt"

baran_results <- read.table(infile, sep = "", header = TRUE)

universe_genes <- baran_results %>% 
    filter(snpn >= 2, indn >= 5,   
           both == "+", hetlr < 1) 

message(paste("\tConsidering",nrow(universe_genes),"genes for imprinting..."))

imprinting_candidates <- universe_genes %>%
    filter(optz3p > 0.8, impglr > 40)

message(paste("\tStarting with",nrow(imprinting_candidates),"imprinted candidate genes..."))

outfile <- "../data/output/imprinting/1KG.baran_imprinted.tsv"
write.table(imprinting_candidates, outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

message("REMOVING eQTL GENES")

message("\tLoading variants with ASE data...")
infile <- "../data/output/imprinting/all_samples.allelic_counts.cov8.tsv"
all_samples_ac <- read.table(infile, sep = "\t", header = TRUE)

### eQTL data from NephQTL is not provided in resources
# The code block must be edited by the user to refer to their
# eQTL dataset

message("\tLoading eQTL variants...")
message("\t\tReading all tested variants...")
g <- fread(cmd= 'zcat ../resources/matrixeqtl/eQTLs_chrAll_40_peers_cis1000kb_Glom.txt.gz')
message("\t\tConverting hg19 to hg38...")
g_hg38 <- fread('../resources/SNP.map.Glom.TORUS.buildCrosswalk.tsv')
g <- g %>% dplyr::select(SNP, beta, gene, `t-stat`,  `p-value`, FDR)
g <- merge(g, g_hg38, by.x='SNP', by.y='hg19') %>%
    mutate(SNP_hg19 = paste0('chr', SNP)) %>%
    mutate(SNP = paste0('chr', hg38)) 
sig_eQTLs <- g %>%
    filter(`p-value` < 1e-6) %>%
    dplyr::select(SNP, gene, `p-value`)
    
###
message("\tIdentifying eQTL genotypes for sample...")
write.table(sig_eQTLs %>% dplyr::select(SNP) %>% distinct() %>% separate(SNP, into = c("CHROM", "POS", "REF", "ALT") , sep = "_") %>% 
            mutate(POS = as.integer(POS)) %>% filter(!is.na(POS)),
            "../data/output/imprinting/1KG.eQTL_GTs.for_query.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

eqtl_hets_list <- list()
i=1
subject_ids <- readLines("../resource/subject_ids.tsv")
for (sample in subject_ids) {    
    message(paste("\t\t", sample, "allelic data..."))
    system(paste0("bcftools query -R ../data/output/imprinting/1KG.eQTL_GTs.for_query.txt -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' ../data/", 
                sample,  "/", sample, ".het.phased.vcf.gz -o ../data/output/imprinting/", sample, ".eQTL_GTs.queried.txt"))

    file_path <- paste0("../data/output/imprinting/", sample, ".eQTL_GTs.queried.txt")
    if(file.exists(file_path)) {
        sample_eqtl_hets <- read.table(file_path, sep = "\t", header = FALSE) %>% mutate(SUBJECT_ID = sample)

        eqtl_hets_list[[i]] <- sample_eqtl_hets
        i <- i+1
    } else {
        message(paste("\t\tFile not found:", file_path))
    }
}
het_eQTLs <- do.call("rbind", eqtl_hets_list)

gene_eQTL_analysis <- function(ensg_gene) {
    gene_variants <- all_samples_ac %>%
        filter(ENSEMBL == ensg_gene)
    
    options(repr.plot.width = 3, repr.plot.height = 3)

    chr <- gene_variants %>%
        pull(CHROM) %>% 
        unique()

    # define region to search for eQTLs
    min <- gene_variants %>%
        pull(START) %>% # gene start
        min()-1000000

    max <- gene_variants %>%
        pull(STOP) %>% # gene stop
        max()+1000000

    # get the GT of each variant-individual in the region
    i=1
    gene_vars_list <- list()
    subject_ids <- readLines("../resource/subject_ids.tsv")
    for (sample in subject_ids) {    
        output <- system(paste0("bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\n]' -r ", chr,":",min,"-",max," ../data/", sample, "/", sample, ".phased.vcf.gz"), intern = TRUE)
        sample_gts <- read.table(text = output, header = FALSE, col.names = c("SUBJECT_ID", "CHROM", "POS", "REF", "ALT", "GT"))
        
        gene_vars_list[[i]] <- sample_gts
        i <- i+1
    }
    gene_vars <- do.call("rbind", gene_vars_list)

    sig_eqtl_list <- sig_eQTLs %>%
        filter(gene == ensg_gene) %>%
        arrange(`p-value`) %>%
        pull(SNP)
    message("\tAnalyzing ", length(sig_eqtl_list), " significant eQTLs for ", ensg_gene, " expression.")
    
    # Look at the most significant eQTL even if none significant.
    if (length(sig_eqtl_list) == 0) {
        sig_eqtl_list <- g %>%
            filter(gene == ensg_gene) %>%
            arrange(`p-value`) %>%
            head(5) %>%
            pull(SNP)
        message("\t\tNon-Significant eQTL:")
    } else {
        message("\t\tSignificant eQTL:")
    }
    
    for (SNP in sig_eqtl_list) {
        exit_code <- single_eqtl_analysis(ensg_gene, SNP)
        if (exit_code == 1) {
            message(paste("\tSUCCESSFUL eQTL FOR", ensg_gene, ":", SNP))
            return(1)
        }
    }
    return(0)
}

single_eqtl_analysis <- function(ensg_gene, eQTL) {
    options(repr.plot.width= 2, repr.plot.height=2)
    
    gene_variants <- all_samples_ac %>%
        filter(ENSEMBL == ensg_gene)
    
    parts <- unlist(strsplit(eQTL, "_"))
    eCHR <- parts[1]
    ePOS <- parts[2]
    eREF <- parts[3]
    eALT <- parts[4]

    chr <- gene_variants %>%
        pull(CHROM) %>% 
        unique()

    min <- gene_variants %>%
        pull(START) %>% # gene start
        min()-1000000

    max <- gene_variants %>%
        pull(STOP) %>% # gene stop
        max()+1000000

    i=1
    gene_vars_list <- list()
    subject_ids <- readLines("../resource/subject_ids.tsv")
    for (sample in subject_ids) {    
        output <- system(paste0("bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\n]' -r ", chr,":",min,"-",max," ../data/", sample, "/", sample, ".phased.vcf.gz"), intern = TRUE)
        sample_gts <- read.table(text = output, header = FALSE, col.names = c("SUBJECT_ID", "CHROM", "POS", "REF", "ALT", "GT"))
        
        gene_vars_list[[i]] <- sample_gts
        i <- i+1
    }
    gene_vars <- do.call("rbind", gene_vars_list)

    gene_vars_subj <- gene_vars %>%
        filter(CHROM == eCHR, POS == ePOS, REF == eREF, ALT == eALT) %>%
        distinct()

    plt_df <- gene_variants %>%
        filter(RNA_TOTAL_COUNT >= 8) %>%
        left_join(gene_vars_subj %>% dplyr::select(SUBJECT_ID, GT), by = c("SUBJECT_ID")) %>%
        left_join(gene_vars, by = c("CHROM", "POS", "REF", "ALT", "SUBJECT_ID")) %>%
        mutate(IN_PHASE = case_when(GT.x == "0|0" ~ "Hom. Ref.",
                                    GT.x == "1|1" ~ "Hom. Alt.",
                                   GT.x == GT.y ~ "Alt.",
                                   TRUE ~ "Ref."),
              eQTL_GT = case_when(GT.x %in% c("0|1", "1|0") ~ "Het.",
                                 GT.x %in% c("0|0", "1|1") ~ "Hom."))
    
    # set alternative hypothesis
    eQTL <- paste(eCHR, ePOS, eREF, eALT, sep = "_")
    
    eQTL_beta <- g %>%
        filter(gene == ensg_gene, SNP == eQTL) %>%
        pull(beta)
    eQTL_beta <- eQTL_beta[[1]]
    
    if (eQTL_beta > 0) { # alt increases
        ref_mono_alt <- "greater"
        alt_mono_alt <- "less"
        }
    else { # alt decreases
        ref_mono_alt <- "less"
        alt_mono_alt <- "greater"
    }
    
    eQTL_p <- g %>%
        filter(gene == ensg_gene, SNP == eQTL) %>%
        pull(`p-value`)
    
    # format eQTL p for printing
    
    if (eQTL_p < 2.2e-16) {
        eQTL_p <- " (p<2.2e-16)"
    } else if (eQTL_p < 0.001) {
        eQTL_p <- format(eQTL_p, scientific = TRUE, digits = 3)
        eQTL_p <- paste0(" (p=", eQTL_p, ")")
    } else {
        eQTL_p <- format(round(eQTL_p, 3), nsmall = 3, scientific = FALSE)
        eQTL_p <- paste0(" (p=", eQTL_p, ")")
    }

    message(paste("\t\t", eQTL, "beta:", eQTL_beta, "p:", eQTL_p))
  
    het_ref_binom_df <- plt_df %>%
        mutate(eQTL_GT_PHASE = paste(eQTL_GT, IN_PHASE)) %>% 
        filter(eQTL_GT_PHASE == "Het. Ref.") %>%
        mutate(REF_AXIS = (RNA_REF_COUNT / (RNA_REF_COUNT + RNA_ALT_COUNT) > 0.5),
              ALT_AXIS = (RNA_ALT_COUNT / (RNA_REF_COUNT + RNA_ALT_COUNT) > 0.5)) %>%
        summarise(N_REF_MONO = sum(REF_AXIS), TOTAL = n())
        
    het_alt_binom_df <- plt_df %>%
        mutate(eQTL_GT_PHASE = paste(eQTL_GT, IN_PHASE)) %>% 
        filter(eQTL_GT_PHASE == "Het. Alt.") %>%
        mutate(REF_AXIS = (RNA_REF_COUNT / (RNA_REF_COUNT + RNA_ALT_COUNT) > 0.5),
              ALT_AXIS = (RNA_ALT_COUNT / (RNA_REF_COUNT + RNA_ALT_COUNT) > 0.5)) %>%
        summarise(N_REF_MONO = sum(REF_AXIS), TOTAL = n())
        
    hom_eQTL_binom_df <- plt_df %>%
        filter(IN_PHASE %in% c("Hom. Ref.", "Hom. Alt.")) %>%
        mutate(REF_AXIS = (RNA_REF_COUNT / (RNA_REF_COUNT + RNA_ALT_COUNT) > 0.9),
              ALT_AXIS = (RNA_ALT_COUNT / (RNA_REF_COUNT + RNA_ALT_COUNT) > 0.9)) %>%
        summarise(N_MONO = sum(REF_AXIS | ALT_AXIS), TOTAL = n())
    
    p_het_ref <- NA
    p_het_alt <- NA
    p_hom <- NA
    if(het_ref_binom_df$TOTAL != 0) {
        het_ref_binom <- binom.test(het_ref_binom_df$N_REF_MONO, het_ref_binom_df$TOTAL, p = 0.5, alternative = ref_mono_alt)
        p_het_ref <- het_ref_binom$p.value
        if (p_het_ref < 2.2e-16) {
            p_het_ref <- "p<2.2e-16"
        } else if (p_het_ref < 0.001) {
            p_het_ref <- format(p_het_ref, scientific = TRUE, digits = 3)
            p_het_ref <- paste0("p=", p_het_ref)
        } else if (is.na(p_het_ref)) {
            p_het_ref <- paste0("p=NA")
            
        } else {
            p_het_ref <- format(round(p_het_ref, 3), nsmall = 3, scientific = FALSE)
            p_het_ref <- paste0("p=", p_het_ref)
        }
    }
    
    if(het_alt_binom_df$TOTAL != 0) {
        het_alt_binom <- binom.test(het_alt_binom_df$N_REF_MONO, het_alt_binom_df$TOTAL, p = 0.5, alternative = alt_mono_alt)
        p_het_alt <- het_alt_binom$p.value
        if (p_het_alt < 2.2e-16) {
            p_het_alt <- "p<2.2e-16"
        } else if (p_het_alt < 0.001) {
            p_het_alt <- format(p_het_alt, scientific = TRUE, digits = 3)
            p_het_alt <- paste0("p=", p_het_alt)
        } else if (is.na(p_het_alt)) {
            p_het_alt <- paste0("p=NA")
            
        } else {
            p_het_alt <- format(round(p_het_alt, 3), nsmall = 3, scientific = FALSE)
            p_het_alt <- paste0("p=", p_het_alt)
        }
    }
    
    if(hom_eQTL_binom_df$TOTAL != 0) {
        # Hom. eQTL binomial test expects mostly balanced signal SNPs, <=10% monoallelic, alternative (imprinted) is >10%
        hom_eQTL_binom <- binom.test(hom_eQTL_binom_df$N_MONO, hom_eQTL_binom_df$TOTAL, p = 0.9, alternative = "less")
        p_hom <- round(hom_eQTL_binom$p.value, 3)
        if (p_hom < 2.2e-16) {
            p_hom <- "p<2.2e-16"
        } else if (p_hom < 0.001) {
            p_hom <- format(p_hom, scientific = TRUE, digits = 3)
            p_hom <- paste0("p=", p_hom)
        } else if (is.na(p_hom)) {
            p_hom <- "p=NA"
        } else {
            p_hom <- format(round(p_hom, 3), nsmall = 3, scientific = FALSE)
            p_hom <- paste0("p=", p_hom)
        }
    }

    ax_lim <- max(plt_df$RNA_REF_COUNT, plt_df$RNA_ALT_COUNT, na.rm = TRUE)

    custom_palette <- c("Alt." = "#984EA3",
                "Ref." = "#E41A1C",    
                "Hom. Ref." = "#377EB8",
                "Hom. Alt." = "#4DAF4A")
    
    PLCG2_glom_counts_plot <- ggplot(data = plt_df) +
        geom_point(aes(x = RNA_REF_COUNT, y = RNA_ALT_COUNT, 
                       shape = eQTL_GT, color = IN_PHASE)) + 
        geom_abline(slope = 1/9, linetype = "dashed") +
        geom_abline(slope = 9, linetype = "dashed") +
        annotate("rect", xmin = 0, xmax = ax_lim, 
                 ymin = ax_lim - ax_lim / 25*2, ymax = ax_lim,
                 fill = "white", alpha = 1) +
        annotate("text", x = ax_lim / 2, y = ax_lim, label = paste0(eCHR,":", ePOS, "_", eREF, "/", eALT, eQTL_p), size = 7/.pt, vjust = 1) + 
        theme_classic() +
        theme(axis.text.x = element_text(size = 7),
              axis.text.y = element_text(size = 7),
              axis.title = element_blank(),
              legend.position = "none") + 
            guides(color = guide_legend(order = 1, ncol = 2, override.aes = list(shape = 15)),
                   shape = guide_legend(order = 2, ncol = 2)) +
        labs(x = "RNA Ref Count", y = "RNA Alt Count", 
             color = "Phased\nwith eQTL Alt:", shape = "GT:") +
        scale_color_manual(values = custom_palette) +
        scale_shape_manual(values=c("Het." = 16,"Hom." = 17)) + 
        xlim(0, ax_lim) +
        ylim(0, ax_lim) +
        annotate("text", x = ax_lim, y = ax_lim - ax_lim / 25*2, label = p_het_ref,
               hjust = 1, vjust = 1, color = "#E41A1C", size = 7/.pt) +
        annotate("text", x = ax_lim, y = ax_lim - ax_lim / 25 * 4, label = p_het_alt,
               hjust = 1, vjust = 1, color = "#984EA3", size = 7/.pt) +
        annotate("text", x = ax_lim, y = ax_lim - ax_lim / 25 * 6, label = p_hom,
               hjust = 1, vjust = 1, color = "black", size = 7/.pt) + 
        facet_wrap(~gene_symbol)

    outfile <- paste0("../data/output/imprinting/figures/", ensg_gene, "/", ensg_gene, ".", eQTL ,".l_plot.pdf")
    dir.create(dirname(outfile), recursive = TRUE)
    pdf(outfile, width=2, height=2)
    print(PLCG2_glom_counts_plot)
    dev.off()
    
    true_count <- 0
    if (exists("het_ref_binom") && !is.null(het_ref_binom$p.value) && het_ref_binom$p.value < 0.05) {
      true_count <- true_count + 1
    }
    if (exists("het_alt_binom") && !is.null(het_alt_binom$p.value) && het_alt_binom$p.value < 0.05) {
      true_count <- true_count + 1
    }
    if (exists("hom_eQTL_binom") && !is.null(hom_eQTL_binom$p.value) && hom_eQTL_binom$p.value < 0.05) {
      true_count <- true_count + 1
    }

    if (true_count >= 2) {
        return(1)
    } else {
        return(0)
    }
}

infile <- "../data/output/imprinting/1KG.baran_imprinted.tsv"
imprinting_candidates <- read.table(infile, sep = "\t", header = TRUE)

true_imprinted <- c()
for (gene in imprinting_candidates$ensg) {
    exit_code <- gene_eQTL_analysis(gene)
    if (exit_code == 1) {
        true_imprinted <- true_imprinted
    } else if (exit_code == 0) {
        true_imprinted <- c(true_imprinted, gene)
    }
}

imprinted <- imprinting_candidates %>% filter(ensg %in% true_imprinted)
outfile <- "../data/output/imprinting/1KG.final_imprinted.tsv"
write.table(imprinted, outfile, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)