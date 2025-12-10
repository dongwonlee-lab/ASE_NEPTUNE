library(tidyverse)

#' Identify DNA-RNA discordant variants
#'
#' Identifies variants where the RNA reads don't
#' match the DNA, generally reads with the alt. allele
#' mapping to a site that is homozygous reference, or 
#' vice versa. Briefly, a discordant variant is one where:
#' - Avg. RNA Alt. Fraction among all Homozygous ref.
#'   individuals is >0.1
#'   OR
#' - Avg. RNA Alt. Fraction among all Homozygous alt.
#'   individuals is <0.9
#'
#' Output: 
#' - all_samples.dna_rna_discordant_vars.bed:
#'   BED file with all discordant variants to be blacklisted
#'   in phASER run.
#'
#' Usage:
#'   Rscript b_find_dna_rna_discord_vars.R
#'
#' Author:
#'   Eric Sakkas
#'   erotokritos.sakkas@childrens.harvard.edu
#'
#' Date:
#'   2025-09-10

subject_ids <- readLines("../resource/subject_ids.tsv")
print(subject_ids)

for (subject_id in subject_ids) {
    print(paste("Reading files for", subject_id))

    file_path <- paste0("../data/output/phaser/", subject_id,"/",subject_id, ".discord_snps_gts.tsv")
    gts <- read.table(file_path, header = FALSE, sep = "\t")

    file_path <- paste0("../data/output/phaser/", subject_id,"/",subject_id, ".discord_snps_pileup.tsv")
    pileup <- read.table(file_path, header = FALSE, sep = "\t")
    pileup$SAFID <- subject_id

    print("Combined Pileup...")

    gts_pileup <- gts %>%
        inner_join(pileup, by = c("V1", "V2", "V5" = "SAFID"))

    rm(gts, pileup)

    print("Combined GT and Pileup...")

    gts_pileup <- gts_pileup %>%
        rename("CHROM"       = "V1",
            "POS"         = "V2",
            "REF"         = "V3.x",
            "ALT"         = "V4.x",
            "TOTAL_COUNT" = "V4.y",
            "SAFID"       = "V5",
            "GT"          = "V6.x",
            "PILEUP"      = "V5.y",
            "BASEQ"       = "V6.y"
        ) %>%
        dplyr::select(-V3.y)

    print("Filtered RNA cov8 and indels...")
    gts_pileup <- gts_pileup %>%
        filter(TOTAL_COUNT >= 8) %>%
        filter(nchar(REF) == 1 | nchar(ALT) == 1)

    count_ref_alt <- function(row) {
        pileup <- tolower(row["PILEUP"])
        ref <- tolower(row["REF"])
        alt <- tolower(row["ALT"])
        ref_count <- sum(strsplit(pileup, "")[[1]] == ref)
        alt_count <- sum(strsplit(pileup, "")[[1]] == alt)
        return(c(ref_count, alt_count))
    }

    results <- lapply(1:nrow(gts_pileup), function(i) count_ref_alt(gts_pileup[i, ]))
    results <- do.call(rbind, results)
    gts_pileup$REF_COUNT <- results[, 1]
    gts_pileup$ALT_COUNT <- results[, 2]
    gts_pileup$PILEUP_TOTAL_COUNT <- gts_pileup$REF_COUNT + gts_pileup$ALT_COUNT

    print("Counted Pileup...")

    gts_pileup <- gts_pileup %>%
        filter(PILEUP_TOTAL_COUNT >=8) %>%
        dplyr::select(-PILEUP, -BASEQ)

    gts_pileup <- gts_pileup %>%
    mutate(
        GT_LAB = case_when(
            GT == "1|0" ~ "HET",
            GT == "0|1" ~ "HET",
            GT == "0/1" ~ "HET",
            GT == "0|0" ~ "HOM_REF",
            GT == "0/0" ~ "HOM_REF",
            GT == "1|1" ~ "HOM_ALT",
            GT == "1/1" ~ "HOM_ALT",
            GT == "./." ~ "MISSING"
        ))
                    
    output_file <- paste0("../data/output/phaser/",subject_id,"/",subject_id,".pileup_counted.tsv")
                    
    write.table(gts_pileup, output_file, row.names = FALSE, sep = "\t", quote = FALSE)
}

all_subj_gts_pileup <- data.frame()

for (subject_id in subject_ids) {
    file_path <- paste0("../data/output/phaser/",subject_id,"/",subject_id,".pileup_counted.tsv")
    pileup_analysis <- read.table(file_path, header = TRUE, sep = "\t")

    pileup_analysis <- pileup_analysis %>% dplyr::select(-GT, -TOTAL_COUNT)

    all_subj_gts_pileup <- rbind(all_subj_gts_pileup, pileup_analysis)
}

print("Summarizing missing GT data...")

snps_miss_frac <- all_subj_gts_pileup %>%
    group_by(CHROM, POS, REF, ALT) %>%
    summarize(N_INDIV_COV8 = n(),
                N_MISSING = sum(GT_LAB=="MISSING"),
                FRAC_MISSING = N_MISSING / N_INDIV_COV8) %>%
    dplyr::select(CHROM, POS, REF, ALT, FRAC_MISSING)

print("Summarizing pileup data across individuals...")

all_subj_gts_pileup <- all_subj_gts_pileup %>%
    mutate(PILEUP_RNA_ALT_FRACTION = ALT_COUNT / PILEUP_TOTAL_COUNT) %>%
    group_by(CHROM, POS, REF, ALT, GT_LAB) %>%
    summarise(N_INDIV_COV8 = n(),
                AVG_RNA_ALT_FRAC = mean(PILEUP_RNA_ALT_FRACTION),
                RNA_ALT_FRAC_LIST = paste(round(PILEUP_RNA_ALT_FRACTION,4), collapse = ","))

all_subj_gts_pileup <- all_subj_gts_pileup %>%
                                left_join(snps_miss_frac, by = c("CHROM", "POS", "REF", "ALT"))

rm(snps_miss_frac)

all_subj_gts_pileup_wider <- all_subj_gts_pileup %>%
    pivot_wider(id_cols = c(CHROM,POS,REF,ALT,FRAC_MISSING),
                names_from = GT_LAB,
                values_from = c(AVG_RNA_ALT_FRAC, N_INDIV_COV8, RNA_ALT_FRAC_LIST))

output_file <- paste0("../data/output/phaser/all_samples.pileup_analysis.tsv")

write.table(all_subj_gts_pileup_wider, output_file, row.names = FALSE, sep = "\t", quote = FALSE)
all_subj_gts_pileup_wider <- read.table(output_file, sep = "\t", header = TRUE)

print("Identifying DNA-RNA discordant variants...")

discord_snps <- all_subj_gts_pileup_wider %>%
    filter(nchar(REF)==1 & nchar(ALT)==1) %>% # remove indels
    filter(AVG_RNA_ALT_FRAC_HOM_REF > 0.1 |
          AVG_RNA_ALT_FRAC_HOM_ALT < 0.9) 

n_discord_vars <- discord_snps %>% nrow()

print(paste("Found",n_discord_vars,"DNA-RNA discordant variants"))

discord_snps_bed <- discord_snps %>%
    dplyr::select(CHROM, POS) %>%
    mutate(POS_BED = POS-1) %>%
    dplyr::select(CHROM, POS_BED, POS) %>%
    mutate(REASON = "DISCORD")

output_file <- paste0("../data/output/phaser/all_samples.dna_rna_discordant_vars.bed")

write.table(discord_snps_bed, output_file, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

print(paste("DNA-RNA discordant variants written to", output_file))
