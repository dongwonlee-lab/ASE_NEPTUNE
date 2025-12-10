#' Generate Sample- and Gene-level TOGA data
#'
#' This module processes TOGA results to generate
#' COV, ASE, and MAE gene sets. It produces both per-sample
#' and per-gene results, along with cohort-level summary statistics.
#' Designed as a post hoc utility, it provides a simple way to
#' analyze TOGA outputs.
#'
#' Outputs:
#' - Per-sample:
#'   - {sample}_cov_genes.txt: List of COV genes
#'   - {sample}_ase_genes.txt: List of ASE genes
#'   - {sample}_mae_genes.txt: List of MAE genes
#'   - {sample}.summary_stats.txt: Summary statistics (COV, ASE, MAE counts and proportions)
#' - Cohort-level:
#'   - all_samples.summary_stats.txt: Summary of COV, ASE, MAE across all samples
#' - Per-gene:
#'   - all_samples.{cov/ase/mae}_genes.txt: Combined per-gene results
#'   - all_genes.samples_with_{cov/ase/mae}.txt: Gene-level sample counts
#'   - all_genes.samples_with_{ase/mae}.prop_sort.txt: Same, sorted by proportion of samples with ASE/MAE
#'
#' Usage:
#'   Rscript ASE_TOGA.R \
#'     -i ./toy_example/TOY_INPUT \
#'     -o ./toy_example/TOY_OUTPUT 
#' 
#' Dependencies:
#'   R ≥ 3.5.x with packages: optparse, dplyr, stringr, tidyr
#' 
#' Author:
#'   Junmo Sung
#'   jsung5@bwh.harvard.edu
#' 
#' Date:
#'   2025-10-01

##### Pre-process #####
# Load libraries
library(optparse)
library(dplyr)
library(stringr)
library(tidyr)

arrange_contig_order <- function(df) {
    #' Rearrange TOGA file by chromosome
    #'
    #' This function reorders a TOGA dataframe by genomic contigs
    #' with autosomes placed in numeric order, followed by
    #' chrX, chrY, and mitochondrial contigs (chrM or chrMT). Any contigs
    #' not matching these patterns are pushed to the end.
    #'
    #' @param df A dataframe containing at least the columns:
    #'   \code{contig}, \code{start}, \code{name}, \code{gene_symbol},
    #'   and \code{sample}.
    #'
    #' @return A dataframe with rows rearranged by chromosome and genomic order.
    #'
    #' @examples
    #' \dontrun{
    #' arranged_df <- arrange_contig_order(toga_df)
    #' }
    #'
    #' @export
    df %>%
        mutate(
        contig_order = suppressWarnings(case_when(
            grepl("^chr[0-9]+$", contig) ~ as.numeric(gsub("chr", "", contig)),
            contig == "chrX" ~ 23,
            contig == "chrY" ~ 24,
            contig %in% c("chrM", "chrMT") ~ 25,
            TRUE ~ 99
        ))
        ) %>%
        arrange(contig_order, start, name, gene_symbol, sample) %>%
        select(-contig_order)
}

# Parse options
option_list <- list(
    make_option(c("-i", "--input"), type="character", default="../data/output/toga/",
                help="Input directory", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output directory. WARNING: if it exists, script will delete its contents.", metavar="character")
    )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input_dir <- opt$input
out_dir <- opt$output

if (dir.exists(out_dir)) {
    message(paste(out_dir, "exists. Deleting contents..."))
    cat("\tDATE, TIME :", format(Sys.time(), "%Y-%m-%d, %H:%M:%S"), "\n")
    unlink(out_dir, recursive = TRUE)
}

dir.create(out_dir)

######################################################
subject_ids <- readLines("../resource/subject_ids.tsv")
######################################################

message("PER INDIVIDUAL RESULT")
cat("\tDATE, TIME :", format(Sys.time(), "%Y-%m-%d, %H:%M:%S"), "\n")
for (sample in subject_ids) {    
    message(paste("Loading", sample, "TOGA data..."))

    file_path <- file.path(input_dir, sample, paste0(sample, "_TOGA_light.txt"))
    toga <- read.delim(file_path)
    
    num_toga<-nrow(toga)
    message(paste0("\t", num_toga, " genes found in TOGA file."))

    ### COV_N genes ###
    cov_genes <-toga %>% filter(ASE_COV %in% c('ASE', 'COV'))
    
    num_cov<-nrow(cov_genes)
    message(paste0("\t", num_cov, " genes with COV"))
    
    # Summed total counts of COV_N genes
    num_totalCount_cov<-sum(cov_genes$total_hap_count)
    
    # Summed total variants of COV_N genes
    num_totalvar_cov<-sum(cov_genes$n_variants)
    
    cov_reorder <- arrange_contig_order(cov_genes)

    out_dir <- sub("/+$", "", out_dir)
    sample_outdir <- file.path(out_dir, sample)
    dir.create(sample_outdir, recursive = TRUE, showWarnings = FALSE)

    outpath <- paste0(sample_outdir,"/",sample,"_cov_genes.txt")
    write.table(cov_reorder,
        file = outpath,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE)
    message(paste("\tCOV genes written to", outpath))

    ### ASE genes ###
    ase<-toga %>% filter(ASE_COV == 'ASE')
    
    # Summed total counts of ASE genes
    num_totalCount_ase<-sum(ase$total_hap_count)
    
    # Summed total variants of ASE genes
    num_totalvar_ase<-sum(ase$n_variants)

    num_ase <- nrow(ase)
    message(paste0("\t", num_ase, " ASE genes found."))

    outpath <- paste0(sample_outdir,"/",sample,"_ase_genes.txt")
    write.table(ase,
        file = outpath,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE)
    message(paste("\tASE genes written to",outpath))
    
    ### MAE genes ###
    mae<-ase %>% filter((hapA_count == 0 | hapB_count == 0))

    num_mae<-nrow(mae)
    message(paste0("\t", num_mae, " monoallelic ASE genes found."))
    
    # Summed total counts of MAE genes
    num_totalCount_mae<-sum(mae$total_hap_count)
    
    # Summed total variants of MAE genes
    num_totalvar_mae<-sum(mae$n_variants)
    
    outpath <- paste0(sample_outdir,"/",sample,"_mae_genes.txt")
    write.table(mae,
        file = outpath,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE)
    message(paste("\tMonoallelic ASE genes written to",outpath))

    ### Create statistics COV_N/ASE/MAE genes per individual ###
    stat <- data.frame(
        n_genes = num_toga,
        n_cov_genes = num_cov,
        n_ase_genes = num_ase,
        n_mae_genes = num_mae,
        total_reads_cov = num_totalCount_cov,
        total_vars_cov = num_totalvar_cov,
        total_reads_ase = num_totalCount_ase,
        total_vars_ase = num_totalvar_ase,
        total_reads_mae = num_totalCount_mae,
        total_vars_mae = num_totalvar_mae,
        prop_cov = num_cov/num_toga,
        prop_ase = num_ase/num_toga,
        prop_ase_cov = num_ase/num_cov,
        prop_mae = num_mae/num_toga,
        prop_mae_cov = num_mae/num_cov,
        prop_mae_ase = num_mae/num_ase)
    
    write.table(stat,
        file = paste0(sample_outdir,"/",sample,".summary_stats.txt"),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE)
}

message("Generating summary statistics for all samples...")
cat("\tDATE, TIME :", format(Sys.time(), "%Y-%m-%d, %H:%M:%S"), "\n")
out_subdir_list <- list.dirs(path = out_dir, full.names = TRUE, recursive = FALSE)

all_stat <- data.frame()
for (out_subdir in out_subdir_list) {
    sample <- basename(out_subdir)
    file_path <- file.path(out_subdir, paste0(sample, ".summary_stats.txt"))
    stat <- read.table(file_path, header = TRUE, sep = "\t")
    stat$sample <- sample
    all_stat <- rbind(all_stat, stat)
}

# Reordering by the sample's name
all_stat <- all_stat[, c("sample", names(all_stat)[!names(all_stat) %in% "sample"])]

outfile <- paste0(out_dir,"/","all_samples.summary_stats.txt")
write.table(all_stat,
    file = outfile,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE)
message(paste("Summary statistics written to", outfile))

message("PER GENE RESULT")
message("Merging COV_N/ASE/MAE across all samples...")
cat("\tDATE, TIME :", format(Sys.time(), "%Y-%m-%d, %H:%M:%S"), "\n")
cov_list <- list()
ase_list <- list()
mae_list <- list()

for (out_subdir in out_subdir_list) {
    out_sample <- basename(out_subdir)

    cov_out_file_path <- file.path(out_subdir, paste0(out_sample, "_cov_genes.txt"))
    ase_out_file_path <- file.path(out_subdir, paste0(out_sample, "_ase_genes.txt"))
    mae_out_file_path <- file.path(out_subdir, paste0(out_sample, "_mae_genes.txt"))
    cov_gene <- read.delim(cov_out_file_path)
    ase_gene <- read.delim(ase_out_file_path)
    mae_gene <- read.delim(mae_out_file_path)
    
    if (nrow(cov_gene) > 0) {
        cov_gene$sample <- out_sample
    }
    if (nrow(ase_gene) > 0) {
        ase_gene$sample <- out_sample
    }
    if (nrow(mae_gene) > 0) {
        mae_gene$sample <- out_sample
    }
    
    # Convert to dataframe format
    cov_list[[out_sample]] <- cov_gene
    ase_list[[out_sample]] <- ase_gene
    mae_list[[out_sample]] <- mae_gene
}

# Merge all dataframes into one dataframe
all_cov <- do.call(rbind, cov_list)
all_ase <- do.call(rbind, ase_list)
all_mae <- do.call(rbind, mae_list)

# Concatenated COV_N, ASE, MAE files across all samples
if (nrow(all_cov) == 0) {
  all_cov_sort <- all_cov
  all_cov_reorder <- all_cov
} else {
  all_cov_sort <- arrange_contig_order(all_cov)
  all_cov_reorder <- all_cov_sort[c("sample", setdiff(names(all_cov_sort), "sample"))]
}

if (nrow(all_ase) == 0) {
  all_ase_sort <- all_ase
  all_ase_reorder <- all_ase
} else {
  all_ase_sort <- arrange_contig_order(all_ase)
  all_ase_reorder <- all_ase_sort[c("sample", setdiff(names(all_ase_sort), "sample"))]
}

if (nrow(all_mae) == 0) {
  all_mae_sort <- all_mae
  all_mae_reorder <- all_mae
} else {
  all_mae_sort <- arrange_contig_order(all_mae)
  all_mae_reorder <- all_mae_sort[c("sample", setdiff(names(all_mae_sort), "sample"))]
}

cov_outfile <- paste0(out_dir,"/","all_samples.cov_genes.txt")
write.table(all_cov_reorder,
    file = cov_outfile,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE)
message(paste("All RNA cov genes written to", cov_outfile))

ase_outfile <- paste0(out_dir,"/","all_samples.ase_genes.txt")
write.table(all_ase_reorder,
    file = ase_outfile,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE)
message(paste("All ASE genes written to", ase_outfile))

mae_outfile <- paste0(out_dir,"/","all_samples.mae_genes.txt")
write.table(all_mae_reorder,
    file = mae_outfile,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE)
message(paste("All monoallelic expression genes written to", mae_outfile))

message("Getting number of samples with COV_N/ASE/MAE gene")

if (nrow(all_cov) == 0) {
  res_cov <- all_cov
} else {
    res_cov <- all_cov %>% group_by(name, gene_symbol) %>%
    summarize(n_observed = n(),
         samples = paste(unique(sample), collapse = ", "),
         .groups = 'drop')
}
if (nrow(all_ase) == 0) {
  res_ase <- all_ase
} else {
    res_ase <- all_ase %>% group_by(name, gene_symbol) %>%
    summarize(n_observed = n(),
         samples = paste(unique(sample), collapse = ", "),
         .groups = 'drop')
}
if (nrow(all_mae) == 0) {
  res_mae <- all_mae
} else {
    res_mae <- all_mae %>% group_by(name, gene_symbol) %>%
    summarize(n_observed = n(),
         samples = paste(unique(sample), collapse = ", "),
         .groups = 'drop')
}

# Sort by number of samples in a descending order
if (nrow(all_cov) == 0) {
  res_cov <- all_cov
} else {
    res_cov<-res_cov[order(-res_cov$n_observed),]
}
if (nrow(all_ase) == 0) {
  res_ase <- all_ase
} else {
    res_ase<-res_ase[order(-res_ase$n_observed),]
}
if (nrow(all_mae) == 0) {
  res_mae <- all_mae
} else {
    res_mae<-res_mae[order(-res_mae$n_observed),]
}

# collapse samples to list of samples
res_cov$samples <- str_split(res_cov$samples, ",", simplify = FALSE)
res_ase$samples <- str_split(res_ase$samples, ",", simplify = FALSE)
res_mae$samples <- str_split(res_mae$samples, ",", simplify = FALSE)

# Count number of samples with RNA coverage >= threshold
if (nrow(all_ase) == 0) {
    merged_ase <- res_ase
} else {
    merged_ase <- left_join(res_ase, res_cov, by = "name", suffix = c("", "_cov"))
}
if (nrow(all_mae) == 0) {
    merged_mae <- res_mae
} else {
    merged_mae <- left_join(res_mae, res_cov, by = "name", suffix = c("", "_cov"))
}
    
# Samples with COV_N but without ASE
merged_ase$samples_cov <- mapply(setdiff, merged_ase$samples_cov, merged_ase$samples)
merged_mae$samples_cov <- mapply(setdiff, merged_mae$samples_cov, merged_mae$samples)
# Convert 'samples_cov' column from list to comma-separated string
merged_ase$samples_cov <- sapply(merged_ase$samples_cov, function(x) paste(x, collapse = ","))
merged_mae$samples_cov <- sapply(merged_mae$samples_cov, function(x) paste(x, collapse = ","))                                                                     

if (nrow(all_ase) == 0) {                                      
    merged_ase <- merged_ase
} else {
    merged_ase <- merged_ase %>% mutate(n_observed_cov = n_observed_cov)
}
if (nrow(all_mae) == 0) {                                      
    merged_mae <- merged_mae
} else {
    merged_mae <- merged_mae %>% mutate(n_observed_cov = n_observed_cov)
}
                                      
# Proportion of samples with ASE (Number of samples with ASE / Number of samples with COV_N)         
if (nrow(all_ase) == 0) {                                      
    merged_ase <- merged_ase
} else {
    merged_ase$pct_observed <- with(merged_ase, n_observed / n_observed_cov) 
}                             
if (nrow(all_mae) == 0) {                                      
    merged_mae <- merged_mae
} else {
    merged_mae$pct_observed <- with(merged_mae, n_observed / n_observed_cov) 
}
                                      
## Subset the dataframe
if (nrow(all_ase) == 0) {                                      
    final_res_ase <- merged_ase
} else {
    final_res_ase <- merged_ase %>% select(name, gene_symbol, n_observed, n_observed_cov, pct_observed, samples,samples_cov)
}    
if (nrow(all_mae) == 0) {                                      
    final_res_mae <- merged_mae
} else {
    final_res_mae <- merged_mae %>% select(name, gene_symbol, n_observed, n_observed_cov, pct_observed, samples,samples_cov)
}
                                     
# Convert 'samples' column back from list to comma-separated string
res_cov$samples <- sapply(res_cov$samples, function(x) paste(x, collapse = ","))   
final_res_ase$samples <- sapply(final_res_ase$samples, function(x) paste(x, collapse = ","))
final_res_mae$samples <- sapply(final_res_mae$samples, function(x) paste(x, collapse = ","))
                                
# Sort by pct_observed in a descending order                                 
if (nrow(all_ase) == 0) {                                      
    final_res_ase_pct_sort <- final_res_ase
} else {
    final_res_ase_pct_sort <- final_res_ase[order(-final_res_ase$pct_observed),]                              
}   
if (nrow(all_mae) == 0) {                                      
    final_res_mae_pct_sort <- final_res_mae
} else {
    final_res_mae_pct_sort <- final_res_mae[order(-final_res_mae$pct_observed),]                              
}

fix_sample_duplicates <- function(data_frames) {
    #' Fix duplicated samples in ASE/MAE result data frames
    #'
    #' This function processes a named list of result data frames and ensures
    #' that the \code{samples_cov} column does not contain any sample IDs
    #' already present in the \code{samples} column. It removes duplicates,
    #' cleans up intermediate columns, and returns the corrected list of
    #' data frames.
    #'
    #' @param data_frames A named list of data frames. Each data frame must
    #'   contain at least the columns \code{samples} and \code{samples_cov}.
    #'
    #' @return A named list of data frames with corrected \code{samples_cov}
    #'   values. Within each data frame:
    #'   \itemize{
    #'     \item \code{samples_cov} will no longer contain sample IDs that
    #'       are present in \code{samples}.
    #'     \item The output retains the same structure as the input, with
    #'       no intermediate list columns.
    #'   }
    #'
    #' @examples
    #' \dontrun{
    #' data_frames <- list(
    #'   final_res_ase = final_res_ase,
    #'   final_res_ase_pct_sort = final_res_ase_pct_sort,
    #'   final_res_mae = final_res_mae,
    #'   final_res_mae_pct_sort = final_res_mae_pct_sort
    #' )
    #'
    #' cleaned <- fix_sample_duplicates(data_frames)
    #' final_res_ase <- cleaned$final_res_ase
    #' }
    #'
    #' @export
    for (name in names(data_frames)) {
        df <- data_frames[[name]]
        if (nrow(df) == 0) next

        df$subdir_list <- strsplit(df$samples, ",\\s*")
        df$samples_cov_list <- strsplit(df$samples_cov, ",\\s*")

        # samples_cov - samples
        df$samples_cov_list <- mapply(function(non_ase, samples) {
        setdiff(non_ase, samples)
        }, df$samples_cov_list, df$subdir_list, SIMPLIFY = FALSE)

        df$samples_cov <- sapply(df$samples_cov_list, paste, collapse = ", ")

        # remove intermediates
        df$subdir_list <- NULL
        df$samples_cov_list <- NULL

        data_frames[[name]] <- df
    }
    return(data_frames)
}

cleaned <- fix_sample_duplicates(list(
    final_res_ase = final_res_ase,
    final_res_ase_pct_sort = final_res_ase_pct_sort,
    final_res_mae = final_res_mae,
    final_res_mae_pct_sort = final_res_mae_pct_sort
))

final_res_ase <- cleaned$final_res_ase
final_res_ase_pct_sort <- cleaned$final_res_ase_pct_sort
final_res_mae <- cleaned$final_res_mae
final_res_mae_pct_sort <- cleaned$final_res_mae_pct_sort

res_cov <- res_cov %>%
    rename("ensembl" = "name", 
           "n_samples_cov" = "n_observed")

final_res_ase <- final_res_ase %>%
    rename("ensembl" = "name", 
           "n_samples_ase" = "n_observed",
           "n_samples_cov" = "n_observed_cov",
           "prop_ase_cov" = "pct_observed")

final_res_ase_pct_sort <- final_res_ase_pct_sort %>%
    rename("ensembl" = "name", 
           "n_samples_ase" = "n_observed",
           "n_samples_cov" = "n_observed_cov",
           "prop_ase_cov" = "pct_observed")

final_res_mae <- final_res_mae %>%
    rename("ensembl" = "name", 
           "n_samples_mae" = "n_observed",
           "n_samples_cov" = "n_observed_cov",
           "prop_mae_cov" = "pct_observed")

final_res_mae_pct_sort <- final_res_mae_pct_sort %>%
    rename("ensembl" = "name", 
           "n_samples_mae" = "n_observed",
           "n_samples_cov" = "n_observed_cov",
           "prop_mae_cov" = "pct_observed")

write.table(res_cov,
    file = paste0(out_dir,"/","all_genes.samples_with_cov.txt"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE)
write.table(final_res_ase,
    file = paste0(out_dir,"/","all_genes.samples_with_ase.txt"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE)
write.table(final_res_ase_pct_sort,
    file = paste0(out_dir,"/","all_genes.samples_with_ase.prop_sort.txt"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE)                                 
write.table(final_res_mae,
    file = paste0(out_dir,"/","all_genes.samples_with_mae.txt"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE)
write.table(final_res_mae_pct_sort,
    file = paste0(out_dir,"/","all_genes.samples_with_mae.prop_sort.txt"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE)                                 
