suppressMessages(library(tidyverse))

#' Appends allelic counts to haplotype counts 
#' 
#' Takes the phASER haplotypic counts output and appends
#' allele counts within the haplotype.
#' The output files are optional inpputs to TOGA.
#'
#' Outputs: 
#' - var_hap_cnt.txt:
#'      For each haplotype, contains:
#' varHapA: Alleles in haplotype A
#' varHapB: Alleles in haplotype B
#' varCntA: Allelic counts in haplotype A
#' varCntB: Allelic counts in haplotype B
#'
#' Usage, do not put paths to `_counts`` files:
#' Rscript extract_var_hap_cnt.R \
#'   -i /path/to/directory_with_phaser_outputs \
#'   -h haplotypic_counts.txt \
#'   -a allelic_counts.txt \
#'   -g gene_ae.txt \
#'   -o /path/to/var_hap_cnt.txt
#'
#' Author:
#'   Junmo Sung
#'   jsung5@bwh.harvard.edu
#'
#' Date:
#'   2025-09-25

# ---- small helper: parse -i/-h/-a/-g/-o flags without extra packages ----
parse_flags <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # if user ran with positional args (legacy), keep supporting it
  if (length(args) == 5 && !any(startsWith(args, "-"))) {
    return(list(
      sampdir = args[1],
      hcfn    = args[2],
      acfn    = args[3],
      gaefn   = args[4],
      outfn   = args[5]
    ))
  }

  # otherwise, expect flagged inputs
  want <- c("-i", "-h", "-a", "-g", "-o")
  if (!all(want %in% args)) {
    stop(paste0(
      "Usage:\n",
      "  Rscript extract_var_hap_cnt.R \\\n",
      "    -i <phaser_dir> \\\n",
      "    -h <haplotypic_counts.txt> \\\n",
      "    -a <allelic_counts.txt> \\\n",
      "    -g <gene_ae.txt> \\\n",
      "    -o <out.txt>\n",
      "\n(Alternatively, legacy positional form is also supported.)\n"
    ))
  }

  # map flags to values
  get_val <- function(flag) {
    idx <- which(args == flag)
    if (length(idx) != 1 || idx == length(args)) {
      stop(sprintf("Flag %s is missing a value.", flag))
    }
    val <- args[idx + 1]
    if (startsWith(val, "-")) stop(sprintf("Flag %s is missing a value.", flag))
    val
  }

  list(
    sampdir = get_val("-i"),
    hcfn    = get_val("-h"),
    acfn    = get_val("-a"),
    gaefn   = get_val("-g"),
    outfn   = get_val("-o")
  )
}

opt <- parse_flags()

sampdir <- opt$sampdir                    # phASER directory
hcfn    <- file.path(sampdir, opt$hcfn)   # haplotypic_counts.txt
acfn    <- file.path(sampdir, opt$acfn)   # allelic_counts.txt
gaefn   <- file.path(sampdir, opt$gaefn)  # gene_ae.txt
outfn   <- opt$outfn                      # output file

# ensure output dir exists
outdir <- dirname(outfn)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- read inputs ----
df_hc  <- readr::read_tsv(hcfn,  show_col_types = FALSE)
df_ac  <- readr::read_tsv(acfn,  show_col_types = FALSE)
df_gae <- readr::read_tsv(gaefn, show_col_types = FALSE)

df_gae1 <- df_gae %>% filter(!is.na(variants))

max_nvar <- max(unlist(lapply(regmatches(df_gae1$variants, gregexpr(',', df_gae1$variants)), length)))

df_gae1_w <- df_gae1 %>% 
        separate_wider_delim(variants, ',',
                             names=paste("variant", 1:(max_nvar+1), sep="_"),
                             too_few = "align_start")

df_gae1_wl <- df_gae1_w %>% pivot_longer(10:(ncol(df_gae1_w)-2),
                      names_to = "variant_no",
                      values_to = "variantID") %>%
        filter(!is.na(variantID)) %>%
        select(-bam)

df_gae1_wlc <- df_gae1_wl %>% inner_join(df_ac, by=c("contig", "variantID"))

max_nvar_hc <- max(unlist(lapply(regmatches(df_hc$variants, gregexpr(',', df_hc$variants)), length)))

df_hc_w <- df_hc %>% 
    select(contig, start, stop, variants, haplotypeA, haplotypeB) %>%
    separate_wider_delim(variants, ',',
                         names=paste("variant", 1:(max_nvar_hc+1), sep="_"),
                         too_few = "align_start") %>%
    separate_wider_delim(haplotypeA, ',',
                         names=paste("blockHapA", 1:(max_nvar_hc+1), sep="_"),
                         too_few = "align_start") %>%
    separate_wider_delim(haplotypeB, ',',
                         names=paste("blockHapB", 1:(max_nvar_hc+1), sep="_"),
                         too_few = "align_start")

df_hc_wl <- df_hc_w %>% pivot_longer(4:ncol(df_hc_w), 
                      names_to = c(".value", "haplo_set"),
                      names_pattern = "(.*)_([0-9]+)") %>%
        filter(!is.na(variant))

df_hc_wl2 <- df_hc_wl %>%
    inner_join(df_hc) %>%
    select(contig, start, stop, haplo_set, variant, blockHapA, blockHapB, blockGWPhase, gwStat) %>%
    rename(blockStart = start,
           blockStop = stop,
           variantID = variant)

df_gae1_wlc_hc <- df_gae1_wlc %>% inner_join(df_hc_wl2)


df_gae1_wlc_hc2 <- df_gae1_wlc_hc %>%
    mutate(blockHapAcnt = ifelse(blockHapA == refAllele, refCount, altCount),
           blockHapBcnt = ifelse(blockHapB == refAllele, refCount, altCount))

df_gae1_wlc_hc3 <- df_gae1_wlc_hc2 %>%
    mutate(hapAcnt = ifelse(gwStat < 0.9, blockHapAcnt,
                            ifelse(blockGWPhase == "0|1", blockHapAcnt, blockHapBcnt)),
           hapBcnt = ifelse(gwStat < 0.9, blockHapBcnt,
                            ifelse(blockGWPhase == "0|1", blockHapBcnt, blockHapAcnt)),
           hapA = ifelse(gwStat < 0.9, blockHapA,
                         ifelse(blockGWPhase == "0|1", blockHapA, blockHapB)),
           hapB = ifelse(gwStat < 0.9, blockHapB,
                         ifelse(blockGWPhase == "0|1", blockHapB, blockHapA)))

df_gae1_wlc_final_long <- df_gae1_wlc_hc3 %>%
                           arrange(contig, position) %>%
                           select(contig, start, stop, name, 
                           variantID, hapA, hapB, hapAcnt, hapBcnt)

df_gae1_wlc_final_collapsed <- df_gae1_wlc_final_long %>%
    group_by(contig, start, stop, name) %>%
    summarise(variants=paste0(variantID, collapse=','),
              varHapA=paste0(hapA, collapse=','),
              varHapB=paste0(hapB, collapse=','),
              varCntA=paste0(hapAcnt, collapse=','),
              varCntB=paste0(hapBcnt, collapse=','))

write.table(df_gae1_wlc_final_collapsed, outfn, row.names=FALSE, quote=FALSE, sep="\t")
message(paste("Outputted to", outfn))