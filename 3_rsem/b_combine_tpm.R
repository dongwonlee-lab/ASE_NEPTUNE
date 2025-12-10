#' Gather TPM data for all samples
#'
#' WARNING: This can only be used once all samples have been processed
#' by `a_rsem_rna_quantification.sh`.
#' 
#' Gathers the transcripts per million (TPM) data for each sample
#' into a single file.
#'
#' Outputs:
#' - all_samples.tpm.tsv:
#'   - gene_id: Identity of the transcriptid_geneid
#'   - SAMPLE_X: TPM for given gene_id in sample
#'
#' Usage:
#'   Rscript b_combine_tpm.R
#'
#' Dependencies:
#'   R ≥ 3.5.x with packages: dplyr, tidyr
#'
#' Author:
#'   Junmo Sung
#'   jsung5@bwh.harvard.edu
#'
#' Date:
#'   2025-09-25

library(dplyr)
library(tidyr)

subject_ids <- readLines("../resource/subject_ids.tsv")
print(subject_ids)

all_subj_rsem <- data.frame()

for (subject_id in subject_ids) {
    print(paste("Reading files for", subject_id))

    file_path <- paste0("../data/output/rsem/", subject_id, ".genes.results")

    sample_rsem <- read.table(file_path, header = TRUE, sep = "\t") %>%
        mutate(SUBJECT_ID = subject_id) %>%
        dplyr::select(gene_id, SUBJECT_ID, TPM)

    all_subj_rsem <- rbind(all_subj_rsem, sample_rsem)
}

all_subj_tpm <- all_subj_rsem %>% pivot_wider(names_from = SUBJECT_ID, values_from = TPM)

outfile <- "../data/output/rsem/all_samples.tpm.tsv"
write.table(all_subj_tpm, outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

print(paste("TPM for all samples written to", outfile))