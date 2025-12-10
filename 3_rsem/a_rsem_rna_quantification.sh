#!/bin/bash
#SBATCH --job-name=rsem
#SBATCH --output=logs_slurm/%j.txt   
#SBATCH --mail-type=FAIL,END   
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=10
#SBATCH --time=23:59:59
#SBATCH --mem=5G

#' Runs RSEM to create normalized expression output
#' 
#' This module uses STAR-WASP aligned BAMs to determine
#' the proportion of certain gene isoforms that each sample
#' expressed. This is used to resolve ASE in overlapping gene
#' regions, and appropriately assign the reads between genes.
#'
#' Outputs: (more information in `rsem-calculate-expression --help`).
#' - sample_name.isoforms.results:
#'      tab-separated file of isoform-level expression estimates.
#' - sample_name.genes.results:
#'      tab-separated file of gene-level expression estimates.
#' - sample_name.stat:
#'      folder of model-related statistics.
#'
#' Usage:
#'   ./a_rsem_rna_quantification.sh
#'
#' Author:
#'   Junmo Sung
#'   jsung5@bwh.harvard.edu
#'
#' Date:
#'   2025-09-25

####### INPUTS #######
REFERENCE=../resource/GRCh38.primary_assembly.genome.fa
GTF=../resource/gencode.v45.primary_assembly.annotation.gtf
REFERENCE_RSEM=../resource/ref/hg38
# [user]$ ls ../resource/ref/
# hg38.chrlist  hg38.grp  hg38.idx.fa  hg38.n2g.idx.fa  hg38.seq  hg38.ti  hg38.transcripts.fa

SAMPLE=HG00513
BAM=../data/output/rna_seq/${SAMPLE}/${SAMPLE}.starwasp2.7.11b.quantmode.Aligned.toTranscriptome.out.bam
OUTDIR=../data/output/rsem
######################

if [ ! -d "$OUTDIR" ]; then
    echo "no directory: creating "$OUTDIR
    mkdir -p $OUTDIR
fi

echo 'Building RSEM reference...'
rsem-prepare-reference \
    --gtf $GTF \
    $REFERENCE \
    $REFERENCE_RSEM

echo 'Estimating isoform expression...'
rsem-calculate-expression --paired-end \
    --alignments \
    -p 10 \
    --estimate-rspd \
    --append-names \
    --no-bam-output \
    $BAM \
    $REFERENCE_RSEM \
    ${OUTDIR}/${SAMPLE}

echo 'Normalized gene expression file written to' ${OUTDIR}/${SAMPLE}.tsv
