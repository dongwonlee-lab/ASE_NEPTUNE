#!/bin/bash
#SBATCH --job-name=phaser
#SBATCH --output=logs_slurm/%j.txt   
#SBATCH --mail-type=FAIL,END  
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=10
#SBATCH --time=23:59:59
#SBATCH --mem=50G

#' Run phASER and phASER Gene AE
#'
#' This module runs phASER and phASER Gene AE using the WGS 
#' and aligned RNA-seq, assigning RNA reads to the appropriate haplotype, 
#' and outputting allele- and haplotype-level counts for all
#' het variants. The program uses the aligned RNA-seq and individual-level
#' blacklists as input.
#'
#' Output:
#' - ../data/output/phaser/phaser_run/${SUBJECT_ID}*:
#'      phASER output, allele and haplotype counts
#' - ../data/output/phaser/phaser_run/${SUBJECT_ID}.gene_ae.txt:
#'      Gene AE output, gene-level haplotype counts
#'
#' Usage:
#'   ./d_phASER_to_geneAE.sh
#'
#' Documentation:
#'   https://github.com/secastel/phaser/tree/master/phaser
#'
#' Author:
#'   Junmo Sung and Eric Sakkas
#'   erotokritos.sakkas@childrens.harvard.edu
#'
#' Date:
#'   2025-09-16

# Environment for phASER; dependencies Cython, scipy, pysam, pandas, intervaltree
# phASER setup.py ran on this python environment
source activate /lab-share/Neph-Lee-e2/Public/esakkas/conda/envs/phaser

### phASER ###
# Assumes local copy of phASER repository in the same directory as our 
# ASE repository. phASER: https://github.com/secastel/phaser/
PHASER_DIR=../../phaser
PHASER_SCRIPT=${PHASER_DIR}/phaser/phaser.py

## phASER inputs ##
SAMPLE=HG00513
VCF=../data/${SAMPLE}/${SAMPLE}.het.phased.vcf.gz

# Aligned RNA, output of X
BAM=../data/output/rna_seq/${SAMPLE}/${SAMPLE}.starwasp2.7.11b.quantmode.Aligned.sortedByCoord.out.samtools_markdup.wasp_vW1.bam

# Regions blacklisted from RNA phasing and haplotype counting, output of X
GENE_BLACKLIST=../data/output/phaser/${SAMPLE}/${SAMPLE}.final_blacklist.bed

# Regions blacklisted from haplotype counting.
# Provided by phASER (https://www.dropbox.com/s/9cn9477bcutvuc7/hg38_haplo_count_blacklist.chr.bed.gz)
HAPLO_BLACKLIST=../resource/hg38_haplo_count_blacklist.chr.bed

## phASER outputs ##
OUTDIR=../data/output/phaser/${SAMPLE}/phaser_run
if [ ! -d "$OUTDIR" ]; then
    echo "no outdir: creating "$OUTDIR
    mkdir -p $OUTDIR
fi
OUTPREFIX=${OUTDIR}/${SAMPLE}

## Running phASER ##
python3 $PHASER_SCRIPT \
    --bam $BAM \
    --vcf $VCF \
    --sample $SAMPLE \
    --baseq 10 \
    --mapq 255 \
    --paired_end 1 \
    --o $OUTPREFIX \
    --haplo_count_blacklist $HAPLO_BLACKLIST \
    --blacklist $GENE_BLACKLIST \
    --write_vcf 1 \
    --remove_dups 1 \
    --pass_only 0 \
    --unphased_vars 1 \
    --gw_phase_method 1 \
    --gw_phase_vcf 1 \
    --gw_af_field 'AF' \
    --gw_phase_vcf_min_confidence 0.90 \
    --threads 10 \
    --show_warning 1 \
    --debug 1 

echo 'phASER results outputed to' $OUTDIR

### phASER GENE AE ###
GENE_AE_SCRIPT=${PHASER_DIR}/phaser_gene_ae/phaser_gene_ae.py

## Gene AE Inputs ##
# Haplotypic counts from phASER output
HAP_COUNTS=${OUTPREFIX}.haplotypic_counts.txt
# Gene feature used for estimating gene level counts
FEATURES=../resource/gencode.v.45.gene_coordinates.bed
OUTPUT=${OUTPREFIX}.gene_ae.txt

python3 $GENE_AE_SCRIPT \
    --haplotypic_counts $HAP_COUNTS \
    --features $FEATURES \
    --o $OUTPUT

echo 'Gene AE results outputed to' $OUTDIR
echo 'phASER and Gene AE done'