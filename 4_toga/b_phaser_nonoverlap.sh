#!/bin/bash
#SBATCH --job-name=phaser
#SBATCH --output=logs_slurm/%j.txt   
#SBATCH --mail-type=FAIL,END  
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=10
#SBATCH --time=23:59:59
#SBATCH --mem=20G

source activate /lab-share/Neph-Lee-e2/Public/esakkas/conda/envs/phaser

#' Runs phASER with gene overlapping regions blackisted
#' 
#' This script creates a new phASER blacklist that removes all
#' overlapping gene regions in the individual and reruns phASER.
#' The output files are used by TOGA algorithm.
#'
#' Outputs: 
#' - In ../data/output/toga/${SAMPLE}/phaser_nonoverlap:
#'      Same as 2_phaser/d_phASER_to_geneAE.sh
#'
#' Usage:
#'   ./b_phaser_nonoverlap.sh
#'
#' Author:
#'   Junmo Sung
#'   jsung5@bwh.harvard.edu
#'
#' Date:
#'   2025-09-25

####### INPUTS #######
SAMPLE=HG00513
######################

echo 'Generating merged haplo-blacklist for' $SAMPLE

HAPLO_BLACKLIST="../resource/hg38_haplo_count_blacklist.chr.bed"  # Update if a different file was used in Step 1
SHARED_VAR_BLACKLIST=../data/output/toga/${SAMPLE}/${SAMPLE}.cov.overlapping_region.bed
MERGED_HAPLO_BLACKLIST=../data/output/toga/${SAMPLE}/${SAMPLE}.merged_output_haplo_blacklist.bed

cat $HAPLO_BLACKLIST $SHARED_VAR_BLACKLIST | sort -k1,1 -k2,2n | bedtools merge > $MERGED_HAPLO_BLACKLIST

echo 'Merged blacklist written to' $MERGED_HAPLO_BLACKLIST

echo 'Running phASER...'
./phaser_to_geneAE.sh \
  --sample ${SAMPLE} \
  --bam ../data/output/rna_seq/${SAMPLE}/${SAMPLE}.starwasp2.7.11b.quantmode.Aligned.sortedByCoord.out.samtools_markdup.wasp_vW1.bam \
  --vcf ../data/${SAMPLE}/${SAMPLE}.het.phased.vcf.gz \
  --o ../data/output/toga/phaser_nonoverlap \
  --phaser ../../phaser/phaser/phaser.py \
  --phaser_geneAE ../../phaser/phaser_gene_ae/phaser_gene_ae.py \
  --haplo_count_blacklist $MERGED_HAPLO_BLACKLIST
echo 'Done!'