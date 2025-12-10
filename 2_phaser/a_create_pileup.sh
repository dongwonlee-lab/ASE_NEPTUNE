#!/bin/bash
#SBATCH --job-name=dna_rna_discord            
#SBATCH --output=logs_slurm/%j.txt
#SBATCH --mail-type=FAIL,END            
#SBATCH --time=23:59:59
#SBATCH --ntasks=1 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

#' Create RNA Pileup data
#'
#' For each variant in the individualized VCF, gathers
#' the alleles that mapped to that location in RNA-seq
#' using samtools mpileup. That data is used to determine
#' loci with discordant DNA and RNA reads, which
#' could impact ASE results in all individuals.
#' WARNING: This must be run for all samples/individuals before 
#'          moving to step (b)
#'
#' Output: 
#' - ${SUBJECT_ID}.discord_snps_gts.tsv:
#'   GTs for all variants in individual's unphased VCF
#' - ${SUBJECT_ID}.discord_snps_pileup.tsv
#'   samtools mpileup for all variants in above
#'
#' Usage:
#'   ./a_create_pileup.sh
#'
#' Author:
#'   Eric Sakkas
#'   erotokritos.sakkas@childrens.harvard.edu
#'
#' Date:
#'   2025-09-10

### Gather GTs for all variants in individual VCF — 1CPU thread and 1GB ###
SNPS_FILE=../resource/gencode.v.45.gene_coordinates.bed
SUBJECT_ID=HG00513
BAM_DIR=../data/output/rna_seq/${SUBJECT_ID}
UNPHASED_VCF=../data/${SUBJECT_ID}/${SUBJECT_ID}.unphased.vcf.gz

OUTDIR=../data/output/phaser/${SUBJECT_ID}
mkdir -p $OUTDIR
if [ ! -d "$OUTDIR" ]; then
    echo "no outdir: creating "$OUTDIR
    mkdir -p $OUTDIR
fi
GTS_FILE=${OUTDIR}/${SUBJECT_ID}.discord_snps_gts.tsv

echo "Writing variant GTs to" $GTS_FILE
echo -e "\tDATE, TIME : $(date '+%Y-%m-%d, %H:%M:%S')"

bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' -o $GTS_FILE -s $SUBJECT_ID -R $SNPS_FILE $UNPHASED_VCF

### samtools mpileup at all vars — 1CPU thread and 1GB ###

BAM_FILE=${BAM_DIR}/${SUBJECT_ID}.starwasp2.7.11b.quantmode.Aligned.sortedByCoord.out.bam
PILE_OUTFILE=${OUTDIR}/${SUBJECT_ID}.discord_snps_pileup.tsv

echo "Writing pileup output to" $PILE_OUTFILE
echo -e "\tDATE, TIME : $(date '+%Y-%m-%d, %H:%M:%S')"

samtools mpileup -Q 10 -q 255 -l "$GTS_FILE" "$BAM_FILE" -o "$PILE_OUTFILE"