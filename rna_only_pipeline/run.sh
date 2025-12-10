#!/bin/bash
#SBATCH --job-name=star_wasp
#SBATCH --output=logs_slurm/%j.txt   
#SBATCH --mail-type=FAIL,END  
#SBATCH --ntasks=1  
#SBATCH --ntasks-per-node=1  
#SBATCH --cpus-per-task=10  
#SBATCH --time=23:59:59   
#SBATCH --mem=50G

#' Run RNA-only STAR pipeline
#'
#' Uses only RNA-seq files (FASTQ) to generate BAM files
#' of aligned RNA reads with STAR and a VCF of the predicted
#' genome using GATK.
#'
#' Outputs:
#' - ${SAMPLE}.starwasp.Aligned.sortedByCoord.out.samtools_markdup.vW1.bam:
#'      BAM of aligned RNA reads
#' - ${SAMPLE}.dp10.woINDEL.woMA.autosomal.het.variant_filtered.vcf.gz:
#'      VCF of predicted genome, filtered to het variants on chr1-22,
#'      no multiallelic variants, only variants supported by >10 reads.
#'
#' Usage:
#'   ./run.sh
#'
#' Author:
#'   Junmo Sung
#'   jsung5@bwh.harvard.edu
#'
#' Date:
#'   2025-09-24

set -euo pipefail

####### INPUTS #######
# To Generate Reference Genome
REFERENCE_FASTA=../resource/GRCh38.primary_assembly.genome.fa
REFERENCE_GTF=../resource/gencode.v45.primary_assembly.annotation.gtf
GENOME_INDEX=../resource/reference_genome_index

# To run STAR
SAMPLE=HG00171
FASTQ_R1=../data/${SAMPLE}/ERR188190_1.fastq.gz
FASTQ_R2=../data/${SAMPLE}/ERR188190_2.fastq.gz
OUTDIR=../data/output/rna_only
######################

### Pre-process: Generate reference genome index before running STAR (If required) ###
echo "Preprocess: Generating reference genome index before running STAR... "

STAR --runThreadN 10 \
    --runMode genomeGenerate \
    --genomeDir $GENOME_INDEX \
    --genomeFastaFiles $REFERENCE_FASTA \
    --sjdbGTFfile $REFERENCE_GTF \
    --sjdbOverhang 100 \
    --genomeSAsparseD 3

echo "######## Step 1: Running STAR without BAM ########"
echo -e "\tDATE, TIME : $(date '+%Y-%m-%d, %H:%M:%S')"

echo 'Running... '
./star.sh \
    ${SAMPLE} \
    ${FASTQ_R1} \
    ${FASTQ_R2} \
    ${GENOME_INDEX} \
    ${OUTDIR}/star # will automatically
echo 'Done!'

echo "######## Step2: Run GATK Pipeline to call variants from RNA-seq reads ########"
echo -e "\tDATE, TIME : $(date '+%Y-%m-%d, %H:%M:%S')"
INPUT_BAM=${OUTDIR}/star/${SAMPLE}/${SAMPLE}.star.Aligned.sortedByCoord.out.samtools_markdup.bam
BAM_OUTDIR=${OUTDIR}/gatk/bam
VCF_OUTDIR=${OUTDIR}/gatk/vcf

TEMP_DIR=${OUTDIR}/tmp
if [ ! -d "$TEMP_DIR" ]; then
    echo "no tempdir: creating "$TEMP_DIR
    mkdir -p $TEMP_DIR
fi

dbSNP_INDEL_IL_DIR=../resource/dbsnps_indel_interval
# [user]$ ls ../resource/dbsnps_indel_interval
#
# Homo_sapiens_assembly38.dbsnp138.vcf.gz      Homo_sapiens_assembly38.known_indels.vcf.gz      Mills_and_1000G_gold_standard.indels.hg38.vcf.gz      wgs_calling_regions.hg38.interval_list
# Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi  Homo_sapiens_assembly38.known_indels.vcf.gz.tbi  Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

echo 'Running... '
./rna_seq.variant_call.sh \
    ${SAMPLE} \
    ${INPUT_BAM} \
    ${BAM_OUTDIR} \
    ${VCF_OUTDIR} \
    ${REFERENCE_FASTA} \
    ${dbSNP_INDEL_IL_DIR} \
    ${TEMP_DIR}
echo 'Done!'

echo "######## Step3: Run STAR-WASP using RNA-seq based vcf file ########"
echo -e "\tDATE, TIME : $(date '+%Y-%m-%d, %H:%M:%S')"
VCF=${VCF_OUTDIR}/${SAMPLE}/${SAMPLE}.autosomal.het.variant_filtered.vcf.gz
STARWASP_OUTDIR=${OUTDIR}/starwasp

echo 'Running... '
./star_wasp.sh \
    ${SAMPLE} \
    ${FASTQ_R1} \
    ${FASTQ_R2} \
    ${VCF} \
    ${GENOME_INDEX} \
    ${STARWASP_OUTDIR}
echo 'Done!'

echo "######## Re-run GATK Pipeline to call variants from wasp-filtered RNA-seq reads ########"
INPUT_BAM=${STARWASP_OUTDIR}/${SAMPLE}/${SAMPLE}.starwasp.Aligned.sortedByCoord.out.samtools_markdup.vW1.bam
BAM_OUTDIR=${OUTDIR}/gatk/bam_vw1
VCF_OUTDIR=${OUTDIR}/gatk/vcf_vw1

echo 'Running... '
./rna_seq.variant_call.sh \
    ${SAMPLE} \
    ${INPUT_BAM} \
    ${BAM_OUTDIR} \
    ${VCF_OUTDIR} \
    ${REFERENCE_FASTA} \
    ${dbSNP_INDEL_IL_DIR} \
    ${TEMP_DIR}
echo 'Done!'

echo "######## Step5: Remove variants with DP less than 10 and those that are multiallelic ########"
AUTOSOMAL_HET_VCF=${VCF_OUTDIR}/${SAMPLE}/${SAMPLE}.autosomal.het.variant_filtered.vcf.gz
FINAL_VCF=${VCF_OUTDIR}/${SAMPLE}/${SAMPLE}.dp10.woINDEL.woMA.autosomal.het.variant_filtered.vcf.gz

# Remove variants that are multiallelic, DP < 10, indels 
echo 'Removing variants with DP less than 10 and those that are multiallelic and indels '
bcftools view -m2 -M2 -v snps -i 'FORMAT/DP>=10' $AUTOSOMAL_HET_VCF -Oz -o $FINAL_VCF

tabix -p vcf $FINAL_VCF
echo 'Done!'

# We recommend the following phASER run:
#
# python3 $PHASER_SCRIPT \
#     --bam $BAM \ <- ${SAMPLE}.starwasp.Aligned.sortedByCoord.out.samtools_markdup.vW1.bam
#     --vcf $VCF \ <- ${SAMPLE}.dp10.woINDEL.woMA.autosomal.het.variant_filtered.vcf.gz
#     --sample $SAMPLE \
#     --baseq 10 \
#     --mapq 255 \
#     --paired_end 1 \
#     --o $OUTPREFIX \
#     --haplo_count_blacklist $HAPLO_BLACKLIST \
#     --blacklist $GENE_BLACKLIST \
#     --write_vcf 1 \
#     --remove_dups 1 \
#     --pass_only 0 \
#     --unphased_vars 1 \
#     --gw_phase_method 0 
#
# We used the default `gencode45_hg38_hla.chr.bed` from phASER for $GENE_BLACKLIST
# and for $HAPLO_BLACKLIST, combined phASER's `hg38_haplo_count_blacklist.chr.bed`
# with all monoallelic variants and variants not found in 1000 Genomes