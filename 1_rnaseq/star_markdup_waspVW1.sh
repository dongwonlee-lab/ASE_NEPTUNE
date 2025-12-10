#!/bin/bash
#SBATCH --job-name=star_wasp
#SBATCH --output=logs_slurm/%j.txt   
#SBATCH --mail-type=FAIL,END   
#SBATCH --ntasks=1  
#SBATCH --ntasks-per-node=1  
#SBATCH --cpus-per-task=10  
#SBATCH --time=23:59:59   
#SBATCH --mem=20G

#' Run STAR (v.2.7.11b) to align RNA reads and WASP-correct
#'
#' This module executes STAR with WASP correction to align
#' the individual's RNA reads from 2 FASTQ files to the reference
#' sequence and output a sorted, WASP-corrected BAM file with 
#' duplicates marked.
#'
#' Usage:
#'   ./star_markdup_waspVW1.sh
#'
#' Output:
#'   ${SAMPLE}.starwasp2.7.11b.quantmode.Aligned.sortedByCoord.out.samtools_markdup.wasp_vW1.bam:
#'       BAM file containing aligned RNA-seq reads that pass WASP-filtering (vW1). Duplicates marked.
#'
#' Author:
#'   Junmo Sung and Eric Sakkas
#'   jsung5@bwh.harvard.edu
#'
#' Date:
#'   2025-09-11

REFERENCE_FASTA=../resource/GRCh38.primary_assembly.genome.fa
REFERENCE_GTF=../resource/gencode.v45.primary_assembly.annotation.gtf

GENOME_INDEX=../resource/reference_genome_index

# Generate REFERENCE index — worked 10 CPU threads and 30 GB
STAR --runThreadN 10 \
    --runMode genomeGenerate \
    --genomeDir $GENOME_INDEX \
    --genomeFastaFiles $REFERENCE_FASTA \
    --sjdbGTFfile $REFERENCE_GTF \
    --sjdbOverhang 100 \
    --genomeSAsparseD 3

### Run STAR(v.2.7.11b) with WASP ###

# Sample ID
SAMPLE=HG00513
DATA_DIR=../data/${SAMPLE}

# Pair-end FASTQ files
FASTQ_R1=${DATA_DIR}/SRR19762919_1.fastq.gz
FASTQ_R2=${DATA_DIR}/SRR19762919_2.fastq.gz

## Create individual-level VCFs-het and non-het filtered
# 1000 Genomes provides file 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz
# Unphased VCF files for 1KG individual HG00513 were created with code:
## for chr in {1..22}; do
##         bcftools view -s HG00513 \
##                 -Oz \
##                 -o HG00513.chr${chr}.unphased.vcf.gz \
##                 --write-index=tbi \
##                 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.recalibrated_variants.vcf.gz
##
##         bcftools view -g het \
##                 -Oz \
##                 -o HG00513.chr${chr}.het.unphased.vcf.gz \
##                 --write-index=tbi \
##                 HG00513.chr${chr}.unphased.vcf.gz
## done
## bcftools concat \
##        $(for i in $(seq 1 22); do echo HG00513.chr${i}.unphased.vcf.gz; done) \
##        -Oz -o HG00513.unphased.vcf.gz 
## bcftools concat \
##        $(for i in $(seq 1 22); do echo HG00513.chr${i}.het.unphased.vcf.gz; done) \
##        -Oz -o HG00513.het.unphased.vcf.gz 
## bcftools index -t HG00513.unphased.vcf.gz
## bcftools index -t HG00513.het.unphased.vcf.gz

# Repeat the above process to create VCF files for other samples, and with the 1KG
# file CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
# to create the population phased het VCFs: ${SAMPLE}.het.phased.vcf.gz

# Phased het VCF
VCF=${DATA_DIR}/${SAMPLE}.het.phased.vcf.gz

OUTDIR=../data/output/rna_seq/${SAMPLE}
if [ ! -d "$OUTDIR" ]; then
    echo "no outdir: creating "$OUTDIR
    mkdir -p $OUTDIR
fi

# Run STAR-WASP — worked 10 CPU threads and 20 GB
STAR --runThreadN 10 \
--genomeDir $GENOME_INDEX \
--readFilesCommand gunzip -c \
--readFilesIn $FASTQ_R1 $FASTQ_R2 \
--outFileNamePrefix $OUTDIR/${SAMPLE}.starwasp2.7.11b.quantmode. \
--outSAMattributes All \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 9 \
--outBAMcompression 10 \
--quantMode TranscriptomeSAM \
--twopassMode Basic \
--waspOutputMode SAMtag \
--varVCFfile <(zcat $VCF) 

# Index output BAM file
samtools index -b $OUTDIR/${SAMPLE}.starwasp2.7.11b.quantmode.Aligned.sortedByCoord.out.bam

### Mark duplicate alignments ###

# Set input/ouput files
BAM=${OUTDIR}/${SAMPLE}.starwasp2.7.11b.quantmode.Aligned.sortedByCoord.out.bam
BAM_COLLATE=${OUTDIR}/${SAMPLE}.samtools_tempcollate.bam
BAM_FIXMATE=${OUTDIR}/${SAMPLE}.samtools_tempfixmate.bam
BAM_SORT=${OUTDIR}/${SAMPLE}.samtools_tempsort.bam
BAM_MARKDUP=${OUTDIR}/${SAMPLE}.starwasp2.7.11b.quantmode.Aligned.sortedByCoord.out.samtools_markdup.bam
STAT_FILE=${OUTDIR}/${SAMPLE}.starwasp2.7.11b.quantmode.Aligned.sortedByCoord.out.samtools_markdup_stats.txt

# samtools collate – shuffles and groups reads together by their names
samtools collate -o $BAM_COLLATE --threads 8 $BAM

# samtools fixmate – fills in mate coordinates and insert size fields.
samtools fixmate -m $BAM_COLLATE --threads 8 $BAM_FIXMATE

# samtools sort – sorts SAM/BAM/CRAM files
samtools sort -o $BAM_SORT --threads 8 $BAM_FIXMATE

# samtools markdup – mark duplicate alignments in a coordinate sorted file
samtools markdup -s -f $STAT_FILE $BAM_SORT --threads 8 $BAM_MARKDUP

# Index the output
samtools index -b $BAM_MARKDUP

###  Select only reads with WASP vW1 tags ###

BAM_MARKDUP_VW1=${OUTDIR}/${SAMPLE}.starwasp2.7.11b.quantmode.Aligned.sortedByCoord.out.samtools_markdup.wasp_vW1.bam 
samtools view -d vW:1 -o $BAM_MARKDUP_VW1 $BAM_MARKDUP

# Index the output
samtools index -b $BAM_MARKDUP_VW1

# Delete intermediate files
rm -f "$BAM_COLLATE" "$BAM_FIXMATE" "$BAM_SORT" "$BAM_MARKDUP"