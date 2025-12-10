#!/bin/bash

# Set-up
SAMPLE=$1
FASTQ_R1=$2
FASTQ_R2=$3
VCF=$4
GENOME_INDEX=$5
OUTDIR=$6

# Create output directory
if [ ! -d "$OUTDIR" ]; then
    echo "no outdir: creating "$OUTDIR
    mkdir -p $OUTDIR
fi
OUTDIR_SAMPLE=${OUTDIR}/${SAMPLE}
if [ ! -d "$OUTDIR_SAMPLE" ]; then
    echo "no outdir: creating "$OUTDIR_SAMPLE
    mkdir -p $OUTDIR_SAMPLE
fi

echo "SAMPLE_ID: " $SAMPLE
echo "FASTQ: " ${FASTQ_R1} ", " ${FASTQ_R2} 
echo "VCF: " $VCF
echo "GENOME_INDEX: " $GENOME_INDEX
echo "OUTDIR: "$OUTDIR
echo "OUTDIR_SAMPLE: "$OUTDIR_SAMPLE
echo " "

## STAR with WASP ##
echo "Runnning STAR-WASP..."
STAR --runThreadN 10 \
--genomeDir $GENOME_INDEX \
--readFilesCommand gunzip -c \
--readFilesIn $FASTQ_R1 $FASTQ_R2 \
--outFileNamePrefix $OUTDIR_SAMPLE/${SAMPLE}.starwasp. \
--outSAMattributes All \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 9 \
--outBAMcompression 10 \
--twopassMode Basic \
--waspOutputMode SAMtag \
--varVCFfile <(zcat $VCF) 

# Index output BAM file
BAM=${OUTDIR_SAMPLE}/${SAMPLE}.starwasp.Aligned.sortedByCoord.out.bam
samtools index -b $BAM
echo "BAM created and indexed: " $BAM
echo " "

## Samtools markdup ##
echo "Running samtools markdup... "
BAM_COLLATE=${OUTDIR_SAMPLE}/${SAMPLE}.starwasp.Aligned.sortedByCoord.out.samtools_tempcollate.bam
BAM_FIXMATE=${OUTDIR_SAMPLE}/${SAMPLE}.starwasp.Aligned.sortedByCoord.out.samtools_tempfixmate.bam
BAM_SORT=${OUTDIR_SAMPLE}/${SAMPLE}.starwasp.Aligned.sortedByCoord.out.samtools_tempsort.bam
BAM_MARKDUP=${OUTDIR_SAMPLE}/${SAMPLE}.starwasp.Aligned.sortedByCoord.out.samtools_markdup.bam
STAT_FILE=${OUTDIR_SAMPLE}/${SAMPLE}.starwasp.Aligned.sortedByCoord.out.samtools_markdup_stats.txt

echo " "
echo "BAM is "$BAM
echo "BAM_COLLATE is "$BAM_COLLATE
echo "BAM_FIXMATE is "$BAM_FIXMATE
echo "BAM_SORT is "$BAM_SORT
echo "BAM_MARKDUP is "$BAM_MARKDUP
echo "STAT_FILE is "$STAT_FILE
echo " "

# samtools collate – shuffles and groups reads together by their names
samtools collate -o $BAM_COLLATE --threads 8 $BAM

# samtools fixmate – fills in mate coordinates and insert size fields.
# -m: Add ms (mate score) tags. These are used by markdup to select the best reads to keep.
samtools fixmate -m $BAM_COLLATE --threads 8 $BAM_FIXMATE

# samtools sort – sorts SAM/BAM/CRAM files
samtools sort -o $BAM_SORT --threads 8 $BAM_FIXMATE

# samtools markdup – mark duplicate alignments in a coordinate sorted file
samtools markdup -s -f $STAT_FILE $BAM_SORT --threads 8 $BAM_MARKDUP

# Index $BAM_MARKDUP
samtools index -b $BAM_MARKDUP
echo "Markdup BAM created and indexed: "$BAM_MARKDUP
echo " "

# Delete intermediate files
rm -f "$BAM_COLLATE" "$BAM_FIXMATE" "$BAM_SORT"

## vW1 filtered ##
echo "vW:1 BAM file generating... "
BAM_VW1=${OUTDIR_SAMPLE}/${SAMPLE}.starwasp.Aligned.sortedByCoord.out.samtools_markdup.vW1.bam
samtools view -d vW:1 -o $BAM_VW1 $BAM_MARKDUP

# Index $BAM_VW1
samtools index -b $BAM_VW1
echo "vW1 filtered BAM created and indexed: "$BAM_VW1
echo " "

# Delete intermediate files
rm -f "$BAM_COLLATE" "$BAM_FIXMATE" "$BAM_SORT" 