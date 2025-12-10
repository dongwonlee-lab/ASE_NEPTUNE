#!/bin/bash

SAMPLE=$1
INPUT_BAM=$2
BAM_OUTDIR=$3
VCF_OUTDIR=$4
REFERENCE_FASTA=$5
dbSNP_INDEL_IL_DIR=$6
TEMP_DIR=$7

# Known snps, indels, interval lists to consider when calling variants
dbSNP=$dbSNP_INDEL_IL_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.gz
dbSNP_indels=$dbSNP_INDEL_IL_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz
dbSNP_indels_2=$dbSNP_INDEL_IL_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
INTERVAL_LIST=$dbSNP_INDEL_IL_DIR/wgs_calling_regions.hg38.interval_list

# Create output directories for bam and vcf files
SAMPLE_BAM_OUTDIR=${BAM_OUTDIR}/${SAMPLE}
if [ ! -d "$SAMPLE_BAM_OUTDIR" ]; then
    mkdir -p $SAMPLE_BAM_OUTDIR
fi

SAMPLE_VCF_OUTDIR=${VCF_OUTDIR}/${SAMPLE}
if [ ! -d "$SAMPLE_VCF_OUTDIR" ]; then
    mkdir -p $SAMPLE_VCF_OUTDIR
fi

# Variables or Directory required for each step
SPLITNCIGARREADS_BAM=${SAMPLE_BAM_OUTDIR}/${SAMPLE}.splitncigar.bam 
FLAG_FILTERED_BAM=${SAMPLE_BAM_OUTDIR}/${SAMPLE}.splitncigar.flag_filtered.bam 
ADDREADGROUPS_BAM=${SAMPLE_BAM_OUTDIR}/${SAMPLE}.add_read_groups.bam 
RECAL=${SAMPLE_BAM_OUTDIR}/${SAMPLE}.recal.csv
APPLYBQSR_BAM=${SAMPLE_BAM_OUTDIR}/${SAMPLE}.apply_bqsr.bam
HAPLOTYPECALLER_VCF=${SAMPLE_VCF_OUTDIR}/${SAMPLE}.vcf.gz
VARIANTFILTRATION_VCF=${SAMPLE_VCF_OUTDIR}/${SAMPLE}.variant_filtered.vcf.gz
AUTOSOMAL_HET_VCF=${SAMPLE_VCF_OUTDIR}/${SAMPLE}.autosomal.het.variant_filtered.vcf.gz

# SplitNCigarReads
gatk SplitNCigarReads \
  -I $INPUT_BAM \
  -R $REFERENCE_FASTA \
  -O $SPLITNCIGARREADS_BAM \
  --tmp-dir $TEMP_DIR
samtools index $SPLITNCIGARREADS_BAM

# Filtered out unpaired reads that have mate unmmaped flag
samtools view -h $SPLITNCIGARREADS_BAM | grep -v -F -f <(samtools view -f 8 -F 1 $SPLITNCIGARREADS_BAM | cut -f1) | samtools view -b > $FLAG_FILTERED_BAM
samtools index $FLAG_FILTERED_BAM

# AddOrReplaceReadGroups
picard AddOrReplaceReadGroups \
    I=$FLAG_FILTERED_BAM \
    O=$ADDREADGROUPS_BAM \
    RGID=1 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=$SAMPLE
samtools index $ADDREADGROUPS_BAM

# BaseRecalibraotr
gatk BaseRecalibrator \
    -I $ADDREADGROUPS_BAM \
    -O $RECAL \
    -R $REFERENCE_FASTA \
    --use-original-qualities \
    -known-sites $dbSNP \
    -known-sites $dbSNP_indels \
    -known-sites $dbSNP_indels_2

# ApplyBQSR
gatk ApplyBQSR \
    -bqsr $RECAL \
    -I $ADDREADGROUPS_BAM \
    -O $APPLYBQSR_BAM \
    -R $REFERENCE_FASTA \
    --add-output-sam-program-record \
    --use-original-qualities
samtools index $APPLYBQSR_BAM

# HaplotypeCaller
gatk HaplotypeCaller \
    -I $APPLYBQSR_BAM \
    -O $HAPLOTYPECALLER_VCF \
    -R $REFERENCE_FASTA \
    -L $INTERVAL_LIST \
    -dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20 \
    --dbsnp $dbSNP \
    --tmp-dir $TEMP_DIR
tabix -p vcf $HAPLOTYPECALLER_VCF

# VariantFiltration
gatk VariantFiltration \
    -V $HAPLOTYPECALLER_VCF \
    -O $VARIANTFILTRATION_VCF \
    -R $REFERENCE_FASTA \
    --window 35 \
    --cluster 3 \
    --filter-name "FS" \
    --filter "FS > 30.0" \
    --filter-name "QD" \
    --filter "QD < 2.0" \
    --tmp-dir $TEMP_DIR
tabix -p vcf $VARIANTFILTRATION_VCF

# Het and chrXY_filter
bcftools view -i 'GT="het"' -t ^chrX,chrY $VARIANTFILTRATION_VCF -Oz -o $AUTOSOMAL_HET_VCF
tabix -p vcf $AUTOSOMAL_HET_VCF

# Delete intermediate files
rm -f "${SPLITNCIGARREADS_BAM}"{,.bai} "$FLAG_FILTERED_BAM"{,.bai} "$ADDREADGROUPS_BAM"{,.bai}

