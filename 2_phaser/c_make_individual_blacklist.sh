#!/bin/bash
#SBATCH --job-name=make_variant_blacklist 
#SBATCH --output=logs_slurm/%j.txt   
#SBATCH --mail-type=FAIL,END  
#SBATCH --time=23:59:59   
#SBATCH --mem=1G

#' Make individual blacklists for phASER 
#'
#' This module executes WGS quality control measures to create
#' individualized variant blacklists for input to phASER.
#' While the parameters can be changed, this executes our
#' recommendations, keeping variants with the following criteria,
#' blacklisting the rest:
#' - Variant was het pre-pop. phasing (not imputed as het)
#' - 10 < number of DNA reads that align to variant <= 80
#' - 40 < avg. DNA Mapping Quality of reads aligned to variant <= 60 (not relevant to NEPTUNE)
#' - 0.2 < fraction of aligned DNA with alternative allele at variant < 0.8 
#' - Variant is a site of DNA-RNA discordance (see `find_dna_rna_discord_vars.sh`)
#'
#' Output:
#' - ${SUBJECT_ID}.final_blacklist.bed:
#'      Individual-level blacklist for phASER
#'
#' Usage:
#'   ./c_make_individual_blacklist.sh 40 60 10 80
#'
#' Author:
#'   Eric Sakkas
#'   erotokritos.sakkas@childrens.harvard.edu
#'
#' Date:
#'   2025-09-10

if [ $# -lt 4 ]; then 
echo "Please enter the correct number of commandline arguments"
echo "sample input:"
echo './c_make_individual_blacklist.sh 40 60 10 80'

echo '$1 = minimum mapping quality (exclusive)'
echo '$2 = maximum mapping quality (inclusive)'
echo '$3 = minimum WGS depth (exclusive)'
echo '$4 = maximum WGS depth (inclusive)'
exit
fi

min_mq=$1
max_mq=$2
min_dp=$3
max_dp=$4

## 1)  Identify VCF files
SUBJECT_ID=HG00171
DATA_DIR=../data/${SUBJECT_ID}

SAMPLE_VCF_UNPHASED=${DATA_DIR}/${SUBJECT_ID}.het.unphased.vcf.gz
SAMPLE_VCF_PHASED=${DATA_DIR}/${SUBJECT_ID}.het.phased.vcf.gz

## 2) Identify revived variants:
echo 'Identifying revived vars...'
echo -e "\tDATE, TIME : $(date '+%Y-%m-%d, %H:%M:%S')"
OUTDIR=../data/output/phaser/${SUBJECT_ID}
OUT_REVIVED_BED=${OUTDIR}/${SUBJECT_ID}.revived_hets.bed

## NEED TO MAKE THESE ACTUAL HET FILES ##

bcftools isec -n~10 -c none $SAMPLE_VCF_PHASED $SAMPLE_VCF_UNPHASED | awk '{print $1"\t"$2-1"\t"$2"\tREVIVED_SNP"}' > $OUT_REVIVED_BED

## 3) Create list of IDs that fail WGS filters in VCF
echo 'Creating Variant Blacklist...'
echo -e "\tDATE, TIME : $(date '+%Y-%m-%d, %H:%M:%S')"
VARIANT_BLACKLIST=${OUTDIR}/${SUBJECT_ID}.filtered_SNPs.bed
bcftools query \
        -f "%CHROM\_%POS\n" \
        -i 'INFO/MQ<='"$min_mq"' | INFO/MQ>'"$max_mq"' | FORMAT/DP<='"$min_dp"' | FORMAT/DP>'"$max_dp"' | FORMAT/AD[:1]/(FORMAT/AD[:0]+FORMAT/AD[:1])<=0.2 | FORMAT/AD[:1]/(FORMAT/AD[:0]+FORMAT/AD[:1])>=0.8' \
        $SAMPLE_VCF_UNPHASED | awk -F _ '{print $1"\t"$2-1"\t"$2}' | bedtools merge | awk '{print $0 "\tWGS_FILTER"}'> $VARIANT_BLACKLIST

## 4) Combine created blacklist with HLA blacklist provided by phASER and discordant variants 
echo 'Combining with HLA Blacklist...'
echo -e "\tDATE, TIME : $(date '+%Y-%m-%d, %H:%M:%S')"
HLA_BLACKLIST=../resource/gencode45_hg38_hla.chr.bed
BAD_SNPS_BED=../data/output/phaser/all_samples.dna_rna_discordant_vars.bed

FINAL_BLACKLIST=${OUTDIR}/${SUBJECT_ID}.final_blacklist.bed
cat $HLA_BLACKLIST $VARIANT_BLACKLIST $OUT_REVIVED_BED $BAD_SNPS_BED | sort -k 1,1 -k2,2n | bedtools merge -c 4 -o distinct > $FINAL_BLACKLIST

echo ''	
echo 'Done! Final blacklist written to' $FINAL_BLACKLIST
