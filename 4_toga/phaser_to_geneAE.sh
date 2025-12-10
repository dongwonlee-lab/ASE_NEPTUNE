#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: phaser_to_geneAE.sh --sample SAMPLE_ID --bam BAM --vcf VCF --o OUTPUT_DIR \
            --phaser PATH/TO/phaser.py --phaser_geneAE PATH/TO/phaser_gene_ae.py
Options:
  # Required
  --sample           Sample ID (e.g., SAMPLE1)
  --bam              Path to BAM
  --vcf              Path to VCF
  --o                Output directory (created if missing)
  --phaser           phASER script (phaser.py)
  --phaser_geneAE    phASER geneAE script (phaser_gene_ae.py)

  # Optional: Defaults used in this script
  --threads (10)
  --features (../resource/gencode.v.45.gene_coordinates.bed)
  --blacklist (../resource/gencode45_hg38_hla.chr.bed)
  --haplo_count_blacklist (../resource/hg38_haplo_count_blacklist.chr.bed)
  --baseq (10)  
  --mapq (255)  
  --paired_end (1)  
  --remove_dups (1)
  --write_vcf (0)  
  --pass_only (0)
  --gw_phase_method (1)  
  --gw_phase_vcf (1)
  --show_warning (1)  
  --debug (1)

  # Optional: Same defaults as in phASER
  --python_string (python3)
  --haplo_count_bam_exclude ()
  --cc_threshold (0.01)
  --isize (0)
  --as_q_cutoff (0.05)
  --include_indels (0)
  --output_read_ids (0)
  --unphased_vars (1)
  --chr_prefix ()
  --gw_af_field ('AF')
  --gw_phase_vcf_min_confidence (0.90)
  --max_block_size (15)
  --temp_dir ()
  --max_items_per_thread (100000)
  --process_slow (0)
  --chr ()
  --unique_ids (0)
  --id_separator ('_')
  --output_network ()
  --help

Example:
  ./phaser_to_geneAE.sh \
    --sample SAMPLE1 \
    --bam /path/SAMPLE1.bam \
    --vcf /path/SAMPLE1.vcf \
    --o /path/to/out \
    --phaser /path/to/phaser.py \
    --phaser_geneAE /path/to/phaser_gene_ae.py
USAGE
}

# -------------------- Defaults --------------------
THREADS=10
FEATURES="../resource/gencode.v.45.gene_coordinates.bed"
BLACKLIST="../resource/gencode45_hg38_hla.chr.bed"
HAPLO_COUNT_BLACKLIST="../resource/hg38_haplo_count_blacklist.chr.bed"
BASEQ=10
MAPQ=255
PAIRED_END=1
WRITE_VCF=0
PASS_ONLY=0
GW_PHASE_METHOD=1
GW_PHASE_VCF=1
SHOW_WARNING=1
DEBUG=1

PYTHON_STRING=python3
HAPLO_COUNT_BAM_EXCLUDE=''
CC_THRESHOLD=0.01
ISIZE=0
AS_Q_CUTOFF=0.05
INCLUDE_INDELS=0
OUTPUT_READ_IDS=0
REMOVE_DUPS=1
UNPHASED_VARS=1
CHR_PREFIX=''
GW_AF_FIELD='AF'
GW_PHASE_VCF_MIN_CONFIDENCE=0.90
MAX_BLOCK_SIZE=15
TEMP_DIR=''
MAX_ITEMS_PER_THREAD=100000
PROCESS_SLOW=0
CHR=''
UNIQUE_IDS=0
ID_SEPARATOR='_'
OUTPUT_NETWORK=''

# -------------------- Parse args (GNU getopt required) --------------------
PARSED_ARGS=$(
  getopt -o h --long \
    sample:,bam:,vcf:,o:,phaser:,phaser_geneAE:,threads:,features:,blacklist:,haplo_count_blacklist:,\
baseq:,mapq:,paired_end:,write_vcf:,pass_only:,gw_phase_method:,gw_phase_vcf:,show_warning:,debug:,\
python_string:,haplo_count_bam_exclude:,cc_threshold:,isize:,as_q_cutoff:,include_indels:,\
output_read_ids:,remove_dups:,unphased_vars:,chr_prefix:,gw_af_field:,gw_phase_vcf_min_confidence:,\
max_block_size:,temp_dir:,max_items_per_thread:,process_slow:,chr:,unique_ids:,id_separator:,\
output_network:,help \
    -n 20250909_phaser_to_geneAE.sh -- "$@"
) || { usage; exit 2; }
eval set -- "$PARSED_ARGS"

SAMPLE_ID=""
BAM=""
VCF=""
OUTPUT_DIR=""
PHASER=""
PHASER_GENE_AE=""

while true; do
  case "$1" in
    --sample)                        SAMPLE_ID="$2"; shift 2;;
    --bam)                           BAM="$2"; shift 2;;
    --vcf)                           VCF="$2"; shift 2;;
    --o)                             OUTPUT_DIR="$2"; shift 2;;
    --phaser)                        PHASER="$2"; shift 2;;
    --phaser_geneAE)                 PHASER_GENE_AE="$2"; shift 2;;
    --threads)                       THREADS="$2"; shift 2;;
    --features)                      FEATURES="$2"; shift 2;;
    --blacklist)                     BLACKLIST="$2"; shift 2;;
    --haplo_count_blacklist)         HAPLO_COUNT_BLACKLIST="$2"; shift 2;;
    --baseq)                         BASEQ="$2"; shift 2;;
    --mapq)                          MAPQ="$2"; shift 2;;
    --paired_end)                    PAIRED_END="$2"; shift 2;;
    --write_vcf)                     WRITE_VCF="$2"; shift 2;;
    --pass_only)                     PASS_ONLY="$2"; shift 2;;
    --gw_phase_method)               GW_PHASE_METHOD="$2"; shift 2;;
    --gw_phase_vcf)                  GW_PHASE_VCF="$2"; shift 2;;
    --show_warning)                  SHOW_WARNING="$2"; shift 2;;
    --debug)                         DEBUG="$2"; shift 2;;

    --python_string)                 PYTHON_STRING="$2"; shift 2;;
    --haplo_count_bam_exclude)       HAPLO_COUNT_BAM_EXCLUDE="$2"; shift 2;;
    --cc_threshold)                  CC_THRESHOLD="$2"; shift 2;;
    --isize)                         ISIZE="$2"; shift 2;;
    --as_q_cutoff)                   AS_Q_CUTOFF="$2"; shift 2;;
    --include_indels)                INCLUDE_INDELS="$2"; shift 2;;
    --output_read_ids)               OUTPUT_READ_IDS="$2"; shift 2;;
    --remove_dups)                   REMOVE_DUPS="$2"; shift 2;;
    --unphased_vars)                 UNPHASED_VARS="$2"; shift 2;;
    --chr_prefix)                    CHR_PREFIX="$2"; shift 2;;
    --gw_af_field)                   GW_AF_FIELD="$2"; shift 2;;
    --gw_phase_vcf_min_confidence)   GW_PHASE_VCF_MIN_CONFIDENCE="$2"; shift 2;;
    --max_block_size)                MAX_BLOCK_SIZE="$2"; shift 2;;
    --temp_dir)                      TEMP_DIR="$2"; shift 2;;
    --max_items_per_thread)          MAX_ITEMS_PER_THREAD="$2"; shift 2;;
    --process_slow)                  PROCESS_SLOW="$2"; shift 2;;
    --chr)                           CHR="$2"; shift 2;;
    --unique_ids)                    UNIQUE_IDS="$2"; shift 2;;
    --id_separator)                  ID_SEPARATOR="$2"; shift 2;;
    --output_network)                OUTPUT_NETWORK="$2"; shift 2;;
    -h|--help)                       usage; exit 0;;
    --) shift; break;;
    *) echo "Unknown option: $1"; usage; exit 2;;
  esac
done

# -------------------- Validate required arguments --------------------
: "${SAMPLE_ID:?Missing --sample}"
: "${BAM:?Missing --bam}"
: "${VCF:?Missing --vcf}"
: "${OUTPUT_DIR:?Missing --o}"
: "${PHASER:?Missing --phaser}"
: "${PHASER_GENE_AE:?Missing --phaser_geneAE}"

# Existence checks (nice to fail early)
[[ -f "$BAM" ]] || { echo "BAM not found: $BAM" >&2; exit 1; }
[[ -f "$VCF" ]] || { echo "VCF not found: $VCF" >&2; exit 1; }
[[ -f "$PHASER" ]] || { echo "phASER script not found: $PHASER" >&2; exit 1; }
[[ -f "$PHASER_GENE_AE" ]] || { echo "phASER geneAE script not found: $PHASER_GENE_AE" >&2; exit 1; }

# -------------------- Setup --------------------
OUTDIR="${OUTPUT_DIR%/}/${SAMPLE_ID}"
mkdir -p "$OUTDIR"
OUTPREFIX="${OUTDIR}/${SAMPLE_ID}"

# -------------------- Run phASER --------------------
echo "[phASER] Running on sample ${SAMPLE_ID}"
python3 "$PHASER" \
  --sample "$SAMPLE_ID" \
  --bam "$BAM" \
  --vcf "$VCF" \
  --o "$OUTPREFIX" \
  --threads "$THREADS" \
  --blacklist "$BLACKLIST" \
  --haplo_count_blacklist "$HAPLO_COUNT_BLACKLIST" \
  --baseq "$BASEQ" \
  --mapq "$MAPQ" \
  --paired_end "$PAIRED_END" \
  --write_vcf "$WRITE_VCF" \
  --pass_only "$PASS_ONLY" \
  --gw_phase_method "$GW_PHASE_METHOD" \
  --gw_phase_vcf "$GW_PHASE_VCF" \
  --show_warning "$SHOW_WARNING" \
  --debug "$DEBUG" \
  --python_string "$PYTHON_STRING" \
  --haplo_count_bam_exclude "$HAPLO_COUNT_BAM_EXCLUDE" \
  --cc_threshold "$CC_THRESHOLD" \
  --isize "$ISIZE" \
  --as_q_cutoff "$AS_Q_CUTOFF" \
  --include_indels "$INCLUDE_INDELS" \
  --output_read_ids "$OUTPUT_READ_IDS" \
  --remove_dups "$REMOVE_DUPS" \
  --unphased_vars "$UNPHASED_VARS" \
  --chr_prefix "$CHR_PREFIX" \
  --gw_af_field "$GW_AF_FIELD" \
  --gw_phase_vcf_min_confidence "$GW_PHASE_VCF_MIN_CONFIDENCE" \
  --max_block_size "$MAX_BLOCK_SIZE" \
  --temp_dir "$TEMP_DIR" \
  --max_items_per_thread "$MAX_ITEMS_PER_THREAD" \
  --process_slow "$PROCESS_SLOW" \
  --chr "$CHR" \
  --unique_ids "$UNIQUE_IDS" \
  --id_separator "$ID_SEPARATOR" \
  --output_network "$OUTPUT_NETWORK"

# -------------------- Run phASER geneAE --------------------
HAPLOTYPIC_COUNTS="${OUTDIR}/${SAMPLE_ID}.haplotypic_counts.txt"
OUTPUT="${OUTDIR}/${SAMPLE_ID}.gene_ae.txt"

echo "[geneAE] Using haplotypic counts ${HAPLOTYPIC_COUNTS}"
[[ -f "$HAPLOTYPIC_COUNTS" ]] || { echo "Missing haplotypic counts: $HAPLOTYPIC_COUNTS" >&2; exit 1; }

python3 "$PHASER_GENE_AE" \
  --haplotypic_counts "$HAPLOTYPIC_COUNTS" \
  --features "$FEATURES" \
  --o "$OUTPUT"

echo "Done. Outputs in: $OUTDIR"
