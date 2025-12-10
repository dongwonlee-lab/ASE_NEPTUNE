import pandas as pd
import numpy as np
import argparse
import os
import datetime

#' Identify overlapping genes among sample's COV genes
#' 
#' Gathers COV genes from Gene AE output using an RNA coverage 
#' coverage threshold, and determines which COV genes overlap
#' within the individual. Outputs this information as a txt 
#' file and a summarized BED file.
#'
#' Outputs: 
#' - {sample_id}.cov.overlapping_info.txt
#'      List of COV genes from phASER geneAE with overlapping information 
#'      and normalized RNA expression added. More info in README
#' - {sample_id}.cov.overlapping_region.bed
#'      BED file containing overlapping genomic intervals. Used
#'      to filter out reads aligned to shared variants between overlapping genes.
#'
#' Usage:
#'   python3 a_create_overlapping_info.py \
#'     -s sample_id \
#'     -i /path to/phaser geneAE file \
#'     -o /path to/output directory \
#'     -e /path to/normalized gene expression file \
#'     -g /path to/gene coordinate file with gene type and gene symbol \ # Optional. Default as `../resource/gencode.v.45.gene_coordinates.type.gene_symbol.bed`
#'     -c coverage threshold # Optional. Default as 20
#'
#' Author:
#'   Junmo Sung
#'   jsung5@bwh.harvard.edu
#'
#' Date:
#'   2025-09-25

# Argument parsing
parser = argparse.ArgumentParser(description="Add overlapping info to phASER gene AE result")
parser.add_argument('-s', '--sample', required=True, help='Sample name (column in expression file)')
parser.add_argument('-i', '--input', required=True, help='Input phASER gene AE file')
parser.add_argument('-o', '--outdir', required=True, help='Output directory')
parser.add_argument('-e', '--exp', required=True, help='Normalized expression file')
parser.add_argument('-g', '--coordinate', default='../resource/gencode.v.45.gene_coordinates.type.gene_symbol.bed', help='Gene coordinate BED file')
parser.add_argument('-c', '--cov_cutoff', type=int, default=20, help='Coverage cutoff')
args = parser.parse_args()

def flatten_region_list(region_string):
    """ Flatten overlapping region for easy data parsing """
    cleaned = region_string.replace("),", ":").replace("(", "").replace(")", "").replace(" ", "")
    return cleaned

def natural_sort_key(chr_label):
    """ Extract numeric value from chromosome label for natural sorting. """
    return int(chr_label.replace('chr', '')) if chr_label.replace('chr', '').isdigit() else chr_label

print(f"Creating overlap info for sample: {args.sample}")
print(("\tDATE, TIME : %s" % (datetime.datetime.now().strftime('%Y-%m-%d, %H:%M:%S'))))

print(f"Creating output directory {args.outdir}/{args.sample}")
out_cov_dir = os.path.join(args.outdir, args.sample)
os.makedirs(args.outdir, exist_ok=True)        # main output directory
os.makedirs(out_cov_dir, exist_ok=True)        # nested cov.overlapping_info directory

phaser_geneAE = pd.read_csv(args.input, sep='\t')
cov_genes = phaser_geneAE[phaser_geneAE['totalCount'] >= args.cov_cutoff].copy()

norm_rna_exp = pd.read_csv(args.exp, sep='\t')
gene_coordinate = pd.read_csv(args.coordinate, sep='\t', header=None, names=['chr', 'start', 'end', 'gene', 'type', 'gene_symbol'])

print(f"Adding overlapping info to COV genes...")
print(("\tDATE, TIME : %s" % (datetime.datetime.now().strftime('%Y-%m-%d, %H:%M:%S'))))

# Pre-process: Split gene_id into name/gene_symbol
if 'gene_id' in norm_rna_exp.columns:
    split_id = norm_rna_exp['gene_id'].str.split('_', n=1, expand=True)
    norm_rna_exp['name'] = split_id[0]
    norm_rna_exp['gene_symbol'] = split_id[1] if split_id.shape[1] > 1 else np.nan

# Merge gene type and gene_symbol into cov_genes
cov_genes = pd.merge(
    cov_genes,
    gene_coordinate[['gene', 'type', 'gene_symbol']],
    how='left',
    left_on='name',
    right_on='gene'
).drop(columns=['bam', 'gene', 'gw_phased'], errors='ignore')

# Add overlapping info
overlapped_flags = []
overlapped_name = []
overlapped_gene_symbol = []
overlapped_gene_type = []
overlapped_length = []
overlapped_region = []

for idx, row in cov_genes.iterrows():
    s_df = cov_genes[cov_genes.index != idx]
    overlaps = s_df[
        (s_df['contig'] == row['contig']) &
        (s_df['start'] <= row['stop']) &
        (s_df['stop'] >= row['start'])
    ]
    if not overlaps.empty:
        overlapped_flags.append(True)
        overlapped_name.append(overlaps['name'].tolist())
        overlapped_gene_symbol.append(overlaps['gene_symbol'].tolist())
        overlapped_gene_type.append(overlaps['type'].tolist())
        overlapped_length.append([
            min(row['stop'], overlap['stop']) - max(row['start'], overlap['start'])
            for _, overlap in overlaps.iterrows()
        ])
        overlapped_region.append([
            (max(row['start'], overlap['start']), min(row['stop'], overlap['stop']))
            for _, overlap in overlaps.iterrows()
        ])
    else:
        overlapped_flags.append(False)
        overlapped_name.append(np.nan)
        overlapped_gene_symbol.append(np.nan)
        overlapped_gene_type.append(np.nan)
        overlapped_length.append(np.nan)
        overlapped_region.append(np.nan)

cov_genes['overlapped_flags'] = overlapped_flags
cov_genes['overlapped_gene_name'] = overlapped_name
cov_genes['overlapped_gene_symbol'] = overlapped_gene_symbol
cov_genes['overlapped_gene_type'] = overlapped_gene_type
cov_genes['overlapped_length'] = overlapped_length
cov_genes['overlapped_region'] = overlapped_region

# Flatten lists for output
for col in ['overlapped_gene_name', 'overlapped_gene_symbol', 'overlapped_gene_type', 'overlapped_length', 'overlapped_region']:
    cov_genes[col] = cov_genes[col].apply(lambda x: ', '.join(map(str, x)) if isinstance(x, list) else x)

# Add normalized expression
if args.sample in norm_rna_exp.columns:
    rsem_s_df = norm_rna_exp[['name', args.sample]].rename(columns={args.sample: 'norm_rna_exp'})
    cov_genes = pd.merge(cov_genes, rsem_s_df, on='name', how='left')

print(f"Saving COV genes with overlapping info to {out_cov_dir}/{args.sample}.cov.overlapping_info.txt")
output_file = os.path.join(out_cov_dir, f"{args.sample}.cov.overlapping_info.txt")
cov_genes.to_csv(output_file, sep='\t', index=False)

print("Create bed file that contains overlapping region")
print(("\tDATE, TIME : %s" % (datetime.datetime.now().strftime('%Y-%m-%d, %H:%M:%S'))))

# Select genes with overlapping region
overlapped_cov_genes = cov_genes[cov_genes['overlapped_flags']].copy()

# Flatten overlapping region
overlapped_cov_genes['overlapped_region'] = overlapped_cov_genes['overlapped_region'].apply(flatten_region_list)

# Explode by overlapping region
overlapped_cov_genes['overlapped_region'] = overlapped_cov_genes['overlapped_region'].str.split(':')
overlapped_cov_genes = overlapped_cov_genes.explode('overlapped_region')

# Split 'overlapped_region' into 'bed_start' and 'bed_end'
overlapped_cov_genes[['bed_start', 'bed_end']] = overlapped_cov_genes['overlapped_region'].str.split(',', expand=True)

# Convert start and end to integers
overlapped_cov_genes['bed_start'] = overlapped_cov_genes['bed_start'].astype(int) 
overlapped_cov_genes['bed_end'] = overlapped_cov_genes['bed_end'].astype(int)

# Create bed file
bed = overlapped_cov_genes[['contig', 'bed_start', 'bed_end']].copy()
bed.columns = ['chr', 'start', 'end']

# Ensure 'chr' column is sorted naturally and sort by numeric chromosome order first, then by start position
bed['chr_num'] = bed['chr'].apply(natural_sort_key)
bed_sorted = bed.sort_values(by=['chr_num', 'start']).drop(columns=['chr_num']).reset_index(drop=True)

# Remove duplicates
bed_flattened = bed_sorted.drop_duplicates(subset=['chr', 'start', 'end']).reset_index(drop=True)

print(f"Saving the bed file to {out_cov_dir}/{args.sample}.cov.overlapping_region.bed")
bed_file = os.path.join(out_cov_dir, f"{args.sample}.cov.overlapping_region.bed")
bed_flattened.to_csv(bed_file, sep='\t', index=False, header=False)
print("Done!")