### Load library ###
import pandas as pd
import numpy as np
import sys
from statsmodels.stats.multitest import multipletests
import os
from scipy.stats import binomtest
import argparse

#' Run TOGA Algorithm
#' 
#' Runs TOGA on allelic and gene ASE data 
#' appropriately assigning ASE counts in overlapping
#' gene regions and outputting 
#' file and a summarized BED file.
#'
#' Inputs:
#'      See parameters detailed in README
#'
#' Outputs: 
#' - HG00171_TOGA_light.txt
#'      List of COV genes from phASER geneAE with ASE 
#'      based on corrected overlapping attributes, see README
#' - HG00171_TOGA_heavy.txt
#'      light.txt file with additional columns, see README
#'
#' Usage:
#'   python3 d_TOGA.py \
#'      * see options in README
#'
#' Author:
#'   Junmo Sung
#'   jsung5@bwh.harvard.edu
#'
#' Date:
#'   2025-09-26

### WELCOME ###
print('######   ######   ########   ###### ')
print('  ##    ##    ##  ##        ##    ##')
print('  ##    ##    ##  ##        ##    ##')
print('  ##    ##    ##  ##  ####  ##    ##')
print('  ##    ##    ##  ##    ##  ########')
print('  ##    ##    ##  ##    ##  ##    ##')
print('  ##     ######   ########  ##    ##')

### Set arguments ###
parser = argparse.ArgumentParser(description="Run TOGA analysis pipeline")

# Required
parser.add_argument('-s', '--sample', required=True, help="Sample name")
parser.add_argument('-i', '--cov_overlapping_info_path', required=True, help="COV genes with overlapping info")
parser.add_argument('-a', '--var_cnt_path', required=True, help="Allelic count file path")
parser.add_argument('-p', '--phase_gene_ae_wo_overlapping_region_path', required=True, help="phASER gene AE results w/o overlapping region")
parser.add_argument('-o', '--output_dir', required=True, help="Output directory")

# Optional
parser.add_argument('-e', '--exon_coordinate_path', default='../resource/gencode.v45.exon.basic.partitioned.bed', help="Exonic coordinates file path")
parser.add_argument('-c', '--cov_cutoff', type=int, default=20, help="Coverage cutoff for non-overlapping variants")
parser.add_argument('-ex', '--exp_cutoff', type=float, default=0, help="RNA expression cutoff")
parser.add_argument('-we', '--w_exon', type=float, default=0.99, help="Weight for exonic variant")
parser.add_argument('-wi', '--w_intron', type=float, default=0.01, help="Weight for intronic variant")
parser.add_argument('-cp', '--cpov_cutoff', type=float, default=0.5, help="CPOV cutoff")
parser.add_argument('-ir', '--irpg_cutoff', type=float, default=0.1, help="IRPG cutoff")
parser.add_argument('-ac', '--adj_cov_cutoff', type=int, default=20, help="Adjusted total coverage cutoff")
parser.add_argument('-r', '--haplotypic_ratio_cutoff', type=float, default=1.5, help="Haplotypic ratio cutoff for ASE")
parser.add_argument('-f', '--binom_fdr_cutoff', type=float, default=0.05, help="Binomial FDR cutoff")
parser.add_argument('-w', '--wrt_heavy', type=int, default=1, help="Save heavy file (1=True, 0=False)")
parser.add_argument('-v', '--var_hap_cnt_path', default="", help="Allelic counts by haplotype file")
parser.add_argument('-vw', '--woOR_var_hap_cnt_path', default="", help="Allelic counts by haplotype w/o overlapping regions file")

args = parser.parse_args()

### Print arguments ###
print('\nArguments:')
for arg in vars(args):
    print(f'  {arg}: {getattr(args, arg)}')

### Creating output direcotry ###
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
    print("No output directory found. Creating one...")

### Load data ### 
cov_overlapping_info=pd.read_csv(args.cov_overlapping_info_path ,sep='\t') # COV genes with overlapping info
phase_gene_ae_wo_overlapping_region = pd.read_csv(args.phase_gene_ae_wo_overlapping_region_path, sep='\t')
var_cnt=pd.read_csv(args.var_cnt_path, sep='\t') # Allelic count (from the one of phASER outputs)
exon_coordinate=pd.read_csv(args.exon_coordinate_path, sep='\t', names=['contig','start','end','gene']) # Exon coordinates

### Functions to use ###
# Cacalculate binomial p-value using aCount and bCount
def get_binom_p(row):
    a = row['aCount']
    b = row['bCount']
    if pd.isnull(a) or pd.isnull(b):
        return np.nan
    n = a + b
    if n < 1:
        return np.nan
    return binomtest(a, n=n, p=0.5).pvalue
# Sort variants by position
def sort_variants(variant_str):
    if isinstance(variant_str, str): 
        variants = variant_str.split(',')
        sorted_variants = sorted(variants, key=lambda x: int(x.split('_')[1])) 
        return ','.join(sorted_variants)
    return variant_str  

# Find overlapping variant
def find_variant_intersection(variants_list):
    variant_sets = [set(variants.split(",")) for variants in variants_list]
    return ",".join(set.intersection(*variant_sets)) if variant_sets else pd.NA

# Exonic status for the union of overlapping varaints
def get_exon_intron_status(row, exon_df):
    gene = row['name']
    df = exon_df[exon_df['gene'] == gene]  
    exon_list = []
    
    for var in row['union_overlapping_variants'].split(','):
        parts = var.split('_')
        contig, pos = parts[0], int(parts[1])  

        # Default to "null" unless found in exon/intron
        variant_status = "null"

        variant_set = set(row['variants'].split(',')) 

        if var in variant_set:
            # Default to intron if the variant is listed in 'variants'
            variant_status = "intron"
            for _, exon in df.iterrows():
                if exon['contig'] == contig and exon['start'] <= pos <= exon['end']:
                    variant_status = "exon"
                    break  # Stop checking once it's confirmed in an exon

        exon_list.append(variant_status)

    return ','.join(exon_list)

# Allelic count for the union of overlapping variants
def get_snp_count(row, s_allelic_count):
    snp_count_list = []

    for var in row['union_overlapping_variants'].split(','):
        variant_set = set(row['variants'].split(',')) 

        if var not in variant_set:
            snp_count_list.append("0")  # Append "0" if variant not in variant_set
        else:
            snp_count = s_allelic_count.loc[s_allelic_count['variantID'] == var, 'totalCount']

            # Check if the variant exists in the allele count data
            if snp_count.empty:
                print(f"Warning: Variant {var} not found in allele count dataset. Row being analyzed: {row.name}")
                snp_count_list.append("Missing")
            else:
                snp_count_list.append(str(snp_count.iloc[0]))  # Convert to string for consistent output

    # If all counts are missing, return NA
    return pd.NA if all(x == "Missing" for x in snp_count_list) else ','.join(snp_count_list)

# Calculate rna expression x weight
def exp_weight(row):
    exp = row['norm_rna_exp']
    var_exon_ls = row['union_variants_overlapping_exon'].split(',')
    var_exon_weight_ls = [
        args.w_exon if exon == "exon" else 
        args.w_intron if exon == "intron" else 
        w_null 
        for exon in var_exon_ls
    ]
    # Convert to NumPy array for element-wise multiplication
    exp_weighted = np.array(var_exon_weight_ls) * exp
    return ",".join(map(str, exp_weighted))

# Calculate adj_snp_count
def multiply_snp_counts(row):
    snp_counts = list(map(float, row['union_overlapping_variants_snp_count'].split(','))) if isinstance(row['union_overlapping_variants_snp_count'], str) else [row['union_overlapping_variants_snp_count']]
    weights = list(map(float, row['union_exp_weight_proportion'].split(','))) if isinstance(row['union_exp_weight_proportion'], str) else [row['union_exp_weight_proportion']]
    
    # Perform element-wise multiplication
    adj_var_count = [snp * weight for snp, weight in zip(snp_counts, weights)]
    return ",".join(map(str, adj_var_count))

# Sort variants by position. Used for var_hap_cnt
def sort_variant_for_hap_counts(row):
    variants = row['variants'].split(',')
    sorted_variants = sorted(variants, key=lambda v: int(v.split('_')[1]))  # Sort by position
    sorted_indices = [variants.index(v) for v in sorted_variants]
    
    # Reorder other columns accordingly
    for col in columns_to_modify:
        values = row[col].split(',')
        row[col] = ','.join([values[i] for i in sorted_indices])
    
    row['variants'] = ','.join(sorted_variants)
    return row

# Get mathing rows
def get_matching_rows(df, df_cols, row, row_cols):
    mask = True
    for df_col, row_col in zip(df_cols, row_cols):
        mask &= (df[df_col] == row[row_col])
    return df[mask]

### Pre-process ###
print('Pre-processing... ')
# Do binomial test on w/wo phaser result
cov_overlapping_info['binom_p'] = cov_overlapping_info.apply(get_binom_p, axis=1)
phase_gene_ae_wo_overlapping_region['binom_p'] = phase_gene_ae_wo_overlapping_region.apply(get_binom_p, axis=1)

# Starting from the cov genes with overlapping info
df_summary = cov_overlapping_info.copy()

# Add non overlapping phaser info to the df_summary
merge_columns = ['contig', 'start', 'stop', 'name']
df_summary = pd.merge(df_summary, phase_gene_ae_wo_overlapping_region[merge_columns + ['aCount', 'bCount', 'totalCount', 'log2_aFC', 'n_variants', 'variants', 'binom_p']], on=merge_columns, suffixes=('', '_blOR'))

# totalCount - totalCount_blOR
df_summary['totalCount_diff'] = df_summary['totalCount'] - df_summary['totalCount_blOR']

# Sort variants by postion
df_summary['variants'] = df_summary['variants'].apply(sort_variants)
df_summary['variants_blOR'] = df_summary['variants_blOR'].apply(sort_variants)
print('Done.')

### STEP1: Greater than RNA expression cutoff ### 
print('Applying RNA expression cutoff... ')
df_summary['GT.EXP_CUTOFF'] = False
df_summary.loc[df_summary['norm_rna_exp'] > args.exp_cutoff, 'GT.EXP_CUTOFF'] = True # NA value in norm_rna_exp will be regareded as GT.EXP_CUTOFF=FALSE
print('Done.')

### STEP2: Genes sharing genomic position ###
print('Finding genes sharing genomic position... ')
df_summary_gp=df_summary[(df_summary['GT.EXP_CUTOFF'] == True)& (df_summary['overlapped_flags'] == True)]
print('Done.')

### STEP3: Genes sharing variants ###
print('Finding genes sharing variants... ')
## First, explode overlapping genes by gene pair 
# Columns to explode
df_summary_tmp = df_summary_gp.copy()
cols_to_explode = ['overlapped_gene_name', 'overlapped_gene_symbol', 'overlapped_gene_type', 'overlapped_length']

for col in cols_to_explode:
    df_summary_tmp[col] = df_summary_tmp[col].str.split(', ')

# Parse overlapped_region
df_summary_tmp['overlapped_region'] = df_summary_tmp['overlapped_region'].str.replace(r'\)\s*,\s*\(', ':', regex=True)
df_summary_tmp['overlapped_region'] = df_summary_tmp['overlapped_region'].str.replace(r'[\(\)]', '', regex=True)
df_summary_tmp['overlapped_region'] = df_summary_tmp['overlapped_region'].str.split(':')  # Convert to list

# Explode by overlapped info(cols_to_explode + 'overlapped_region')
exploded_df_summary_tmp = df_summary_tmp.explode(cols_to_explode + ['overlapped_region'], ignore_index=True)

## Second, select genes with paired overlapping genes because some of paired genes were resolved by rna expression already in the dataframe.
paired_df_summary_tmp = [] 
for idx, row in exploded_df_summary_tmp.iterrows():
    z = exploded_df_summary_tmp.drop(index=idx) 
    gene = row['name']
    overlapped_gene = row['overlapped_gene_name']
    # Check if both 'gene' and 'overlapped_gene' exist in the remaining dataframe
    if ((z['overlapped_gene_name'] == gene) & (z['name'] == overlapped_gene)).any():
        paired_df_summary_tmp.append(row)  
paired_df_summary_tmp = pd.DataFrame(paired_df_summary_tmp)
paired_df_summary_tmp.reset_index(inplace=True)
paired_df_summary_tmp.drop(columns=['index'], inplace=True)

## Third,assign the gene cluster number grouped by overlapping region
paired_df_summary_tmp['gene_group_by_overlapping_region'] = np.nan
gene_group_by_overlapping_region = 0
region_dict = {}

for idx, row in paired_df_summary_tmp.iterrows():
    if not np.isnan(row['gene_group_by_overlapping_region']):
        continue  # Skip if already assigned
    overlapped_region = row['overlapped_region']
    
    # Check if this gene pairing was already assigned
    if overlapped_region in region_dict:
        assigned_group = region_dict.get(overlapped_region)
        paired_df_summary_tmp.at[idx, 'gene_group_by_overlapping_region'] = assigned_group
    else:
        # Find the matching row
        paired_idx = paired_df_summary_tmp[
            paired_df_summary_tmp['overlapped_region'] == overlapped_region
        ].index

        if len(paired_idx) > 0:
            paired_df_summary_tmp.at[idx, 'gene_group_by_overlapping_region'] = gene_group_by_overlapping_region
            paired_df_summary_tmp.at[paired_idx[0], 'gene_group_by_overlapping_region'] = gene_group_by_overlapping_region
            region_dict[overlapped_region] = gene_group_by_overlapping_region
            gene_group_by_overlapping_region += 1

## Fourth, only select genes that have shared variants within overlapping region
no_variants_group=[]
for group in paired_df_summary_tmp['gene_group_by_overlapping_region'].unique():
    df = paired_df_summary_tmp[paired_df_summary_tmp['gene_group_by_overlapping_region'] == group]
    
    # Convert 'variants' column into sets of variants
    variant_sets = df['variants'].dropna().apply(lambda x: set(x.split(',')))
    
    # Find the intersection of all variant sets in this group
    common_variants = set.intersection(*variant_sets) if not variant_sets.empty else set()
    
    # If there are no common variants across all rows, add the group to the list
    if not common_variants:
        no_variants_group.append(int(group))

shared_var_paired_df_summary_tmp = paired_df_summary_tmp[
    ~paired_df_summary_tmp['gene_group_by_overlapping_region'].isin(no_variants_group)
]

## Fifth, find overlapping variant and sort them by position
grouped_variants = shared_var_paired_df_summary_tmp.groupby("gene_group_by_overlapping_region")["variants"].apply(find_variant_intersection).reset_index()
shared_var_paired_df_summary_tmp = shared_var_paired_df_summary_tmp.merge(grouped_variants, on="gene_group_by_overlapping_region", suffixes=("", "_overlapping"))

# Sort overlapping variants by position
shared_var_paired_df_summary_tmp['variants_overlapping'] = shared_var_paired_df_summary_tmp['variants_overlapping'].apply(sort_variants)

## Sixth, assign the gene cluster number grouped by overlapping variant
shared_var_paired_df_summary_tmp['gene_group_by_overlapping_variant'] = np.nan
gene_group_by_overlapping_variant = 0
variant_dict = {}

for idx, row in shared_var_paired_df_summary_tmp.iterrows():
    if not np.isnan(row['gene_group_by_overlapping_variant']):
        continue  # Skip if already assigned
    overlapped_variant = row['variants_overlapping']
    
    # Check if this gene pairing was already assigned
    if overlapped_variant in variant_dict:
        assigned_group = variant_dict.get(overlapped_variant)
        shared_var_paired_df_summary_tmp.at[idx, 'gene_group_by_overlapping_variant'] = assigned_group
    else:
        # Find the matching row
        paired_idx = shared_var_paired_df_summary_tmp[
            shared_var_paired_df_summary_tmp['variants_overlapping'] == overlapped_variant
        ].index

        if len(paired_idx) > 0:
            shared_var_paired_df_summary_tmp.at[idx, 'gene_group_by_overlapping_variant'] = gene_group_by_overlapping_variant
            shared_var_paired_df_summary_tmp.at[paired_idx[0], 'gene_group_by_overlapping_variant'] = gene_group_by_overlapping_variant
            variant_dict[overlapped_variant] = gene_group_by_overlapping_variant
            gene_group_by_overlapping_variant += 1
shared_var_paired_df_summary_tmp['gene_group_by_overlapping_variant'] = shared_var_paired_df_summary_tmp['gene_group_by_overlapping_variant'].astype(int)
print('Done.')

### STEP4: Pre-process before calculating CPOV and IRPG ###
print('Preparing before calculating CPOV... ')
## First, assign gene group by the union of overlapping variants
overlapping_gene_clusters = {}
n = 0

for idx, row in shared_var_paired_df_summary_tmp.iterrows():
    gene = row['name']
    overlapping_gene = row['overlapped_gene_name']
    
    # Ensure overlapping_gene is treated as a set
    if isinstance(overlapping_gene, str):  # Check if it's a string
        overlapping_gene = {overlapping_gene}  # Convert to a set
    else:
        overlapping_gene = set() if pd.isna(overlapping_gene) else {overlapping_gene}

    # Create a new cluster if gene is not assigned to any existing cluster
    if not any(gene in cluster for cluster in overlapping_gene_clusters.values()):
        overlapping_gene_clusters[n] = {gene} | overlapping_gene
        n += 1
    else:
        for key, cluster in overlapping_gene_clusters.items():
            if gene in cluster:
                cluster.update(overlapping_gene)  # Properly update with a set

# Create a mapping of gene to cluster ID
gene_to_cluster = {gene: cluster_id for cluster_id, genes in overlapping_gene_clusters.items() for gene in genes}

# Assign cluster ID to 'df_summary' based on 'name'
shared_var_paired_df_summary_tmp['gene_group_by_union_overlapping_variants'] = shared_var_paired_df_summary_tmp['name'].map(gene_to_cluster)

## Second, get variants from the set of union
for idx in set(shared_var_paired_df_summary_tmp['gene_group_by_union_overlapping_variants']):
    df = shared_var_paired_df_summary_tmp[
        shared_var_paired_df_summary_tmp['gene_group_by_union_overlapping_variants'] == idx
    ]
    # Flatten the set into a comma-separated string
    flattened_variants = ",".join(set(','.join(df['variants_overlapping']).split(',')))
    
    # Assign the flattened string to all rows of the matching index
    shared_var_paired_df_summary_tmp.loc[
        shared_var_paired_df_summary_tmp['gene_group_by_union_overlapping_variants'] == idx, 
        'union_overlapping_variants'
    ] = flattened_variants

# Sort union_overlapping_variants
shared_var_paired_df_summary_tmp['union_overlapping_variants'] = shared_var_paired_df_summary_tmp['union_overlapping_variants'].apply(sort_variants)

## Third, get all the variants(including overlapping variants) from the union set
for idx in set(shared_var_paired_df_summary_tmp['gene_group_by_union_overlapping_variants']):
    df = shared_var_paired_df_summary_tmp[
        shared_var_paired_df_summary_tmp['gene_group_by_union_overlapping_variants'] == idx
    ]
    # Flatten the set into a comma-separated string
    flattened_variants = ",".join(set(','.join(df['variants']).split(',')))
    
    # Assign the flattened string to all rows of the matching index
    shared_var_paired_df_summary_tmp.loc[
        shared_var_paired_df_summary_tmp['gene_group_by_union_overlapping_variants'] == idx, 
        'union_variants'
    ] = flattened_variants

# Sort union_variants
shared_var_paired_df_summary_tmp['union_variants'] = shared_var_paired_df_summary_tmp['union_variants'].apply(sort_variants)

## Fourth, collapse by 'gene_group_by_union_variants_overlapping' to get the same gene name
exclude_columns = ['overlapped_gene_name', 'overlapped_gene_symbol', 'overlapped_gene_type', 'overlapped_length', 'overlapped_region', 'gene_group_by_overlapping_region', 'variants_overlapping', 'gene_group_by_overlapping_variant']

shared_var_paired_df_summary_tmp_prior_cpov_ls = []
for group in shared_var_paired_df_summary_tmp['gene_group_by_union_overlapping_variants'].unique():
    df = shared_var_paired_df_summary_tmp[
        shared_var_paired_df_summary_tmp['gene_group_by_union_overlapping_variants'] == group
    ].copy()
    
    # Aggregate overlapped gene information
    df_overlapped_info = df[['name'] + exclude_columns].drop_duplicates().groupby('name').agg({
    'overlapped_gene_name': lambda x: ':'.join(map(str, x.dropna())), 
    'overlapped_gene_symbol': lambda x: ':'.join(map(str, x.dropna())),
    'overlapped_gene_type': lambda x: ':'.join(map(str, x.dropna())),
    'overlapped_length': lambda x: ':'.join(map(str, x.dropna())),
    'overlapped_region': lambda x: ':'.join(map(str, x.dropna())),
    'gene_group_by_overlapping_region': lambda x: ':'.join(map(str, x.dropna())),
    'variants_overlapping': lambda x: ':'.join(map(str, x.dropna())),
    'gene_group_by_overlapping_variant': lambda x: ':'.join(map(str, x.dropna())),
    }).reset_index()
    
    # Drop excluded columns and remove duplicates
    df = df.drop(columns=exclude_columns).drop_duplicates()
    
    # Merge with aggregated overlapped information
    df = pd.merge(df, df_overlapped_info, on='name', how='left')
    
    # Store the processed dataframe
    shared_var_paired_df_summary_tmp_prior_cpov_ls.append(df)

# Concatenate all processed dataframes
shared_var_paired_df_summary_tmp_prior_cpov = pd.concat(shared_var_paired_df_summary_tmp_prior_cpov_ls, ignore_index=True)

## Fifth, get exonic status for union of overlapping variants
shared_var_paired_df_summary_tmp_prior_cpov['union_variants_overlapping_exon'] = shared_var_paired_df_summary_tmp_prior_cpov.apply(
    lambda row: get_exon_intron_status(row, exon_coordinate), axis=1
)

## Sixth, get allelic count union of overlapping variants
shared_var_paired_df_summary_tmp_prior_cpov['union_overlapping_variants_snp_count'] = shared_var_paired_df_summary_tmp_prior_cpov.apply(
    lambda row: get_snp_count(row, var_cnt), axis=1
)

## Seventh, calculate rna expression x weight
w_null = 0
shared_var_paired_df_summary_tmp_prior_cpov['union_exp_weight'] = shared_var_paired_df_summary_tmp_prior_cpov.apply(exp_weight, axis=1)

## Eighth, calculate exp x weight in proportion (W(vi,gj))
# Convert 'exp_weight' to lists of floats
shared_var_paired_df_summary_tmp_prior_cpov['union_exp_weight'] = shared_var_paired_df_summary_tmp_prior_cpov['union_exp_weight'].apply(lambda x: list(map(float, x.split(','))))

shared_var_paired_df_summary_tmp_cpov = pd.DataFrame()
for group in shared_var_paired_df_summary_tmp_prior_cpov['gene_group_by_union_overlapping_variants'].unique():
    k = shared_var_paired_df_summary_tmp_prior_cpov[shared_var_paired_df_summary_tmp_prior_cpov['gene_group_by_union_overlapping_variants'] == group].copy()
    
    # Find max list length in this group(max list length should be same since used union_overlapping_variants)
    max_len = max(k['union_exp_weight'].apply(len))
    
    # Compute the sum at each index across all rows
    total_exp_weight = [sum(row[i] if i < len(row) else 0 for row in k['union_exp_weight']) for i in range(max_len)]
    
    # Store the total sum in all rows of the group
    k['union_total_exp_weight'] = [total_exp_weight] * len(k)
    
    # Concatenate the results
    shared_var_paired_df_summary_tmp_cpov = pd.concat([shared_var_paired_df_summary_tmp_cpov, k], ignore_index=True)
# Calculate exp x weight in proportion (W(vi,gj))
shared_var_paired_df_summary_tmp_cpov['union_exp_weight_proportion'] = shared_var_paired_df_summary_tmp_cpov.apply(
    lambda row: [t / total if total != 0 else 0 for t, total in zip(row['union_exp_weight'], row['union_total_exp_weight'])], axis=1
)

# Revert list to comma separated string fomrat
shared_var_paired_df_summary_tmp_cpov['union_exp_weight_proportion'] = shared_var_paired_df_summary_tmp_cpov['union_exp_weight_proportion'].apply(lambda x: ','.join(map(str, x)))
shared_var_paired_df_summary_tmp_cpov['union_exp_weight'] = shared_var_paired_df_summary_tmp_cpov['union_exp_weight'].apply(lambda x: ','.join(map(str, x)))
shared_var_paired_df_summary_tmp_cpov['union_total_exp_weight'] = shared_var_paired_df_summary_tmp_cpov['union_total_exp_weight'].apply(lambda x: ','.join(map(str, x)))

## Nienth, calculate CPOV(f(gj)), adjusted snp count, and adjusted total snp count(this total count contains reads that align to more than one variant)
shared_var_paired_df_summary_tmp_cpov['union_adj_snp_count'] = shared_var_paired_df_summary_tmp_cpov.apply(multiply_snp_counts, axis=1)

# Calculate total count of adj_snp_count(f`(gj))
shared_var_paired_df_summary_tmp_cpov['union_total_adj_snp_count'] = shared_var_paired_df_summary_tmp_cpov['union_adj_snp_count'].apply(
    lambda x: sum(map(float, x.split(','))) if isinstance(x, str) else 0
)

print('Calculating CPOV...')
# Calculate CPOV(f(gj))
shared_var_paired_df_summary_tmp_cpov['CPOV'] = shared_var_paired_df_summary_tmp_cpov['union_total_adj_snp_count'] / shared_var_paired_df_summary_tmp_cpov.groupby('gene_group_by_union_overlapping_variants')['union_total_adj_snp_count'].transform('sum')
print('Done.')

### STEP5: Calculate IRPG(Incorrect Read Proportion in Gene) ###
print('Calculating IRPG... ')
shared_var_paired_df_summary_tmp_cpov['IRPG'] = ((1 - shared_var_paired_df_summary_tmp_cpov['CPOV']) *  shared_var_paired_df_summary_tmp_cpov['totalCount_diff']) / shared_var_paired_df_summary_tmp_cpov['totalCount']
print('Done.')

### STEP6: Calculate the total count adjusted by TPM(adj_COV) ###
print('Calculating adjusted total count by CPOV... ')
shared_var_paired_df_summary_tmp_cpov['overlappedCount_by_CPOV'] = shared_var_paired_df_summary_tmp_cpov['totalCount_diff']*shared_var_paired_df_summary_tmp_cpov['CPOV']
shared_var_paired_df_summary_tmp_cpov['adj_totalCount'] = shared_var_paired_df_summary_tmp_cpov['totalCount_blOR'] + shared_var_paired_df_summary_tmp_cpov['overlappedCount_by_CPOV']
print('Done.')

### STEP7: Filtering Step ###
# Merge it to original df_summary
print('Filtering step before tagging each gene... ')
columns = list(shared_var_paired_df_summary_tmp_cpov.loc[:, 'gene_group_by_union_overlapping_variants':].columns) + ['name']
to_remove = [
    'overlapped_gene_name',
    'overlapped_gene_symbol',
    'overlapped_gene_type',
    'overlapped_length',
    'overlapped_region'
]
columns = [col for col in columns if col not in to_remove]
cpov_to_merge = shared_var_paired_df_summary_tmp_cpov[columns]
df_summary = pd.merge(df_summary, cpov_to_merge, on='name', how='left')

# Genes with irpg < irpg_cutoff 
shared_var_paired_df_summary_tmp_cpov_irpg=shared_var_paired_df_summary_tmp_cpov[shared_var_paired_df_summary_tmp_cpov['IRPG'] >= args.irpg_cutoff]

# Genes with adj_cov <  adj_cov_cutoff
adj_cov_below = shared_var_paired_df_summary_tmp_cpov_irpg[
    shared_var_paired_df_summary_tmp_cpov_irpg['adj_totalCount'] < args.adj_cov_cutoff
]

# Genes with adj_cov >=  adj_cov_cutoff 
adj_cov_above_equal = shared_var_paired_df_summary_tmp_cpov_irpg[
    shared_var_paired_df_summary_tmp_cpov_irpg['adj_totalCount'] >= args.adj_cov_cutoff
]

# Genes with cpov < cpov_cutoff is regarded as resovled
adj_cov_above_equal_cpov_unresolved=adj_cov_above_equal[adj_cov_above_equal['CPOV'] >= args.cpov_cutoff]
adj_cov_above_equal_cpov_resolved_tmp=adj_cov_above_equal[adj_cov_above_equal['CPOV'] < args.cpov_cutoff]

# We regard genes as resolved if the reads from non overlapping region are equal or above cov_cutoff and use only reads from non_overlapping regions
adj_cov_above_equal_cpov_excluded = adj_cov_above_equal_cpov_resolved_tmp[adj_cov_above_equal_cpov_resolved_tmp['totalCount_blOR']<args.cov_cutoff]
adj_cov_above_equal_cpov_resolved = adj_cov_above_equal_cpov_resolved_tmp[adj_cov_above_equal_cpov_resolved_tmp['totalCount_blOR']>=args.cov_cutoff]
print('Done.')

### STEP8: Assign gene cluster(=genes that are mutually overlapping, universe is cov20) ###
print('Grouping genes that are are mutually overlapping by regions... ')
overlapped = df_summary[df_summary['overlapped_flags'] == True]

overlapping_gene_clusters = {}
n = 0

# Set of primary gene names
set_comp = set(overlapped['name'])

for idx, row in overlapped.iterrows():
    gene = row['name']
    overlapping_genes = row['overlapped_gene_name'].split(", ")  # Ensure proper handling of multiple genes

    # Filter only overlapping genes that exist in set_comp
    overlapping_genes = set(overlapping_genes).intersection(set_comp)

    if not overlapping_genes:
        continue  # Skip if no valid overlapping genes

    if not any(gene in cluster for cluster in overlapping_gene_clusters.values()):
        # Create a new cluster
        overlapping_gene_clusters[n] = set([gene]) | overlapping_genes
        n += 1
    else:
        # Append to existing cluster if gene already exists in a cluster
        for key, cluster in overlapping_gene_clusters.items():
            if gene in cluster:
                cluster.update(overlapping_genes)

# Create a mapping of gene to cluster ID
gene_to_cluster = {gene: cluster_id for cluster_id, genes in overlapping_gene_clusters.items() for gene in genes}

# Assign cluster ID to 'df_summary' based on 'name'
df_summary['gene_cluster'] = df_summary['name'].map(gene_to_cluster)
print('Done.')

### STEP9: Assign TAGs ###
# R: Resolved
# RC: Resolved by CPOV
# E: Excluded
# U: Unresolved
print('Assigning tags to each gene(R,RC,E,U)... ')
df_summary['TAG'] = np.where(
    (df_summary['GT.EXP_CUTOFF'] == False) | 
    (df_summary['name'].isin(adj_cov_below['name'])) |
    (df_summary['name'].isin(adj_cov_above_equal_cpov_excluded['name'])), 
    'E',
    np.where(
        df_summary['name'].isin(adj_cov_above_equal_cpov_unresolved['name']), 
        'U',
        np.where(
            df_summary['name'].isin(adj_cov_above_equal_cpov_resolved['name']),
            'RC',
            'R'
        )
    )
)

# Refine columns by ':' separated for mulitple overlapping genes
col1 = ['overlapped_gene_name', 'overlapped_gene_symbol', 'overlapped_gene_type', 'overlapped_length']

for col in col1:
    df_summary[col] = df_summary[col].astype(str).str.replace(', ', ':', regex=False)
df_summary['overlapped_region'] = df_summary['overlapped_region'].astype(str)
df_summary['overlapped_region'] = df_summary['overlapped_region'].str.replace(r'\), ', ':', regex=True)
df_summary['overlapped_region'] = df_summary['overlapped_region'].str.replace(r'[\(\)]', '', regex=True)
df_summary['overlapped_region'] = df_summary['overlapped_region'].str.replace(', ', ',', regex=False)
print('Done.')

### STEP10: Perform ASE call ###
# R,U: Use reads from all variants
# RC: Use reads from only non-overlapping variants
# E: Exclude from new ASE call. Will be filled with NA
print('Performing ASE call... ')
toga_columns = {
    'aCount':('aCount_blOR', 'aCount'),
    'bCount':('bCount_blOR', 'bCount'),
    'totalCount':('totalCount_blOR',   'totalCount'),
    'log2_aFC':('log2_aFC_blOR', 'log2_aFC'),
    'n_variants':('n_variants_blOR', 'n_variants'),
    'variants':('variants_blOR', 'variants'),
    'binom_p':('binom_p_blOR', 'binom_p'),
}

for col, (blor_col, reg_col) in toga_columns.items():
    toga_col = f'TOGA_{col}'
    df_summary[toga_col] = np.select(
        [
            df_summary['TAG'] == 'RC',
            df_summary['TAG'].isin(['R', 'U'])
        ],
        [
            df_summary[blor_col],
            df_summary[reg_col]
        ],
        default=np.nan
    )

# Calculate FDR and call ASE genes
valid_rows = df_summary['TOGA_binom_p'].notna()
p_values = df_summary.loc[valid_rows, 'TOGA_binom_p']

fdr_values = multipletests(p_values, method='fdr_bh')[1]

df_summary['TOGA_FDR'] = np.nan
df_summary.loc[valid_rows, 'TOGA_FDR'] = fdr_values

ase_cov = np.where(
    df_summary['TOGA_FDR'].isna(),
    None,  # Use None instead of np.nan for missing value (works with strings in object dtype)
    np.where(
        (df_summary['TOGA_log2_aFC'].abs() >= np.log2(args.haplotypic_ratio_cutoff)) & (df_summary['TOGA_FDR'] < args.binom_fdr_cutoff), 
        'ASE', 
        'COV'
    )
)
df_summary['TOGA_ASE_COV'] = ase_cov

# Add MHF(Major Haplotypic Fraction) ##
df_summary['TOGA_MHF'] = df_summary[['TOGA_aCount', 'TOGA_bCount']].max(axis=1) / df_summary['TOGA_totalCount']

# Add variant counts aligned to each haplotype
if args.var_hap_cnt_path and args.woOR_var_hap_cnt_path:
    print("Found files of variants counts for each haplotype:")
    print('var_hap_cnt_path', args.var_hap_cnt_path)
    print('woOR_var_hap_cnt_path', args.woOR_var_hap_cnt_path)
    
    # Load data
    var_hap_cnt=pd.read_csv(args.var_hap_cnt_path, sep='\t') # Extracted allelic counts for each haplotype
    woOR_var_hap_cnt=pd.read_csv(args.woOR_var_hap_cnt_path, sep='\t') # Extracted allelic counts for each haplotype(excluding variants from overlapping regions)

    columns_to_modify = ['varHapA', 'varHapB', 'varCntA', 'varCntB']
    var_hap_cnt = var_hap_cnt.apply(sort_variant_for_hap_counts, axis=1)
    woOR_var_hap_cnt = woOR_var_hap_cnt.apply(sort_variant_for_hap_counts, axis=1)
    
    df_summary_var_tmp = df_summary.copy()
    columns_to_include = ['varHapA', 'varHapB', 'varCntA', 'varCntB']
    # Ensure the columns are of object dtype, so strings can be assigned
    for col in columns_to_include:
        df_summary_var_tmp[col] = pd.Series([np.nan] * len(df_summary_var_tmp), dtype=object)
    
    matching_columns_1 = ['contig', 'start', 'stop', 'name', 'variants']
    matching_columns_2 = ['contig', 'start', 'stop', 'name', 'TOGA_variants']

    for index, row in df_summary_var_tmp.iterrows():
        tag = row['TAG']
        if tag == 'E':
            continue
        elif tag == 'RC':
            match = get_matching_rows(woOR_var_hap_cnt, matching_columns_1, row, matching_columns_2)
        else:
            match = get_matching_rows(var_hap_cnt, matching_columns_1, row, matching_columns_2)
        if not match.empty:
            df_summary_var_tmp.loc[index, columns_to_include] = match[columns_to_include].values[0]
    
    ## Rename columns ##
    df_summary = df_summary_var_tmp.copy()
    
    rename_dict_1 = {
        'varHapA': 'hapA_allele',
        'varHapB': 'hapB_allele',
        'varCntA': 'hapA_allele_count',
        'varCntB': 'hapB_allele_count'
    }
    df_summary = df_summary.rename(columns=rename_dict_1)

# Save the heavy version of df_summary
if args.wrt_heavy == 1:
    heavy_path = f'{args.output_dir}/{args.sample}_TOGA_heavy.txt'
    df_summary.to_csv(heavy_path, sep='\t', index=False)
    print(f"Saved heavy TOGA file: {heavy_path}")
print('Done.')

### STEP11: Create light version of df_summary ###
print('Saving the final result... ')
df_summary_light = df_summary.copy()

# Columns to keep
columns_to_keep = ['contig', 'start', 'stop', 'name', 'gene_symbol', 'type', 'TOGA_aCount', 'TOGA_bCount', 'TOGA_totalCount', 
                   'TOGA_log2_aFC', 'TOGA_MHF', 'hapA_allele', 'hapB_allele', 'hapA_allele_count', 'hapB_allele_count',
                  'TOGA_n_variants', 'TOGA_variants', 'TOGA_binom_p', 'TOGA_FDR', 
                   'overlapped_flags', 'overlapped_gene_name', 'overlapped_gene_symbol', 'overlapped_gene_type', 'overlapped_length', 'overlapped_region', 'variants_overlapping',
                   'gene_group_by_union_overlapping_variants', 'gene_cluster', 'norm_rna_exp', 'CPOV', 'IRPG', 'adj_totalCount', 'TAG', 'TOGA_ASE_COV'
                  ]
columns_to_keep_existing = [col for col in columns_to_keep if col in df_summary_light.columns]
df_summary_light = df_summary_light[columns_to_keep_existing]

# Rename the columns
rename_dict_2 = {
    'TOGA_aCount': 'hapA_count',
    'TOGA_bCount': 'hapB_count',
    'TOGA_totalCount': 'total_hap_count',
    'TOGA_log2_aFC': 'log2_aFC',
    'TOGA_MHF': 'MHF',
    'TOGA_n_variants': 'n_variants',
    'TOGA_variants': 'variants',
    'TOGA_binom_p': 'binom_p',
    'TOGA_FDR': 'binom_fdr',
    'TOGA_ASE_COV': 'ASE_COV',
    'gene_group_by_union_overlapping_variants': 'gene_group_by_variants',
    'gene_cluster': 'gene_group_by_regions',
    'adj_totalCount': 'adj_total_hap_count'
}
df_summary_light = df_summary_light.rename(columns=rename_dict_2)

# float to int
cols_float_to_int = ['n_variants', 'hapA_count', 'hapB_count', 'total_hap_count','gene_group_by_variants', 'gene_group_by_regions']
df_summary_light[cols_float_to_int] = df_summary_light[cols_float_to_int].apply(pd.to_numeric, errors='coerce')
df_summary_light[cols_float_to_int] = df_summary_light[cols_float_to_int].astype('Int64')

# Keep decimal to 6
col_to_round = ["CPOV", "adj_total_hap_count", "IRPG"]
df_summary_light[col_to_round] = df_summary_light[col_to_round].astype(float).round(6)

# nan -> NaN
df_summary_light[['overlapped_gene_name', 'overlapped_gene_symbol', 'overlapped_gene_type', 'overlapped_length', 'overlapped_region', "CPOV", "adj_total_hap_count", "IRPG"]] = (
    df_summary_light[['overlapped_gene_name', 'overlapped_gene_symbol', 'overlapped_gene_type', 'overlapped_length', 'overlapped_region', "CPOV", "adj_total_hap_count", "IRPG"]]
    .replace({'nan': np.nan})
)

# Save the light version of df_summary
light_path = f'{args.output_dir}/{args.sample}_TOGA_light.txt'
df_summary_light.to_csv(light_path, sep='\t', index=False)
print(f"Saved TOGA file: {light_path}")
print('Done.')
