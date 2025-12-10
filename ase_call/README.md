# Calling ASE Genes
**Last updated:** 09/03/2025

<img src="../image/ASE-Call.jpg" alt="screenshot" width="1000"/>

This script processes phASER geneAE outputs to generate COV, ASE, and MAE data. It produces per-sample and per-gene results, along with cohort-level summary statistics. Designed as a post hoc utility, it provides users with a simple and convenient way to analyze phASER results.

Reference: https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/

# Terminology
The pipeline uses the following definitions:
1. **COV**: Genes with total haplotypic read coverage ≥ threshold (default: 20). Coverage is the sum of reads mapped to both haplotypes.
2. **ASE**: Allelic-Specific Expressed genes. Subset of COV genes filtered by:
   * Exclusion of blacklist variants (via phASER processing)
   * Haplotypic count ratio ≥ threshold (default: 1.5), i.e. one haplotype has ≥ 1.5x more reads than the other haplotype.
   * Binomial test with FDR < threshold (default: 0.05), using:
        - Binomial p-values assuming a success probability of 0.5
        - Benjamini-Hochberg correction across all COV genes
3. **MAE**: Monoallelic Expressed genes. COV genes with all reads mapping to one haplotype (i.e., aCount = 0 or bCount = 0).

# Dependency
Runs on R≥3.5.x with the following packages: `optparse`, `dplyr`, `stringr`, `tidyr`

# Usage
```
Rscript ASE.R \
    -i ./example_in \
    -o ./example_out \
    -g ../resource/gencode.v.45.gene_coordinates.gene_symbol.bed \
    -c 20 \
    -r 1.5 \
    -f 0.05
```
> Note: The output of `2_phaser/` is written to `../data/output/phaser/` and works as an input to the `-i` parameter.

# Arguments
## Required
* -i, --input: Input directory containing sample subdirectories. Each must include phASER `*.gene_ae.txt` results.
* -o, --output: Output directory where results will be written. WARNING: If this directory exists, it will be deleted and overwritten.

## Optional 
* -g, --gene-coords (`../resource/gencode.v.45.gene_coordinates.gene_symbol.bed`): BED file of gene coordinates with gene symbols.
* -c, --cov-threshold (20): Minimum total haplotypic coverage for COV & ASE & MAE genes.
* -r, --hap-ratio (1.5): Minimum haplotypic count ratio (e.g., aCount/bCount or bCount/aCount) for identifying ASE genes.
* -f, --binom-fdr (0.05): FDR cutoff for binomial test to define ASE genes.
* -h, --help: Show help message with available arguments.

# Outputs
## Per-Sample Results
Gene-level haplotypic counts are generated for each individual. These files follow the structure of phASER's geneAE output, with two additional columns:
- **binom_p**: p-value from the binomial test  
- **binom_q**: FDR-adjusted p-value  

1. `/{sample}/{sample}_cov_genes.txt`  
   List of COV genes.  

2. `/{sample}/{sample}_ase_genes.txt`  
   List of ASE genes.  

3. `/{sample}/{sample}_mae_genes.txt`  
   List of MAE genes.  

4. `/{sample}/{sample}.summary_stats.txt`  
   Summary statistics of COV, ASE, and MAE genes in the individual:  
   - **n_genes**: Total number of genes from the coordinate file  
   - **n_cov_genes**: Number of COV genes  
   - **n_ase_genes**: Number of ASE genes  
   - **n_mae_genes**: Number of MAE genes  
   - **total_reads_cov**: Total RNA read counts across COV genes  
   - **total_vars_cov**: Total number of variants across COV genes  
   - **total_reads_ase**: Total read counts across ASE genes  
   - **total_vars_ase**: Total variant counts across ASE genes  
   - **total_reads_mae**: Total read counts across MAE genes  
   - **total_vars_mae**: Total variant counts across MAE genes  
   - **prop_cov**: proportion of all genes that met COV threshold. cov_gene ÷ total_gene  
   - **prop_ase**: ase_gene ÷ total_gene  
   - **prop_ase_cov**: ase_gene ÷ cov_gene  
   - **prop_mae**: mae_gene ÷ total_gene  
   - **prop_mae_cov**: mae_gene ÷ cov_gene  
   - **prop_mae_ase**: mae_gene ÷ ase_gene  

5. `all_samples.summary_stats.txt`  
   Cohort-level summary of COV, ASE, and MAE genes across all samples.  

## Per-Gene Results

1. `all_samples.{cov/ase/mae}_genes.txt`  
   Merged COV/ASE/MAE data across all samples. Each row corresponds to a single gene in a sample.  

2. `all_genes.samples_with_{cov/ase/mae}.txt`  
   Gene-level summary of sample counts:  
   - **ensembl**: ENSEMBL ID for gene
   - **gene_symbol**: HUGO ID for gene
   - **n_samples_{cov/ase/mae}**: Number of samples where the gene is COV/ASE/MAE  
   - **prop_{ase/mae}_cov**: Proportion of samples where the gene is ASE or MAE (`n_observed ÷ n_observed_cov`)  
   - **samples**: Sample IDs with ASE or MAE  
   - **samples_cov**: Sample IDs with COV but not ASE/MAE  

3. `all_genes.samples_with_{cov/ase/mae}.prop_sort.txt` 
   Same as `all_genes.samples_with_{cov/ase/mae}.txt`, but sorted by **prop_{ase/mae}_cov** in descending order.
