# The landscape of allele-specific expression in human kidneys

This is the GitHub repository for the code used to analyze allele-specific expression in human kidneys.


# Codes in this Github repository

`1_rnaseq`: Aligning RNA-seq reads to the reference genome using [STAR](https://pmc.ncbi.nlm.nih.gov/articles/PMC3530905/), removing duplicates with [samtools markdup](https://www.htslib.org/doc/samtools-markdup.html) and removing reference bias with [WASP](https://pmc.ncbi.nlm.nih.gov/articles/PMC4626402/)

`2_phaser`: Haplotype phasing using [phASER](https://www.nature.com/articles/ncomms12817) and estimaing gene-level haplotypic counts using phASER geneAE

`3_phaser`: Transcript quantification using [RSEM](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)

`4_toga`: Resolving overlapping ASE genes and calling ASE Genes from phASER data.

`5_imprinting`: Use population-level ASE data to identify genomic imprinting candidates among genes.

`ase_call`: Alternative to `4_toga` for calling ASE genes from gene-level haplotypic counts without our overlapping genes analysis.

`rna_only_pipeline`: Alternative to `1_rnaseq` if only RNA-seq data present. Calling heterozygous variants from RNA-seq reads, optimized for ASE analysis

`DEG_GSEA`: Performing differential gene expression using [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) and gene set enrichment analysis using [clusterProfiler](https://www.liebertpub.com/doi/10.1089/omi.2011.0118)

`Cell-Fraction`: Estimating cell-type fraction of ASE genes using [KPMP single cell RNA-seq data](https://atlas.kpmp.org/repository/?size=n_20_n&filters%5B0%5D%5Bfield%5D=dois&filters%5B0%5D%5Bvalues%5D%5B0%5D=10.48698%2Fmgd0-gz70&filters%5B0%5D%5Btype%5D=any)

# Dependencies
This repository relies on R libraries, Python libraries, command line executables, and cloning github repositories to the local system. 
## 1. Local Repository Dependencies:

In the same directory as this repository, clone [phaser](https://github.com/secastel/phaser) and [RNA_Imprinting](https://github.com/zaitlenLab/RNA_Imprinting/).

Repository structure:
```
[user ASE_NEPTUNE]$ ls ../
ASE_NEPTUNE  phaser  RNA_Imprinting
```

## 2. Command-Line Executables:

The following software packages must be installed on your system and accessible from the command line: [STAR](https://github.com/alexdobin/STAR), [bcftools](https://samtools.github.io/bcftools/bcftools.html), [samtools](https://www.htslib.org/doc/samtools.html), [bedtools](https://bedtools.readthedocs.io/en/latest/), [rsem](https://github.com/deweylab/RSEM), [tabix](https://www.htslib.org/doc/tabix.html), [gatk](https://gatk.broadinstitute.org/hc/en-us), [picard](https://github.com/broadinstitute/picard)

## 3. R Libraries:

The R libraries should be installed within your R environment: [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [org.Hs.eg.dB](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html), [Seurat](https://satijalab.org/seurat/), [SeuratDisk](https://github.com/mojaveazure/seurat-disk), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [tidyverse](https://tidyverse.org/), [optparse](https://cran.r-project.org/web/packages/optparse/index.html), [data.table](https://r-datatable.com/), [R.matlab](https://cran.r-project.org/web/packages/R.matlab/index.html), [hdf5r](https://cran.r-project.org/web/packages/hdf5r/index.html), [rhdf5](https://bioconductor.org/packages/release/bioc/html/rhdf5.html), [arrow](https://arrow.apache.org/docs/r/)

## 4. Python Libraries:

These Python libraries are required: [pandas](https://pypi.org/project/pandas/), [numpy](https://numpy.org/), [statsmodels](https://www.statsmodels.org/stable/index.html), [scipy](https://scipy.org/)

# Reference

Ana C. Onuchic-Whitford†, Junmo Sung†, Eric D. Sakkas†, Michelle T. McNulty, Christopher L. O’Connor, Anya Greenberg, Jihoon G. Yoon, Sowmya Badina, NEPTUNE (Nephrotic Syndrome Study Network), Laura Mariani, Markus Bitzer, Matthew G. Sampson*, Dongwon Lee*. [The landscape of allele-specific expression in human kidneys](https://www.biorxiv.org/content/10.64898/2025.12.03.692205v1), BioRxiv, 2025

† These authors contributed equally to this work.

*Corresponding authors. Email: dongwon.lee@childrens.harvard.edu, matthew.sampson@childrens.harvard.edu
