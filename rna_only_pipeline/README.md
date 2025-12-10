# Variant Calling for ASE analysis from RNA-seq
**Last updated:** 09/09/2025

<img src="../image/RNA_ONLY_PIPELINE.jpg" alt="screenshot" width="600"/>

Accurate identification of true heterozygous variants is essential for reliable allele-specific expression (ASE) quantification. However, variant calling from RNA-seq data alone can introduce artifacts due to poorly mapped reads or sequencing error. For example, in our study, we observed an average of 6 monoallelically expressed genes per individual where heterozygous variants in these genes were misclassified as homozygous. To overcome these limitations, we developed a workflow that integrates **STAR-WASP** with the **GATK germline SNP calling pipeline**, providing a more robust strategy to identify true heterozygous variants directly from RNA-seq data.

The entire pipeline is implemented in the **`run.sh`** script.

# References
### STAR-WASP
STAR: [Dobin et al., Bioinformatics, 2013](https://academic.oup.com/bioinformatics/article/29/1/15/272537)

WASP: [van de Geijn et al., Nature Methods, 2015](https://www.nature.com/articles/nmeth.3582)

STAR Manual (w/ WASP as an option): [STAR 2.7.11b Documentation (GitHub)](https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf)

### GATK SNP Calling from RNA-seq
RNA-seq Short Variant Discovery (SNPs & Indels): [GATK Documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels)

# Useful Files
[Gencode v45 Reference Genome FASTQ & GTF file](https://www.gencodegenes.org/human/release_45.html)

[Known snps, indels, interval lists used when calling variants](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/)
