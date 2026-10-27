# Genetic regulation of the sperm methylome in cattle

## Project Overview

This repository provides analysis scripts for investigating the genetic regulation of DNA methylation in bovine sperm. Using whole-genome sequencing (WGS) and whole-genome bisulfite sequencing (WGBS) data from 125 bulls, the study examines associations between DNA methylation blocks (meBlocks) and four classes of genetic variants: single-nucleotide polymorphisms (SNPs), insertions and deletions (InDels), short tandem repeats (STRs), and structural variants (SVs).

The analyses include variant identification, DNA methylation quantification, meQTL mapping, fine-mapping, and integration with GWAS of complex traits.

## Analysis Scripts

### 1. Genotype identification

- **`1. SNP and Indel Genotype identification.sh`**: SNP and InDel calling and filtering from WGS data.
- **`1. STR Genotype identification.sh`**: STR genotyping and integration using GangSTR, ExpansionHunter, HipSTR, and EnsembleTR.
- **`1. SV Genotype identification.sh`**: SV discovery and genotyping using multiple tools, including BayesTyper.

### 2. DNA methylation quantification

- **`2. Mask CpG SNPs.R`**: Exclusion of CpG sites overlapping SNPs at the cytosine or adjacent guanine position.
- **`3. DNA methylation quantification.sh`**: WGBS processing, methylation quantification, and meBlock construction.

### 3. Genetic association and fine-mapping

- **`4. meGene identification.py`**: cis-meQTL mapping for DNA methylation blocks using tensorQTL.
- **`5. Fine-mapping.r`**: Fine-mapping of meQTL signals using SuSiE to identify credible sets and posterior inclusion probabilities.

### 4. Chromatin annotation

- **`6. MACS3 peak calling.sh`**: ATAC-seq peak calling using MACS3.
- **`7. TOBIAS footprinting.sh`**: Transcription factor footprinting and binding-site analysis using TOBIAS.

### 5. Integration with complex traits

- **`8. Cololocalization.r`**: meQTL–GWAS colocalization analysis using coloc.abf.
- **`9. Methylation-wide association study.sh`**: SPrediXcan-based methylation-wide association analysis using GWAS summary statistics.

## Notes

Scripts illustrate the study-specific analysis workflow and may require adjustment of file paths, software environments, input formats, and computational resources before use. Genomic coordinates are based on the ARS-UCD1.2 bovine reference assembly.

## Contact Information

For questions or suggestions, please contact:
- **Ying Yu**: yuying@cau.edu.cn
- **Lingzhao Fang**: lingzhao.fang@qgg.au.dk
