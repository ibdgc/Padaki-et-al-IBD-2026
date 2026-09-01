# Fine-Mapping

This directory contains scripts used for locus-specific fine-mapping of the PTGER4, LINC01266, and MHC association signals using SuSiE RSS.

## Directory structure

- `PTGER4/`  
  Fine-mapping workflow for the PTGER4 locus using a selected window around the lead association signal.

- `LINC01266/`  
  Fine-mapping workflow for the LINC01266 locus using a selected kb window around the lead association signal.

- `MHC/`  
  Fine-mapping workflow for the broader MHC region on chromosome 6 (`chr6:29,000,000-33,500,000`).

## General workflow

For PTGER4 and LINC01266, the workflow:

- Extracts locus-specific additive GWAS summary statistics
- Cleans variants with invalid or missing summary statistics
- Extracts matching variants from the African-ancestry analysis dataset
- Computes in-sample LD after variant filtering
- Reorders summary statistics to match the LD matrix
- Runs SuSiE RSS and generates posterior inclusion probabilities and 95% credible sets

The MHC workflow follows the same general fine-mapping framework but uses a predefined regional BED interval and applies genotype missingness and minor allele count filters before LD calculation.

## Requirements

The workflows require:

- PLINK2
- PLINK 1.9
- R
- R packages `data.table` and `susieR`
- `ggplot2` for MHC fine-mapping output plots

The broader software environment is documented in the repository-level `environment.yml`.
