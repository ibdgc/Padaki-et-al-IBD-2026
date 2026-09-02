# African Ancestry IBD Genetic Analysis

This directory contains the analysis workflows, parameters, and commands used for the genetic analyses reported in the African ancestry IBD study. The repository is intended to document the analytical workflow and provide methodological details.

## Analysis workflows

- `01_data_harmonization/`  
  Harmonization and merging of whole-genome sequencing (WGS) and Blended Genome-Exome (BGE) datasets.

- `02_gwas/`  
  PLINK2 workflow for genome-wide association analyses of IBD, Crohn's disease (CD), and ulcerative colitis (UC).

- `03_rare_variant_analysis/`  
  Preparation of rare-variant gene sets and SKAT-O gene-based association testing.

- `04_finemapping/`  
  SuSiE fine-mapping workflows for the PTGER4, LINC01266, and MHC association signals using in-sample African-ancestry LD.

- `05_prs/`  
  Polygenic risk score construction using published GWAS effect sizes, ancestry-specific allele frequencies, and liability-scale conversion.

## Reproducibility

Each analysis directory contains a README describing the relevant workflow, parameters, and required inputs.

The software environment used for these analyses is documented in `environment.yml`.

