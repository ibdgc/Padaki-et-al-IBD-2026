# Rare Variant Analysis

Scripts in this directory were used to prepare rare-variant gene sets and run gene-based SKAT-O association analyses.

## Scripts

- `generate_cadd_missense_variants.sh`  
  Prepares the CADD≥20 rare-variant subset used for SKAT-O analysis. The workflow applies a 1% frequency threshold, generates gene–variant sets, extracts the corresponding genotypes with PLINK2, and reformats variant IDs to match the exported genotype matrix.

- `make_skat_lists_plink2.pl`  
  Filters annotated variants by allele frequency, CADD score, and functional consequence, assigns qualifying variants to genes, and generates gene–variant lists using variants present in the PLINK2 dataset.

- `run_skato.R`  
  Runs SKAT-O gene-based rare-variant association tests using PLINK2-exported genotype dosages and gene-based variant sets. The null model adjusts for the first five ancestry principal components (PC0–PC4), and genes with fewer than two qualifying variants are skipped.

## Workflow

- Filter and annotate qualifying rare variants
- Generate gene–variant lists
- Extract variant genotypes with PLINK2
- Format gene sets to match PLINK2 genotype column names
- Run SKAT-O association tests for CD, UC, or IBD

## Requirements

The workflow requires:

- Perl
- PLINK2
- R
- R packages `data.table` and `SKAT`

The broader software environment is documented in the repository-level `environment.yml`.

