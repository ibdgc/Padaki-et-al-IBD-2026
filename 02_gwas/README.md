# Genome-Wide Association Analysis

This directory contains the PLINK2 workflow used to perform genome-wide association analyses in the African-ancestry cohort.

## Script

- `00_run_gwas.sh`  
  Subsets the analysis dataset to African-ancestry samples, removes duplicate samples, performs a control-control association test for variant-level QC, applies variant filtering, and runs association tests for IBD, Crohn's disease (CD), and ulcerative colitis (UC).

## Workflow

- Subset to African-ancestry samples
- Remove duplicate samples
- Run control-control association testing
- Exclude variants identified by the control-control QC step
- Apply genotype missingness and Hardy-Weinberg equilibrium filters
- Run additive logistic regression for IBD, CD, and UC
- Adjust association models for PC0-PC4

## Requirements

The workflow requires PLINK2.

The broader software environment is documented in the repository-level `environment.yml`.

