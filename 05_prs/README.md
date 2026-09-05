# Polygenic Risk Score Analysis

This directory contains the core files used for the polygenic risk score (PRS) workflow.

## Files

- `calc_prs.pl`  
  Main PRS script. It uses ancestry information, sample metadata, GWAS effect sizes, and variant annotations to construct PRS values. The script also performs odds-ratio to liability-scale conversion and adjusts score expectations using ancestry-specific allele frequencies.

- `all_liu_effects.txt`  
  Harmonized GWAS effect-size input file containing effect sizes and allele frequencies reported relative to the reference allele.

## Additional required inputs

The workflow also requires:

- A Bystro ancestry-score file for the analyzed samples
- A PLINK `.fam` or `.psam` sample metadata file
- A Bystro annotation file containing the variants represented in the GWAS effect-size file

These individual-level or consortium-derived inputs are not included in this repository.

## Example usage

```bash
./calc_prs.pl   ancestry_scores.csv   samples.fam   all_liu_effects.txt   annotations.txt   0.0066   1e-8   500000
```

Arguments:

1. Bystro ancestry-score file
2. PLINK `.fam` or `.psam` sample file
3. Harmonized GWAS effect-size file
4. Bystro annotation file
5. Disease prevalence parameter used for liability-scale conversion
6. P-value threshold for variant inclusion
7. Window size, in base pairs, used for simple P+T-style thinning

## Variant selection

Variant selection uses a simple P+T approach. Variants passing the specified p-value threshold and located within the specified genomic window are grouped, and the variant with the smallest p-value in each group is retained.

A very small window can be used when pruning is not desired and all supplied variants should be retained.

## Notes

The example above uses a prevalence parameter of `0.0066`, a p-value threshold of `1e-8`, and a thinning window of `500000` bp. These values should be adjusted as appropriate for the phenotype and input variant set being analyzed.

The broader software environment is documented in the repository-level `environment.yml`.
