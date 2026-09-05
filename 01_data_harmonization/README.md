# Data Harmonization

Scripts in this directory were used to harmonize and merge the whole-genome sequencing (WGS) and Blended Genome-Exome (BGE) datasets prior to downstream analyses.

## Scripts

- `wgs_bge_harmonization_pipeline.sh`  
  Main workflow script coordinating dataset preparation, variant harmonization, subsetting, and VCF merging.

- `add_gnomad_to_pvar.pl`  
  Rewrites PLINK2 variant IDs using a `CHROM-POS-REF-ALT` format to enable consistent matching between datasets.

- `intersect_variants.pl`  
  Generates variant keep/remove lists for the WGS and BGE datasets. Autosomal variants are retained if they are shared between datasets or annotated as exonic; variants containing `*` alleles are excluded.

- `merge_split_vcf.pl`  
  Merges the harmonized WGS and BGE VCFs into a combined multi-sample VCF.

## Workflow

The main workflow performs:

- Input preparation
- PASS variant filtering
- Multiallelic splitting
- Variant ID harmonization
- Variant selection
- Dataset subsetting
- VCF merge

Run with:

```bash
sh wgs_bge_harmonization_pipeline.sh
```

The workflow requires PLINK2, bcftools, bgzip, pigz, and Perl. The software environment is documented in the repository-level `environment.yml`.

