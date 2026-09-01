#!/bin/sh

### STEP 1 — Unpack the tarball and clean up
# The tarball likely contains the BGE (blended exome + imputation) VCF and annotation.
tar -xvf v1_merge_niddk_imputed_phased.tar
rm v1_merge_niddk_imputed_phased.tar

### STEP 2 — Expand the annotation file
# Decompresses the large annotation TSV so tools like cut/grep run faster.
gzip -d merged_output_vcf-20250213020152097.annotation.tsv.gz

### STEP 3 — Extract variant ID–related columns from the annotation
# From each annotation file, grab certain fields (chromosome, position, ref, alt, + key annotations).
# Output: text files with only those relevant columns for both datasets.
cut -f 1,2,4,5,16,17,20,166 big_daly_vcf-20240330122349397.annotation.tsv > big_daly.all_variants.txt
cut -f 1,2,4,5,16,17,20,166 merged_output_vcf-20250213020152097.annotation.tsv > merged_bge.all_variants.txt

### STEP 4 — Convert BGE VCF into PLINK2 format
# Produces .pgen (genotypes), .pvar (variants), .psam (samples).
plink2 --vcf merged_output.vcf.gz --make-pgen --out bge_merge_v1

### STEP 5 — Prepare big_daly for harmonization
# big_daly = large WGS-based dataset.
# Keep only PASS variants (bcftools has issues with non-PASS records).
plink2 --pfile big_daly_v1 --var-filters --make-pgen --out big_daly_v1_pass_only

### Export to VCF for bcftools processing
plink2 --pfile big_daly_v1_pass_only --export vcf bgz --out big_daly_v1_pass_only

### Split multiallelics into biallelic records
bcftools norm -m- big_daly_v1_pass_only.vcf.gz -Oz -o big_daly_v1_pass_split.vcf.gz

### Convert back to PLINK2 format after splitting
plink2 --make-pgen --vcf big_daly_v1_pass_split.vcf.gz --out big_daly_v1_pass_split

### STEP 6 — Harmonize variant IDs to gnomAD style
# Ensures both datasets use the same CHR-POS-REF-ALT naming.
./add_gnomad_to_pvar.pl bge_merge_v1.pvar bge_merge_v1_gnomad.pvar
mv bge_merge_v1.pvar bge_merge_v1.pvar.save
mv bge_merge_v1_gnomad.pvar bge_merge_v1.pvar

./add_gnomad_to_pvar.pl big_daly_v1_pass_split.pvar big_daly_v1_pass_split_gnomad.pvar
mv big_daly_v1_pass_split.pvar big_daly_v1_pass_split.pvar.save
mv big_daly_v1_pass_split_gnomad.pvar big_daly_v1_pass_split.pvar

### STEP 7 — Define keep/remove variant lists
# Uses annotation + custom Perl script to keep:
# - chr1–22 only
# - no variants with "*" alleles
# - all exonic variants
# - all variants present in both datasets
./intersect_variants.pl big_daly_v1_pass_split big_daly.all_variants.txt bge_merge_v1 merged_bge.all_variants.txt

# At this point, summary counts show how many variants were kept/removed.

### STEP 8 — Subset both datasets to the keep lists
plink2 --pfile big_daly_v1_pass_split --out big_daly_v1_subset --make-pgen --extract big_daly_v1_pass_split.keep.txt
plink2 --pfile bge_merge_v1          --out bge_merge_v1_subset  --make-pgen --extract bge_merge_v1.keep.txt

### STEP 9 — Export subsets to VCF for merging
plink2 --pfile big_daly_v1_subset --export vcf bgz --out big_daly_v1_subset
plink2 --pfile bge_merge_v1_subset --export vcf bgz --out bge_merge_v1_subset


### STEP 10 — Try standard bcftools merge
bgzip -r big_daly_v1_subset.vcf.gz
bgzip -r bge_merge_v1_subset.vcf.gz
bcftools index -t big_daly_v1_subset.vcf.gz
bcftools index -t bge_merge_v1_subset.vcf.gz

# If bcftools merge attempt fails in case headers/allele definitions are too different, STEP 11.

### STEP 11 — Fall back to custom Perl merge
# Convert bgz to  plain gz (Perl script expects standard gzip).
mv bge_merge_v1_subset.vcf.gz bge_merge_v1_subset.vcf.bgz
mv big_daly_v1_subset.vcf.gz   big_daly_v1_subset.vcf.bgz
bgzip -@8 -d bge_merge_v1_subset.vcf.bgz  -c | pigz -p 56 > bge_merge_v1_subset.vcf.gz
bgzip -@8 -d big_daly_v1_subset.vcf.bgz   -c | pigz -p 56 > big_daly_v1_subset.vcf.gz

# Custom merge using bge_header.txt as canonical header.
./merge_split_vcf.pl bge_header.txt big_daly_v1_subset.vcf.gz bge_merge_v1_subset.vcf.gz | pigz -p 56 > mega_v1.vcf.gz

