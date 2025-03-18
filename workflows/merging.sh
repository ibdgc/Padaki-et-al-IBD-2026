#!/bin/sh

#at this point the plink files should contain proper variant ids and they should be consistant across the
# the two datasets.  
# Next we make a list of overlap variants that are either in bge or exonic in big_daly

./make_bge_exonic_list.pl bge_merge_v1.pvar big_daly.all_variants.txt big_daly_bge_subset.variants.txt

# big_daly_bge_subset.variants.txt contains a list of exonic variants in big_daly, and any non-coding
# variant found in both big_daly and bge

#create plink and vcf files with just the these overlap variants

plink2 --extract big_daly_bge_subset.variants.txt --pfile big_daly_v1 --make-pgen --out big_daly_v1_bgesubset
plink2 --pfile big_daly_v1_bgesubset --export vcf --out big_daly_v1_bgesubset
gzip big_daly_v1_bgesubset.vcf

#upload the big_daly_v1_bgesubset.vcf.gz to bystro and have it annotated.
#download the resutls from bystro

#Now do the susbseting the other way around.   Here we are creating a subset of merged_bge that
#includes all the exon variants from bge together with any non-coding variants found in big_daly_v1

./make_bge_exonic_list.pl big_daly_v1_bgesubset.pvar merged_bge.all_variants.txt  bge_subseted_bigdaly.txt 

# bge_subseted_bigdaly.txt is a list of exonic variants found in bge, together with all the non-coding
# variants that are found in both big_daly and bge.

plink2 --extract bge_subseted_bigdaly.txt  --pfile bge_merge_v1 --make-pgen --out bge_merge_v1_bigdalysubset
plink2 --pfile bge_merge_v1_bigdalysubset --export vcf --out bge_merge_v1_bigdalysubset
gzip bge_merge_v1_bigdalysubset.vcf

#upload the bge_merge_v1_bigdalysubset.vcf.gz  to bystro and have it annotated.
#download the resutls from bystro

# bge_merge_v1_bigdalysubset and big_daly_v1_bgesubset should be fully "compatible" and able to be merged
# together into a single dataset as bge_merge_v1_bigdalysubset contains only exonic variants called in bge,
# together with non-coding variants called in both sets, and big_daly_v1_bgesubset is the converse:   
# Exonic variants from big_daly and non-coding called in both big_daly and bge.
# We will check several of the broad scale qc metrics, and if nothing looks disastrously wrong, we
# will merge these too datasets, and start QC'ing the combined data.

# If this merge is going to work it will entirely come down to what happens at sites found in 
# one dataset not the other, and sites biallelic in one dataset, but multi-allelic in the other.
# plink2, currently, wants no part in trying that merger.  
# bcftools is probably the right choice for this.

# Time to convert from vcf.gz to vcf.bgz format

# All of this shit failed badly... none of it helped.
# gzip -c -d bge_merge_v1_bigdalysubset.vcf.gz |  bgzip -c > bge_merge_v1_bigdalysubset.vcf.bgz
# gzip -c -d big_daly_v1_bgesubset.vcf.gz |  bgzip -c > big_daly_v1_bgesubset.vcf.bgz
# forgot to make the index
# bgzip -r big_daly_v1_bgesubset.vcf.bgz
# bgzip -r bge_merge_v1_bigdalysubset.vcf.bgz
# create tabix files
# bcftools index big_daly_v1_bgesubset.vcf.bgz
# bcftools index bge_merge_v1_bigdalysubset.vcf.bgz
# bcftools merge --use-header bge_merge_v1_bigdalysubset.vcf.bgz  --missing-to-ref big_daly_v1_bgesubset.vcf.bgz bge_merge_v1_bigdalysubset.vcf.bgz  -o mega_bge_big_daly_v1.vcf.gz
# That was a waste of time.

# Use plink2 make the header it's good for nothing else at the moment.
plink2 --make-pgen --multiallelics-already-joined --out mega_bge_big_daly_v1 --pfile big_daly_v1_bgesubset --pmerge bge_merge_v1_bigdalysubset

# this merger step was a complete failure.   Does not work, but a nice .pvar header started.  We will use as the header 
# for our VCF file.
# had to write a custom script to complete the merge.   Neither plink2 nor bcftools wanted to help.
# using the header from plink2

./merge_vcf.pl mega_bge_big_daly_v1-merge.pvar big_daly_v1_bgesubset.vcf.gz bge_merge_v1_bigdalysubset.vcf.gz > mega_bge_big_daly_v1.vcf

# import into plink2 to test correctness of formating.  In this step we 
# A) Get rid of all sites not "PASS" quality
# B) Get rid of the X chromosome because only WGS has X data at the moment.  

plink2 --make-pgen --chr 1-22 --make-pgen --out mega_bge_big_daly_v1 --var-filter --vcf mega_bge_big_daly_v1.vcf

# Now let's get rid of of all variants that have 0 minor allele count, have only 1 allele, and drop 2 samples with a ton of missing data.
# output format is plink2

plink2 --pfile mega_bge_big_daly_v1 --mind 0.05 --mac 1 --min-alleles 2 --make-pgen --out mega_bge_big_daly_v1_1

# Now make a vcf file 
plink2 --pfile mega_bge_big_daly_v1_1 --export vcf --out mega_bge_big_daly_v1_1
gzip mega_bge_big_daly_v1_1.vcf.gz






