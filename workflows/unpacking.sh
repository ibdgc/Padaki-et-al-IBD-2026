#!/bin/sh

#unpack the tarball and cleanup a bit
tar -xvf v1_merge_niddk_imputed_phased.tar
rm v1_merge_niddk_imputed_phased.tar

#unzip the annotation.  It takes a ton of room, but facilities command line tool efficiency
gzip -d merged_output_vcf-20250213020152097.annotation.tsv.gz

#there is SOOO much in the annotation, but for now all I really want are the variant IDs
cut -f 1,2,4,5,16,17,20,166 big_daly_vcf-20240330122349397.annotation.tsv > big_daly.all_variants.txt  
cut -f 1,2,4,5,16,17,20,166 merged_output_vcf-20250213020152097.annotation.tsv > merged_bge.all_variants.txt 

#make a plink2 format file out of the BGE file.  I already have a version for big_daly handy
plink2 --vcf merged_output.vcf.gz --make-pgen --out bge_merge_v1

#add rs numbers to the plink2 file so that variant ids are unique and can be referenced uniquely
add_rs_to_pvar.pl merged_bge.all_variants.txt bge_merge_v1.pvar bge_try1.pvar

#my version of big daly had been annotated years ago, so replace the variant ids from a newer version of bystro
#to ensure we are using the same dbsnp identifiers.

add_rs_to_pvar.pl big_daly.all_variants.txt big_daly_v1.pvar big_daly_v1.try1.pvar 

#replace stuff
mv bge_merge_v1.pvar bge_merge_v1.pvar.save
mv bge_try1.pvar bge_merge_v1.pvar
mv big_daly_v1.pvar big_daly_v1.pvar.save
mv big_daly_v1.try1.pvar big_daly_v1.pvar

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

