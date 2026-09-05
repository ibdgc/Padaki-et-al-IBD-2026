#!/bin/bash

################# 1) Generate the CADD20 missense gene–variant list #################
$(which perl) ./make_skat_lists_plink2.pl \
  ANNOTATION_FILE \
  mega_v5_AFR_filt.pvar \
  0.01 20 y n > cadd20_miss_variants_skat.txt

#Note: This should again produce a 2-column file:
#GENE   VARID
#GENE   VARID
#...



############### (2) Create ids + setid files #############
# variant list
cut -f2 cadd20_miss_variants_skat.txt | sort -u > cadd20_miss_variants.ids

# gene set file (gene followed by all its variants)
awk '
{
  g=$1; v=$2;
  if(!(g in seen)){ order[++n]=g; seen[g]=1 }
  sets[g]=sets[g] " " v
}
END{
  for(i=1;i<=n;i++){
    g=order[i];
    print g sets[g]
  }
}' cadd20_miss_variants_skat.txt > cadd20_miss_sets.setid


################ (3) Export genotypes for these variants (PLINK2) ################
plink2 --pfile FDR_filtered/mega_v5_AFR_filt \
  --extract cadd20_miss_variants.ids \
  --export A \
  --out mega_v5_AFR_filt_CADD20


############## (4) Match the set IDs to the .raw header ###############
#.raw columns will  look like: VARID_A1 (e.g., 1-962412-C-T_C). Format variant list in the same manner.

awk '{
  gene=$1;
  printf "%s", gene;

  for (i=2; i<=NF; i++) {
    v=$i;
    n=split(v, a, "-");
    ref=a[3];
    printf " %s_%s", v, ref;
  }
  printf "\n"
}' cadd20_miss_sets.setid > cadd20_miss_sets.rawcols.setid

