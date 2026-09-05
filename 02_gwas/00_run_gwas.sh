#!/bin/bash
#1. Subset to African ancestry
#2. Remove duplicates 
#3. Control–control comparison (QC step)
#     - snp exclusion list FDR < 1.0 
#4. Variant-level filtering (ctrl–ctrl FDR, HWE, missingness)
#5. Association testing for phenotypes (IBD, CD, UC)

plink2 --pfile mega_v2 \
  --keep-if ANCESTRY == AFR \
  --make-pgen --out mega_v3_AFR

plink2 --pfile mega_v3_AFR \
  --remove remove_plink.txt \
  --make-pgen --out mega_v3_AFR_nodup

plink2 --pfile mega_v3_AFR_nodup \
  --glm \
  --pheno-name CTRLTEST \
  --covar-name PC0,PC1,PC2,PC3,PC4 \
  --out mega_AFR_controltest

#Also filter for hwe and missingness
plink2 \
  --pfile mega_v3_AFR_nodup \
  --exclude afr.ctrltest_exclude.snplist \
  --geno 0.01 \
  --hwe 1e-7 midp \
  --make-pgen \
  --out mega_v3_filt_AFR

#Run GLMs on phenotypes 
plink2 \
  --pfile mega_v3_filt_AFR \
  --glm \
  --pheno-name IBD \
  --covar-name PC0,PC1,PC2,PC3,PC4 \
  --out afr_IBD

plink2 \
  --pfile mega_v3_filt_AFR \
  --glm \
  --pheno-name CD \
  --covar-name PC0,PC1,PC2,PC3,PC4 \
  --out afr_CD

plink2 \
  --pfile mega_v3_filt_AFR \
  --glm \
  --pheno-name UC\ 
  --covar-name PC0,PC1,PC2,PC3,PC4 \
  --out afr_UC






