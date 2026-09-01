#!/bin/bash 

#4) Make locus-only BED from PLINK2 (for PLINK1 LD)
plink2 \
  --pfile FDR_filtered/mega_v5_AFR_filt \
  --extract PTGER4.clean.snps.txt \
  --rm-dup force-first \
  --make-bed \
  --out PTGER4_final_tmp

#B) Modified LD step (the fix that worked)
#5) Compute LD with a MAF filter (prevents NaNs)
plink \
  --bfile PTGER4_final_tmp \
  --maf 0.001 \
  --geno 0.01 \
  --r square \
  --allow-no-sex \
  --threads 8 \
  --out PTGER4_LD_maf

# Sanity Check 
grep -ci nan PTGER4_LD_maf.ld
# 0

#6) Create a “post-MAF” SNP list matching the LD matrix order
plink \
  --bfile PTGER4_final_tmp \
  --maf 0.001 \
  --geno 0.01 \
  --make-bed \
  --allow-no-sex \
  --out PTGER4_maf_tmp

cut -f2 PTGER4_maf_tmp.bim > PTGER4.maf.snps.txt

#7) Reorder sumstats to match the post-MAF SNP order
awk -F'\t' '
NR==FNR {order[++i]=$1; next}
NR==1 {hdr=$0; next}
{row[$3]=$0}
END{
  print hdr
  for (j=1; j<=i; j++) if (order[j] in row) print row[order[j]]
}
' PTGER4.maf.snps.txt PTGER4.sumstats.clean.tsv > PTGER4.final.sumstats.maf.ordered.tsv


#8) Fix the “blank header line” issue (if this occurs)
awk 'NF>0{p=1} p' PTGER4.final.sumstats.maf.ordered.tsv > PTGER4.final.sumstats.maf.ordered.noblank.tsv

{
  echo -e "chr\tpos\tsnp\tbeta\tse\tz\tp\tn"
  cat PTGER4.final.sumstats.maf.ordered.noblank.tsv
} > PTGER4.final.sumstats.maf.ordered.withheader.tsv

