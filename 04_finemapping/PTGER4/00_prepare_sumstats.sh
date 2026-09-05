#!/bin/bash
#1) Extract PTGER4 window SNP list from full GWAS (ADD only)
awk -F'\t' -v start=40176318 -v end=40676318 '
NR>1 && $1==5 && $2>=start && $2<=end && $8=="ADD" {print $3}
' afr_CD.CD.glm.logistic.hybrid > PTGER4.window250kb.snps.txt

# (optional but recommended)
sort -u PTGER4.window250kb.snps.txt > PTGER4.window250kb.snps.uniq.txt
mv PTGER4.window250kb.snps.uniq.txt PTGER4.window250kb.snps.txt

#2) Extract matching summary stats for SuSiE from full GWAS (ADD only)
awk -F'\t' -v start=40176318 -v end=40676318 '
NR==1 { print "chr","pos","snp","beta","se","z","p","n"; next }
$1==5 && $2>=start && $2<=end && $8=="ADD" {
  beta = log($10);
  print $1,$2,$3,beta,$11,$12,$13,$9
}' OFS='\t' afr_CD.CD.glm.logistic.hybrid > PTGER4.window250kb.sumstats.tsv

#3) Clean sumstats (drop inf/NA/SE==0), and derive “clean snps”
awk -F'\t' '
NR==1 {print; next}
($4=="inf" || $4=="-inf" || $5=="NA" || $6=="NA" || $7=="NA" || $5==0) {next}
{print}
' PTGER4.window250kb.sumstats.tsv > PTGER4.sumstats.clean.tsv

cut -f3 PTGER4.sumstats.clean.tsv | tail -n +2 > PTGER4.clean.snps.txt


