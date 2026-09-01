#!/bin/bash
set -euo pipefail

# -------------------------
# USER EDITS
# -------------------------
PREFIX="LINC01266"
CHR=3              # e.g. 5
LEAD_POS=616521         # e.g. 123456789  (hg38 position of your lead SNP)
#WINDOW=250000
WINDOW=200000
GWAS="afr_CD.CD.glm.logistic.hybrid"   # change if using UC/IBD, etc.
# -------------------------

START=$((LEAD_POS - WINDOW))
END=$((LEAD_POS + WINDOW))

echo "[INFO] ${PREFIX} window: chr${CHR}:${START}-${END}"

# 1) Extract window SNP list from full GWAS (ADD only)
awk -F'\t' -v chr="${CHR}" -v start="${START}" -v end="${END}" '
NR>1 && $1==chr && $2>=start && $2<=end && $8=="ADD" {print $3}
' "${GWAS}" > "${PREFIX}.window200kb.snps.txt"

sort -u "${PREFIX}.window200kb.snps.txt" -o "${PREFIX}.window200kb.snps.txt"

# 2) Extract matching summary stats for SuSiE (ADD only)
# NOTE: PTGER4 code uses beta = log(OR) from column $10 and SE from $11, Z from $12, P from $13, N from $9.
awk -F'\t' -v chr="${CHR}" -v start="${START}" -v end="${END}" '
NR==1 { print "chr","pos","snp","beta","se","z","p","n"; next }
$1==chr && $2>=start && $2<=end && $8=="ADD" {
  beta = log($10);
  print $1,$2,$3,beta,$11,$12,$13,$9
}' OFS='\t' "${GWAS}" > "${PREFIX}.window200kb.sumstats.tsv"

# 3) Clean sumstats (drop inf/NA/SE==0), and derive clean SNP list
awk -F'\t' '
NR==1 {print; next}
($4=="inf" || $4=="-inf" || $5=="NA" || $6=="NA" || $7=="NA" || $5==0) {next}
{print}
' "${PREFIX}.window200kb.sumstats.tsv" > "${PREFIX}.sumstats.clean.tsv"

cut -f3 "${PREFIX}.sumstats.clean.tsv" | tail -n +2 > "${PREFIX}.clean.snps.txt"

echo "[INFO] Clean SNPs: $(wc -l < ${PREFIX}.clean.snps.txt)"

