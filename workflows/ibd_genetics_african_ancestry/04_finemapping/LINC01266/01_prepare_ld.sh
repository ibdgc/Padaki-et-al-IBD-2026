#!/bin/bash
set -euo pipefail

# -------------------------
# USER EDITS
# -------------------------
PREFIX="LINC01266"
PFILE="../FDR_filtered/mega_v5_AFR_filt"   # plink2 prefix
THREADS=8
MAF=0.001
GENO=0.01
# -------------------------

# 4) Make locus-only BED from PLINK2 (for PLINK1 LD)
plink2 \
  --pfile "${PFILE}" \
  --extract "${PREFIX}.clean.snps.txt" \
  --rm-dup force-first \
  --make-bed \
  --out "${PREFIX}_final_tmp"

# 5) Compute LD with a MAF filter (prevents NaNs)
plink \
  --bfile "${PREFIX}_final_tmp" \
  --maf "${MAF}" \
  --geno "${GENO}" \
  --r square \
  --allow-no-sex \
  --threads "${THREADS}" \
  --out "${PREFIX}_LD_maf"

# Sanity check
echo "[INFO] NaN count in LD:"
grep -ci nan "${PREFIX}_LD_maf.ld" || true

# 6) Create a post-filter SNP list matching LD matrix order
plink \
  --bfile "${PREFIX}_final_tmp" \
  --maf "${MAF}" \
  --geno "${GENO}" \
  --make-bed \
  --allow-no-sex \
  --out "${PREFIX}_maf_tmp"

cut -f2 "${PREFIX}_maf_tmp.bim" > "${PREFIX}.maf.snps.txt"

# 7) Reorder sumstats to match post-filter SNP order
awk -F'\t' '
NR==FNR {order[++i]=$1; next}
NR==1 {hdr=$0; next}
{row[$3]=$0}
END{
  print hdr
  for (j=1; j<=i; j++) if (order[j] in row) print row[order[j]]
}
' "${PREFIX}.maf.snps.txt" "${PREFIX}.sumstats.clean.tsv" > "${PREFIX}.final.sumstats.maf.ordered.tsv"

# 8) Fix blank header line issue (if present)
awk 'NF>0{p=1} p' "${PREFIX}.final.sumstats.maf.ordered.tsv" > "${PREFIX}.final.sumstats.maf.ordered.noblank.tsv"

{
  echo -e "chr\tpos\tsnp\tbeta\tse\tz\tp\tn"
  cat "${PREFIX}.final.sumstats.maf.ordered.noblank.tsv"
} > "${PREFIX}.final.sumstats.maf.ordered.withheader.tsv"

echo "[INFO] Final sumstats rows: $(($(wc -l < ${PREFIX}.final.sumstats.maf.ordered.withheader.tsv) - 1))"
echo "[INFO] LD dim should be: $(wc -l < ${PREFIX}_LD_maf.ld)"

