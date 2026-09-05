#!/usr/bin/env bash
set -euo pipefail

########################################
# USER EDITS
########################################
GWAS="afr_UC.UC.glm.logistic.hybrid"          # full plink2 --glm output
PFILE_PREFIX="mega_v5_AFR_filt"  # AFR plink2 prefix (has .pgen/.pvar/.psam)
CLASS_BED="HLA_CLASS_I_II.merged_loci.bed" # combine both Loci for rerun 
OUTROOT="hla_finemap_out"
THREADS=16

# Filters (for LD stability + compute)
GENO=0.01      # variant missingness filter for LD dataset
MAC=20         # minimum minor allele count
L_SIGNALS=10   # SuSiE max number of effects

# R runner
SUSIE_R_SCRIPT="03_run_susie_block.R"
########################################

mkdir -p "${OUTROOT}"

# sanity checks
[[ -f "${GWAS}" ]] || { echo "[ERROR] Missing GWAS: ${GWAS}"; exit 1; }
[[ -f "${PFILE_PREFIX}.pgen" ]] || { echo "[ERROR] Missing pgen: ${PFILE_PREFIX}.pgen"; exit 1; }
[[ -f "${CLASS_BED}" ]] || { echo "[ERROR] Missing loci bed: ${CLASS_BED}"; exit 1; }
[[ -f "${SUSIE_R_SCRIPT}" ]] || { echo "[ERROR] Missing R script: ${SUSIE_R_SCRIPT}"; exit 1; }

command -v plink  >/dev/null || { echo "[ERROR] plink (1.9) not found in PATH"; exit 1; }
command -v plink2 >/dev/null || { echo "[ERROR] plink2 not found in PATH"; exit 1; }
command -v Rscript >/dev/null || { echo "[ERROR] Rscript not found in PATH"; exit 1; }

echo "[INFO] GWAS: ${GWAS}"
echo "[INFO] Genotypes: ${PFILE_PREFIX}"
echo "[INFO] Loci: ${CLASS_BED}"
echo "[INFO] Filters: GENO=${GENO} MAC=${MAC} ; SuSiE L=${L_SIGNALS}"
echo "[INFO] Threads: ${THREADS}"

i=0
while IFS=$'\t' read -r CHR START END; do
  i=$((i+1))
  TAG="$(basename "${CLASS_BED}" .merged_loci.bed).locus${i}.chr6_${START}_${END}"
  OUTDIR="${OUTROOT}/${TAG}"
  mkdir -p "${OUTDIR}"

  echo "=================================================="
  echo "[INFO] ${TAG}"
  echo "[INFO] Interval: ${CHR}:${START}-${END}"
  echo "=================================================="

  ########################################
  # 1) Extract sumstats (ADD only) for interval
  #    GWAS columns: 1 CHROM,2 POS,3 ID,8 TEST,9 N,10 OR,11 LOG(OR)_SE,12 Z,13 P
  ########################################
  awk -F'\t' -v start="${START}" -v end="${END}" '
    NR==1 { next }
    $1==6 && $2>=start && $2<=end && $8=="ADD" {
      beta = log($10);
      print $1,$2,$3,beta,$11,$12,$13,$9
    }' OFS='\t' "${GWAS}" > "${OUTDIR}/sumstats.raw.tsv"

  ########################################
  # 2) Clean sumstats and WRITE A GUARANTEED HEADER (prevents blank-header bug)
  ########################################
  {
    echo -e "chr\tpos\tsnp\tbeta\tse\tz\tp\tn"
    awk -F'\t' '
      ($4=="inf" || $4=="-inf" || $5=="NA" || $5==0 || $6=="NA" || $7=="NA") {next}
      {print}
    ' OFS='\t' "${OUTDIR}/sumstats.raw.tsv"
  } > "${OUTDIR}/sumstats.clean.tsv"

  NSUM=$(( $(wc -l < "${OUTDIR}/sumstats.clean.tsv") - 1 ))
  echo "[INFO] Clean sumstats SNPs: ${NSUM}"
  if [[ "${NSUM}" -lt 50 ]]; then
    echo "[WARN] Too few SNPs after cleaning; skipping locus."
    continue
  fi

  ########################################
  # 3) SNP list for genotype extraction
  ########################################
  cut -f3 "${OUTDIR}/sumstats.clean.tsv" | tail -n +2 > "${OUTDIR}/snps.clean.txt"

  ########################################
  # 4) Extract genotypes to BED (locus-only)
  ########################################
  plink2 \
    --pfile "${PFILE_PREFIX}" \
    --extract "${OUTDIR}/snps.clean.txt" \
    --rm-dup force-first \
    --make-bed \
    --threads "${THREADS}" \
    --out "${OUTDIR}/geno_tmp" >/dev/null

  ########################################
  # 5) Apply GENO + MAC filters for stable LD
  ########################################
  plink \
    --bfile "${OUTDIR}/geno_tmp" \
    --geno "${GENO}" \
    --mac "${MAC}" \
    --make-bed \
    --allow-no-sex \
    --out "${OUTDIR}/geno_mac" >/dev/null

  # Extract final SNP order from BIM (whitespace-safe)
  awk '{print $2}' "${OUTDIR}/geno_mac.bim" > "${OUTDIR}/snps.final.txt"
  NGENO=$(wc -l < "${OUTDIR}/snps.final.txt")
  echo "[INFO] SNPs after GENO+MAC: ${NGENO}"
  if [[ "${NGENO}" -lt 50 ]]; then
    echo "[WARN] Too few SNPs after GENO+MAC; skipping locus."
    continue
  fi

  ########################################
  # 6) Compute LD (r square)
  ########################################
  plink \
    --bfile "${OUTDIR}/geno_mac" \
    --r square \
    --allow-no-sex \
    --threads "${THREADS}" \
    --out "${OUTDIR}/LD" >/dev/null

  NLD=$(wc -l < "${OUTDIR}/LD.ld")
  echo "[INFO] LD rows: ${NLD}"
  if [[ "${NLD}" -ne "${NGENO}" ]]; then
    echo "[ERROR] LD rows (${NLD}) != SNPs (${NGENO})."
    exit 1
  fi

  ########################################
  # 7) Reorder sumstats to match BIM/LD order
  #    sumstats.clean.tsv is guaranteed to have the correct header now
  ########################################
  awk -F'\t' '
    NR==FNR {order[++i]=$1; next}
    NR==1 {hdr=$0; next}
    {row[$3]=$0}
    END{
      print hdr
      for (j=1; j<=i; j++) if (order[j] in row) print row[order[j]]
    }' "${OUTDIR}/snps.final.txt" "${OUTDIR}/sumstats.clean.tsv" > "${OUTDIR}/sumstats.ordered.tsv"

  # Defensive: if first line is blank for any reason, rebuild it
  if [[ "$(awk -F'\t' 'NR==1{print NF; exit}' "${OUTDIR}/sumstats.ordered.tsv")" -eq 0 ]]; then
    echo "[WARN] Blank header detected in sumstats.ordered.tsv; repairing..."
    { echo -e "chr\tpos\tsnp\tbeta\tse\tz\tp\tn"; awk 'NF>0' "${OUTDIR}/sumstats.ordered.tsv"; } \
      > "${OUTDIR}/sumstats.ordered.fixed.tsv"
    mv "${OUTDIR}/sumstats.ordered.fixed.tsv" "${OUTDIR}/sumstats.ordered.tsv"
  fi

  NORD=$(( $(wc -l < "${OUTDIR}/sumstats.ordered.tsv") - 1 ))
  echo "[INFO] Ordered sumstats rows: ${NORD}"
  if [[ "${NORD}" -ne "${NGENO}" ]]; then
    echo "[ERROR] Ordered sumstats rows (${NORD}) != SNPs (${NGENO})."
    echo "[DEBUG] First few snps.final:"
    head -5 "${OUTDIR}/snps.final.txt" | sed 's/^/[snps.final] /'
    echo "[DEBUG] First few sumstats.ordered:"
    head -6 "${OUTDIR}/sumstats.ordered.tsv" | sed 's/^/[sumstats.ordered] /'
    exit 1
  fi

  ########################################
  # 8) Run SuSiE
  ########################################
  Rscript "${SUSIE_R_SCRIPT}" \
    "${OUTDIR}/sumstats.ordered.tsv" \
    "${OUTDIR}/LD.ld" \
    "${OUTDIR}" \
    "${L_SIGNALS}"

done < "${CLASS_BED}"

echo "[DONE] Fine-mapping outputs under ${OUTROOT}/"

