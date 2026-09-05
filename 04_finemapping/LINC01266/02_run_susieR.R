library(data.table)
library(susieR)

# -------------------------
# USER EDITS
# -------------------------
PREFIX <- "LINC01266"
L <- 10
# -------------------------

sumstats_file <- paste0(PREFIX, ".final.sumstats.maf.ordered.withheader.tsv")
ld_file       <- paste0(PREFIX, "_LD_maf.ld")

ss <- fread(sumstats_file)
R  <- as.matrix(fread(ld_file, header = FALSE))

stopifnot(nrow(ss) == nrow(R), ncol(R) == nrow(R))

z <- ss$beta / ss$se
fit <- susie_rss(z = z, R = R, n = median(ss$n), L = L)

out <- data.table(
  snp = ss$snp,
  chr = ss$chr,
  pos = ss$pos,
  pip = fit$pip,
  p   = ss$p,
  z   = z
)
setorder(out, -pip)
fwrite(out, "LINC01266.susie.pip.tsv", sep = "\t")


cat("\nTop 20 SNPs by PIP:\n")
print(out[1:20])

cs <- susie_get_cs(fit, coverage = 0.95)

cat("\nCredible set sizes:\n")
if (!is.null(cs$cs)) {
  print(sapply(cs$cs, length))
} else {
  cat("No credible sets returned.\n")
}

# CS1 table
if (!is.null(cs$cs) && length(cs$cs) >= 1) {
  cs1_idx <- cs$cs[[1]]
  cs1 <- ss[cs1_idx, .(chr, pos, snp, beta, se, z, p)]
  cs1[, pip := fit$pip[cs1_idx]]
  setorder(cs1, -pip)

  fwrite(cs1, paste0(PREFIX, ".CS1.variants.tsv"), sep = "\t")

  cat("\nCS1 variants (written to ", paste0(PREFIX, ".CS1.variants.tsv"), "):\n", sep="")
  print(cs1)
}

