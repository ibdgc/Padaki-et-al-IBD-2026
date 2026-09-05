#!/usr/bin/env Rscript
library(data.table)
library(susieR)
library(ggplot2)
library(Rfast)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 03_run_susie_block.R <sumstats.ordered.tsv> <LD.ld> <outdir> <L>")
}
sumstats_file <- args[1]
ld_file       <- args[2]
outdir        <- args[3]
L             <- as.integer(args[4])

ss <- fread(sumstats_file)
R  <- as.matrix(fread(ld_file, header=FALSE))

stopifnot(nrow(R) == ncol(R))
stopifnot(nrow(ss) == nrow(R))

z <- ss$beta / ss$se
fit <- susie_rss(z=z, R=R, n=median(ss$n), L=L)
cs  <- susie_get_cs(fit, coverage=0.95)

# PIP table
pip_tab <- data.table(
  snp=ss$snp, chr=ss$chr, pos=ss$pos,
  pip=fit$pip, p=ss$p, z=z
)
setorder(pip_tab, -pip)
fwrite(pip_tab, file.path(outdir, "susie.pip.tsv"), sep="\t")

# CS membership table
cs_list <- cs$cs
cs_dt <- rbindlist(lapply(seq_along(cs_list), function(k){
  idx <- cs_list[[k]]
  if (is.null(idx) || length(idx)==0) return(NULL)
  data.table(cs_id=paste0("L",k),
             snp=ss$snp[idx], chr=ss$chr[idx], pos=ss$pos[idx],
             pip=fit$pip[idx], p=ss$p[idx])
}), fill=TRUE)

if (!is.null(cs_dt) && nrow(cs_dt)>0) {
  setorder(cs_dt, cs_id, -pip)
  fwrite(cs_dt, file.path(outdir, "susie.credible_sets.tsv"), sep="\t")
}

# Purity table
pur <- cs$purity
if (!is.null(pur)) {
  pur_dt <- as.data.table(pur, keep.rownames="cs_id")
  sizes <- sapply(cs_list, length)
  pur_dt[, size := sizes[cs_id]]
  fwrite(pur_dt, file.path(outdir, "susie.cs_purity.tsv"), sep="\t")
}

# Summary text
sink(file.path(outdir, "susie.summary.txt"))
cat("n_snps:", nrow(ss), "\n")
cat("L:", fit$L, "\n")
cat("n_cs:", ifelse(is.null(cs_list), 0, length(cs_list)), "\n")
if (!is.null(cs_list)) cat("cs_sizes:", paste(sapply(cs_list, length), collapse=","), "\n")
if (!is.null(pur)) { cat("\nPurity:\n"); print(pur) }
cat("\nTop 20 by PIP:\n")
print(pip_tab[1:min(20,.N)])
sink()

# Simple plot: PIP vs position, highlight CS1 and lead SNP
lead_snp <- pip_tab[1]$snp
pip_tab[, category := "Other"]
if (!is.null(cs_list) && length(cs_list) >= 1 && !is.null(cs_list[[1]])) {
  pip_tab[snp %in% ss$snp[cs_list[[1]]], category := "CS1"]
}
pip_tab[snp == lead_snp, category := "Lead SNP"]
pip_tab[, category := factor(category, levels=c("Other","CS1","Lead SNP"))]

p <- ggplot(pip_tab, aes(x=pos, y=pip, color=category)) +
  geom_point(alpha=0.8, size=1.6) +
  theme_bw(base_size=12) +
  labs(x=paste0("Position (chr", unique(pip_tab$chr), ")"),
       y="PIP", title="SuSiE fine-mapping (PIP)") +
  scale_color_manual(values=c("Other"="grey70","CS1"="#d73027","Lead SNP"="#fdae61")) +
  guides(color=guide_legend(title=""))

ggsave(file.path(outdir, "pip_by_position.png"), p, width=8, height=4.5, dpi=300)

