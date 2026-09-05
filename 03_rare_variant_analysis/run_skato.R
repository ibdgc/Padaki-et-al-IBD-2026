library(data.table)
library(SKAT)

#------- ############# Input files #####################

# RAW_FILE:
#   PLINK2 --export A output restricted to variant set being tested.
#   Rows = individuals, columns = variant dosages.
# PSAM_FILE:
#   Sample metadata from PLINK2, including phenotypes and ancestry PCs.
# SET_FILE:
#   Gene : variant mapping, with variant IDs matching the .raw column names
#   (i.e. PLINK-style VARID_A1, such as 1-962412-C-T_C).

RAW_FILE  <- "mega_v5_AFR_filt_CADD20.raw"
PSAM_FILE <- "mega_v5_AFR_filt.psam"
SET_FILE  <- "cadd20_miss_sets.rawcols.setid"

#------- ################## Load genotype matrix ###################3

# Read the PLINK2 .raw file. The first 6 columns are metadata
# (FID, IID, PAT, MAT, SEX, PHENOTYPE); remaining columns are
# variant dosage columns.
geno <- fread(RAW_FILE)

# Save sample IDs in the exact order used by the genotype matrix
iid <- geno$IID

# Extract only genotype dosage columns
X <- as.matrix(geno[, -(1:6), with = FALSE])
colnames(X) <- colnames(geno)[-(1:6)]

cat("Genotype matrix:", nrow(X), "samples x", ncol(X), "variants\n")


#--------- ###################### Load PSAM and align samples ##################

# Load PLINK2 psam file containing phenotypes and covariates.
# Reorder rows to exactly match the IID order in the genotype matrix.
psam <- fread(PSAM_FILE)
psam <- psam[match(iid, psam$IID), ]

# Sanity check: all genotype samples must be present in psam
if (any(is.na(psam$IID))) {
  stop("IID mismatch between .raw and .psam")
}


#-------- ######### Construct phenotype and covariates (CD/UC/IBD vs controls) ############
# PLINK2 phenotype coding:
#   1 = control, 2 = case, 0/NA = missing
# Select status for analyses: CD/UC/IBD
pheno_raw <- psam$UC

# Keep only individuals with valid case/control labels
keep <- which(pheno_raw %in% c(1, 2))

# Recode phenotype to SKAT format:
#   0 = control, 1 = case
y <- ifelse(pheno_raw[keep] == 2, 1, 0)

# Subset genotype matrix and covariates to analyzed samples
Xk <- X[keep, , drop = FALSE]
Z  <- as.matrix(psam[keep, .(PC0, PC1, PC2, PC3, PC4)])

cat("Phenotype samples kept:", length(y),
    "(cases =", sum(y == 1),
    ", controls =", sum(y == 0), ")\n")


#-------- ############### Fit SKAT null model #####################

# Logistic regression of case/control  status on ancestry PCs only.
# This defines the baseline model; SKAT tests gene effects
# on top of this.
null_model <- SKAT_Null_Model(y ~ Z, out_type = "D")


#--------- ################ Load gene–variant sets ##############

# Each row of SET_FILE is:
#   GENE  var1  var2  var3  ...
# fill=Inf allows variable numbers of variants per gene.
sets <- fread(SET_FILE, header = FALSE, fill = Inf)
cat("Gene sets loaded:", nrow(sets), "\n")

# Variant IDs present in the genotype matrix (PLINK .raw columns)
geno_ids <- colnames(Xk)
cat("Example geno IDs:", paste(head(geno_ids, 5), collapse=", "), "\n")


#---------- ############# Run SKAT-O gene-based tests #############

# Results table:
#   GENE   = gene identifier
#   NVAR   = number of variants used in the test
#   P_SKATO = SKAT-O p-value (optimal.adj)
results <- data.table(
  GENE = character(),
  NVAR = integer(),
  P_SKATO = numeric()
)

for (i in 1:nrow(sets)) {

  # Gene label (symbol or RefSeq ID)
  gene <- sets[i, V1]

  # Extract variant IDs listed for this gene
  vids <- unlist(sets[i, -1, with = FALSE])
  vids <- vids[!is.na(vids)]

  # Retain only variants that exist in the genotype matrix
  vids <- vids[vids %in% geno_ids]

  # Require at least two variants to perform a set-based test
  # (single-variant genes are skipped)
  if (length(vids) < 2) next

  # Construct gene-specific genotype matrix
  G <- Xk[, vids, drop = FALSE]

  # Run SKAT-O (small-sample adjusted)
  p <- SKAT(G, null_model, method = "optimal.adj")$p.value

  # Store result
  results <- rbind(results, list(gene, length(vids), p))

  # Progress indicator
  if (i %% 200 == 0) {
    cat("Processed", i, "genes\n")
  }
}


#-------- ############### Write output ####################3

# Output file contains one row per tested gene
fwrite(results, "skat_one_percent_cadd20_UC.tsv", sep = "\t")

cat("Genes tested:", nrow(results), "\n")
cat("Done.\n")

