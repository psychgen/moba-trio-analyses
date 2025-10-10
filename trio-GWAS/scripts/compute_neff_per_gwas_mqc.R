#!/usr/bin/env Rscript

# Calculate effective sample size per snp 

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript compute_neff_per_gwas.R <gwas_file> <pheno_file> <pheno>")
}
gwas_file     <- args[1]
pheno_file    <- args[2]
pheno_name    <- args[3]

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

# ---- Load GWAS (expects SNP, SE) ----
gwas <- fread(gwas_file)
stopifnot(all(c("SNP","SE") %in% names(gwas)))

# ---- Read snp.info (expects SNP, MAF) and merge ----
if (!file.exists("scratch/snp.info")) stop("snp.info not found at scratch/snp.info")
snpinfo <- fread("scratch/snp.info", select = c("SNP","MAF"))
stopifnot(all(c("SNP","MAF") %in% names(snpinfo)))
snpinfo[, MAF := as.numeric(MAF)]
gwas <- merge(gwas, snpinfo, by = "SNP", all.x = TRUE)

# ---- Get sd_resid ----
pheno <- fread(pheno_file)
if (!"phen" %in% names(pheno)) setnames(pheno, 1, "phen")
sd_resid <- pheno[phen == pheno_name, as.numeric(sd_resid)]
if (length(sd_resid) == 0) stop("sd_resid for phenotype '", pheno_name, "' not found in ", pheno_file)

# ---- Compute Neff ----
gwas[, `:=`(SE = as.numeric(SE), MAF = as.numeric(MAF))]
gwas[, Neff := NA_real_]
valid_idx <- with(gwas, which(!is.na(SE) & !is.na(MAF) & MAF > 0.01 & MAF < 1))
if (length(valid_idx) > 0) {
  gwas[valid_idx, Neff := round(1 / (SE[valid_idx]^2) * sd_resid^2 / (2 * MAF[valid_idx] * (1 - MAF[valid_idx])), 6)]
}

# ---- Overwrite the original GWAS file with Neff col included ----
fwrite(gwas, gwas_file, sep = "\t")
cat(sprintf("Calculated effective N and included as column 'Neff' in: %s (%d SNPs)\n", gwas_file, nrow(gwas)))