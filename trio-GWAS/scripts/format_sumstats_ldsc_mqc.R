# Prepare trios gwas summary stats files for LDSC 
# IB aug 25 


rm(list=ls())

#---------------- command-line parsing ----------------#
args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage: Rscript format_sumstats_ldsc_mqc.R <POP_TABLE> <COND_TABLE> <PREFIX>",
    "[BIM_FILE] [OUT_DIR]\n", file = stderr())
  quit(status = 1)
}

if (length(args) < 3 || length(args) > 5 ) usage()

popfile     <- args[1]
confile     <- args[2]
prefix      <- args[3]
bim_file <- if (length(args) >= 4) args[4] else "scratch/MoBaPsychGen_v1-ec-eur-batch-basic-qc-exc.bim"
out_dir     <- if (length(args) == 5) args[5] else "results/ldsc_aug25"

stopifnot(file.exists(popfile), file.exists(confile), file.exists(bim_file))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Packages
library(data.table)
library(dplyr)

#----------------------- Format ldsc ----------------------------------------#

format_ldsc <- function(dt, snpref) {
  
  # standardise column names
  setnames(dt, intersect(names(dt), c("Wald_P", "Offspring_P")), "P")
  if ("Offspring_Effect" %in% names(dt)) setnames(dt, "Offspring_Effect", "Effect")
  setnames(dt, intersect(names(dt), c("Wald_Stat", "Offspring_Z")), "Z")
  setnames(dt, intersect(names(dt), c("SD", "Offspring_SD")), "SE")
  setnames(dt, intersect(names(dt), "Chromosome"), "chr")
  setnames(dt, intersect(names(dt), "Basepair"), "pos")
  
  #dt[, c("chr", "pos") := tstrsplit(MarkerName, ":", keep = 1:2)]
  #dt[, chr := sub("^chr", "", chr)]
  dt[, pos := as.integer(pos)]
  
  setkey(dt, chr, pos)
  rsid_dt <- snpref[dt]
  rsid_dt <- rsid_dt[!duplicated(Predictor)]
  
  #rsid_dt <- rsid_dt[Freq1 >= 0.1 & Freq1 <= 0.4]
  
  rsid_dt %>%
    rename(SNP = Predictor) %>%
    select(chr, pos, SNP, rsID, Z, Effect, SE, P, A1, A2, Num_Samples)
}

#----------------------- Load rsids ref -------------------------------------#
bim <- fread(bim_file,
                select    = c(1, 2, 4),
                col.names = c("chr", "rsID", "pos"))

setkey(bim, chr, pos)

#-------------------- Format pop and trios sumstats -------------------------#

pop_ldsc <- format_ldsc(fread(popfile), bim)
con_ldsc <- format_ldsc(fread(confile), bim)

#---------------------- Write outputs----------------------------------------#
pop_out <- file.path(out_dir, paste0(prefix, "_pop_ldsc.txt"))
con_out <- file.path(out_dir, paste0(prefix, "_trio_ldsc.txt"))

fwrite(pop_ldsc, pop_out, sep = "\t")
fwrite(con_ldsc, con_out, sep = "\t")

cat("Written:\n  ", normalizePath(pop_out), "\n  ", normalizePath(con_out), "\n")