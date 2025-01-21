# MoBa QC paper: filter summary statitics by INFO > 0.95

rm(list=ls())

# Packages
library(data.table)
library(dplyr)

# Set wd
setwd("N:/durable/projects/qc_paper/trio-GWAS")
dir ="N:/durable/projects/trio_gwas/"

# File paths
height_pop_file <- paste(dir, "results/03/height/model06_trios_height.basic", sep="")
height_trio_file <- paste(dir, "results/03/height/model06_trios_height.trios", sep="")
exam_pop_file <- paste(dir, "results/03/exam/model06_trios_exam.basic", sep="")
exam_trio_file <- paste(dir, "results/03/exam/model06_trios_exam.trios", sep="")
snpinfo_file <- paste(dir, "scratch/snp.info", sep="")

# INFO scores
snpinfo <- fread(snpinfo_file)

head(snpinfo)
summary(snpinfo$INFO_release4)

snpinfo <- snpinfo %>%
  rename(Predictor = "SNP")

# Height 
height_pop <- fread(height_pop_file)
height_trio <- fread(height_trio_file)

# Exam scores
exam_pop <- fread(exam_pop_file)
exam_trio <- fread(exam_trio_file)

# Merge and merge

height_pop_trim <- height_pop %>%
  left_join(snpinfo, by="Predictor") %>%
  filter(INFO_release4 >= 0.95)

height_trio_trim <- height_trio %>%
  left_join(snpinfo, by="Predictor") %>%
  filter(INFO_release4 >= 0.95)

exam_pop_trim <- exam_pop %>%
  left_join(snpinfo, by="Predictor") %>%
  filter(INFO_release4 >= 0.95)

exam_trio_trim <- exam_trio %>%
  left_join(snpinfo, by="Predictor") %>%
  filter(INFO_release4 >= 0.95)

# Save filtered summary data

write.table(height_pop_trim, "N:/durable/projects/qc_paper/trio-GWAS/sumstats/info_0.95/height_pop_trim.txt", col.names=T, row.names=F, quote=F, sep='\t')
write.table(height_trio_trim, "N:/durable/projects/qc_paper/trio-GWAS/sumstats/info_0.95/height_trio_trim.txt", col.names=T, row.names=F, quote=F, sep='\t')

write.table(exam_pop_trim, "N:/durable/projects/qc_paper/trio-GWAS/sumstats/info_0.95/exam_pop_trim.txt", col.names=T, row.names=F, quote=F, sep='\t')
write.table(exam_trio_trim, "N:/durable/projects/qc_paper/trio-GWAS/sumstats/info_0.95/exam_trio_trim.txt", col.names=T, row.names=F, quote=F, sep='\t')
