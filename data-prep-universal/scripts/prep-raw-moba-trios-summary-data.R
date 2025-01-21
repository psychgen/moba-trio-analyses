# Prepare trios GWAS summary stats files for LDSC

# Required columns:
# SNP: The rsID or identifier of the SNP.
# CHR: Chromosome number.
# BP: Base-pair position.
# A1: Effect allele.
# A2: Non-effect allele.
# Z, OR, or BETA: Z-score, odds ratio, or effect size of the association.
# P: P-value of the association.
# N: Sample size for each SNP.

rm(list=ls())

library(data.table)
library(dplyr)

setwd("N:/durable/projects/qc_paper/trio-GWAS")

#ldsc_snps_file <-paste0("N:/durable/common/software/ldsc-master/eur_w_ld_chr/w_hm3.noMHC.snplist")
#ldsc_snps <- fread(ldsc_snps_file,h=T)

bimfile <- paste0("N:/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.bim") # for rsIDs
bim<-fread(bimfile)

head(bim)
bim <- bim %>%
  rename(Chromosome = V1, rsID = V2, Basepair = V4, A1 = V5, A2 = V6)

# Height

height_pop_n <- fread("sumstats/raw/model06_trios_height.details") # n = 21261
#height_pop<-fread("sumstats/raw/model06_trios_height.basic") # INFO>0.80
#height_con<-fread("sumstats/raw/model06_trios_height.trios")
height_pop<-fread("sumstats/info_0.95/height_pop_trim.txt") # INFO>0.95
height_con<-fread("sumstats/info_0.95/height_trio_trim.txt")

head(height_pop)

height_pop_trim <- height_pop %>%
  select(Chromosome, Basepair, A1, A2, Wald_Stat, Wald_P)

height_pop_trim <- height_pop_trim %>%
  merge(bim, by = c("Chromosome", "Basepair", "A1", "A2")) %>%
  rename(CHR = Chromosome, BP = Basepair, Z = Wald_Stat, P = Wald_P, SNP = rsID) %>%
  select(-V3)

head(height_con)
height_con_trim <- height_con %>%
  select(Chromosome, Basepair, A1, A2, Offspring_Z, Offspring_P) %>%
  merge(bim, by = c("Chromosome", "Basepair", "A1", "A2")) %>%
  rename(CHR = Chromosome, BP = Basepair, Z = Offspring_Z, P = Offspring_P, SNP = rsID) %>%
  select(-V3)

#write.table(height_pop_trim, file = "height_pop_ldsc.txt", quote = F, row.names = F, col.names = T, sep = "\t")
#write.table(height_con_trim, file = "height_con_ldsc.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(height_pop_trim, file = "height_pop_filtered_ldsc.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(height_con_trim, file = "height_con_filtered_ldsc.txt", quote = F, row.names = F, col.names = T, sep = "\t")

# EA 

ea_n <- fread("sumstats/raw/model06_trios_exam.details") # 42442
#ea_pop <- fread("sumstats/raw/model06_trios_exam.basic") # INFO>0.80
#ea_con <- fread("sumstats/raw/model06_trios_exam.trios")
ea_pop<-fread("sumstats/info_0.95/exam_pop_trim.txt") # INFO>0.95
ea_con<-fread("sumstats/info_0.95/exam_trio_trim.txt")

head(ea_pop)
ea_pop_trim <- ea_pop %>%
  select(Chromosome, Basepair, A1, A2, Wald_Stat, Wald_P) %>%
  merge(bim, by = c("Chromosome", "Basepair", "A1", "A2")) %>%
  rename(CHR = Chromosome, BP = Basepair, Z = Wald_Stat, P = Wald_P, SNP = rsID) %>%
  select(-V3)

head(ea_con)
ea_con_trim <- ea_con %>%
  select(Chromosome, Basepair, A1, A2, Offspring_Z, Offspring_P) %>%
  merge(bim, by = c("Chromosome", "Basepair", "A1", "A2")) %>%
  rename(CHR = Chromosome, BP = Basepair, Z = Offspring_Z, P = Offspring_P, SNP = rsID) %>%
  select(-V3)

#write.table(ea_pop_trim, file = "exam_pop_ldsc.txt", quote = F, row.names = F, col.names = T, sep = "\t")
#write.table(ea_con_trim, file = "exam_con_ldsc.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(ea_pop_trim, file = "exam_pop_filtered_ldsc.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(ea_con_trim, file = "exam_con_filtered_ldsc.txt", quote = F, row.names = F, col.names = T, sep = "\t")
