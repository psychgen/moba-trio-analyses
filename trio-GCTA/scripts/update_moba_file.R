library(tidyverse)
library(haven)

pheno <- read.table("//ess01/p471/data/durable/projects/qc_paper/data-prep-universal/scratch/phenotypes_outcomes_mqc.pheno", header = TRUE)
covs <- read.table("//ess01/p471/data/durable/projects/qc_paper/data-prep-universal/scratch/covariates_noyob_mqc.cov", header = TRUE)

pheno_old <- read.table("//ess01/p471/data/durable/projects/qc_paper/data-prep-universal/scratch/old/phenotypes_outcomes.pheno", header = TRUE) %>%
  rename("fid" = FID, 
         "iid" = IID)
covs_old <- read.table("//ess01/p471/data/durable/projects/qc_paper/data-prep-universal/scratch/old/covariates.cov", header = TRUE)  %>%
  rename("fid" = FID, 
         "iid" = IID) %>%
  select(fid,iid, starts_with("Age"))

Moba_dat <- read.table("N:/durable/projects/trio_gcta_qc/data/MoBaPsychGen_v1-ec-eur-batch-basic-qc_105210.moba", sep = ",", header = TRUE) 

Moba_dat2_child <- left_join(Moba_dat, pheno_old, by = c("fid", "iid")) %>%
  left_join(covs_old, by = c("fid", "iid") ) %>%
  slice_tail(n = 35070)# %>%
  #filter(!is.na(Age_height.x) & is.na(Age_height.y))
# seems that people with extreme age values were removed last time and not this time 

Moba_dat2_child %>%
  filter(height.x != height.y) %>%
  nrow()

Moba_dat2_child %>%
  filter(exam.x != exam.y) %>%
  nrow()

table(covs_old$Age_height)
table(covs_old$Age_exam)

Moba_dat2 <- Moba_dat %>%
  mutate(Age_exam = if_else(Age_exam != 10, NA, Age_exam), 
         Age_height = if_else(Age_height > 8, NA, Age_height), 
         Age_dep = if_else(Age_dep > 9 | Age_dep < 8, NA, Age_dep), 
         Age_sleep = if_else(Age_sleep > 8, NA, Age_sleep), 
         exam = if_else(is.na(Age_exam) & !is.na(exam), NA_integer_, exam, missing = NA_integer_),
         height = if_else(is.na(Age_height), NA_integer_, height, missing = NA_integer_), 
         smfq_dep = if_else(is.na(Age_dep), NA_integer_, smfq_dep, missing = NA_integer_), 
         sleep_hrs = if_else(is.na(Age_sleep), NA_integer_, sleep_hrs, missing = NA_integer_))



write.csv(Moba_dat2, paste0("//ess01/p471/data/durable/projects/trio_gcta_qc/data/",genefile,"v2_105210", ".moba"), row.names = F)


  