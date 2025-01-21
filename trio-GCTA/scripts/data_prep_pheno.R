#author: Espen Eilertsen, Laura Hegemann 
#prep phenotypic data / select sample for trio-GCTA

library(tidyverse)
library(haven)

# data

pheno <- read.table("//ess01/p471/data/durable/projects/qc_paper/data-prep-universal/scratch/phenotypes_outcomes.pheno", header = TRUE)
covs <- read.table("//ess01/p471/data/durable/projects/qc_paper/data-prep-universal/scratch/covariates.cov", header = TRUE)

dup_ids <- read.table("//ess01/p471/data/durable/projects/qc_paper/data-prep-universal/scratch/duplicated_ids.txt", header = TRUE)

pheno_fid <- pheno$FID

#trio

pheno <- pheno %>%
  filter(!(IID %in% dup_ids$IID))

covs <- covs %>%
  filter(!(IID %in% dup_ids$IID))

# We only want those from moba that have been genotyped\valid genotype data
genedir = "N:/durable/data/genetic/MoBaPsychGen_v1"
genefile = "MoBaPsychGen_v1-ec-eur-batch-basic-qc"
fam = read.table(paste0(genedir, "/", genefile, ".fam"), sep = " ", na.strings = "0",
                 col.names =  c("fid", "iid", "father", "mother", "sex", "phenotype"))

# Get keys from sentrix ids to moba ids
# Mothers and fathers can have multiple pregnancies, and each pregnancy can have multiple children
ce_child = read.delim ('//ess01/p471/data/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt') %>%
  select(IID, ID_2306, Role) %>% 
  filter (Role == "Child") %>%
  select(-Role)

ce_mother = read.delim ('//ess01/p471/data/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt') %>%
  select(IID, ID_2306, Role) %>% 
  filter (Role == "Mother")%>%
  select(-Role)

ce_father = read.delim ('//ess01/p471/data/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt') %>%
  select(IID, ID_2306, Role) %>% 
  filter (Role == "Father")%>%
  select(-Role)


f = inner_join(fam, ce_father, by = c("iid" = "IID")) %>%
  select(-c(father,mother))
m = inner_join(fam, ce_mother, by = c("iid" = "IID")) %>%
  select(-c(father,mother))
c = inner_join(fam, ce_child, by = c("iid" = "IID"))


length(unique(f$iid))
length(unique(f$ID_2306))
length(unique(m$iid))
length(unique(m$ID_2306))
length(unique(c$iid))
length(unique(c$ID_2306))



cm = inner_join(c, m, by = c("mother" = "iid"), suffix = c("", "_m"))

cmf = inner_join(cm, f, by = c("father" = "iid"), suffix = c("_c", "_f"))

# When mother and fathers have several children they are duplicated, 
# take the first record

#data
pheno <-  pheno %>%
  rename(fid_c = FID,
         iid = IID) 
cmf_pheno <- left_join(cmf, pheno) %>%
  mutate(height = as.numeric(height))

cmf_filt_data1 <- cmf_pheno %>%
  filter(if_any(c("height", "exam"), ~!is.na(.)))

cmf_filt_data2 = cmf_filt_data1 %>%
  group_by(ID_2306_m) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(ID_2306_f) %>%
  slice_head(n = 1) %>%
  ungroup()


covs_s <- covs %>%
  rename(fid = FID, 
         iid = IID) %>% 
  select(iid, Age_exam, Age_height) 

cmf_filt_data = cmf_filt_data2 %>%
  left_join(., covs_s) %>%
  filter((!is.na(Age_height)| !is.na(exam))) 


## Rearrange according to the .fam file
cmf_filt_data$famid = 1:nrow(cmf_filt_data)

f_arr = f %>%
  filter(iid %in% cmf_filt_data$father) %>%
  arrange(match(iid, cmf_filt_data$father)) %>%
  mutate(famid = cmf_filt_data$famid) %>%
  select(-ID_2306)
m_arr = m %>%
  filter(iid %in% cmf_filt_data$mother) %>%
  arrange(match(iid, cmf_filt_data$mother)) %>%
  mutate(famid = cmf_filt_data$famid) %>%
  select(-ID_2306)

cmf_filt_data <- cmf_filt_data %>%
  rename(ID_2306 = ID_2306_c)



c_arr = c %>%
  filter(iid %in% cmf_filt_data$iid) %>%
  arrange(match(iid, cmf_filt_data$iid)) %>%
  mutate(famid = cmf_filt_data$famid) %>%
  left_join(., select(cmf_filt_data, exam, height, ID_2306, Age_height, Age_exam, ID_2306))  

dat_ready = bind_rows(m_arr, f_arr, c_arr) %>%
  mutate(Age_exam = if_else(is.na(exam), NA_integer_, Age_exam, missing = NA_integer_))

 

write.csv(dat_ready, paste0("//ess01/p471/data/durable/projects/trio_gcta_qc/data/",genefile,"_107088", ".moba"), row.names = F)

