#author: Espen Eilertsen, Laura Hegemann 
#prep phenotypic data / select sample for trio-GCTA

library(tidyverse)
library(haven)

# data

pheno <- read.table("//ess01/p471/data/durable/projects/qc_paper/data-prep-universal/scratch/phenotypes_outcomes_mqc.pheno", header = TRUE)
covs <- read.table("//ess01/p471/data/durable/projects/qc_paper/data-prep-universal/scratch/covariates_noyob_mqc.cov", header = TRUE)
covs_yob <- read.table("//ess01/p471/data/durable/projects/qc_paper/data-prep-universal/scratch/covariates_yob_mqc.cov", header = TRUE)


dup_ids <- read.table("//ess01/p471/data/durable/projects/qc_paper/data-prep-universal/scratch/duplicated_ids.txt", header = TRUE)

pheno_fid <- pheno$FID

#check withdrawls removed. 
EC_lk <- read_table("//ess01/p471/data/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt") %>%
  select(ID_2306, FID,IID, Role) %>%
  filter(Role == "Child") %>%
  separate(ID_2306,c("PREG_ID_2306", "remove"), "_") %>%
  select(!remove)

sv_info <- read_spss("N:/durable/data/MoBaPhenoData/PDB2306_MoBa_v12/SPSS/PDB2306_SV_INFO_V12.sav") %>%
  mutate(PREG_ID_2306 = as.character(PREG_ID_2306)) %>%
  right_join(EC_lk) %>%
  right_join(covs)

sum(is.na(sv_info$PREG_ID_2306))
sum(pheno$PREG_ID_2306 %in% sv_info$PREG_ID_2306 == FALSE)

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


covs_s <- covs %>%
  rename(fid = FID, 
         iid = IID) %>% 
  select(iid, Age_exam, Age_height, Age_sleep, Age_dep) 

cmf_filt_data1 =cmf_pheno %>%
  left_join(., covs_s) %>%
  mutate(exam = if_else(is.na(Age_exam) & !is.na(exam), NA_integer_, exam, missing = NA_integer_),
         height = if_else(is.na(Age_height), NA_integer_, height, missing = NA_integer_), 
         smfq_dep = if_else(is.na(Age_dep), NA_integer_, smfq_dep, missing = NA_integer_), 
         sleep_hrs = if_else(is.na(Age_sleep), NA_integer_, sleep_hrs, missing = NA_integer_),
         Age_exam = if_else(is.na(exam), NA_integer_, Age_exam, missing = NA_integer_),
         Age_height = if_else(is.na(height), NA_integer_, Age_height, missing = NA_integer_), 
         Age_dep = if_else(is.na(smfq_dep), NA_integer_, Age_dep, missing = NA_integer_), 
         Age_sleep = if_else(is.na(sleep_hrs), NA_integer_, Age_sleep, missing = NA_integer_)) %>%
  filter(!(is.na(exam) & is.na(height) & is.na(smfq_dep) & is.na(sleep_hrs))) 

cmf_filt_data = cmf_filt_data1 %>%
  group_by(ID_2306_m) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(ID_2306_f) %>%
  slice_head(n = 1) %>%
  ungroup()







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
  left_join(., cmf_filt_data) %>%
  select(fid, iid, sex, phenotype, famid,exam, height, smfq_dep, sleep_hrs, Age_height, Age_exam, Age_dep, Age_sleep)  

dat_ready = bind_rows(m_arr, f_arr, c_arr)


dat_ready[c(1,1+nrow(c_arr),1+(2*nrow(c_arr))),]
dat_ready[c(5000,5000+nrow(c_arr),5000+(2*nrow(c_arr))),]




sum(!is.na(dat_ready$exam))
sum(!is.na(dat_ready$Age_exam))
sum(!is.na(dat_ready$height))
sum(!is.na(dat_ready$Age_height))
sum(!is.na(dat_ready$sleep_hrs))
sum(!is.na(dat_ready$Age_sleep))
sum(!is.na(dat_ready$smfq_dep))
sum(!is.na(dat_ready$Age_dep))

sum(is.na(dat_ready$sex))


write.csv(dat_ready, paste0("//ess01/p471/data/durable/projects/trio_gcta_qc/data/",genefile,"_105210", ".moba"), row.names = F)


#test subset for gcta

test_sub <- cmf_filt_data %>%
  slice_sample(n =5000)

c_arr_test = c %>%
  filter(iid %in% test_sub$iid) %>%
  arrange(match(iid, test_sub$iid)) %>%
  mutate(famid = test_sub$famid) %>%
  left_join(., test_sub) %>%
  select(fid, iid, sex, phenotype, famid,exam, height, smfq_dep, sleep_hrs, Age_height, Age_exam, Age_dep, Age_sleep)  

f_arr_test = f %>%
  filter(iid %in% test_sub$father) %>%
  arrange(match(iid, test_sub$father)) %>%
  mutate(famid = test_sub$famid) %>%
  select(-ID_2306)
m_arr_test = m %>%
  filter(iid %in% test_sub$mother) %>%
  arrange(match(iid, test_sub$mother)) %>%
  mutate(famid = test_sub$famid) %>%
  select(-ID_2306)

dat_ready_test = bind_rows(m_arr_test, f_arr_test, c_arr_test)


dat_ready_test[c(1,1+nrow(c_arr_test),1+(2*nrow(c_arr_test))),]
dat_ready_test[c(5000,5000+nrow(c_arr_test),5000+(2*nrow(c_arr_test))),]

ids_keep <- dat_ready_test %>%
  select(fid,iid)

write.table(ids_keep, "//ess01/p471/data/durable/people/Laura_H/Cluster_backup/trio_gcta_qc/test_grm_keepiid.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
