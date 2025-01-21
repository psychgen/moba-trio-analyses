
# 01_data_prep.R

# In this script, we read in outcome and covariate datasets created in the trio-GWAS
# projects, attach to genotools-curated PGS datasets, and add a few necessary control
# variables (parental height and edu) using phenotools, ready to perform trio
# PGS and pTDT analyses

library(phenotools)
library(genotools)
library(tidyverse)


# Using Neil/Isabella's data prep script to ensure consistency with trio-GWAS
# //ess01/P471/data/durable/projects/qc_paper/data-prep-universal/00_main-phenos-height-ea


covs <- read.table(
  "//ess01/P471/data/durable/projects/qc_paper/data-prep-universal/scratch/covariates.cov", header=T)

phenos <- read.table(
  "//ess01/P471/data/durable/projects/qc_paper/data-prep-universal/scratch/phenotypes_outcomes.pheno", header=T)

# also make a phenotools "backbone" with parental covariates

available_variables("moba") %>% variable_search("bmi")

moba <- curate_dataset(c("bmi_derived_f_far","bmi_derived_m_q1","bmi_derived_f_q1"),
                       out_format = "merged_df")

# Here we process a dataset of polygenic scores with which to perform trio-PGS and
# ptdt analyses of the two outcomes:

# For validation purposes, we get the PGS as created in PRSice and ldpred2

#PRSice

pgs_names <- c("height2",
               "bmi2",
               "ea3",
               "smoke2" )

pgs <- fetch_pgs(pgs_names,
                 pgs_software="prsice2", 
                 pgs_directory = "//ess01/P471/data/durable/common/pgs_directory/pgs/",
                 maf="0.01",
                 clump="250_1_0.1" ) 


pgs_prcsd <-pgs %>% 
  process_pgs() 

save(pgs_prcsd, file= "./data/pgs_prcsd.RData")

# Join PGS to pheno: one row per child (from one row per family member)

pheno_pgs <- pgs_prcsd %>% 
  filter(Role=="Child") %>% 
  select(preg_id,IID, BARN_NR,`height2_p<5e-08`:`smoke2_p<1_res`) %>% 
  rename_with(~ paste0("child_",.), `height2_p<5e-08`:`smoke2_p<1_res`) %>% 
  mutate(BARN_NR = as.double(BARN_NR)) %>% 
  left_join(phenos %>% as_tibble()) %>% 
  left_join(moba %>% 
              select(preg_id:birth_yr, matches("height"))) %>% 
  select(IID,FID,preg_id,BARN_NR, birth_yr,m_id,f_id,exam,height,matches("heightcm"),everything()) %>% 
  left_join(pgs_prcsd %>% 
              filter(Role=="Mother") %>% 
              select(m_id,`height2_p<5e-08`:`smoke2_p<1_res`) %>% 
              rename_with(~ paste0("mother_",.), `height2_p<5e-08`:`smoke2_p<1_res`)) %>% 
  left_join(pgs_prcsd %>% 
              filter(Role=="Father") %>% 
              select(f_id,`height2_p<5e-08`:`smoke2_p<1_res`) %>% 
              rename_with(~ paste0("father_",.), `height2_p<5e-08`:`smoke2_p<1_res`)) %>%
  left_join(covs %>%
              select(FID,IID,SEX, Age_height,Age_exam))

# Perform PGS-PCA within Role to reduce n PGS per trait down to 1

dat <-  pgs_pca(pheno_pgs,indid = "IID",pgs_var_stem = c(paste0("child_", pgs_names),
                                                         paste0("mother_", pgs_names),
                                                         paste0("father_", pgs_names)))

# Reduce down to final dataset

dat_final <- dat %>% 
  select(IID:height, "mheight"=heightcm_derived_m_q1, heightcm_derived_f_far,heightcm_derived_f_q1, SEX,Age_height,Age_exam, matches(".pgs.pc") ) %>% 
  mutate(fheight= ifelse(!is.na(heightcm_derived_f_far),heightcm_derived_f_far,heightcm_derived_f_q1 )) %>% 
  select(IID:mheight, fheight,SEX,Age_height,Age_exam,  matches(".pgs.pc") )

## Add in ldpred2 scores

ldpgs <- fetch_pgs(paste0(pgs_names,"_0.95")) 

ldpgs_prcsd <-ldpgs %>% 
  process_pgs() 

# Join PGS to pheno: one row per child (from one row per family member)

dat_final <- dat_final %>%
  left_join(ldpgs_prcsd %>% 
              filter(Role=="Child") %>% 
              select(preg_id,IID, BARN_NR,matches("_res")) %>% 
              rename_with(~ paste0("child_",.), matches("_res")) %>% 
              mutate(BARN_NR = as.double(BARN_NR))) %>% 
  left_join(ldpgs_prcsd %>% 
              filter(Role=="Mother") %>% 
              select(m_id,matches("_res")) %>% 
              rename_with(~ paste0("mother_",.), matches("_res"))) %>% 
  left_join(ldpgs_prcsd %>% 
              filter(Role=="Father") %>% 
              select(f_id,matches("_res")) %>% 
              rename_with(~ paste0("father_",.),matches("_res")))

dat_final |> select(matches("pgs.pc"),matches("0.95")) |> cor(use="pairwise.complete.obs")

#correlations between ldpred and prsice scores range from 0.77 to 0.82 as expected

# Calculate and add in parental years of education from SSB

#Code based on Neil Davies' stata code in the parental_education_fullsample project
#The first digit indicates the level of education
#This corresponds to the ISCED definitions

#   Level name - coding based on Fartein's suggestions
# 0 No education and pre-school education 1
# 1 Primary education 7
# 2 Lower secondary education 10
# 3 Upper secondary, basic 11
# 4 Upper secondary, final year 13
# 5 Post-secondary not higher education 14
# 6 First stage of higher education, undergraduate level 16
# 7 First stage of higher education, graduate level 18
# 8 Second stage of higher education (postgraduate education) 21
# 9 Unspecified


ssb = read_csv("N:/durable/data/SSB/Utdanning/W21_5323_TAB_PERSON.csv")

ssb <- ssb %>% 
  mutate(father_eduyears =case_when(str_sub(as.character(NUS2000_FAR_16),start=1,end=1)=="0" ~ 1,
                                    str_sub(as.character(NUS2000_FAR_16),start=1,end=1)=="1" ~ 7,
                                    str_sub(as.character(NUS2000_FAR_16),start=1,end=1)=="2" ~ 9,
                                    str_sub(as.character(NUS2000_FAR_16),start=1,end=1)=="3" ~ 11,
                                    str_sub(as.character(NUS2000_FAR_16),start=1,end=1)=="4" ~ 13,
                                    str_sub(as.character(NUS2000_FAR_16),start=1,end=1)=="5" ~ 14,
                                    str_sub(as.character(NUS2000_FAR_16),start=1,end=1)=="6" ~ 16,
                                    str_sub(as.character(NUS2000_FAR_16),start=1,end=1)=="7" ~ 18,
                                    str_sub(as.character(NUS2000_FAR_16),start=1,end=1)=="8" ~ 21),
         mother_eduyears =case_when(str_sub(as.character(NUS2000_MOR_16),start=1,end=1)=="0" ~ 1,
                                    str_sub(as.character(NUS2000_MOR_16),start=1,end=1)=="1" ~ 7,
                                    str_sub(as.character(NUS2000_MOR_16),start=1,end=1)=="2" ~ 9,
                                    str_sub(as.character(NUS2000_MOR_16),start=1,end=1)=="3" ~ 11,
                                    str_sub(as.character(NUS2000_MOR_16),start=1,end=1)=="4" ~ 13,
                                    str_sub(as.character(NUS2000_MOR_16),start=1,end=1)=="5" ~ 14,
                                    str_sub(as.character(NUS2000_MOR_16),start=1,end=1)=="6" ~ 16,
                                    str_sub(as.character(NUS2000_MOR_16),start=1,end=1)=="7" ~ 18,
                                    str_sub(as.character(NUS2000_MOR_16),start=1,end=1)=="8" ~ 21)) %>% 
  select(PID,mother_eduyears,father_eduyears)

ssb_link <- haven::read_sav("N:/durable/data/Linkage_files/SSB_link/PDB2306_kobling_SSB&KUHR_Combined_20220112.sav") 

ssb = ssb %>% 
  left_join(ssb_link %>% 
              select("preg_id"=PREG_ID_2306,BARN_NR,"PID"=PID_NPR_2306))

dat_final= dat_final %>% 
  left_join(ssb %>% 
              mutate(preg_id=as.character(preg_id)) %>% 
              select(preg_id,BARN_NR, matches("eduyears")))


save(dat_final,pgs_names, file= "./data/dat_final.RData")
