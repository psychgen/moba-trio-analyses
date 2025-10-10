
# 01_data_prep.R

# In this script, we read in outcome and covariate datasets created in the trio-GWAS
# projects, attach to genotools-curated PGS datasets, and add a few necessary control
# variables (parental height and edu) using phenotools, ready to perform trio
# PGS and pTDT analyses

library(phenotools)
library(genotools)
library(tidyverse)


# Using Neil/Isabella's data prep script to ensure consistency with trio-GWAS
# //ess01/P471/data/durable/projects/trio_gwas/scripts/00_main-phenos-height-ea


covs <- read.table(
  "//ess01/P471/data/durable/projects/qc_paper/data-prep-universal/scratch/covariates_mqc.cov", header=T)

phenos <- read.table(
  "//ess01/P471/data/durable/projects/qc_paper/data-prep-universal/scratch/phenotypes_outcomes_mqc.pheno", header=T)

# also make a phenotools "backbone" with parental covariates

available_variables("moba") %>% variable_search("bmi")

moba <- curate_dataset(c("bmi_derived_f_far","bmi_derived_m_q1","bmi_derived_f_q1",
                         "scl_dep_f_far","scl_dep_m_q3",
                         "UM311", #How long awake at night (mothers, age 14 Q)
                         "G_58_2"), #How often awake at night (fathers, 2015 Q)
                       out_format = "merged_df")

# Here we process a dataset of polygenic scores with which to perform trio-PGS and
# ptdt analyses of the two outcomes:

pgs_names <- c("height2018",
               "bmi2018",
               "EA2018",
               "smok2019",
               "mdd2025")


## Add in ldpred2 scores

ldpgs <- fetch_pgs(pgs_names, within_fam=T) 

ldpgs_prcsd <-ldpgs %>% 
  process_pgs() 


save(ldpgs_prcsd,file="./data/pgs_prcsd.RData")

# Join PGS to pheno,covs,moba: one row per child (from one row per family member)

dat_final <- phenos %>%
  left_join(ldpgs_prcsd %>% 
              filter(Role=="Child") %>% 
              select(preg_id,IID, BARN_NR,matches("_res")) %>% 
              rename_with(~ paste0("child_",.), matches("_res")) %>% 
              mutate(BARN_NR = as.double(BARN_NR))) %>% 
  left_join(moba %>%
              select(preg_id:birth_yr, matches("height|scl|UM|G_"))) %>%
  left_join(ldpgs_prcsd %>% 
              filter(Role=="Mother") %>% 
              select(m_id,matches("_res")) %>% 
              rename_with(~ paste0("mother_",.), matches("_res"))) %>% 
  left_join(ldpgs_prcsd %>% 
              filter(Role=="Father") %>% 
              select(f_id,matches("_res")) %>% 
              rename_with(~ paste0("father_",.),matches("_res"))) %>%
  left_join(covs %>%
              select(IID,SEX,YOB,Age_height,Age_exam,Age_sleep,Age_dep)) %>%
    mutate(fheight= ifelse(!is.na(heightcm_derived_f_far),heightcm_derived_f_far,heightcm_derived_f_q1 )) %>% 
  select(FID,IID,preg_id, BARN_NR, m_id,f_id,birth_yr,YOB,SEX,'mheight'=heightcm_derived_m_q1, fheight,matches("Age_"), everything() )
  
  

# This section involves processing of parental phenotypes used only in the PTDT analyses
# to exclude children for whom the expectation of over/undertransmission based on phenotype
# selection does not hold (i.e., those whose parents are also similarly extreme for the phenotype-
# or in this case, often a related phenotype)

# Since these phenotypes are only used in PGS analyses, they are not in the data-prep-universal pipeline

##PARENTAL EDUCATION

# Calculate and add in parental years of education from SSB

#Code based on Neil Davies' stata code in the parental_education_fullsample project
#The first digit indicates the level of education
#This corresponds to the ISCED definitions

#   Level name - years of education
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


ssb = read_csv("//ess01/P471/data/durable/data/SSB/Utdanning/W21_5323_TAB_PERSON.csv")

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

ssb_link <- haven::read_sav("//ess01/P471/data/durable/data/Linkage_files/SSB_link/PDB2306_kobling_SSB&KUHR_Combined_20220112.sav") 

ssb = ssb %>% 
  left_join(ssb_link %>% 
              select("preg_id"=PREG_ID_2306,BARN_NR,"PID"=PID_NPR_2306))

dat_final= dat_final %>% 
  left_join(ssb %>% 
              mutate(preg_id=as.character(preg_id)) %>% 
              select(preg_id,BARN_NR, matches("eduyears")))





save(dat_final,pgs_names, file= "./data/dat_final.RData")
