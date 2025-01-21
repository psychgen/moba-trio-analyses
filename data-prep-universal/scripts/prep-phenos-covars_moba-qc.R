# Phenotypes and covariates for MoBa QC paper 2024 
# IB 

rm(list=ls())

library(phenotools)
library(data.table)
library(dplyr)
library(haven)
library(readr)
library(stringr)

dir(path="N:/durable/projects/trio_gwas")
#File path for the phenotypic data
moba_data="N:/durable/data/"

# Input files used
sv_info_file <- paste(moba_data, "MoBaPhenoData/PDB2306_MoBa_v12/SPSS/PDB2306_SV_INFO_V12.sav",sep="")
q7_file <- paste(moba_data, "MoBaPhenoData/PDB2306_MoBa_v12/SPSS/PDB2306_Q7yrs_v12.sav",sep="") 
utd05_file <- paste(moba_data,"SSB/Utdanning/W21_5323_NASJONALE_PROVER.csv",sep="") 
stdexam_file <- "N:/durable/projects/height_bmi_ea/scratch/exam_dta_for_merge.dta"  # previously cleaned education data by ND (scripts/n/001_exam_data.do) and PR (scripts/exam_score_merge_prep_poppy.do)
genLink_c_file <- paste(moba_data,"Linkage_files/Genetics_link/2022_01_06_MoBaGeneticsTot_Child_PDB2306.sav",sep="")
ssb_link_file <- paste(moba_data,"Linkage_files/SSB_link/PDB2306_kobling_SSB&KUHR_Combined_20220112.sav",sep="") 
covar_file <- paste(moba_data,"genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt",sep="") 

# Output files used
covariates_file <- "N:/durable/projects/qc_paper/trio-GWAS/scratch/covariates.cov"
phenotypes_file <- "N:/durable/projects/qc_paper/trio-GWAS/scratch/phenotypes_outcomes.pheno"

# Function for outliers
findouts <- function(data, upper_sd = 5, lower_sd = 5) {
  mean_data <- mean(data, na.rm = TRUE)
  sd_data <- sd(data, na.rm = TRUE)
  upper <- upper_sd * sd_data + mean_data 
  lower <- mean_data - lower_sd * sd_data
  replace(data, data > upper | data < lower, NA)
}

# Read in SV INFO file 
sv_info <- read_sav(sv_info_file) # updated list of consenting MoBa participants; check if need to remove participants from cov file 

## Prepare covariates ## 

# Read in covariates file
cov <- fread(covar_file)

cov <- cov %>%
  filter(Role == "Child") %>%
  select(SENTRIXID, FID, IID, SEX, YOB, PC1:PC20, starts_with("batch"))

summary(cov)

cov <- cov %>%
  select(-c(batch_norment_jan2015, batch_norment_jun2015, batch_norment_mar2021_1532)) # remove batches with no variance + one batch to avoid over-fitting 

table(cov$SEX, useNA = "ifany")
cov <- cov %>%
  mutate(SEX = ifelse(SEX == 0, NA, SEX)) %>%
  mutate(SEX = ifelse(SEX == 2, 0, SEX)) # recode sex so that 1 = M and 0 = F

# Read in genetics-moba linkage file 
gen_link <- read_sav(genLink_c_file)

gen_link <- gen_link %>%
  select(1:3) %>%
  rename(SENTRIXID = "SENTRIX_ID")

# Join covariates file with gen_link file 
cov <- cov %>%
  left_join(gen_link, by = "SENTRIXID") %>%
  filter(PREG_ID_2306 %in% sv_info$PREG_ID_2306) %>%
  select(-SENTRIXID) %>%
  relocate(c(PREG_ID_2306, BARN_NR), .before = FID)

head(cov)

# Offspring age at measurement covariates  
# Phenotypes: 
# Height at 7 yrs 
# 5th grade reading (10-11yrs; NB: National tests are in the fall (sept/oct) of 5th grade so most children will be 11 - see which of the tests has the least missing on assessment year variable 
# 5th grade english (10-11yrs)
# 5th grade maths (10-11yrs)

# Height 7 years - age at measurement: 

q7 <- read_sav(q7_file) %>% # Read in Q7 MoBA questionnaire 
  select(PREG_ID_2306, BARN_NR, AGE_RETURN_MTHS_Q7) # age in months when questionnaire was returned 

summary(q7)

q7 <- q7 %>% 
  mutate(Age_height = round((AGE_RETURN_MTHS_Q7/12), digits=2)) # do I need to do this? 

summary(q7$Age_height) # lots of missing from non-response at q7

cov <- cov %>%
  left_join(q7, by=c("PREG_ID_2306", "BARN_NR")) %>%
  select(-AGE_RETURN_MTHS_Q7)

# 5th grade exams - age at measurement:

# Read in moba-ssb linkage file 
ssb_link <- read_sav(ssb_link_file) %>% # Linkage file for LNR to PREG_ID
  filter(SUP_Type == "child") 

# National exam data; filter for 5th grade results
utd05 <- read_csv(utd05_file) %>% 
  filter(str_detect(PROVE, "05"))  

# Link moba IDs 
ssb_link <- ssb_link %>%
  select(PREG_ID_2306, BARN_NR, PID_NPR_2306)

utd05 <- utd05 %>%
  #left_join(ssb_link, by=c("PID" = "PID_NPR_2306")) %>%
  select(PREG_ID_2306, BARN_NR, PROVE, AARGANG) # select needed variables from exam data

cov_tmp <- cov %>%
  select(PREG_ID_2306, BARN_NR, YOB) # we need YOB variable from covariates file to calculate age at measurement - create temp reduced df

npLES <- utd05 %>%
  filter(PROVE == "NPLES05") %>%
  left_join(cov_tmp, by=c("PREG_ID_2306", "BARN_NR")) %>%
  select(PREG_ID_2306, BARN_NR, YOB, AARGANG)
head(npLES)

npLES$year <- substr(npLES$AARGANG, 1, 4) # Remove AARGANG month (month that test was taken) as no month of birth variable for MoBa children anyway  
npLES$year <- as.numeric(npLES$year)

npLES <- npLES %>%
  mutate(Age_NPLES05 = year - YOB) %>% 
  select(PREG_ID_2306, BARN_NR, Age_NPLES05) 

summary(npLES$Age_NPLES05)
sum(!is.na(npLES$Age_NPLES05))

# 5th grade maths:

npREG <- utd05 %>%
  filter(PROVE == "NPREG05") %>%
  left_join(cov_tmp, by=c("PREG_ID_2306", "BARN_NR")) %>%
  select(PREG_ID_2306, BARN_NR, YOB, AARGANG)
head(npREG)

npREG$year <- substr(npREG$AARGANG, 1, 4)
npREG$year <- as.numeric(npREG$year)

npREG <- npREG %>%
  mutate(Age_NPREG05 = year - YOB) %>% 
  select(PREG_ID_2306, BARN_NR, Age_NPREG05) 

summary(npREG$Age_NPREG05)
sum(!is.na(npREG$Age_NPREG05))

#5th grade english: 

npENG <- utd05 %>%
  filter(PROVE == "NPENG05") %>%
  left_join(cov_tmp, by=c("PREG_ID_2306", "BARN_NR")) %>%
  select(PREG_ID_2306, BARN_NR, YOB, AARGANG)
head(npENG)

npENG$year <- substr(npENG$AARGANG, 1, 4)
npENG$year <- as.numeric(npENG$year)

npENG <- npENG %>%
  mutate(Age_NPENG05 = year - YOB) %>% 
  select(PREG_ID_2306, BARN_NR, Age_NPENG05) 

summary(npENG$Age_NPENG05)
sum(!is.na(npENG$Age_NPENG05))

# Join datasets together 
AGE_NP05 <- npREG %>%
  merge(npLES, by=c("PREG_ID_2306","BARN_NR")) %>%
  merge(npENG, by=c("PREG_ID_2306","BARN_NR"))

# Actually we don't need separate ages for the tests as a) children were the same age when taking these tests and b) will create average score on all 3 tests
summary(AGE_NP05) # Age_NPREG has least missing - use other age variables to fill in missing

AGE_NP05 <- AGE_NP05 %>%
  mutate(Age_exam = coalesce(Age_NPREG05, Age_NPLES05, Age_NPENG05)) %>%
  select(PREG_ID_2306, BARN_NR, Age_exam)

summary(AGE_NP05)

AGE_NP05 <- AGE_NP05 %>%
  filter(!is.na(PREG_ID_2306)) %>%
  filter(!is.na(BARN_NR))
  
table(AGE_NP05$Age_exam, useNA = "ifany")

cov <- cov %>%
  left_join(AGE_NP05, by=c("PREG_ID_2306", "BARN_NR"))

# Finalise covariates file

head(cov)
cov <- cov %>%
  select(FID, IID, SEX, Age_height, Age_exam, PC1:PC20, starts_with("batch"))

## Prepare phenotypes ##

# Height

df <- curate_dataset(variables_required = c("bmi_derived_c_7yr","bmi_derived_c_8yr"), return_items = T, out_format = "merged_df")

df <- df %>%
  rename(PREG_ID_2306 = preg_id) %>%
  select(PREG_ID_2306, BARN_NR, heightcm_derived_c_7yr, heightcm_derived_c_8yr)

colSums(!is.na(df[c('heightcm_derived_c_7yr', 'heightcm_derived_c_8yr')])) # larger N at 7yr 

# Combine height at ages 7 and 8 years (use non-missing value between 7yrs and 8yrs)
df <- df %>%
  mutate(height = coalesce(heightcm_derived_c_7yr, heightcm_derived_c_8yr))

summary(df$height) #n=55294

df <- df %>%
  mutate(height_trim = findouts(height, upper_sd = 5, lower_sd = 5)) # no outliers +- 5SD from mean 

# Prepare data on exam scores

stdexam <- read_dta(stdexam_file) %>%
  select(PREG_ID_2306, BARN_NR, zNPENG05, zNPLES05, zNPREG05) # 5th grade exam scores
head(stdexam)

stdexam <- stdexam %>%
  mutate(exam = rowMeans(across(zNPENG05:zNPREG05), na.rm = TRUE))

# Merge with height data

df <- df %>%
  mutate(PREG_ID_2306 = as.numeric(PREG_ID_2306)) %>%
  left_join(stdexam, by=c("PREG_ID_2306", "BARN_NR")) %>%
  filter(PREG_ID_2306 %in% sv_info$PREG_ID_2306)

head(df)
df<-as.data.frame(df)
df[is.na(df)] <- NA

# Link FID and IID 

df <- df %>%
  left_join(gen_link, by=c("PREG_ID_2306","BARN_NR")) %>% #link with sentrix ID
  rename(IID = "SENTRIXID")
   
cov_tmp <- cov %>%
  select(FID, IID) 

df <- cov_tmp %>%
  left_join(df, by="IID") %>%
  select(FID, IID, height, exam)

# Save 
write.table(cov, covariates_file, col.names=T, row.names=F, quote=F, sep='\t')
write.table(df, phenotypes_file, col.names=T, row.names=F, quote=F, sep='\t')
