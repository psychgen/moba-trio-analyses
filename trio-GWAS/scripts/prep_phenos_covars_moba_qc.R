# Phenotypes and covariates for MoBa QC paper 2024 
# Isabella Badini
# updated august 25 

rm(list=ls())

library(haven)
library(readr)
library(stringr)
library(dplyr)
library(data.table)
library(phenotools)

dir(path="N:/durable/projects/trio_gwas")
#File path for the phenotypic data
moba_data="N:/durable/data/"

# Input files used
sv_info_file <- paste(moba_data, "MoBaPhenoData/PDB2306_MoBa_v12/SPSS/PDB2306_SV_INFO_V12.sav",sep="")
q7_file <- paste(moba_data, "MoBaPhenoData/PDB2306_MoBa_v12/SPSS/PDB2306_Q7yrs_v12.sav",sep="") 
q8_file <- paste(moba_data, "MoBaPhenoData/PDB2306_MoBa_v12/SPSS/PDB2306_Q8yrs_v12.sav",sep="") 
utd05_file <- paste(moba_data,"SSB/Utdanning/W21_5323_NASJONALE_PROVER.csv",sep="") 
stdexam_file <- "N:/durable/projects/height_bmi_ea/scratch/exam_dta_for_merge.dta"  # previously cleaned education data by ND (scripts/n/001_exam_data.do) and PR (scripts/exam_score_merge_prep_poppy.do)
genLink_c_file <- paste(moba_data,"Linkage_files/Genetics_link/2022_01_06_MoBaGeneticsTot_Child_PDB2306.sav",sep="")
ssb_link_file <- paste(moba_data,"Linkage_files/SSB_link/PDB2306_kobling_SSB&KUHR_Combined_20220112.sav",sep="") 
covar_file <- paste(moba_data,"genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt",sep="") 

# Output files used
covariates_file <- "N:/durable/projects/qc_paper/trio-GWAS/scratch/covariates_mqc.cov"
phenotypes_file <- "N:/durable/projects/qc_paper/trio-GWAS/scratch/phenotypes_outcomes_mqc.pheno"

# Function for outliers
findouts <- function(data, upper_sd = 4, lower_sd = 4) {
  mean_data <- mean(data, na.rm = TRUE)
  sd_data <- sd(data, na.rm = TRUE)
  upper <- upper_sd * sd_data + mean_data 
  lower <- mean_data - lower_sd * sd_data
  replace(data, data > upper | data < lower, NA)
}

# Read in SV INFO file 
sv_info <- read_sav(sv_info_file) # updated list of consenting MoBa participants; check if need to remove participants from cov file 

#---------------- Prep covariates file -----------------------#

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
  mutate(SEX = ifelse(SEX == 0, NA, SEX)) #%>%
  #mutate(SEX = ifelse(SEX == 2, 0, SEX)) # recode sex so that 1 = M and 0 = F

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
# sleep duration at 7 yrs
# smfq at 8 yrs
# 5th grade reading, english, maths (10-11yrs; NB: National tests are in the fall (sept/oct) of 5th grade so most children will be 11 - see which of the tests has the least missing on assessment year variable 
# !might be worth using age at 7/8yrs to fill in age at 8/7yrs respectively (have not done this!)

# Height 7 years - age at measurement: 
q7 <- read_sav(q7_file) %>% # Read in Q7 MoBA questionnaire 
  select(PREG_ID_2306, BARN_NR, AGE_RETURN_MTHS_Q7) # age in months when questionnaire was returned 

summary(q7)

q7 <- q7 %>% 
  mutate(Age_height = round((AGE_RETURN_MTHS_Q7/12), digits=2))  

summary(q7$Age_height)

#cov <- cov %>%
#  left_join(q7, by=c("PREG_ID_2306", "BARN_NR")) %>%
#  select(-AGE_RETURN_MTHS_Q7)

# Sleep duration 7 years - age at measurement:
q7 <- q7 %>% 
  mutate(Age_sleep = round((AGE_RETURN_MTHS_Q7/12), digits=2))  

summary(q7$Age_sleep)

cov <- cov %>%
  left_join(q7, by=c("PREG_ID_2306", "BARN_NR")) %>%
  select(-AGE_RETURN_MTHS_Q7)

# Depressive symptoms (smfq) 8 years - age at measurement:
q8 <- read_sav(q8_file) %>% # Read in Q7 MoBA questionnaire 
  select(PREG_ID_2306, BARN_NR, AGE_RETURN_MTHS_Q8AAR) # age in months when questionnaire was returned 
summary(q8)

q8 <- q8 %>% 
  mutate(Age_dep = round((AGE_RETURN_MTHS_Q8AAR/12), digits=2))  
summary(q8$Age_dep)

cov <- cov %>%
  left_join(q8, by=c("PREG_ID_2306", "BARN_NR")) %>%
  select(-AGE_RETURN_MTHS_Q8AAR)

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
  left_join(ssb_link, by=c("PID" = "PID_NPR_2306")) %>%
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

# Remove outliers from age covariates
cov <- cov %>%
  mutate(Age_height = findouts(Age_height, upper_sd = 4, lower_sd = 4)) %>% #335 outliers removed
  mutate(Age_sleep = findouts(Age_sleep, upper_sd = 4, lower_sd = 4)) %>% #335
  mutate(Age_dep = findouts(Age_dep, upper_sd = 4, lower_sd = 4)) %>% #221
  mutate(Age_exam = findouts(Age_exam, upper_sd = 4, lower_sd = 4)) #390 outliers removed


# Finalise covariates file

head(cov)
cov <- cov %>%
  select(FID, IID, SEX, YOB, Age_height, Age_exam, Age_sleep, Age_dep, PC1:PC20, starts_with("batch"))

################### Prepare phenotypes #########################################

#---------------- Height -----------------------#

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
  mutate(height = findouts(height, upper_sd = 4, lower_sd = 4)) 

#---------------- Exam scores -----------------------#

stdexam <- read_dta(stdexam_file) %>%
  select(PREG_ID_2306, BARN_NR, zNPENG05, zNPLES05, zNPREG05) # 5th grade exam scores
head(stdexam)

stdexam <- stdexam %>%
  mutate(exam = rowMeans(across(zNPENG05:zNPREG05), na.rm = TRUE))

stdexam <- stdexam %>%
  mutate(exam = findouts(exam, upper_sd = 4, lower_sd = 4)) 

stdexam <- stdexam %>%
  mutate(exam_norm = scale(exam, center = TRUE, scale = TRUE)[,1])

# Merge with height data

df <- df %>%
  mutate(PREG_ID_2306 = as.numeric(PREG_ID_2306)) %>%
  left_join(stdexam, by=c("PREG_ID_2306", "BARN_NR")) %>%
  filter(PREG_ID_2306 %in% sv_info$PREG_ID_2306)

head(df)
df<-as.data.frame(df)
df[is.na(df)] <- NA

#---------------- Sleep duration -----------------------#
tmp <- curate_dataset(variables_required = c("JJ418", "NN387"), return_items = T, out_format = "merged_df")
colnames(tmp) <- gsub("_raw$", "", colnames(tmp)) 

#describe(tmp$LL507) #5yrs n=26381
describe(tmp$JJ418) #7yrs n=52032
describe(tmp$NN387) #8yrs n=35047

# Check correlations between items at different waves
#corvars <- tmp[,c("LL507", "JJ418", "NN387")]
#corm <- cor(corvars, use = "complete.obs")
#print(corm)
# correlations not particularly high; don't average across 
# tmp$avg_sleep_duration <- rowMeans(tmp[, c("sleep_duration_5yr", "sleep_duration_7yr", "sleep_duration_8yr")], na.rm = TRUE)

head(tmp$JJ418) # go with 7yr item since it has the largest n 
table(tmp$JJ418, useNA = "ifany")

head(tmp$NN387)
table(tmp$NN387, useNA = "ifany") # use available data at age 8 to fill in missing data at 7yrs (r = 0.41) 

# Convert likert scale into hours
tmp <- tmp %>%  
  mutate(JJ418_hrs = case_when(
    JJ418 == 0 ~ NA_real_,
    JJ418 == 1 ~ 8,
    JJ418 == 2 ~ 9,
    JJ418 == 3 ~ 10,
    JJ418 == 4 ~ 11,
    JJ418 == 5 ~ 12,
    JJ418 == 6 ~ 8.5,
    JJ418 == 7 ~ 9.5,
    JJ418 == 8 ~ 10.5,
    JJ418 == 9 ~ 11.5,
    TRUE ~ NA_real_
  ),
  NN387_hours = case_when(
    NN387 == 0 ~ NA_real_,  # More than 1 checkbox filled in (treat as NA)
    NN387 == 1 ~ 8,         # 8 hrs or less
    NN387 == 2 ~ 9,         # 9hrs
    NN387 == 3 ~ 10,        # 10hrs
    NN387 == 4 ~ 11,        # 11hrs
    NN387 == 5 ~ 12,        # 12hrs+
    NN387 == 6 ~ 8.5,       # 8hrs + 9hrs
    NN387 == 7 ~ 9.5,       # 9hrs + 10hrs
    NN387 == 8 ~ 10.5,      # 10hrs + 11hrs
    NN387 == 9 ~ 11.5,      # 11hrs + 12hrs or more
    TRUE ~ NA_real_         
  ))

# Fill in missing data for sleep duration at 7 years using data from 8 years
tmp <- tmp %>%
  mutate(sleep_hrs = coalesce(JJ418_hrs, NN387_hours))

#describe(tmp$sleep_hrs) #52028 obs; 61847 NAs (54.3% missing) - 7yrs only
describe(tmp$sleep_hrs) #56717 obs; 56992 NAs (50.1% missing)

tmp <- tmp %>%
  mutate(sleep_hrs = findouts(sleep_hrs, upper_sd = 4, lower_sd = 4)) # no outliers

tmp <- tmp %>%
  mutate(sleep_hrs = scale(sleep_hrs, center = TRUE, scale = TRUE)[,1])

# merge with df 
tmp <- tmp %>% 
  mutate(preg_id = as.numeric(preg_id)) %>%
  rename(PREG_ID_2306 = preg_id) %>%
  select(PREG_ID_2306, BARN_NR, sleep_hrs)

df <- df %>%
  left_join(tmp, by=c("PREG_ID_2306", "BARN_NR")) %>%
  filter(PREG_ID_2306 %in% sv_info$PREG_ID_2306)

#---------------- Depression -----------------------#
tmp <- curate_dataset(variables_required = c("smfq_dep_c_8yr"), return_items = F, out_format = "merged_df")

describe(tmp$smfq_dep_c_8yr) #n=43164

tmp <- tmp %>%
  mutate(smfq_dep = findouts(smfq_dep_c_8yr, upper_sd = 4, lower_sd = 4)) %>% #358 outliers
  mutate(smfq_dep = scale(smfq_dep, center = TRUE, scale = TRUE)[,1]) 

# merge with df 
tmp <- tmp %>% 
  mutate(preg_id = as.numeric(preg_id)) %>%
  rename(PREG_ID_2306 = preg_id) %>%
  select(PREG_ID_2306, BARN_NR, smfq_dep)

df <- df %>%
  left_join(tmp, by=c("PREG_ID_2306", "BARN_NR")) %>%
  filter(PREG_ID_2306 %in% sv_info$PREG_ID_2306)

#---------------- Final preperations -----------------------#
names(df)

df <- df %>%
  select(PREG_ID_2306, BARN_NR, height, exam, smfq_dep, sleep_hrs)

describe(df)

# Link FID and IID 
df <- df %>%
  left_join(gen_link, by=c("PREG_ID_2306","BARN_NR")) %>% #link with sentrix ID
  rename(IID = "SENTRIXID")
   
cov_tmp <- cov %>%
  select(FID, IID) 

df <- cov_tmp %>%
  left_join(df, by="IID") %>%
  select(FID, IID, height, exam, smfq_dep, sleep_hrs)

head(df)
sum(duplicated(df$IID)) # 73 dups
df_nodups <- df %>% 
  filter(!duplicated(IID))

head(cov)
sum(duplicated(cov$IID)) # 73 dups
cov_nodups <- cov %>% 
  filter(!duplicated(IID)) %>%
  filter(IID %in% df_nodups$IID)

# Save 
write.table(cov_nodups, covariates_file, col.names=T, row.names=F, quote=F, sep='\t')
write.table(df_nodups, phenotypes_file, col.names=T, row.names=F, quote=F, sep='\t')

# I need separate cov files:
cov_noyob <- cov_nodups %>% 
  select(-YOB)
write.table(cov_noyob, "N:/durable/projects/qc_paper/trio-GWAS/scratch/covariates_noyob_mqc.cov", col.names = T, row.names = F, sep='\t', quote = F)

cov_height <- cov_nodups %>%
  select(-Age_exam, -Age_sleep, -Age_dep, -YOB)
write.table(cov_height, "N:/durable/projects/qc_paper/trio-GWAS/scratch/covariates_height_mqc.cov", col.names = T, row.names = F, sep='\t', quote = F)

cov_exam <- cov_nodups %>%
  select(-Age_height, -Age_sleep, -Age_dep, -Age_exam)
write.table(cov_exam, "N:/durable/projects/qc_paper/trio-GWAS/scratch/covariates_exam_mqc.cov", col.names = T, row.names = F, sep='\t', quote = F)

cov_exam_noyob <- cov_nodups %>%
  select(-Age_height, -Age_sleep, -Age_dep, -Age_exam, -YOB)
write.table(cov_exam_noyob, "N:/durable/projects/qc_paper/trio-GWAS/scratch/covariates_exam_noyob_mqc.cov", col.names = T, row.names = F, sep='\t', quote = F)

cov_sleep <- cov_nodups %>%
  select(-Age_height, -Age_exam, -Age_dep, -YOB)
write.table(cov_sleep, "N:/durable/projects/qc_paper/trio-GWAS/scratch/covariates_sleep_mqc.cov", col.names = T, row.names = F, sep='\t',quote = F)

cov_dep <- cov_nodups %>%
  select(-Age_height, -Age_sleep, -Age_exam, -YOB)
write.table(cov_dep, "N:/durable/projects/qc_paper/trio-GWAS/scratch/covariates_smfq_mqc.cov", col.names = T, row.names = F, sep='\t', quote = F)

