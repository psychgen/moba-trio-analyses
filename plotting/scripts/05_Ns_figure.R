
library(tidyverse)
library(phenotools)
library(ggthemes)
library(RColorBrewer)
library(viridis)

# In this script we attemp to visualise the available data for genotyped family
# members in MoBa, by approximating number of responses per genotyped family member,
# per topic area over time

# First, we establish which kids, mothers, and father have data at each timepoint:

available_variables() |>  variable_search("scl")

base = curate_dataset(c("scl_short_m_q1","scl_full_m_q3",
                        "scl_full_m_6m", "scl_full_m_18m","scl_full_m_3yr", "scl_full_m_5yr",
                        "scl_full_m_8yr", "scl_full_m_14m", "bmi_derived_c_14c",
                        "scl_full_f_far","scl_full_f_far2"),
                      out_format="merged_df")

# Ascertain parental sibs to enable counting of aunt/uncle/cousin responses

#Read in a file with genetic relationships to ID sibs of ABC parents

kin <- readr::read_tsv ("//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-rel.kin")

kin = kin %>%
  filter(InfType%in%c("FS","Dup/MZ")) # filter down to probably sibships

#Read in the covariate file to get IIDs for MoBa adults

cov <- readr::read_tsv ("//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt")

cov = cov %>%
  filter(Role %in% c("Mother","Father")) |>
  select(ID_2306, SENTRIXID)

# Join to get sibIDs

sibs = cov |>
  left_join(kin |>
              select("SENTRIXID" = ID1,
                     ID2)) |>
  left_join(kin |>
              select("SENTRIXID" = ID2,
                     ID1)) |>
  pivot_longer(cols = c(-ID_2306,-SENTRIXID), values_to="au_iid") |>
  drop_na(au_iid) |>
  left_join(cov |>
              select("au_iid"=SENTRIXID, "au_mobaid"=ID_2306))

# Join to dataset with PREG_IDs to make child level

au = base |>
  select(preg_id,BARN_NR,m_id,f_id) |>
  pivot_longer(cols = c(-preg_id,-BARN_NR), values_to="p_id") |>
  left_join(sibs |>
              select("p_id" = ID_2306, au_iid, au_mobaid)) |>
  select(-name) |>
  drop_na(au_iid)


#Re-read in the covariate file to get IIDs for MoBa sibs

cov <- readr::read_tsv ("//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt")

cov = cov %>%
  filter(!Role %in% c("Mother","Father")) |>
  select(ID_2306, SENTRIXID)

# Join to get sibIDs

sibs = cov |>
  left_join(kin |>
              select("SENTRIXID" = ID1,
                     ID2)) |>
  left_join(kin |>
              select("SENTRIXID" = ID2,
                     ID1)) |>
  pivot_longer(cols = c(-ID_2306,-SENTRIXID), values_to="sib_iid") |>
  drop_na(sib_iid) |>
  left_join(cov |>
              select("sib_iid"=SENTRIXID, "sib_mobaid"=ID_2306))

# Join to dataset with PREG_IDs to make child level

sib = base |>
  left_join(sibs |>
              separate(ID_2306, into = c("preg_id","BARN_NR"), sep="_") |>
              mutate(BARN_NR=as.double(BARN_NR)) |>
              separate(sib_mobaid, into = c("sib_preg_id","sib_BARN_NR"), sep="_") |>
              select(matches("preg_id|BARN_NR"))) |>
  drop_na(sib_preg_id)

# Summarise to datasets containing counts of valid responses at each age for: 
# - MoBa child (regardless genotyping)
ch_total = base |> 
  mutate(Nq6m = ifelse(!is.na(scl_full_m_6m),1,0),
         Nq18m = ifelse(!is.na(scl_full_m_18m),1,0),
         Nq3yr = ifelse(!is.na(scl_full_m_3yr),1,0),
         Nq5yr = ifelse(!is.na(scl_full_m_5yr),1,0),
         Nq8yr = ifelse(!is.na(scl_full_m_8yr),1,0),
         q14m = ifelse(!is.na(scl_full_m_14m),1,0),
         q14c = ifelse(!is.na(bmi_derived_c_14c),1,0),
         Nq14 = ifelse(q14m==1|q14c==1,1,0)) |> 
  select(matches("Nq")) |> 
  pivot_longer(everything()) |>
  filter(value==1) |> 
  group_by(name) |> 
  summarise(n=n()) |> 
  mutate(Category="Child (all MoBa)")

# - MoBa child (genotyped)
ch_geno = base |> 
  filter(geno_child==T) |> 
  mutate(Nq6m = ifelse(!is.na(scl_full_m_6m),1,0),
         Nq18m = ifelse(!is.na(scl_full_m_18m),1,0),
         Nq3yr = ifelse(!is.na(scl_full_m_3yr),1,0),
         Nq5yr = ifelse(!is.na(scl_full_m_5yr),1,0),
         Nq8yr = ifelse(!is.na(scl_full_m_8yr),1,0),
         q14m = ifelse(!is.na(scl_full_m_14m),1,0),
         q14c = ifelse(!is.na(bmi_derived_c_14c),1,0),
         Nq14 = ifelse(q14m==1|q14c==1,1,0)) |> 
  select(matches("Nq")) |> 
  pivot_longer(everything()) |>
  filter(value==1) |> 
  group_by(name) |> 
  summarise(n=n()) |> 
  mutate(Category="Child")
# - Mothers (genotyped)
set.seed(23523)
m_geno = base |> 
  filter(geno_mother==T) |> 
  mutate(Nq1 = ifelse(!is.na(scl_short_m_q1),1,0),
         Nq3 = ifelse(!is.na(scl_full_m_q3),1,0),
         Nq6m = ifelse(!is.na(scl_full_m_6m),1,0),
         Nq18m = ifelse(!is.na(scl_full_m_18m),1,0),
         Nq3yr = ifelse(!is.na(scl_full_m_3yr),1,0),
         Nq5yr = ifelse(!is.na(scl_full_m_5yr),1,0),
         Nq8yr = ifelse(!is.na(scl_full_m_8yr),1,0),
         Nq14 = ifelse(!is.na(scl_full_m_14m),1,0)) |> 
  select(matches("Nq"), m_id) |>
  group_by(m_id) |> 
  slice_sample(n=1) |> 
  ungroup() |> 
  select(!m_id) |> 
  pivot_longer(everything()) |>
  filter(value==1) |> 
  group_by(name) |> 
  summarise(n=n()) |> 
  mutate(Category="Mother")
# - Fathers (genotyped)
set.seed(23523)

f_geno = base |> 
  filter(geno_father==T) |>
  mutate(f2_year = 2015-birth_yr) |> 
  group_by(f_id) |> 
  slice_sample(n=1) |> 
  ungroup() |> 
  select(!f_id) |> 
  mutate(Nqf = ifelse(!is.na(scl_full_f_far),1,0),
         Nqf2_6 =ifelse(!is.na(scl_full_f_far2) & f2_year==6,1,0),
         Nqf2_7 =ifelse(!is.na(scl_full_f_far2) & f2_year==7,1,0),
         Nqf2_8 =ifelse(!is.na(scl_full_f_far2) & f2_year==8,1,0),
         Nqf2_9 =ifelse(!is.na(scl_full_f_far2) & f2_year==9,1,0),
         Nqf2_10 =ifelse(!is.na(scl_full_f_far2) & f2_year==10,1,0),
         Nqf2_11 =ifelse(!is.na(scl_full_f_far2) & f2_year==11,1,0),
         Nqf2_12 =ifelse(!is.na(scl_full_f_far2) & f2_year==12,1,0),
         Nqf2_13 =ifelse(!is.na(scl_full_f_far2) & f2_year==13,1,0),
         Nqf2_14 =ifelse(!is.na(scl_full_f_far2) & f2_year==14,1,0)) |> 
  select(matches("Nq")) |>
  pivot_longer(everything()) |>
  filter(value==1) |> 
  group_by(name) |> 
  summarise(n=n()) |> 
  mutate(Category="Father")
# - Siblings to a randomly selected MoBa child (genotyped)

sib_geno = sib |>
  select('preg_id'=sib_preg_id, 'BARN_NR'=sib_BARN_NR) |> 
  mutate(BARN_NR=as.double(BARN_NR)) |> 
  left_join(base) |> 
  filter(geno_child==T) |> 
  mutate(Nq6m = ifelse(!is.na(scl_full_m_6m),1,0),
         Nq18m = ifelse(!is.na(scl_full_m_18m),1,0),
         Nq3yr = ifelse(!is.na(scl_full_m_3yr),1,0),
         Nq5yr = ifelse(!is.na(scl_full_m_5yr),1,0),
         Nq8yr = ifelse(!is.na(scl_full_m_8yr),1,0),
         q14m = ifelse(!is.na(scl_full_m_14m),1,0),
         q14c = ifelse(!is.na(bmi_derived_c_14c),1,0),
         Nq14 = ifelse(q14m==1|q14c==1,1,0)) |> 
  select(matches("Nq")) |> 
  pivot_longer(everything()) |>
  filter(value==1) |> 
  group_by(name) |> 
  summarise(n=n()) |> 
  mutate(Category="Sib")

# - Aunts/uncles to a randomly selected MoBa child (genotyped)

aunt_geno = au |> 
  filter(str_detect(au_mobaid,"M")) |> 
  select('m_id'=au_mobaid,'own_preg_id' = preg_id , 'own_BARN_NR' = BARN_NR) |> 
  left_join(base) |> 
  filter(!(preg_id==own_preg_id & BARN_NR ==own_BARN_NR)) |> 
  filter(geno_mother==T) |> 
  mutate(Nq1 = ifelse(!is.na(scl_short_m_q1),1,0),
         Nq3 = ifelse(!is.na(scl_full_m_q3),1,0),
         Nq6m = ifelse(!is.na(scl_full_m_6m),1,0),
         Nq18m = ifelse(!is.na(scl_full_m_18m),1,0),
         Nq3yr = ifelse(!is.na(scl_full_m_3yr),1,0),
         Nq5yr = ifelse(!is.na(scl_full_m_5yr),1,0),
         Nq8yr = ifelse(!is.na(scl_full_m_8yr),1,0),
         Nq14 = ifelse(!is.na(scl_full_m_14m),1,0)) |> 
  select(matches("Nq"), m_id) |>
  group_by(m_id) |> 
  slice_sample(n=1) |> 
  ungroup() |> 
  select(!m_id) |> 
  pivot_longer(everything()) |>
  filter(value==1) |> 
  group_by(name) |> 
  summarise(n=n()) |> 
  mutate(Category="Aunt")

unc_geno = au |> 
  filter(str_detect(au_mobaid,"F")) |> 
  select('f_id'=au_mobaid,'own_preg_id' = preg_id , 'own_BARN_NR' = BARN_NR) |> 
  left_join(base) |> 
  filter(!(preg_id==own_preg_id & BARN_NR ==own_BARN_NR)) |> 
  filter(geno_father==T) |>
  mutate(f2_year = 2015-birth_yr) |> 
  group_by(f_id) |> 
  slice_sample(n=1) |> 
  ungroup() |> 
  select(!f_id) |> 
  mutate(Nqf = ifelse(!is.na(scl_full_f_far),1,0),
         Nqf2_6 =ifelse(!is.na(scl_full_f_far2) & f2_year==6,1,0),
         Nqf2_7 =ifelse(!is.na(scl_full_f_far2) & f2_year==7,1,0),
         Nqf2_8 =ifelse(!is.na(scl_full_f_far2) & f2_year==8,1,0),
         Nqf2_9 =ifelse(!is.na(scl_full_f_far2) & f2_year==9,1,0),
         Nqf2_10 =ifelse(!is.na(scl_full_f_far2) & f2_year==10,1,0),
         Nqf2_11 =ifelse(!is.na(scl_full_f_far2) & f2_year==11,1,0),
         Nqf2_12 =ifelse(!is.na(scl_full_f_far2) & f2_year==12,1,0),
         Nqf2_13 =ifelse(!is.na(scl_full_f_far2) & f2_year==13,1,0),
         Nqf2_14 =ifelse(!is.na(scl_full_f_far2) & f2_year==14,1,0)) |> 
  select(matches("Nq")) |>
  pivot_longer(everything()) |>
  filter(value==1) |> 
  group_by(name) |> 
  summarise(n=n()) |> 
  mutate(Category="Uncle")

### 

all = ch_total |> 
  bind_rows(ch_geno) |> 
  bind_rows(m_geno) |> 
  bind_rows(f_geno) |> 
  bind_rows(sib_geno) |> 
  bind_rows(aunt_geno) |> 
  bind_rows(unc_geno) |> 
  mutate(age = case_when(name=="Nq1"~-1,
                         name=="Nq3"~0,
                         name=="Nq6m"~.5,
                         name=="Nq18m"~1.5,
                         name=="Nq3yr"~3,
                         name=="Nq5yr"~5,
                         name=="Nq8yr"~8,
                         name=="Nq14"~14,
                         name=="Nqf"~-1,
                         name=="Nqf2_6"~6,
                         name=="Nqf2_7"~7,
                         name=="Nqf2_8"~8,
                         name=="Nqf2_9"~9,
                         name=="Nqf2_10"~10,
                         name=="Nqf2_11"~11,
                         name=="Nqf2_12"~12,
                         name=="Nqf2_13"~13,
                         name=="Nqf2_14"~14 ),
         n_neuro = case_when(name=="Nq1"&Category%in%c("Child","Child (all MoBa)","Sib")~0*n,
                             name=="Nq1"&Category%in%c("Mother","Aunt")~0*n,
                             name=="Nq3"&Category%in%c("Child","Child (all MoBa)","Sib")~0*n,
                             name=="Nq3"&Category%in%c("Mother","Aunt")~0*n,
                             name=="Nq6m"&Category%in%c("Child","Child (all MoBa)","Sib")~2*n,
                             name=="Nq6m"&Category%in%c("Mother","Aunt")~0*n,
                             name=="Nq18m"&Category%in%c("Child","Child (all MoBa)","Sib")~6*n,
                             name=="Nq18m"&Category%in%c("Mother","Aunt")~0*n,
                             name=="Nq3yr"&Category%in%c("Child","Child (all MoBa)","Sib")~12*n,
                             name=="Nq3yr"&Category%in%c("Mother","Aunt")~1*n,
                             name=="Nq5yr"&Category%in%c("Child","Child (all MoBa)","Sib")~13*n,
                             name=="Nq5yr"&Category%in%c("Mother","Aunt")~0*n,
                             name=="Nq8yr"&Category%in%c("Child","Child (all MoBa)","Sib")~8*n,
                             name=="Nq8yr"&Category%in%c("Mother","Aunt")~1*n,
                             name=="Nq14"&Category%in%c("Child","Child (all MoBa)","Sib")~6*n,
                             name=="Nq14"&Category%in%c("Mother","Aunt")~1*n,
                             name=="Nqf"~2*n,
                             str_detect(name, "Nqf2")~1*n),
         n_psych = case_when(name=="Nq1"&Category%in%c("Child","Child (all MoBa)","Sib")~0*n,
                             name=="Nq1"&Category%in%c("Mother","Aunt")~6*n,
                             name=="Nq3"&Category%in%c("Child","Child (all MoBa)","Sib")~0*n,
                             name=="Nq3"&Category%in%c("Mother","Aunt")~6*n,
                             name=="Nq6m"&Category%in%c("Child","Child (all MoBa)","Sib")~1*n,
                             name=="Nq6m"&Category%in%c("Mother","Aunt")~6*n,
                             name=="Nq18m"&Category%in%c("Child","Child (all MoBa)","Sib")~2*n,
                             name=="Nq18m"&Category%in%c("Mother","Aunt")~6*n,
                             name=="Nq3yr"&Category%in%c("Child","Child (all MoBa)","Sib")~4*n,
                             name=="Nq3yr"&Category%in%c("Mother","Aunt")~7*n,
                             name=="Nq5yr"&Category%in%c("Child","Child (all MoBa)","Sib")~3*n,
                             name=="Nq5yr"&Category%in%c("Mother","Aunt")~3*n,
                             name=="Nq8yr"&Category%in%c("Child","Child (all MoBa)","Sib")~5*n,
                             name=="Nq8yr"&Category%in%c("Mother","Aunt")~6*n,
                             name=="Nq14"&Category%in%c("Child","Child (all MoBa)","Sib")~20*n,
                             name=="Nq14"&Category%in%c("Mother","Aunt")~6*n,
                             name=="Nqf"~6*n,
                             str_detect(name, "Nqf2")~5*n),
         n_physical = case_when(name=="Nq1"&Category%in%c("Child","Child (all MoBa)","Sib")~0*n,
                                name=="Nq1"&Category%in%c("Mother","Aunt")~7*n,
                                name=="Nq3"&Category%in%c("Child","Child (all MoBa)","Sib")~0*n,
                                name=="Nq3"&Category%in%c("Mother","Aunt")~7*n,
                                name=="Nq6m"&Category%in%c("Child","Child (all MoBa)","Sib")~6*n,
                                name=="Nq6m"&Category%in%c("Mother","Aunt")~4*n,
                                name=="Nq18m"&Category%in%c("Child","Child (all MoBa)","Sib")~6*n,
                                name=="Nq18m"&Category%in%c("Mother","Aunt")~2*n,
                                name=="Nq3yr"&Category%in%c("Child","Child (all MoBa)","Sib")~2*n,
                                name=="Nq3yr"&Category%in%c("Mother","Aunt")~1*n,
                                name=="Nq5yr"&Category%in%c("Child","Child (all MoBa)","Sib")~1*n,
                                name=="Nq5yr"&Category%in%c("Mother","Aunt")~1*n,
                                name=="Nq8yr"&Category%in%c("Child","Child (all MoBa)","Sib")~1*n,
                                name=="Nq8yr"&Category%in%c("Mother","Aunt")~1*n,
                                name=="Nq14"&Category%in%c("Child","Child (all MoBa)","Sib")~3*n,
                                name=="Nq14"&Category%in%c("Mother","Aunt")~2*n,
                                name=="Nqf"~1*n,
                                str_detect(name, "Nqf2")~3*n),
         n_healthb = case_when(name=="Nq1"&Category%in%c("Child","Child (all MoBa)","Sib")~0*n,
                               name=="Nq1"&Category%in%c("Mother","Aunt")~11*n,
                               name=="Nq3"&Category%in%c("Child","Child (all MoBa)","Sib")~0*n,
                               name=="Nq3"&Category%in%c("Mother","Aunt")~9*n,
                               name=="Nq6m"&Category%in%c("Child","Child (all MoBa)","Sib")~6*n,
                               name=="Nq6m"&Category%in%c("Mother","Aunt")~7*n,
                               name=="Nq18m"&Category%in%c("Child","Child (all MoBa)","Sib")~8*n,
                               name=="Nq18m"&Category%in%c("Mother","Aunt")~4*n,
                               name=="Nq3yr"&Category%in%c("Child","Child (all MoBa)","Sib")~8*n,
                               name=="Nq3yr"&Category%in%c("Mother","Aunt")~1*n,
                               name=="Nq5yr"&Category%in%c("Child","Child (all MoBa)","Sib")~3*n,
                               name=="Nq5yr"&Category%in%c("Mother","Aunt")~2*n,
                               name=="Nq8yr"&Category%in%c("Child","Child (all MoBa)","Sib")~2*n,
                               name=="Nq8yr"&Category%in%c("Mother","Aunt")~3*n,
                               name=="Nq14"&Category%in%c("Child","Child (all MoBa)","Sib")~8*n,
                               name=="Nq14"&Category%in%c("Mother","Aunt")~5*n,
                               name=="Nqf"~9*n,
                               str_detect(name, "Nqf2")~8*n),
         n_socenv = case_when(name=="Nq1"&Category%in%c("Child","Child (all MoBa)","Sib")~0*n,
                              name=="Nq1"&Category%in%c("Mother","Aunt")~12*n,                
                              name=="Nq3"&Category%in%c("Child","Child (all MoBa)","Sib")~0*n,
                              name=="Nq3"&Category%in%c("Mother","Aunt")~8*n,
                              name=="Nq6m"&Category%in%c("Child","Child (all MoBa)","Sib")~1*n,
                              name=="Nq6m"&Category%in%c("Mother","Aunt")~4*n,
                              name=="Nq18m"&Category%in%c("Child","Child (all MoBa)","Sib")~3*n,
                              name=="Nq18m"&Category%in%c("Mother","Aunt")~6*n,
                              name=="Nq3yr"&Category%in%c("Child","Child (all MoBa)","Sib")~4*n,
                              name=="Nq3yr"&Category%in%c("Mother","Aunt")~6*n,
                              name=="Nq5yr"&Category%in%c("Child","Child (all MoBa)","Sib")~4*n,
                              name=="Nq5yr"&Category%in%c("Mother","Aunt")~3*n,
                              name=="Nq8yr"&Category%in%c("Child","Child (all MoBa)","Sib")~6*n,
                              name=="Nq8yr"&Category%in%c("Mother","Aunt")~5*n,
                              name=="Nq14"&Category%in%c("Child","Child (all MoBa)","Sib")~11*n,
                              name=="Nq14"&Category%in%c("Mother","Aunt")~6*n,
                              name=="Nqf"~10*n,
                              str_detect(name, "Nqf2")~9*n)) |> 
  rename('n_overall'=n, 'wave'= name) |> 
  pivot_longer(cols=matches("n_")) |> 
  mutate(Role = factor(case_when(Category%in%c("Child","Child (all MoBa)")~"Children",
                                 Category%in%c("Mother","Father","Sib")~"1\u00B0 relatives",
                                 Category%in%c("Aunt","Uncle") ~"2\u00B0 relatives" ),
                       levels= c("Children","1\u00B0 relatives","2\u00B0 relatives")))|> 
  mutate(Category=factor(Category,levels=c("Child (all MoBa)","Child","Mother","Father","Sib","Aunt","Uncle"),
                         labels=c("Child (all MoBa)","Child","Mother","Father","Sibling","Aunt","Uncle")),
         Domain = factor(name, levels=c("n_neuro","n_psych","n_physical","n_healthb","n_socenv","n_overall"),
                         labels=c("Neurodevelopmental","Mental health &\nlife satisfaction",
                                  "Physical\nhealth & pain","Health behaviours\n& medication",
                                  "Social &\nenvironmental", "Overall")))



ggplot(all |> 
         filter(name!= "n_overall",
                Category!="Child (all MoBa)"), aes(x=age, y=value, group=Category,fill= factor(Category) ))+
  geom_smooth(method = "loess",
              se = FALSE,
              formula = 'y ~ x',
              span = 0.8,
              aes(colour=factor(Category))) +
  stat_smooth(
    geom = 'area', method = 'loess', span = 0.8,
    alpha = 1/12)+
  scale_colour_manual("Participant", 
                      values= c("#E31A1C",  "#33A02C", "#B2DF8A","#FB9A99", "#1F78B4", "#A6CEE3") )+
  scale_fill_manual("Participant",                      
                    values= c("#E31A1C",  "#33A02C", "#B2DF8A","#FB9A99", "#1F78B4", "#A6CEE3") )+
  facet_grid(Role~Domain, scales="free_y", space= "free_y") + 
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t =  0, r = 20, b = 0, l = 0), angle=90),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 120, l = 0)),
        #axis.text.y  = element_blank(),
        #axis.ticks.y  = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle=0),
        legend.position = c(0.8,-0.35),
        legend.direction = "horizontal",
        plot.caption = element_text(size=10,face="italic"))+
  scale_y_continuous("N observations",labels=scales::comma, breaks= c(seq(0,800000,100000)))+
  scale_x_continuous("Child age at measurement")

ggsave("./moba/alt1_grid.tiff",  device="tiff", width=26,height=22,units="cm",dpi=300,bg="white")

## Bar plot for relative type counts

Ns = tibble("Role"= factor(c("Child\n(all MoBa)",
                      "Child",
                      "Mother",
                      "Father",
                      "Sibling",
                      "Aunt",
                      "Uncle"), levels=c("Child\n(all MoBa)",
                                         "Child",
                                         "Mother",
                                         "Father",
                                         "Sibling",
                                         "Aunt",
                                         "Uncle")), 
            "N"= c(nrow(base), 
                   nrow(base |> filter(geno_child==T)),
                   nrow(base |> filter(geno_mother==T) |> select(m_id) |> distinct()),
                   nrow(base |> filter(geno_father==T) |> select(f_id) |> distinct()),
                   nrow(sib),
                   nrow(au |>  filter(str_detect(au_mobaid,"M"))),
                   nrow(au |> filter(str_detect(au_mobaid,"F"))))) 

c("#B15928", "#FFFF99", "#6A3D9A", "#CAB2D6", "#FF7F00", "#FDBF6F", 
  "#E31A1C", "#FB9A99", "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3"
)

ggplot( Ns,  aes(x=Role,  y= N, fill= Role, colour=Role))+
  geom_bar(stat = "identity", position=position_dodge(1),alpha=0.6, linewidth=1, width=0.8)+
  scale_colour_manual(values= c("#E31A1C","#E31A1C", "#33A02C", "#B2DF8A", "#FB9A99", "#1F78B4", "#A6CEE3"))+
  scale_fill_manual(values= c("#E31A1C","#E31A1C", "#33A02C", "#B2DF8A", "#FB9A99", "#1F78B4", "#A6CEE3"))+
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t =  0, r = 20, b = 0, l = 0), angle=90),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 5, l = 0)),
        #axis.text.y  = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour="grey10"),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle=0),
        legend.position = "none",
        legend.direction = "horizontal",
        plot.caption = element_text(size=10,face="italic"))+
  scale_x_discrete("")+
  scale_y_continuous("N participants",labels=scales::comma)

ggsave("./moba/basicNs.tiff",  device="tiff", width=12,height=8,units="cm",dpi=300,bg="white")




### ALTERNATIVES

ggplot(all |> 
         filter(name!= "n_overall",
                !Category %in% c("Aunt","Uncle")), aes(x=age, y=value, group=factor(Domain),fill= factor(Domain)))+
  geom_point(aes( colour=factor(Domain))) +
  scale_colour_brewer("Domain", palette = "Set2", direction=-1)+
  scale_fill_brewer("Domain", palette = "Set2", direction=-1)+
  geom_smooth(method = "loess",
              se = FALSE,
              formula = 'y ~ x',
              span = 0.8,
              aes(colour=factor(Domain))) +
  stat_smooth(
    geom = 'area', method = 'loess', span = 0.8,
    alpha = 1/12)+
  facet_wrap(Category~., scales="free") + 
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t =  0, r = 20, b = 0, l = 0), angle=90),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 20, l = 0)),
        #axis.text.y  = element_blank(),
        #axis.ticks.y  = element_blank(),
        strip.text = element_text(face="bold"),
        legend.position = "right",
        legend.direction = "vertical",
        plot.caption = element_text(size=10,face="italic"))+
  scale_y_continuous("N observations*",labels=scales::comma)+
  scale_x_continuous("Child age at measurement")+
  labs(caption="*Calculated as N individuals in category x N measures in domain at measurement occasion")

ggsave("./moba/alt2_grid.tiff",  device="tiff", width=18,height=28,units="cm",dpi=300,bg="white")

ggsave("./moba/alt2_wrap.tiff",  device="tiff", width=24,height=24,units="cm",dpi=300,bg="white")


