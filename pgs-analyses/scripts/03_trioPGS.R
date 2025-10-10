# 03_trioPGS.R

# In this script we run the trio PGS analyses for both outcomes

library(tidyverse)
library(lavaan)
library(patchwork)

load("./data/dat_final.RData")


run_trioPGS <- function(out, covs, scores, data, clust){
  
  
  all_res=tibble()
  
  for(score in scores){
    
    chscore=paste0("child_",score)
    mscore=paste0("mother_",score)
    fscore=paste0("father_",score)
    
    basemod = paste0(out, " ~ ",chscore, " + ", paste0(covs,collapse = " + ") )
    
    basefit= lavaan::sem(basemod, 
                         data=data, 
                         estimator="MLR",
                         se="robust",
                         cluster=clust)
    
    baseline_est <- standardizedsolution(basefit, type="std.nox") %>% 
      filter (op == "~" ) %>% 
      mutate(Type="Child-only")
    
    triomod = paste0(out, " ~ ",chscore, " + ", mscore, " + ", fscore, " + ", paste0(covs,collapse = " + "))
    
    triofit= lavaan::sem(triomod, 
                         data=data, 
                         estimator="MLR",
                         se="robust",
                         cluster=clust) 
    
    trio_ests <- standardizedsolution(triofit, type="std.nox") %>% 
      filter (op == "~" ) %>% 
      mutate(Type="Trio")
    
    all_ests <- baseline_est %>% 
      as_tibble() %>% 
      bind_rows(trio_ests %>% 
                  as_tibble)
    
    all_res = rbind(all_res,all_ests)
  }
  
  return(all_res)
  
} 



pgs_names <- c("height2018",
               "bmi2018",
               "EA2018",
               "smok2019",
               "mdd2025")

pgs_display <- c("Height","BMI", "Edu. attain.", "Ever smoked", "MDD")



trio_ex = run_trioPGS(out="exam", 
                      covs=c("SEX"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_final,
                      clust="m_id")
trio_ht = run_trioPGS(out="height", 
                      covs=c("SEX","Age_height"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_final,
                      clust="m_id")
trio_dp = run_trioPGS(out="smfq_dep", 
                      covs=c("SEX","Age_dep"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_final,
                      clust="m_id")
trio_sl = run_trioPGS(out="sleep_hrs", 
                      covs=c("SEX","Age_sleep"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_final,
                      clust="m_id")
# Do some post-processing of results

trio_res = trio_ex %>% 
  bind_rows(trio_ht) |> 
  bind_rows(trio_dp) |> 
  bind_rows(trio_sl) |> 
  separate(rhs, into = c("Member", "PGS"), sep="_") %>% 
  mutate(PGS = factor(str_remove_all(PGS,"_0.95_pgs_res"), levels= pgs_names, labels = pgs_display )) %>% 
  select('Outcome'=lhs,Type, Member,PGS,est.std,se,pvalue,ci.lower,ci.upper) %>%
  drop_na(PGS)

# Save out

save(trio_res, file = "./output/triopgs_res.RData")

# Sensitivity analyses

#1) unrelated trios only

library(genotools)
ref = unrelate() |> 
  filter(Role=="Child")

dat_unrelated = ref |> 
  ungroup()|> 
  select(preg_id, BARN_NR) |> 
  mutate(BARN_NR= as.double(BARN_NR)) |> 
  left_join(dat_final)

trio_ex = run_trioPGS(out="exam", 
                      covs=c("SEX"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_unrelated,
                      clust=NULL)
trio_ht = run_trioPGS(out="height", 
                      covs=c("SEX","Age_height"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_unrelated,
                      clust=NULL)
trio_dp = run_trioPGS(out="smfq_dep", 
                      covs=c("SEX","Age_dep"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_unrelated,
                      clust=NULL)
trio_sl = run_trioPGS(out="sleep_hrs", 
                      covs=c("SEX","Age_sleep"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_unrelated,
                      clust=NULL)
# Do some post-processing of results

trio_res = trio_ex %>% 
  bind_rows(trio_ht) |> 
  bind_rows(trio_dp) |> 
  bind_rows(trio_sl) |> 
  separate(rhs, into = c("Member", "PGS"), sep="_") %>% 
  mutate(PGS = factor(str_remove_all(PGS,"_0.95_pgs_res"), levels= pgs_names, labels = pgs_display )) %>% 
  select('Outcome'=lhs,Type, Member,PGS,est.std,se,pvalue,ci.lower,ci.upper) %>%
  drop_na(PGS)

# Save out

save(trio_res, file = "./output/triopgs_res_unrelated.RData")

#2) cluster on family id

trio_ex = run_trioPGS(out="exam", 
                      covs=c("SEX"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_final,
                      clust="FID")
trio_ht = run_trioPGS(out="height", 
                      covs=c("SEX","Age_height"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_final,
                      clust="FID")
trio_dp = run_trioPGS(out="smfq_dep", 
                      covs=c("SEX","Age_dep"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_final,
                      clust="FID")
trio_sl = run_trioPGS(out="sleep_hrs", 
                      covs=c("SEX","Age_sleep"), #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_pgs_res"), 
                      data=dat_final,
                      clust="FID")

# Do some post-processing of results

trio_res = trio_ex %>% 
  bind_rows(trio_ht) |> 
  bind_rows(trio_dp) |> 
  bind_rows(trio_sl) |> 
  separate(rhs, into = c("Member", "PGS"), sep="_") %>% 
  mutate(PGS = factor(str_remove_all(PGS,"_0.95_pgs_res"), levels= pgs_names, labels = pgs_display )) %>% 
  select('Outcome'=lhs,Type, Member,PGS,est.std,se,pvalue,ci.lower,ci.upper) %>%
  drop_na(PGS)

# Save out

save(trio_res, file = "./output/triopgs_res_altclust.RData")
