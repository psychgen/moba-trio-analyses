# 01_qc_figure.R

# The purpose of this script is to create a composite figure incorporating all 
# main text display elements for the QC component of the manuscript

source("./set_aes_pars.R")

library(tidyverse)
library(ggthemes)
library(ggsci)
library(patchwork)

# 3*2-panel figure, vertical LHS double-panel for abbreviated pipeline, horizontal 
# double-panel at the bottom for kinship, two single panel RHS plots - suggest
# these are used, respectively, to summarise something about individuals and 
# something about SNPs

# RHS plot 1: Individuals

covs <-readr::read_delim("//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt")

covs_long = covs |> 
  select(FID, Role, imputation_batch) |> 
  mutate(Batch = factor(case_when(imputation_batch=="release1-hce" ~ "R1-HCE",
                                  imputation_batch=="release1-omni" ~ "R1-OMN",
                                  imputation_batch=="release1-gsa" ~ "R1-GSA",
                                  imputation_batch=="release2" ~ "R2",
                                  imputation_batch=="release3" ~ "R3",
                                  imputation_batch=="release4" ~ "R4") ))

individuals_batch = covs_long |> group_by(Batch, Role) |>  summarise(n=n())

families_batch = covs_long |> group_by(Batch) |>  distinct(FID) |> summarise(n=n())

p1.1 = ggplot( )+
  geom_col(data= families_batch, aes(y="Father", x=n ,colour= "grey60"),fill= "grey80", linewidth=0, width=2)+
  geom_col(data= individuals_batch, aes(y=Role, x=n, fill=Batch ),alpha=0.8, colour= "grey60",linewidth=1.1)+
  facet_grid(Batch ~.)+
  scale_fill_npg()+
  scale_color_identity(guide="legend", name=NULL, labels="N unique families")+
  theme_few()+
  guides(fill= "none")+
  xlab("N individuals genotyped in batch")+
  theme(legend.position = c(0.8,0.93),
        axis.title.y = element_blank())


# RHS plot 2: SNPs

snpstats <-readr::read_delim("//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.snpstats")

snps_long = snpstats |> 
  select(matches("INFO")) |> 
  pivot_longer(everything(),values_to="INFO score") |> 
  mutate(Batch = factor(case_when(name=="INFO_release1-hce" ~ "R1-HCE",
                          name=="INFO_release1-omni" ~ "R1-OMN",
                          name=="INFO_release1-gsa" ~ "R1-GSA",
                          name=="INFO_release2" ~ "R2",
                          name=="INFO_release3" ~ "R3",
                          name=="INFO_release4" ~ "R4") ))

# Density plot of info scores by batch
p1.2 = ggplot(snps_long, aes(x= `INFO score`, colour= Batch, fill=Batch))+
  geom_histogram(aes(y=..density..),colour= "grey60",alpha=0.2, position="dodge")+
  geom_density(alpha=0.6, fill=NA, linewidth=1.5)+
  coord_cartesian(xlim=c(0.95,0.99), ylim=c(0,40))+
  scale_color_npg()+
  scale_fill_npg()+
  theme_few()+
  ylab("Density")+
  xlab("INFO score for imputation quality by batch")+
  theme(axis.line = element_line(colour = "black"),
        panel.border =  element_blank())

# Batch-wise imputation summaries

imp_long = snpstats |> 
  select(matches("IMPUTED")) |> 
  pivot_longer(everything(),values_to="IMP") |> 
  group_by(name) |> 
  summarise(imp = sum(IMP, na.rm=T)) |> 
  mutate(prop_imp = imp/nrow(snpstats))|> 
  mutate(Batch = factor(case_when(name=="IMPUTED_release1-hce" ~ "R1-HCE",
                                  name=="IMPUTED_release1-omni" ~ "R1-OMN",
                                  name=="IMPUTED_release1-gsa" ~ "R1-GSA",
                                  name=="IMPUTED_release2" ~ "R2",
                                  name=="IMPUTED_release3" ~ "R3",
                                  name=="IMPUTED_release4" ~ "R4") ))

p1.3=ggplot(imp_long, aes(x= prop_imp, y=fct_rev(Batch), fill=Batch))+
  geom_col(colour= "grey60",linewidth=1.1)+
  coord_cartesian(xlim=c(0.9,0.99))+
  scale_fill_npg()+
  theme_few()+
  xlab("Proportion of SNPs imputed per batch")+
  ylab("Batch")+
  theme(legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.border =  element_blank())
  

  
# Summary of N_genotyped
geno_sum =
  snpstats |> 
  filter(N_genotyped>0) %>% 
  group_by(N_genotyped) |> 
  summarise(n=n(),prop=n()/nrow(snpstats)) |> 
  mutate(`N batches` = as.character(N_genotyped))

p1.4=ggplot(geno_sum, aes(x= `N batches`, y=prop, fill=`N batches`))+
  geom_col()+
  coord_cartesian(ylim=c(0,0.05))+
  scale_fill_few()+
  theme_few()+
  ylab("SNPs directly genotyped (proportional of total)")+
  xlab("N batches in which genotyped")+
  theme(legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.border =  element_blank())


ggsave("./qc-pipeline/p1-3.tiff",p1.3, dpi=200, device= "tiff", width=13, height= 10, units ="cm" )



p1.234= ((p1.4+p1.3)/p1.2) +
  plot_layout(widths=c(2,1,1,2))+
  plot_annotation(theme = theme(plot.background = element_rect(colour= "black", size=1)))

ggsave("./qc-pipeline/p1-234.tiff",p1.234, dpi=200, device= "tiff", width=17, height= 20, units ="cm" )

p1.1 = p1.1 + 
  plot_annotation(theme = theme(plot.background = element_rect(colour= "black", size=1)))

ggsave("./qc-pipeline/p1-1.tiff",p1.1, dpi=200, device= "tiff", width=13, height= 13, units ="cm" )


p1.1/p1.234+
  plot_annotation(theme = theme(plot.background = element_rect(colour= "black", size=1)))
