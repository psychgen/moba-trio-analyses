# 03_gwas_figure.R

# The purpose of this script is to create a composite figure incorporating all 
# main text display elements for the GWAS component of the manuscript

source("./set_aes_pars.R")

#devtools::install_local("//ess01/p471/data/durable/common/software/miamiplot-master.zip")

library(tidyverse)
library(ggthemes)
library(ggsci)
library(miamiplot)
library(patchwork)

# Read in data


# Panel A h2 shrinkage

tmp_height_con <- read_lines("./trio-GWAS/results/ldsc_aug25/height_trio_EA3.log")
tmp_height_pop <- read_lines("./trio-GWAS/results/ldsc_aug25/height_pop_EA3.log") 
tmp_exam_con <- read_lines("./trio-GWAS/results/ldsc_aug25/exam_norm_trio_EA3.log") 
tmp_exam_pop <- read_lines("./trio-GWAS/results/ldsc_aug25/exam_norm_pop_EA3.log") 
tmp_dep_con <- read_lines("./trio-GWAS/results/ldsc_aug25/smfq_dep_trio_EA3.log")
tmp_dep_pop <- read_lines("./trio-GWAS/results/ldsc_aug25/smfq_dep_pop_EA3.log") 
tmp_sleep_con <- read_lines("./trio-GWAS/results/ldsc_aug25/sleep_hrs_trio_EA3.log") 
tmp_sleep_pop <- read_lines("./trio-GWAS/results/ldsc_aug25/sleep_hrs_pop_EA3.log") 

out=data.frame()
for (tmp in list( tmp_height_con,tmp_height_pop,tmp_exam_con,tmp_exam_pop,
                  tmp_dep_con,tmp_dep_pop,tmp_sleep_con,tmp_sleep_pop)){
  
  res =  str_remove_all( str_remove_all(tmp[grep("Total Observed scale h2: ",tmp)][[1]], "Total Observed scale h2: "), "[[)(]]")
  
  out = out |> 
    rbind(tibble("file"=tmp[[11]], "values"= res) |> 
            separate(values, into=c("h2","SE"), sep=" ") |> 
            mutate(h2=as.numeric(h2),
                   SE=as.numeric(SE)))
}
ploth2= out |> 
  mutate(Outcome = factor(rep(c("Height","Educational\nachievment","Depressive\nsymptoms","Sleep\nduration"),each=2),
                          levels=c("Height","Educational\nachievment","Sleep\nduration","Depressive\nsymptoms")),
         Model = factor(rep(c("Trio","Population"),4),
                        levels=c("Population","Trio")))|> 
  rowwise() |> 
  mutate(lci= h2-1.96*SE,
         uci= h2+1.96*SE)



p2.1=ggplot(ploth2, aes(x=Outcome, y=h2, fill=Model))+
  geom_col(position=position_dodge(1), colour="black")+
  geom_errorbar(aes(ymin=lci,ymax=uci),  width=0, linewidth=1, position= position_dodge(1), colour="grey40" )+
  scale_fill_manual(values=c("#FF9933","#99CCFF"))+
  theme_few()+
  xlab("Outcome trait and GWAS model")+
  ylab("Estimated SNP heritability\n(SNP h2)")+
  theme(axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.line = element_line(colour = "black"),
        panel.border =  element_blank())

ggsave(p2.1, filename = "./trio-GWAS/plots/h2_comp.tiff", device = "tiff", width = 14, 
       height = 10, units = "cm", dpi = 300)




# Panel B rG shrinkage

out= data.frame()

for (file in list.files("./trio-GWAS/results/ldsc_aug25/")[str_detect(list.files("./trio-GWAS/results/ldsc_aug25/"),".log")] ){
  
  tmp <- read_lines(paste0("./trio-GWAS/results/ldsc_aug25/", file)) 
  
  res =  str_remove_all( str_remove_all(tmp[grep("Genetic Correlation: ",tmp)], "Genetic Correlation: "), "[[)(]]")
  
  out = out |> 
    rbind(tibble("file"=file, "values"= res) |> 
            separate(values, into=c("rG","SE"), sep=" ") |> 
            mutate(rG=as.numeric(rG),
                   SE=as.numeric(SE)))
  
}


plot_rgs=out |> 
  filter(!str_detect(file,"pop_v_trio")) |> 
  mutate(Model = factor(ifelse(str_detect(file, "trio"), "Trio", "Population")),
         Trait = factor(case_when(str_detect(file,"EA3") ~ "Educ. attain.",
                                  str_detect(file,"smoke") ~ "Smoking",
                                  str_detect(file,"locke") ~ "BMI",
                                  str_detect(file,"wood") ~ "Height",
                                  str_detect(file,"mdd") ~ "MDD")),
         Outcome = factor(case_when(str_detect(file,"exam") ~ "Educational\nachievement",
                                    str_detect(file,"height_") ~ "Height",
                                    str_detect(file,"smfq") ~ "Depressive\nsymptoms",
                                    str_detect(file,"sleep") ~ "Sleep\nduration"),
                          levels=c("Height","Educational\nachievement","Sleep\nduration","Depressive\nsymptoms"))) |> 
  rowwise() |> 
  mutate(lci= rG-1.96*SE,
         uci= rG+1.96*SE)

p2.2 = ggplot(plot_rgs |>  filter(file!="sleep_hrs_trio_lockebmi.log"), aes(x= rG, y=Trait,  fill=Model, shape= Model, group= Model))+
  geom_vline(xintercept=0, linetype= 2, linewidth=1.2, colour="grey80")+
  geom_errorbarh(aes(xmin=lci,xmax=uci), height=0, linewidth=ci_lw, position= position_dodge(pd), colour="grey60" ) +
  geom_point(size=pt_sz, stroke=stroke_val,position= position_dodge(0.8), colour=stroke_col)+
  facet_grid(Outcome~.)+
  scale_fill_manual(values=c("#FF9933","#99CCFF"))+
  scale_shape_manual(values=c(21,22))+
  theme_few()+
  ylab("Externally GWASed trait")+
  xlab("rG with MoBa a) population \nand b) trio GWAS")+
  coord_cartesian(xlim=c(-1,1))+
  guides(fill=guide_legend(reverse = F),
         shape=guide_legend(reverse = F))+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.background = element_rect(colour="grey80"))

ggsave(p2.2, filename = "./trio-GWAS/plots/rg_comp.tiff", device = "tiff", width = 14, 
       height = 18, units = "cm", dpi = 600)


# Supplement X miami plot height

height_con = data.table::fread("./trio-GWAS/sumstats/model06_trios_height_mqc.trios")
height_pop = data.table::fread("./trio-GWAS/sumstats/model06_trios_height_mqc.basic")

height_both = height_con |>
  mutate(model="Trio") |>
  select(model,'CHR' = Chromosome,'BP' = Basepair ,'P'=Offspring_P) |>
  bind_rows(height_pop |>
              mutate(model="Population") |>
              select(model,'CHR' = Chromosome,'BP' = Basepair ,'P'=Wald_P))

#plot
miami_both<- ggmiami(data = height_both,
                     split_by = "model",
                     split_at = "Trio",
                     p = "P", chr = "CHR", pos = "BP",
                     genome_line = 1e-05,
                     genome_line_color = "grey40",
                     suggestive_line = 5e-08,
                     suggestive_line_color = "blue",
                     chr_colors = NULL,
                     upper_chr_colors = c("#99CCFF","#0072B2"),
                     lower_chr_colors = c("#FF9933","#D55E00"),
                     upper_ylab = "Trio",
                     lower_ylab = "Population")

ggsave(miami_both, filename = "./trio-GWAS/plots/miami_height.tiff", device = "tiff", width = 18,
       height = 10, units = "cm", dpi = 300)

rm(height_con,height_pop,height_both)
gc()
# Supplement X (supp) miami plot educational acheivement

exam_con = data.table::fread("./trio-GWAS/sumstats/model06_trios_exam_norm_mqc.trios")
exam_pop = data.table::fread("./trio-GWAS/sumstats/model06_trios_exam_norm_mqc.basic")

exam_both = exam_con  |>
  mutate(model="Trio") |>
  select(model,'CHR' = Chromosome,'BP' = Basepair ,'P'=Offspring_P) |>
  bind_rows(exam_pop |>
              mutate(model="Population") |>
              select(model,'CHR' = Chromosome,'BP' = Basepair ,'P'=Wald_P))


#plot
miami_both2<- ggmiami(data = exam_both, 
                      split_by = "model", 
                      split_at = "Trio", 
                      p = "P", chr = "CHR", pos = "BP",
                      genome_line = 1e-05,
                      genome_line_color = "grey40",
                      suggestive_line = 5e-08,
                      suggestive_line_color = "blue",
                      chr_colors = NULL,
                      upper_chr_colors = c("#99CCFF","#0072B2"),
                      lower_chr_colors = c("#FF9933","#D55E00"),
                      upper_ylab = "Trio",
                      lower_ylab = "Population")

ggsave(miami_both2, filename = "./trio-GWAS/plots/miami_exam.tiff", device = "tiff", width = 18, 
       height = 10, units = "cm", dpi = 300)
rm(exam_con,exam_pop,exam_both)
gc()
# Supplement X (supp) miami plot sleep duration

slp_con = data.table::fread("./trio-GWAS/sumstats/model06_trios_sleep_hrs_mqc.trios")
slp_pop = data.table::fread("./trio-GWAS/sumstats/model06_trios_sleep_hrs_mqc.basic")

slp_both = slp_con  |>
  mutate(model="Trio") |>
  select(model,'CHR' = Chromosome,'BP' = Basepair ,'P'=Offspring_P) |>
  bind_rows(slp_pop |>
              mutate(model="Population") |>
              select(model,'CHR' = Chromosome,'BP' = Basepair ,'P'=Wald_P))


#plot
miami_both3<- ggmiami(data = slp_both, 
                      split_by = "model", 
                      split_at = "Trio", 
                      p = "P", chr = "CHR", pos = "BP",
                      genome_line = 1e-05,
                      genome_line_color = "grey40",
                      suggestive_line = 5e-08,
                      suggestive_line_color = "blue",
                      chr_colors = NULL,
                      upper_chr_colors = c("#99CCFF","#0072B2"),
                      lower_chr_colors = c("#FF9933","#D55E00"),
                      upper_ylab = "Trio",
                      lower_ylab = "Population")

ggsave(miami_both3, filename = "./trio-GWAS/plots/miami_slp.tiff", device = "tiff", width = 18, 
       height = 10, units = "cm", dpi = 300)
rm(slp_con,slp_pop,slp_both)
gc()
# Supplement X (supp) miami plot depressive symptoms 

dep_con = data.table::fread("./trio-GWAS/sumstats/model06_trios_smfq_dep_mqc.trios")
dep_pop = data.table::fread("./trio-GWAS/sumstats/model06_trios_smfq_dep_mqc.basic")

dep_both = dep_con  |>
  mutate(model="Trio") |>
  select(model,'CHR' = Chromosome,'BP' = Basepair ,'P'=Offspring_P) |>
  bind_rows(dep_pop |>
              mutate(model="Population") |>
              select(model,'CHR' = Chromosome,'BP' = Basepair ,'P'=Wald_P))


#plot
miami_both4<- ggmiami(data = dep_both, 
                      split_by = "model", 
                      split_at = "Trio", 
                      p = "P", chr = "CHR", pos = "BP",
                      genome_line = 1e-05,
                      genome_line_color = "grey40",
                      suggestive_line = 5e-08,
                      suggestive_line_color = "blue",
                      chr_colors = NULL,
                      upper_chr_colors = c("#99CCFF","#0072B2"),
                      lower_chr_colors = c("#FF9933","#D55E00"),
                      upper_ylab = "Trio",
                      lower_ylab = "Population")

ggsave(miami_both4, filename = "./trio-GWAS/plots/miami_dep.tiff", device = "tiff", width = 18, 
       height = 10, units = "cm", dpi = 300)
rm(dep_con,dep_pop,dep_both)
gc()
