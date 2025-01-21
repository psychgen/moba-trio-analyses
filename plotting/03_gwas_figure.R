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

tmp_height_con <- read_lines("./trio-GWAS/results/heightcon_EA3.log")
tmp_height_pop <- read_lines("./trio-GWAS/results/heightpop_EA3.log") 
tmp_exam_con <- read_lines("./trio-GWAS/results/examcon_EA3.log") 
tmp_exam_pop <- read_lines("./trio-GWAS/results/exampop_EA3.log") 

out=data.frame()
for (tmp in list( tmp_height_con,tmp_height_pop,tmp_exam_con,tmp_exam_pop)){

res =  str_remove_all( str_remove_all(tmp[grep("Total Observed scale h2: ",tmp)][[1]], "Total Observed scale h2: "), "[[)(]]")

out = out |> 
  rbind(tibble("file"=tmp[[11]], "values"= res) |> 
          separate(values, into=c("h2","SE"), sep=" ") |> 
          mutate(h2=as.numeric(h2),
                 SE=as.numeric(SE)))
}
ploth2= out |> 
  mutate(Outcome = rep(c("Height","Educational achievment"),each=2),
         Model = factor(rep(c("Trio","Population"),2),
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
  ylab("Estimated SNP heritability (SNP h2)")+
  theme(axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.line = element_line(colour = "black"),
        panel.border =  element_blank())

ggsave(p2.1, filename = "./trio-GWAS/plots/h2_comp.tiff", device = "tiff", width = 10, 
       height = 10, units = "cm", dpi = 600)




# Panel B miami plot height

height_con = data.table::fread("./trio-GWAS/sumstats/height_con_filtered_ldsc.txt")
height_pop = data.table::fread("./trio-GWAS/sumstats/height_pop_filtered_ldsc.txt")

height_both = height_con |> 
  mutate(model="Trio") |> 
  select(model,CHR,BP,P) |> 
  bind_rows(height_pop |> 
              mutate(model="Population") |> 
              select(model,CHR,BP,P))

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
       height = 10, units = "cm", dpi = 600)

# Panel C rG shrinkage

out= data.frame()

for (file in list.files("./trio-GWAS/results/") ){

  tmp <- read_lines(paste0("./trio-GWAS/results/", file)) 
  
  res =  str_remove_all( str_remove_all(tmp[grep("Genetic Correlation: ",tmp)], "Genetic Correlation: "), "[[)(]]")

  out = out |> 
    rbind(tibble("file"=file, "values"= res) |> 
                separate(values, into=c("rG","SE"), sep=" ") |> 
                mutate(rG=as.numeric(rG),
                       SE=as.numeric(SE)))
  
}


plot_rgs=out |> 
  mutate(Model = factor(ifelse(str_detect(file, "con"), "Trio", "Population")),
         Trait = factor(case_when(str_detect(file,"EA3") ~ "Educ. attain.",
                           str_detect(file,"smoke") ~ "Smoking",
                           str_detect(file,"locke") ~ "BMI",
                           str_detect(file,"wood") ~ "Height")),
         Outcome = ifelse(str_detect(file,"exam"), "Educational achievement","Height")) |> 
  rowwise() |> 
         mutate(lci= rG-1.96*SE,
                uci= rG+1.96*SE)

p2.2 = ggplot(plot_rgs, aes(x= rG, y=Trait,  fill=Model, shape= Model))+
  geom_vline(xintercept=0, linetype= 2, linewidth=1.2, colour="grey80")+
  geom_errorbarh(aes(xmin=lci,xmax=uci), height=0, linewidth=ci_lw, position= position_dodge(pd), colour="grey60" ) +
  geom_point(size=pt_sz, stroke=stroke_val,position= position_dodge(0.8), colour=stroke_col)+
  facet_grid(Outcome~.)+
  scale_fill_manual(values=c("#FF9933","#99CCFF"))+
  scale_shape_manual(values=c(21,22))+
  theme_few()+
  ylab("Externally GWASed trait")+
  xlab("rG with MoBa a) population \nand b) trio GWAS")+
    guides(fill=guide_legend(reverse = F),
           shape=guide_legend(reverse = F))+
  theme(legend.position = c(0.71,0.93),
        legend.direction = "horizontal",
        legend.background = element_rect(colour="grey80"))

ggsave(p2.2, filename = "./trio-GWAS/plots/rg_comp.tiff", device = "tiff", width = 14, 
       height = 10, units = "cm", dpi = 600)

# Panel X (supp) miami plot educational acheivement

exam_con = data.table::fread("./trio-GWAS/sumstats/exam_con_filtered_ldsc.txt")
exam_pop = data.table::fread("./trio-GWAS/sumstats/exam_pop_filtered_ldsc.txt")

exam_both = exam_con |> 
  mutate(model="Trio") |> 
  select(model,CHR,BP,P) |> 
  bind_rows(exam_pop |> 
              mutate(model="Population") |> 
              select(model,CHR,BP,P))

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
       height = 10, units = "cm", dpi = 600)

