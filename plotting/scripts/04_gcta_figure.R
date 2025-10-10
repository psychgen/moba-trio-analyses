# 04_gcta_figure.R

# The purpose of this script is to create a composite figure incorporating all 
# main text display elements for the trio GCTA component of the manuscript

source("./set_aes_pars.R")

library(tidyverse)
library(ggthemes)
library(ggsci)
library(patchwork)
library(readxl)
library(ggtext)



#load trio GCTA results

ests = read_xlsx("./trio-GCTA/results/estimates_afterreview_forbarplot.xlsx")
fits = read_xlsx("./trio-GCTA/results/fits_afterreview_trioGCTA.xlsx")
lrt = read_xlsx("./trio-GCTA/results/lrt_afterreview_trioGCTA.xlsx")



ests = ests |> 
  mutate(model = factor(model, levels = c("Full", "No covariance", "Direct")),
         variable = str_remove_all(variable, "\n"),
         outcome= factor(case_when( outcome=="Height" ~ "Height",
                                    outcome=="Edu. attain." ~ "Educational\nachievement",
                                    outcome=="Sleep Hrs." ~ "Sleep\nduration",
                                    outcome=="Depression symptoms" ~"Depressive\nsymptoms"),,
                         levels=c("Height","Educational\nachievement","Sleep\nduration","Depressive\nsymptoms")
                         ))

p4.1 <- ggplot(ests, aes(fill=variable, y=Percent_variance, x=model)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~ outcome, switch = "y") +  
  expand_limits(y = c(0,51)) +
  scale_fill_npg() +
  geom_text(aes(x = model, y = 47, label = best_fit_height), color = "black", size = 10) +
  geom_text(aes(x = model, y = 39, label = best_fit_exam), color = "black", size = 10) +
  geom_text(aes(x = model, y = 20, label = best_fit_sleep), color = "black", size = 10) +
  geom_text(aes(x = model, y = 15, label = best_fit_dep), color = "black", size = 10) +
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
  legend.position = "bottom",
  legend.direction = "horizontal") +
  labs(y = "Percent Variance Explained",
       title = "",
       fill = "Genetic Effects",
       x = "Model") +
  guides(fill = guide_legend(override.aes = list(size=10) ))

p4.1

p4.1sp = plot_spacer()/p4.1+
  plot_layout(heights=c(2.5,1.5))

ggsave("./trio-GCTA/plots/p4-1.tiff",p4.1sp, dpi=300, device= "tiff", width=26, height= 24, units ="cm" )

# Fits

fits_tests = fits |> 
  mutate(outcome= str_remove_all(outcome, "_95_|_hrs|_dep")) |> 
  left_join(lrt) |> 
  mutate(log10p = log10(`p..Chisq.`),
         model = factor(model, levels = c("full", "nocov", "direct","null"),
                        labels = c("Full", "No covariance","Direct","No genetic")),
           outcome= factor(case_when( outcome=="height" ~ "Height",
                                    outcome=="exam" ~ "Educational\nachievement",
                                    outcome=="sleep" ~ "Sleep\nduration",
                                    outcome=="smfq" ~"Depressive\nsymptoms"),,
                         levels=c("Height","Educational\nachievement","Sleep\nduration","Depressive\nsymptoms")),
         aic_y = ifelse(!is.na(log10p), log10p-10, log10(0.05))) 

p4.2 = ggplot(fits_tests , aes(x=model, y= log10p, group=outcome, fill=log10p ))+
  geom_line(colour="grey80", linewidth=1)+
  geom_point(shape=21, size= pt_sz*1.5, stroke=stroke_val, colour=stroke_col)+
  geom_label(aes(y=aic_y, label=aic), fill="white", nudge_x=0, size=2.6)+
  geom_label(label= "AIC value", x=1, y=-80, fill="white", size=2.6)+
  scale_fill_distiller("LRT p-value", palette="YlOrRd",breaks= log10(c(0.05,
                                                                       0.000000000000000000000000000005,
                                                                       0.000000000000000000000000000000000000000000000000000000000005,
                                                                       5e-90)),
                       labels= c("0.05",'5 x 10<sup>-30</sup>','5 x 10<sup>-60</sup>','5 x 10<sup>-90</sup>'))+
  scale_y_continuous("P-value of sequential LRT\nvs. next more complex model\n(log10 scale)",
                     breaks= log10(c(0.05,
                                     0.000000000000000000000000000005,
                                     0.000000000000000000000000000000000000000000000000000000000005,
                                     5e-90)),
                     labels= c("0.05",'5 x 10<sup>-30</sup>','5 x 10<sup>-60</sup>','5 x 10<sup>-90</sup>'))+
  facet_grid(.~outcome)+
  coord_cartesian(ylim = c(-120,0))+
  xlab("Model")+
  theme_few()+
  theme(axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        legend.position = "bottom",
        legend.direction = "horizontal")

p4.2sp = plot_spacer()/p4.2+
  plot_layout(heights=c(2.5,1.5))

p4.2sp  
ggsave("./trio-GCTA/plots/p4-2.tiff",p4.2sp, dpi=300, device= "tiff", width=26, height= 24, units ="cm" )

# Alternative fits
p4.2 = ggplot(fits_tests |> 
                mutate(flag=ifelse(model=="Full"&outcome!="Educational\nachievement", 1, 0)) |> 
                filter(!flag==1), aes(x=model, y= log10p,group=outcome,  fill=log10p ))+
  geom_line(colour="grey80", linewidth=1,alpha=1, position=position_dodge(0.2) )+
  geom_point(size= pt_sz, stroke=stroke_val, colour=stroke_col,alpha=1, position=position_dodge(0.2), shape=21)+
  facet_grid(.~outcome, scales="free")+
  # geom_label(aes(y=aic_y, label=aic), fill="white", nudge_x=0)+
  #geom_label(label= "AIC value", x=1, y=-80, fill="white")+
  scale_fill_distiller("LRT p-value", palette="YlOrRd",breaks= log10(c(0.05,
                                                                       5e-90)),
                       labels= c("0.05",'5 x 10<sup>-90</sup>'))+
  scale_y_continuous("P-value of sequential LRT\nvs. next more complex model\n(log10 scale)",
                     breaks= log10(c(0.05,
                                     0.000000000000000000000000000005,
                                     0.000000000000000000000000000000000000000000000000000000000005,
                                     5e-90)),
                     labels= c("0.05",'5 x 10<sup>-30</sup>','5 x 10<sup>-60</sup>','5 x 10<sup>-90</sup>'))+
  coord_cartesian(ylim = c(-100,0))+
  xlab("Model")+
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_markdown(),
        legend.text = element_markdown(),
        legend.position = "bottom",
        legend.direction = "horizontal")+
  guides( shape = "none",
          fill = guide_colorbar(reverse=T))


p4.2sp = plot_spacer()/p4.2+
  plot_layout(heights=c(3,1))

p4.2sp  
ggsave("./trio-GCTA/plots/p4-2.tiff",p4.2sp, dpi=300, device= "tiff", width=26, height= 24, units ="cm" )


# Precision

ests_ses= ests |> 
  filter(!is.na(best_fit_height)|!is.na(best_fit_exam)) |> 
  mutate(se = case_when(value==0.256 ~ 0.016,
                        value==0.039 ~ 0.013,
                        value==0.061 ~ 0.013,
                        value==0.316 ~ 0.026,
                        value==0.044 ~ 0.023,
                        value==0.029 ~ 0.022)) |>  ##Manually adding standard errors for now
  filter(!str_detect(variable, "Covariance")) 


p4.3<- ggplot(ests_ses, aes(x=variable,y=value,fill=variable, colour= variable, group=outcome, shape = outcome)) +
  geom_hline(yintercept=0, linetype=2, colour=null_col, linewidth=null_lw)+  
  geom_errorbar(aes(ymin=value-se,ymax=value+se),position=position_dodge(pd), linewidth=ci_lw, width=0) +
  geom_point(size=pt_sz, position=position_dodge(pd), colour= stroke_col, stroke=stroke_val)+
  scale_shape_manual(values=c(21,22))+
  scale_fill_npg() +
  scale_colour_npg()+
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  scale_y_continuous("Proportion of variance explained \nin best-fitting model")+
  scale_x_discrete("Source of genetic variation")+
  guides(fill="none", colour="none",
         shape = guide_legend(title="",override.aes =  list(size = 5, alpha=1)))

p4.3

p4.23 =(p4.3+p4.2)
p4.23
ggsave("./trio-GCTA/plots/p4-23_new.tiff",p4.23, dpi=dpi, device= "tiff", width=22, height= 10, units ="cm" )


