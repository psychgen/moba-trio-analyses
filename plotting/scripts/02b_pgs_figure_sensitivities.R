# 02_pgs_figure.R

# The purpose of this script is to create a composite figure incorporating all 
# main text display elements for the PGS component of the manuscript

source("./set_aes_pars.R")

library(tidyverse)
library(ggthemes)
library(ggsci)
library(patchwork)


# Read in the results files

load("./pgs-analyses/triopgs_res_unrelated.RData")
trio_res_un = trio_res
load("./pgs-analyses/triopgs_res_altclust.RData")
trio_res_fid = trio_res
rm(trio_res)


# Child ests attenuation

trio_res_un <- trio_res_un |>
  mutate(Outcome=factor(case_when(Outcome =="exam" ~ "Educational\nachievement",
                           Outcome =="height" ~ "Height",
                           Outcome =="smfq_dep" ~ "Depressive\nsymptoms",
                           Outcome =="sleep_hrs" ~ "Sleep\nduration"),
                        levels=c("Height","Educational\nachievement","Sleep\nduration","Depressive\nsymptoms")))

p3.1<- ggplot(trio_res_un %>% 
              filter(Member=="child"), aes(x=est.std,y=Member,fill=Type,colour=Type, shape=Type)) +
  geom_vline(xintercept=0, linetype= 2, linewidth=null_lw, colour=null_col)+
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper), height=0, linewidth=ci_lw, position= position_dodge(pd) ) +
  geom_point(size=pt_sz, stroke=stroke_val,position= position_dodge(pd),  colour=stroke_col)+
  facet_grid(PGS~Outcome,scales = "free",switch="y", space="free")+
  scale_fill_manual(values=c("#D55E00","#0072B2"))+
  scale_colour_manual(values=c("#D55E00","#0072B2"))+
  scale_shape_manual(values=c(21,22))+
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t =  0, r = 5, b = 0, l = 0), angle=90),
        axis.title.x = element_text(margin =  margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        strip.text.y.left = element_text(angle=0),
        legend.position = "top",
        legend.direction = "horizontal") +
  scale_x_continuous("Standardised effect of\nchild PGS on outcome", breaks=c(0,0.1,0.2,0.3,0.4), labels = c("0.0","0.1","0.2","0.3","0.4")  )+
  #coord_cartesian(xlim= c(-0.02, 0.32))+
  scale_y_discrete("PGS trait")+
  guides(fill=guide_legend(title="Model"),
         colour=guide_legend(title="Model"),
         shape=guide_legend(title="Model"))

p3.1


p3.2<- ggplot(trio_res_un %>% 
              filter(Member!="child"), aes(x=est.std,y=PGS,fill=Member, colour=Member)) +
  geom_vline(xintercept=0, linetype= 2, linewidth=null_lw, colour=null_col)+
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper), height=0, linewidth=ci_lw, position= position_dodge(pd) ) +
  geom_point(size=pt_sz, stroke=stroke_val,position= position_dodge(pd), shape=23, colour=stroke_col)+
  facet_grid(Outcome~.,scales = "free",space="free")+
  scale_colour_brewer(type=seq, palette = "Accent", direction=-1)+
  scale_fill_brewer( palette = "Accent", direction=-1)+
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t =  0, r = 0, b = 0, l = 0), angle=90),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        #axis.text.y  = element_blank(),
        #axis.ticks.y  = element_blank(),
        strip.text.y = element_text(angle=0),
        legend.position = "top",
        legend.direction = "horizontal")+
  scale_x_continuous("Standardised effect of\nparental PGS on height\n in trio-PGS model",
                     breaks = c(-0.1,0.0,0.1),labels = c("-0.1","0.0","0.1"))+
  coord_cartesian(xlim= c(-0.1, 0.1))+
  scale_y_discrete("PGS trait")+
  guides(fill=guide_legend(title="Parent"),
         colour=guide_legend(title="Parent"))
p3.2


p3.1sp = plot_spacer()/p3.1+
  plot_layout(heights=c(2.5,1.5))
p3.2sp = plot_spacer()+p3.2+
  plot_layout(widths=c(3.3,0.7))

p3.12= (p3.1+p3.2) +
  plot_layout(widths=c(3,1))

#ggsave("./pgs-analyses/plots/fig3.12.tiff", p3.12, device="tiff", width=26,height=12,units="cm",dpi=600,bg="white")
ggsave("./pgs-analyses/plots/fig3.1_unrelated.tiff", p3.1sp, device="tiff", width=26,height=24,units="cm",dpi=300,bg="white")
ggsave("./pgs-analyses/plots/fig3.2_unrelated.tiff", p3.2sp, device="tiff", width=26,height=24,units="cm",dpi=300,bg="white")




# Child ests attenuation

trio_res_fid <- trio_res_fid |>
  mutate(Outcome=factor(case_when(Outcome =="exam" ~ "Educational\nachievement",
                                  Outcome =="height" ~ "Height",
                                  Outcome =="smfq_dep" ~ "Depressive\nsymptoms",
                                  Outcome =="sleep_hrs" ~ "Sleep\nduration"),
                        levels=c("Height","Educational\nachievement","Sleep\nduration","Depressive\nsymptoms")))

p3.1<- ggplot(trio_res_fid %>% 
                filter(Member=="child"), aes(x=est.std,y=Member,fill=Type,colour=Type, shape=Type)) +
  geom_vline(xintercept=0, linetype= 2, linewidth=null_lw, colour=null_col)+
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper), height=0, linewidth=ci_lw, position= position_dodge(pd) ) +
  geom_point(size=pt_sz, stroke=stroke_val,position= position_dodge(pd),  colour=stroke_col)+
  facet_grid(PGS~Outcome,scales = "free",switch="y", space="free")+
  scale_fill_manual(values=c("#D55E00","#0072B2"))+
  scale_colour_manual(values=c("#D55E00","#0072B2"))+
  scale_shape_manual(values=c(21,22))+
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t =  0, r = 5, b = 0, l = 0), angle=90),
        axis.title.x = element_text(margin =  margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        strip.text.y.left = element_text(angle=0),
        legend.position = "top",
        legend.direction = "horizontal") +
  scale_x_continuous("Standardised effect of\nchild PGS on outcome", breaks=c(0,0.1,0.2,0.3,0.4), labels = c("0.0","0.1","0.2","0.3","0.4")  )+
  #coord_cartesian(xlim= c(-0.02, 0.32))+
  scale_y_discrete("PGS trait")+
  guides(fill=guide_legend(title="Model"),
         colour=guide_legend(title="Model"),
         shape=guide_legend(title="Model"))

p3.1


p3.2<- ggplot(trio_res_fid %>% 
                filter(Member!="child"), aes(x=est.std,y=PGS,fill=Member, colour=Member)) +
  geom_vline(xintercept=0, linetype= 2, linewidth=null_lw, colour=null_col)+
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper), height=0, linewidth=ci_lw, position= position_dodge(pd) ) +
  geom_point(size=pt_sz, stroke=stroke_val,position= position_dodge(pd), shape=23, colour=stroke_col)+
  facet_grid(Outcome~.,scales = "free",space="free")+
  scale_colour_brewer(type=seq, palette = "Accent", direction=-1)+
  scale_fill_brewer( palette = "Accent", direction=-1)+
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t =  0, r = 0, b = 0, l = 0), angle=90),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        #axis.text.y  = element_blank(),
        #axis.ticks.y  = element_blank(),
        strip.text.y = element_text(angle=0),
        legend.position = "top",
        legend.direction = "horizontal")+
  scale_x_continuous("Standardised effect of\nparental PGS on height\n in trio-PGS model",
                     breaks = c(-0.1,0.0,0.1),labels = c("-0.1","0.0","0.1"))+
  coord_cartesian(xlim= c(-0.1, 0.1))+
  scale_y_discrete("PGS trait")+
  guides(fill=guide_legend(title="Parent"),
         colour=guide_legend(title="Parent"))
p3.2


p3.1sp = plot_spacer()/p3.1+
  plot_layout(heights=c(2.5,1.5))
p3.2sp = plot_spacer()+p3.2+
  plot_layout(widths=c(3.3,0.7))

p3.12= (p3.1+p3.2) +
  plot_layout(widths=c(3,1))

#ggsave("./pgs-analyses/plots/fig3.12.tiff", p3.12, device="tiff", width=26,height=12,units="cm",dpi=600,bg="white")
ggsave("./pgs-analyses/plots/fig3.1_fid.tiff", p3.1sp, device="tiff", width=26,height=24,units="cm",dpi=300,bg="white")
ggsave("./pgs-analyses/plots/fig3.2_fid.tiff", p3.2sp, device="tiff", width=26,height=24,units="cm",dpi=300,bg="white")








