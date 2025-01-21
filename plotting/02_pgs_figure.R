# 02_pgs_figure.R

# The purpose of this script is to create a composite figure incorporating all 
# main text display elements for the PGS component of the manuscript

source("./set_aes_pars.R")

library(tidyverse)
library(ggthemes)
library(ggsci)
library(patchwork)


# Read in the results files

load("./pgs-analyses/ptd_edu.RData")
ptd_res_edu = ptdres
load("./pgs-analyses/ptd_height.RData")
ptd_res_ht = ptdres
load("./pgs-analyses/triopgs_edu.RData")
load("./pgs-analyses/triopgs_hei.RData")

# Height

# Child ests attenuation

p3.1<- ggplot(trio_res_ht %>% 
              filter(Member=="child"), aes(x=est.std,y=Member,fill=Type,colour=Type, shape=Type)) +
  geom_vline(xintercept=0, linetype= 2, linewidth=null_lw, colour=null_col)+
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper), height=0, linewidth=ci_lw, position= position_dodge(pd) ) +
  geom_point(size=pt_sz, stroke=stroke_val,position= position_dodge(pd),  colour=stroke_col)+
  facet_grid(PGS~.,scales = "free_y",switch="y")+
  scale_fill_manual(values=c("#D55E00","#0072B2"))+
  scale_colour_manual(values=c("#D55E00","#0072B2"))+
  scale_shape_manual(values=c(21,22))+
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        strip.text.y.left = element_text(angle=0),
        legend.position = "top",
        legend.direction = "horizontal") +
  scale_x_continuous("Standardised effect of\nchild PGS on height" )+
  coord_cartesian(xlim= c(-0.02, 0.42))+
  scale_y_discrete("PGS trait")+
  guides(fill=guide_legend(title="Model"),
         colour=guide_legend(title="Model"),
         shape=guide_legend(title="Model"))

p3.1


p3.2<- ggplot(trio_res_ht %>% 
              filter(Member!="child"), aes(x=est.std,y=Type,fill=Member, colour=Member)) +
  geom_vline(xintercept=0, linetype= 2, linewidth=null_lw, colour=null_col)+
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper), height=0, linewidth=ci_lw, position= position_dodge(pd) ) +
  geom_point(size=pt_sz, stroke=stroke_val,position= position_dodge(pd), shape=23, colour=stroke_col)+
  facet_grid(PGS~.,scales = "free_y",switch="y")+
  scale_colour_brewer(type=seq, palette = "Accent", direction=-1)+
  scale_fill_brewer( palette = "Accent", direction=-1)+
  theme_few()+
  theme(axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        #legend.title = element_text(size =12),
        strip.text.y.left =   element_blank(),
        legend.position = "top",
        legend.direction = "horizontal") +
  scale_x_continuous("Standardised effect of\nparental PGS on height\n in trio-PGS model",
                     breaks = c(0.0,0.1,0.2))+
  coord_cartesian(xlim= c(-0.02, 0.21))+
  scale_y_discrete("")+
  guides(fill=guide_legend(title="Parent"),
         colour=guide_legend(title="Parent"))
p3.2


p3.12= (p3.1+p3.2) +
  plot_layout(widths=c(2,1))

ggsave("./pgs-analyses/plots/fig3.12.tiff", p3.12, device="tiff", width=22,height=12,units="cm",dpi=600,bg="white")


#ptdt

# Create a standard format plot for the ptd results

p3.3<- ggplot(ptd_res_ht, aes(x=Group,y=mean,fill=Group,colour=Group)) +
  geom_hline(yintercept=0, linetype=2, colour=null_col, linewidth=null_lw)+  
  geom_errorbar(aes(ymin=mean-1.96*se,ymax=mean+1.96*se),position=position_dodge(pd), linewidth=ci_lw, width=0) +
  geom_point(size=pt_sz, position=position_dodge(pd), shape= 21, colour= stroke_col, stroke=stroke_val)+
  scale_colour_brewer(type="qual", palette = "Set2")+
  scale_fill_brewer(type="qual", palette = "Set2")+
  facet_grid(.~ptd_pheno, scales = "free_x", switch = "x")+
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x  = element_blank(),
        axis.ticks.x  = element_blank()) +
  scale_y_continuous("Relative over-/under-transmission\n of alleles (mid-parent PGS SDs)")+
  scale_x_discrete("PGS trait and group (probands ascertained\n for extreme height & their siblings)")

p3.3


p3.123= p3.12/p3.3 +
  plot_layout(heights=c(2,1))

ggsave("./pgs-analyses/plots/fig3.123.tiff", p3.123, device="tiff", width=22,height=16,units="cm",dpi=600,bg="white")

ggsave("./pgs-analyses/plots/fig3.3.tiff", p3.3, device="tiff", width=22,height=10,units="cm",dpi=600,bg="white")


# Educational achievement


# Child ests attenuation

p3.4<- ggplot(trio_res_ex %>% 
                filter(Member=="child"), aes(x=est.std,y=Member,fill=Type,colour=Type, shape=Type)) +
  geom_vline(xintercept=0, linetype= 2, linewidth=null_lw, colour=null_col)+
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper), height=0, linewidth=ci_lw, position= position_dodge(pd) ) +
  geom_point(size=pt_sz, stroke=stroke_val,position= position_dodge(pd),  colour=stroke_col)+
  facet_grid(PGS~.,scales = "free_y",switch="y")+
  scale_fill_manual(values=c("#D55E00","#0072B2"))+
  scale_colour_manual(values=c("#D55E00","#0072B2"))+
  scale_shape_manual(values=c(21,22))+
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        strip.text.y.left = element_text(angle=0),
        legend.position = "top",
        legend.direction = "horizontal") +
  scale_x_continuous("Standardised effect of child\nPGS on educational achievement" )+
  coord_cartesian(xlim= c(-0.1, 0.42))+
  scale_y_discrete("PGS trait")+
  guides(fill=guide_legend(title="Model"),
         colour=guide_legend(title="Model"),
         shape=guide_legend(title="Model"))

p3.4


p3.5<- ggplot(trio_res_ex %>% 
                filter(Member!="child"), aes(x=est.std,y=Type,fill=Member, colour=Member)) +
  geom_vline(xintercept=0, linetype= 2, linewidth=null_lw, colour=null_col)+
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper), height=0, linewidth=ci_lw, position= position_dodge(pd) ) +
  geom_point(size=pt_sz, stroke=stroke_val,position= position_dodge(pd), shape=23, colour=stroke_col)+
  facet_grid(PGS~.,scales = "free_y",switch="y")+
  scale_colour_brewer(type=seq, palette = "Accent", direction=-1)+
  scale_fill_brewer( palette = "Accent", direction=-1)+
  theme_few()+
  theme(axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        #legend.title = element_text(size =12),
        strip.text.y.left =   element_blank(),
        legend.position = "top",
        legend.direction = "horizontal") +
  scale_x_continuous("Standardised effect of parental\nPGS on educational achievement\nin trio-PGS model",
                     breaks = c(-0.1,0.0,0.1,0.2))+
  coord_cartesian(xlim= c(-0.1, 0.21))+
  scale_y_discrete("")+
  guides(fill=guide_legend(title="Parent"),
         colour=guide_legend(title="Parent"))
p3.5


p3.45= (p3.4+p3.5) +
  plot_layout(widths=c(2,1))

ggsave("./pgs-analyses/plots/fig3.45.tiff", p3.45, device="tiff", width=22,height=12,units="cm",dpi=600,bg="white")


#ptdt

# Create a standard format plot for the ptd results

ptd_res_edu = ptd_res_edu |> 
  mutate(Group = ifelse(Group==">90th percentile\nexam performance",
                        ">90th percentile\neducational achievement", 
                        Group))

p3.6<- ggplot(ptd_res_edu, aes(x=Group,y=mean,fill=Group,colour=Group)) +
  geom_hline(yintercept=0, linetype=2, colour=null_col, linewidth=null_lw)+  
  geom_errorbar(aes(ymin=mean-1.96*se,ymax=mean+1.96*se),position=position_dodge(pd), linewidth=ci_lw, width=0) +
  geom_point(size=pt_sz, position=position_dodge(pd), shape= 21, colour= stroke_col, stroke=stroke_val)+
  scale_colour_brewer(type="qual", palette = "Set2")+
  scale_fill_brewer(type="qual", palette = "Set2")+
  facet_grid(.~ptd_pheno, scales = "free_x", switch = "x")+
  theme_few()+
  theme(axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x  = element_blank(),
        axis.ticks.x  = element_blank()) +
  scale_y_continuous("Relative over-/under-transmission\n of alleles (mid-parent PGS SDs)")+
  scale_x_discrete("PGS trait and group (probands ascertained\n for high educational achivement & their siblings)")

p3.6


p3.456= p3.45/p3.6 +
  plot_layout(heights=c(2,1))

ggsave("./pgs-analyses/plots/fig3.456.tiff", p3.456, device="tiff", width=22,height=16,units="cm",dpi=600,bg="white")

ggsave("./pgs-analyses/plots/fig3.6.tiff", p3.6, device="tiff", width=22,height=10,units="cm",dpi=600,bg="white")






