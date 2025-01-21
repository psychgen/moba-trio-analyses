# 03_trioPGS.R

# In this script we run the trio PGS analyses for both outcomes

library(tidyverse)
library(lavaan)
library(patchwork)

load("./data/dat_final.RData")


run_trioPGS <- function(out, covs, scores, data){
  
  
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
                     cluster="m_id")

    baseline_est <- standardizedsolution(basefit, type="std.nox") %>% 
        filter (op == "~" ) %>% 
      mutate(Type="Child-only")
    
    triomod = paste0(out, " ~ ",chscore, " + ", mscore, " + ", fscore, " + ", paste0(covs,collapse = " + "))
      
    triofit= lavaan::sem(triomod, 
                         data=data, 
                         estimator="MLR",
                         se="robust",
                         cluster="m_id") 
    
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

pgs_names <- c("height2",
               "bmi2",
               "ea3",
               "smoke2" )

pgs_display <- c("Height","BMI", "Edu. attain.", "Ever smoked")


trio_ex = run_trioPGS(out="exam", 
                      covs="SEX", #Age excluded for exam analyses because of zero variance
                      scores=paste0(pgs_names, "_0.95_pgs_res"), 
                      data=dat_final)

# Do some post-processing of results

trio_res_ex = trio_ex %>% 
  separate(rhs, into = c("Member", "PGS"), sep="_") %>% 
  mutate(PGS = factor(str_remove_all(PGS,"_0.95_pgs_res"), levels= pgs_names, labels = pgs_display )) %>% 
  select(Type, Member,PGS,est.std,se,pvalue,ci.lower,ci.upper) %>%
  drop_na(PGS)

# Save out

save(trio_res_ex, file = "./output/triopgs_edu.RData")
         

# Plot attenuation of child only estimates

p2<- ggplot(trio_res_ex %>% 
              filter(Member=="child"), aes(x=est.std,y=fct_rev(Type),fill=Type,colour=Type, shape=Type)) +
  geom_vline(xintercept=0, linetype=2, colour="grey80", size=1.2)+  
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper),position=position_dodge(0.8), size=0.6, height=0) +
  geom_point(size=3, position=position_dodge(0.8),  colour= "black", stroke=0.6)+
  scale_shape_manual(values = c(22,23))+
  scale_colour_brewer(type=seq, palette = "Paired", direction=-1)+
  scale_fill_brewer( palette = "Paired", direction=-1)+
  theme(text = element_text(size=12),
        axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey10"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_line(colour="grey90"),
        strip.text.y.left =   element_text(size=12, face="bold", colour="white", angle=0),
        legend.title = element_text(size =12),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_rect(fill="grey10")) +
  facet_grid(PGS~., scales = "free_y", switch= "y")+
  scale_x_continuous("Standardised effect of PGS on exam scores" )+
  coord_cartesian(xlim= c(-0.2, 0.3))+
  scale_y_discrete("PGS trait")

p2

ggsave("./output/plots/fig2.tiff", device="tiff", width=14,height=9,units="cm",dpi=320,bg="white")

# Plot effect of parental PGS

trio_res_ex =  trio_res_ex %>% 
  mutate(Member = case_when(Member=="mother"~"Mother",
                            Member=="father"~"Father",
                            TRUE~"child"))

p3<- ggplot(trio_res_ex %>% 
              filter(Member!="child"), aes(x=est.std,y=fct_rev(Member),fill=Member, colour=Member)) +
  geom_vline(xintercept=0, linetype=2, colour="grey80", size=1.2)+  
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper),position=position_dodge(0.8), size=0.6, height=0) +
  geom_point(size=3, position=position_dodge(0.8), shape= 21, colour= "black", stroke=0.6)+
  scale_colour_brewer(type=seq, palette = "Accent", direction=-1)+
  scale_fill_brewer( palette = "Accent", direction=-1)+
  theme(text = element_text(size=12),
        axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey10"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_line(colour="grey90"),
        strip.text.y.left =   element_blank(),
        legend.title = element_text(size =12),
       # legend.position = "bottom",
      #  legend.direction = "horizontal",
        strip.background = element_blank()) +
  facet_grid(PGS~., scales = "free_y", switch= "y")+
  scale_x_continuous("Standardised effect of PGS on exam scores" )+
  coord_cartesian(xlim= c(-0.2, 0.3))+
  scale_y_discrete("PGS trait")
p3

ggsave("./output/plots/fig3.tiff", device="tiff", width=14,height=9,units="cm",dpi=320,bg="white")


save(p2,p3, file= "./scratch/pgsplots.RData")


# Complete & plot the trioPGS for height also


trioh = run_trioPGS(out="height", 
                    covs=c("SEX","Age_height"),
                    scores=paste0(pgs_names, "_0.95_pgs_res"),
                    data=dat_final)

# Do some post-processing of results

trio_res_ht = trioh %>% 
  separate(rhs, into = c("Member", "PGS"), sep="_") %>% 
  mutate(PGS = factor(str_remove_all(PGS,"_0.95_pgs_res"), levels= pgs_names, labels = pgs_display )) %>% 
  select(Type, Member,PGS,est.std,se,pvalue,ci.lower,ci.upper)%>%
  drop_na(PGS)

# Save out

save(trio_res_ht, file = "./output/triopgs_hei.RData")


# Plot attenuation of child only estimates

p5<- ggplot(trio_res_ht %>% 
              filter(Member=="child"), aes(x=est.std,y=fct_rev(Type),fill=Type,colour=Type, shape=Type)) +
  geom_vline(xintercept=0, linetype=2, colour="grey80", size=1.2)+  
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper),position=position_dodge(0.8), size=0.6, height=0) +
  geom_point(size=3, position=position_dodge(0.8),  colour= "black", stroke=0.6)+
  scale_shape_manual(values = c(22,23))+
  scale_colour_brewer(type=seq, palette = "Paired", direction=-1)+
  scale_fill_brewer( palette = "Paired", direction=-1)+
  theme(text = element_text(size=12),
        axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey10"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_line(colour="grey90"),
        strip.text.y.left =   element_text(size=12, face="bold", colour="white", angle=0),
        legend.title = element_text(size =12),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_rect(fill="grey10")) +
  facet_grid(PGS~., scales = "free_y", switch= "y")+
  scale_x_continuous("Standardised effect of PGS on height" )+
  coord_cartesian(xlim= c(-0.2, 0.4))+
  scale_y_discrete("PGS trait")

p5

ggsave("./output/plots/fig5.tiff", device="tiff", width=14,height=9,units="cm",dpi=320,bg="white")

# Plot effect of parental PGS

trio_res_ht = trio_res_ht %>% 
  mutate(Member = case_when(Member=="mother"~"Mother",
                            Member=="father"~"Father",
                            TRUE~"child"))

p6<- ggplot(trio_res_ht %>% 
              filter(Member!="child"), aes(x=est.std,y=fct_rev(Member),fill=Member, colour=Member)) +
  geom_vline(xintercept=0, linetype=2, colour="grey80", size=1.2)+  
  geom_errorbarh(aes(xmin=ci.lower,xmax=ci.upper),position=position_dodge(0.8), size=0.6, height=0) +
  geom_point(size=3, position=position_dodge(0.8), shape= 21, colour= "black", stroke=0.6)+
  scale_colour_brewer(type=seq, palette = "Accent", direction=-1)+
  scale_fill_brewer( palette = "Accent", direction=-1)+
  theme(text = element_text(size=12),
        axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y  = element_blank(),
        axis.ticks.y  = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey10"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_line(colour="grey90"),
        strip.text.y.left =   element_blank(),
        legend.title = element_text(size =12),
        # legend.position = "bottom",
        #  legend.direction = "horizontal",
        strip.background = element_blank()) +
  facet_grid(PGS~., scales = "free_y", switch= "y")+
  scale_x_continuous("Standardised effect of PGS on height" )+
  coord_cartesian(xlim= c(-0.2, 0.4))+
  scale_y_discrete("PGS trait")
p6

ggsave("./output/plots/fig6.tiff", device="tiff", width=14,height=9,units="cm",dpi=320,bg="white")


p5 + p6 + plot_layout(axis_titles = "collect" )

save(p5,p6, file= "./scratch/pgs_supp_plots.RData")

