# 02_ptdt_heightex.R

# In this script we select trios on the basis of "extreme" height among children, 
# dropping those with "extremely" tall parents, attaching siblings as comparators,
# and run a ptdt analysis

library(tidyverse)

# Set seed to keep the sampling consistent

set.seed(15246)

load("./data/dat_final.RData")

# Filter to select the ptdt sample for the height analysis

dat_ptdt <- dat_final %>% 
  filter(!is.na(height),
         if_all(matches("pgs_res"), ~!is.na(.))) %>% 
  group_by(m_id) %>%
  slice_sample(n=1)%>%
  rowwise() %>% 
  mutate(pheight = mean(c(mheight,fheight))) %>% 
  ungroup() %>% 
  mutate(ext_ch = ifelse(height> (mean(height,na.rm=T)+2*sd(height,na.rm=T)), "Yes","No" ),
         ext_p = ifelse(pheight> (mean(pheight,na.rm=T)+2*sd(pheight,na.rm=T)), "Yes","No" )) %>% 
  filter(ext_ch=="Yes",
         ext_p=="No")

# Add in sibling PGS

load("./data/pgs_prcsd.RData")

pgs_sibs <- dat_final %>% 
  filter(!IID %in% dat_ptdt$IID,
         m_id %in% dat_ptdt$m_id,
         if_all(matches("pgs_res"), ~!is.na(.))) %>% 
  select(m_id,f_id,FID,matches("child")) %>% 
  rename_at(vars(matches("child")), ~str_replace_all(.,"child","sib")) %>% 
  group_by(m_id,f_id,FID) %>% 
  slice_sample(n=1) %>% 
  ungroup()

dat_ptdt <- dat_ptdt %>% 
  left_join(pgs_sibs)

# Create pTD score for each invidual, by:
## Calculating mid-parental average (rowwise)
## Calculating deviation (rowwise)
## Scaling by mid-parental SD

make_ptd <- function(scores, data){
  
  
  alldat <- data %>% 
    select(preg_id,BARN_NR)
  
  for(score in scores){
    
    chscore=paste0("child_",score)
    mscore=paste0("mother_",score)
    fscore=paste0("father_",score)
    siscore=paste0("sib_",score)
    
    tmpdat <- data %>% 
      drop_na(.data[[chscore]],.data[[mscore]],.data[[fscore]]) %>% 
      rowwise() %>% 
      mutate(mp.pgs=(.data[[mscore]]+.data[[fscore]])/2 ,
             ptd.unsc= .data[[chscore]]-mp.pgs ) %>% 
      ungroup() %>% 
      mutate(ptd.sc=ptd.unsc/sd(mp.pgs))
    
    tmpdat[[str_remove_all(paste0("ptd_",score), "_0.95_pgs_res")]] <- tmpdat$ptd.sc
    
    tmpdat2 <- data %>% 
      drop_na(.data[[siscore]],.data[[mscore]],.data[[fscore]]) %>% 
      rowwise() %>% 
      mutate(mp.pgs=(.data[[mscore]]+.data[[fscore]])/2 ,
             ptd.unsc= .data[[siscore]]-mp.pgs ) %>% 
      ungroup() %>% 
      mutate(ptd.sib=ptd.unsc/sd(tmpdat$mp.pgs))
    
    tmpdat2[[str_remove_all(paste0("ptd.sib_",score), "_0.95_pgs_res")]] <- tmpdat2$ptd.sib
    
    alldat = alldat %>% 
      left_join(tmpdat %>% 
                  select(preg_id,BARN_NR,matches("ptd_"))) %>% 
      left_join(tmpdat2 %>% 
                  select(preg_id,BARN_NR,matches("ptd.sib_")))
    
    
  }
  
  return(alldat)
  
} 

# Apply the function 

ptdat = dat_ptdt %>% 
  left_join(make_ptd(paste0(pgs_names,"_0.95_pgs_res"),dat_ptdt))


# Get the results, with formatted PGS names


pgs_names <- c("height2",
               "bmi2",
               "ea3",
               "smoke2" )

pgs_display <- c("Height","BMI", "Edu. attain.", "Ever smoked")


# Small function to get means for each ptd column, adjusted for age and sex

lem <- function(x){ 
  tmp <- lm(x~1+SEX+scale(Age_height), data= ptdat)
  return(broom::tidy(tmp))}

# Create the plotting dataset

ptdres <- ptdat%>%  
  summarise(across(matches("ptd_"),.fns= list(reg= lem), .names = "{.col}")) %>% 
  pivot_longer(cols=everything(), names_to = "ptd_pheno") %>% 
  mutate(mean= value$estimate,
         se=value$std.error,
         p=value$p.value,
         fdrp = p.adjust(p, method="fdr"),
         ptd_pheno = factor(ptd_pheno,levels=paste0("ptd_",pgs_names), labels=pgs_display),
         Group = "+2Ds height") %>% 
  bind_rows(ptdat%>%  
              summarise(across(matches("ptd.sib_"),.fns= list(reg= lem), .names = "{.col}")) %>% 
              pivot_longer(cols=everything(), names_to = "ptd_pheno") %>% 
              mutate(mean= value$estimate,
                     se=value$std.error,
                     p=value$p.value,
                     fdrp = p.adjust(p, method="fdr"),
                     ptd_pheno = factor(ptd_pheno,levels=paste0("ptd.sib_",pgs_names), labels=pgs_display),
                     Group= "Siblings"))  %>% 
  filter(value=="(Intercept)")
  

# Save out the plotting dataset

save(ptdres, file = "./output/ptd_height.RData")

# Create a standard format plot for the ptd results

p1<- ggplot(ptdres, aes(x=Group,y=mean,fill=Group,colour=Group)) +
  geom_hline(yintercept=0, linetype=2, colour="grey80", size=1.2)+  
  geom_errorbar(aes(ymin=mean-1.96*se,ymax=mean+1.96*se),position=position_dodge(0.8), linewidth=0.6, width=0) +
  geom_point(size=3, position=position_dodge(0.8), shape= 21, colour= "black", stroke=0.6)+
  scale_colour_brewer(type="qual", palette = "Set2")+
  scale_fill_brewer(type="qual", palette = "Set2")+
  theme(text = element_text(size=12),
        axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        #axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x  = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey10"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_line(colour="grey90"),
        strip.text = element_text(size=12, face="bold", colour="white"),
        legend.title = element_text(size =12),
        strip.background = element_rect(fill="grey10")) +
  facet_grid(.~ptd_pheno, scales = "free_x")+
  scale_y_continuous("Relative over-/under-transmission\n of alleles (mid-parent PGS SDs)")+
  scale_x_discrete("PGS trait and group (probands ascertained\n for extreme height vs. their siblings)")

p1

# Save the pdt plot, both as a Tiff and an RData file

ggsave("./output/plots/fig1.tiff", device="tiff", width=28,height=9,units="cm",dpi=320,bg="white")

save(p1, file="./scratch/ptdt_plot.RData")

