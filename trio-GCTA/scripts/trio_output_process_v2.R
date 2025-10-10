#author: Laura Hegemann 
#process results for downstream use



library(tidyverse)
library(ggsci)
library(xlsx)

setwd("//ess01/p471/data/durable/people/Laura_H/Cluster_backup/trio_gcta_qc/results/") #LJH changed to user-invariant filepath


results <- list.files(pattern = "afterreview2.csv")


results

exam <- read.csv(results[1]) %>%
  mutate(outcome = "exam")
height <- read.csv(results[2]) %>%
  mutate(outcome = "height")
sleep <- read.csv(results[3]) %>%
  mutate(outcome = "sleep")
dep <- read.csv(results[4]) %>%
  mutate(outcome = "depression")


result_df <- rbind(exam,height,sleep,dep) %>%
  mutate(best_fit_height = case_when((outcome == "height") & model == "nocov" ~ "*", 
                              TRUE ~" "),
         best_fit_exam = case_when((outcome == "exam") & model == "full" ~ "*", 
                                     TRUE ~" "), 
         best_fit_sleep = case_when((outcome == "sleep") & model == "nocov" ~ "*", 
                                   TRUE ~" "), 
         best_fit_dep = case_when((outcome == "depression") & model == "nocov" ~ "*", 
                                    TRUE ~" ")) %>%
  #filter(best_fit == "yes") %>%
  filter(model != "null") %>%
  filter(variable != "mp") %>%
  filter(variable != "e") %>%
  mutate(value = round(value,3),
         Percent_variance = value*100, 
         model = case_when(model == "full" ~ "Full", 
                                  model == "nocov" ~ "No covariance", 
                                  model == "direct" ~ "Direct"),
         model = factor(model, levels = c("Full", "No covariance", "Direct")),
         outcome = case_when(outcome == "exam" ~ "Edu. attain.",
                             outcome == "height" ~ "Height", 
                             outcome == "sleep" ~"Sleep Hrs.", 
                             outcome == "depression" ~ "Depression symptoms"),
         variable = case_when(variable == "m" ~ "\n Maternal \n", 
                              variable == "p" ~ "\nPaternal \n", 
                              variable == "o" ~ "\nChild \n",
                              variable == "om" ~ "Maternal/Child\nCovariance\n",
                              variable == "op" ~ "Paternal/Child\nCovariance\n"),
         #label = case_when(best_fit == "yes" ~ "*")
         variable = factor(variable, levels = c("\nChild \n", "Maternal/Child\nCovariance\n", "\nPaternal \n", "Paternal/Child\nCovariance\n","\n Maternal \n" ))
         )
 
write.xlsx(result_df, file = "N:/durable/projects/trio_gcta_qc/outputs/estimates_afterreview_forbarplot.xlsx")


p <- ggplot(result_df, aes(fill=variable, y=Percent_variance, x=model)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~ outcome, switch = "y") +  
  expand_limits(y = c(0,51)) +
  scale_colour_brewer(type=seq, palette = "Paired", direction=1)+
  scale_fill_brewer( palette = "Paired", direction=1) +
  geom_text(aes(x = model, y = 45, label = best_fit_height), color = "black", size = 10) +
  geom_text(aes(x = model, y = 35, label = best_fit_exam), color = "black", size = 10) +
  geom_text(aes(x = model, y = 15, label = best_fit_sleep), color = "black", size = 10) +
  geom_text(aes(x = model, y = 10, label = best_fit_dep), color = "black", size = 10) +
  theme(text = element_text(size=12),
        axis.title.y = element_text(margin =  margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin =  margin(t = 20, r = 0, b = 0, l = 0)),
        panel.background = element_rect(fill = "white", colour = "grey10"),
        panel.spacing = unit(0.2, "lines"),
        panel.grid.major = element_line(colour="grey90"),
        panel.grid.minor = element_line(colour="grey90"),
        strip.text =   element_text(size=12, face="bold", colour="white", angle=0),
        legend.title = element_text(size =12),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_rect(fill="grey10")) +
  labs(y = "Percent Variance Explained",
       title = "",
       fill = "Genetic Effects",
       x = "Model") +
  guides(alpha = "none", 
         color = "none",
        # fill = "none"
  )



tiff("//ess01/p471/data/durable/projects/trio_gcta_qc/outputs/trio_gcta_afterreview.tiff",  width = 20, height = 10, units = 'in', res = 320)
p 
dev.off()

fits_files <-  as_tibble(list.files(pattern = "afterreview2_fit.csv")) %>%
  rename("File" = value) %>%
  mutate(outcome = str_remove(File, "afterreview2_fit.csv$"), 
         data = map(File, read.csv), 
         data = map2(data, outcome, ~ .x %>%
             mutate(model = c("full", "nocov", "direct", "null"), 
                    outcome = .y)
        
         ))

fit_data <- bind_rows(fits_files$data)

lrt_files <- as_tibble(list.files(pattern = "afterreview2_lrtest.csv")) %>%
  rename("File" = value ) %>%
  mutate(outcome = str_remove(File, "_(\\w+)(?:_\\w+)?\\.csv$"), 
         data = map(File, read.csv), 
         data = map2(data, outcome, ~ .x %>%
                       mutate(model = c("full", "nocov", "direct", "null"), 
                              outcome = .y)
                     
         ))

lrt_data <- bind_rows(lrt_files$data)

nonstd_chisq <- function(llratio, j){ #j = number of uncorrelated random effects removed that fall on the boundary
  
  if(is.na(llratio) | is.na(j)){
    return(NA_real_)
  }else{
    j <- abs(j)
    
    s <- 2^-j
    sum <- 0
    
    for (i in 0:j) {
      c <- choose(j, i)
      prob <- pchisq(llratio, df = i, lower.tail = FALSE)
      sum <- sum + (c * prob)
    }
    
    P = s*sum
    
    return(P)
  }
}

lrt_data <- lrt_data %>%
  mutate(p_final = map2_dbl(Î.Deviance, Î.DOF, nonstd_chisq), #I.deviance is the  difference in -2logLik  between models
         p_final = if_else(model == "nocov", p..Chisq., p_final))


cor_files <- as_tibble(list.files(pattern = "afterreview2_cor.csv")) %>%
  rename("File" = value ) %>%
  mutate(outcome = str_remove(File, "_(\\w+)(?:_\\w+)?\\.csv$"), 
         data = map(File, read.csv), 
         data = map2(data, outcome, ~ .x %>%
                       mutate(outcome = .y)
                     
         ))

cor_data <- bind_rows(cor_files$data) %>%
  filter(outcome == "exam")

write.xlsx(fit_data, file = "N:/durable/projects/trio_gcta_qc/outputs/fits_afterreview_trioGCTA.xlsx")
write.xlsx(lrt_data, file = "N:/durable/projects/trio_gcta_qc/outputs/lrt_afterreview_trioGCTA.xlsx")
write.xlsx(cor_data, file = "N:/durable/projects/trio_gcta_qc/outputs/cor_trioGCTA_afterreview.xlsx")

result_df_all <- rbind(exam,height,dep,sleep) %>%
  mutate(value = round(value,3)) %>%
  rename("Genetic Effect" = "variable", 
         "Estimate" = "value") %>%
  pivot_wider(names_from = "Genetic Effect", values_from = Estimate) %>%
  mutate(best_fit = case_when((outcome == "height" | outcome == "depression" | outcome == "sleep" ) & (model == "nocov") ~ "yes", 
                       (outcome == "exam") & model == "full" ~ "yes", 
                       TRUE ~"no")) %>%
  filter(best_fit == "yes")

write.xlsx(result_df_all, file = "N:/durable/projects/trio_gcta_qc/outputs/estimates_95_update_TrioGCTA.xlsx")

