#author: Laura Hegemann 
#process results for downstream use



library(tidyverse)
library(ggsci)
library(xlsx)

setwd("//ess01/p471/data/durable/people/Laura_H/Cluster_backup/trio_gcta_qc/results/") #LJH changed to user-invariant filepath


results <- list.files(pattern = "exam_95_update.csv|height_95_update.csv")

results

exam <- read.csv(results[1]) %>%
  mutate(outcome = "exam")
height <- read.csv(results[2]) %>%
  mutate(outcome = "height")


result_df <- rbind(exam,height) %>%
  mutate(best_fit_height = case_when((outcome == "height") & model == "nocov" ~ "*", 
                              TRUE ~" "),
         best_fit_exam = case_when((outcome == "exam") & model == "nocov" ~ "*", 
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
                                                  outcome == "height" ~ "Height"),
         variable = case_when(variable == "m" ~ "\n Maternal \n", 
                              variable == "p" ~ "\nPaternal \n", 
                              variable == "o" ~ "\nChild \n",
                              variable == "om" ~ "Maternal/Child\nCovariance\n",
                              variable == "op" ~ "Paternal/Child\nCovariance\n"),
         #label = case_when(best_fit == "yes" ~ "*")
         variable = factor(variable, levels = c("\nChild \n", "Maternal/Child\nCovariance\n", "\nPaternal \n", "Paternal/Child\nCovariance\n","\n Maternal \n" ))
         )
 
write.xlsx(result_df, file = "N:/durable/projects/trio_gcta_qc/outputs/estimates_95_update_forbarplot.xlsx")


p <- ggplot(result_df, aes(fill=variable, y=Percent_variance, x=model)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~ outcome, switch = "y") +  
  expand_limits(y = c(0,51)) +
  scale_colour_brewer(type=seq, palette = "Paired", direction=1)+
  scale_fill_brewer( palette = "Paired", direction=1) +
  geom_text(aes(x = model, y = 50, label = best_fit_height), color = "black", size = 10) +
  geom_text(aes(x = model, y = 45, label = best_fit_exam), color = "black", size = 10) +
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



tiff("//ess01/p471/data/durable/projects/trio_gcta_qc/outputs/trio_gcta_95.tiff",  width = 20, height = 10, units = 'in', res = 320)
p 
dev.off()

fits_files <-  as_tibble(list.files(pattern = "95_update_fit.csv")) %>%
  rename("File" = value) %>%
  mutate(outcome = str_remove(File, "_95_update_fit.csv$"), 
         data = map(File, read.csv), 
         data = map2(data, outcome, ~ .x %>%
             mutate(model = c("full", "nocov", "direct", "null"), 
                    outcome = .y)
        
         ))

fit_data <- bind_rows(fits_files$data)

lrt_files <- as_tibble(list.files(pattern = "95_update_lrtest.csv")) %>%
  rename("File" = value ) %>%
  mutate(outcome = str_remove(File, "_(\\w+)(?:_\\w+)?\\.csv$"), 
         data = map(File, read.csv), 
         data = map2(data, outcome, ~ .x %>%
                       mutate(model = c("full", "nocov", "direct", "null"), 
                              outcome = .y)
                     
         ))

lrt_data <- bind_rows(lrt_files$data)

#full <- lrt_data %>% filter(model == "full") %>% mutate(p_adj = NA_real_)
#nocov_p <- lrt_data %>% filter(model == "nocov") %>% mutate(p_adj = p.adjust(.$p..Chisq., method = "fdr", n = 6))
#direct_p <- lrt_data %>% filter(model == "direct") %>% mutate(p_adj = p.adjust(.$p..Chisq., method = "fdr", n = 6))
#null_p <- lrt_data %>% filter(model == "null") %>% mutate(p_adj = p.adjust(.$p..Chisq., method = "fdr", n = 6))

#lrt_data2 <- bind_rows(full, nocov_p, direct_p, null_p)

cor_files <- as_tibble(list.files(pattern = "95_update_cor.csv")) %>%
  rename("File" = value ) %>%
  mutate(outcome = str_remove(File, "_(\\w+)(?:_\\w+)?\\.csv$"), 
         data = map(File, read.csv), 
         data = map2(data, outcome, ~ .x %>%
                       mutate(outcome = .y)
                     
         ))

cor_data <- bind_rows(cor_files$data)

write.xlsx(fit_data, file = "N:/durable/projects/trio_gcta_qc/outputs/fits_95_update_trioGCTA.xlsx")
write.xlsx(lrt_data, file = "N:/durable/projects/trio_gcta_qc/outputs/lrt_95_update_trioGCTA.xlsx")
#write.xlsx(cor_data, file = "N:/durable/projects/trio_gcta_qc/outputs/cor_trioGCTA.xlsx")

result_df_all <- rbind(exam,height) %>%
  mutate(value = round(value,3)) %>%
  rename("Genetic Effect" = "variable", 
         "Estimate" = "value") %>%
  pivot_wider(names_from = "Genetic Effect", values_from = Estimate) %>%
  mutate(best_fit = case_when((outcome == "height") & (model == "direct"  | model == "nocov") ~ "yes", 
                       (outcome == "exam") & model == "nocov" ~ "yes", 
                       TRUE ~"no"))

write.xlsx(result_df_all, file = "N:/durable/projects/trio_gcta_qc/outputs/estimates_95_update_TrioGCTA.xlsx")

