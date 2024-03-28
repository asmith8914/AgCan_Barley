## Comparing PAs for predictions of advanced trials and response to selection
library(tidyverse)

cor_df <- read.csv("Data/gp/gp_pa_AYT_rep_unrep.csv")

round_unrep <- cor_df %>% 
  mutate(cor = as.numeric(sprintf("%.2f", cor))) %>%
  filter(rep_level == "unrep") %>%
  group_by(trait) %>%
  summarise(max = max(cor), min = min(cor))

comp_dif <- cor_df %>% 
  mutate(cor = as.numeric(sprintf("%.2f", cor))) %>%
  filter(rep_level == "unrep" & trait != "yield") %>%
  arrange(cohort, trait) %>%
  group_by(cohort, trait) %>%
  mutate(pa_dif = cor - lag(cor))

comp_dif$pa_dif[!is.na(comp_dif$pa_dif)]

gp_range <- cor_df %>% 
  mutate(cor = as.numeric(sprintf("%.2f", cor))) %>%
  filter(rep_level == "unrep" & trait != "yield" & cor_type == "gp_cor") %>%
  group_by(trait) %>%
  summarise(max = max(cor), min = min(cor))

ps_range <- cor_df %>% 
  mutate(cor = as.numeric(sprintf("%.2f", cor))) %>%
  filter(rep_level == "unrep" & trait != "yield" & cor_type == "pheno_cor") %>%
  group_by(trait) %>%
  summarise(max = max(cor), min = min(cor))

gp_rep <- cor_df %>% 
  mutate(cor = as.numeric(sprintf("%.2f", cor))) %>%
  filter(rep_level == "rep" & trait != "yield" & cor_type == "gp_cor") %>%
  group_by(trait) %>%
  summarise(max = max(cor), min = min(cor))

ps_rep <- cor_df %>% 
  mutate(cor = as.numeric(sprintf("%.2f", cor))) %>%
  filter(rep_level == "rep" & trait != "yield" & cor_type == "pheno_cor") %>%
  group_by(trait) %>%
  summarise(max = max(cor), min = min(cor))

# Response to selection
resp_sel <- read.csv("Data/gp/resptosel_res_trialmodel.csv")

yield_dif <- resp_sel %>% 
  filter(trait == "yield") %>%
  arrange(cohort, method) %>%
  group_by(cohort) %>%
  mutate(dif = change_mean - lag(change_mean))

yield_sum <- resp_sel %>%
  filter(trait == "yield") %>%
  group_by(method) %>%
  summarise(mean = mean(change_mean))

swtw <- resp_sel %>% 
  filter((trait == "sw" & cohort == "2020") | (trait == "tw" & cohort == "2019"))

# Also need to compare for the other values calculated for repsonse to selection
resp_sel_allayt <- read.csv("Data/gp/resptosel_res_allselaytlines.csv")
resp_sel_allpyt <- read.csv("Data/gp/resptosel_res_allsellines.csv")

yield_dif_allayt <- resp_sel_allayt %>% 
  filter(trait == "yield") %>%
  arrange(cohort, method) %>%
  group_by(cohort) %>%
  mutate(dif = change_mean - lag(change_mean))

yield_sum_allayt <- resp_sel_allayt %>%
  filter(trait == "yield") %>%
  group_by(method) %>%
  summarise(mean = mean(change_mean))

yield_dif_allpyt <- resp_sel_allpyt %>% 
  filter(trait == "yield") %>%
  arrange(cohort, method) %>%
  group_by(cohort) %>%
  mutate(dif = change_mean - lag(change_mean))

yield_sum_allpyt <- resp_sel_allpyt %>%
  filter(trait == "yield") %>%
  group_by(method) %>%
  summarise(mean = mean(change_mean))
