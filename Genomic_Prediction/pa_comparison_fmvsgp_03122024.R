## Comparing PAs for family means, genomic predictions, and combined predictions
library(tidyverse)

overall_res <- read.csv("Data/gp/gp_pa_fm_comb_rel2.csv")

max_min <- overall_res %>%
  group_by(trait) %>%
  summarise(max = max(overall_pa), min = min(overall_pa))
max_min %>% mutate(max = sprintf("%.2f", max), min = sprintf("%.2f", min))

rounded <- overall_res %>%
  mutate(overall_pa = sprintf("%.2f", overall_pa), stdev = sprintf("%.3f", stdev))

max_min2 <- rounded %>%
  group_by(trait) %>%
  filter(overall_pa == max(overall_pa) | overall_pa == min(overall_pa)) %>%
  arrange(trait, overall_pa)

comb_highergp <- rounded %>%
  filter(!(trait == "yield" & Cohort %in% c("2019", "2020")) & !(trait == "tw" & Cohort == "2020")) %>%
  filter(method != "Family Means")

diff_comb_gp <- comb_highergp %>%
  arrange(Cohort, trait) %>%
  mutate(overall_pa = as.numeric(overall_pa)) %>%
  group_by(Cohort, trait) %>%
  mutate(pa_dif = overall_pa - lag(overall_pa))

mean(diff_comb_gp$pa_dif, na.rm = TRUE)

sqrt(var(diff_comb_gp$pa_dif, na.rm = TRUE))  

comb_higherfm <- rounded %>%
  filter(!(trait == "ph" & Cohort == "2018") & !(trait == "hd" & Cohort == "2020")) %>%
  filter(method != "Genomic Prediction")

diff_comb_fm <- comb_higherfm %>%
  arrange(Cohort, trait) %>%
  mutate(overall_pa = as.numeric(overall_pa)) %>%
  group_by(Cohort, trait) %>%
  mutate(pa_dif = overall_pa - lag(overall_pa))

mean(diff_comb_fm$pa_dif, na.rm = TRUE)

sqrt(var(diff_comb_fm$pa_dif, na.rm = TRUE))  

fmvsgp <- rounded %>%
  filter(method != "Combined Prediction") %>%
  filter((trait  %in% c("yield", "height", "tw") & Cohort == "2018") | (trait == "hd" & Cohort == "2019") | (trait  %in% c("sw", "hd") & Cohort == "2020"))

diff_fm_gp <- fmvsgp %>%
  arrange(Cohort, trait) %>%
  mutate(overall_pa = as.numeric(overall_pa)) %>%
  group_by(Cohort, trait) %>%
  mutate(pa_dif = overall_pa - lag(overall_pa))

mean(diff_fm_gp$pa_dif, na.rm = TRUE)

sqrt(var(diff_fm_gp$pa_dif, na.rm = TRUE))  
