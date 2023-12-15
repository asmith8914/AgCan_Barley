## Catch up AYT heading date!
# 09/15/2023

# Objectives: Calculate days to heading

library(tidyverse)
library(readxl)

# Read in heading data
ayt_pheno <- read_xlsx("C:/Users/alexa/OneDrive - University of Guelph/AgCan_Barley/AgCan_Barley_Data/Phenodata_adv_hd_cor.xlsx", sheet = "AYT_2019-20-21")

# Filter out 2022 and any remaining missing data (no days to heading, or planting or heading dates)
ayt_pheno2 <- ayt_pheno %>%
  filter(YEAR != 2022)

ayt_pheno_hd <- ayt_pheno2 %>%
  filter(!is.na(HD_DAYS) | (!is.na(HD)) & !is.na(`PLANTING DATE`))

# For materials and methods, which data is missing?
missing <- ayt_pheno2 %>%
  filter(is.na(HD_DAYS) & (is.na(HD) | is.na(`PLANTING DATE`)))

# Confirm locations are the same
ayt_pheno %>% ayt_pheno2
  filter(YEAR != 2022) %>%
  group_by(YEAR, LOCATION) %>%
  summarise()
# There is a St. Rosalie 2020 trial (not described in the excel description but it was in my M&M)

# Calculate the days to heading
ayt_pheno_dth <- ayt_pheno_hd %>%
  mutate(HEADING_DATE = case_when(grepl("[0-9][0-9]", ayt_pheno_hd$HD) == TRUE ~ as.character(HD), 
                                  (grepl("[0-9][0-9]", ayt_pheno_hd$HD) == FALSE & is.na(ayt_pheno_hd$HD) == FALSE) ~ paste0("0", ayt_pheno_hd$HD),
                                  (grepl("[0-9][0-9]", ayt_pheno_hd$HD) == FALSE & is.na(ayt_pheno_hd$HD) == TRUE) ~ NA)) %>%
  mutate(heading_ymd = case_when((is.na(HM) == FALSE & is.na(HEADING_DATE) == FALSE) ~ paste0(YEAR, "-0", HM, "-", HEADING_DATE),
                                 (is.na(HM) == TRUE | is.na(HEADING_DATE) == TRUE) ~ NA))

ayt_pheno_dth2 <- ayt_pheno_dth %>%
  mutate(days_to_heading = case_when(is.na(heading_ymd) == TRUE ~ HD_DAYS,
                                     is.na(heading_ymd) == FALSE ~ as.numeric(as.Date(heading_ymd) - as.Date(`PLANTING DATE`))))

# Join the days to heading with the original (complete) dataset
ayt_pheno_all <- left_join(ayt_pheno2, ayt_pheno_dth2)

# There is one observation with "-9" as days to heading. Replace with NA
ayt_pheno_all2 <- ayt_pheno_all %>%
  mutate(days_to_heading = case_when(days_to_heading > 0 ~ days_to_heading,
                                     days_to_heading < 0 ~ NA))

# Save the days to heading.
phenodata_save <- ayt_pheno_all2[, c(1:10,17)]

# Clean up entry names and column names
adv_fixed <- phenodata_save %>%
  mutate(NAME = gsub("'", "", NAME),
         NAME = gsub("AAC Ling", "AAC-Ling", NAME),
         NAME = gsub("AAC_Ling", "AAC-Ling", NAME),
         NAME = gsub("AAC Synergy", "AAC-Synergy", NAME)) %>%
  rename(yield = `YIELD_KGHA-1`, testweight = `TESTWEIGHT_hLKG-1`) %>%
  mutate(YEAR = as.factor(YEAR), REP = as.factor(REP))

table(adv_fixed$NAME)
  
# Save data
write.csv(adv_fixed, "Data/pheno/Phenodata_adv_hd_cor.csv", row.names = FALSE)
