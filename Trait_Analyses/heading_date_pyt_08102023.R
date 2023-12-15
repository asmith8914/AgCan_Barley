## Obtain days to heading
# 08/10/2023

# Objective: add/fix planting dates to Phenodata, then calculate days to heading

library(tidyverse)

orig_phenodata <- read.csv("C:/Users/alexa/OneDrive - University of Guelph/AgCan_Barley/AgCan_Barley_Data/phenodata/Phenodata_prelim_fixed.csv")

# Make a new column with the correct planting dates. Also removing 2021 prelims.
phenodata_pl<- orig_phenodata %>%
  filter(YEAR != 2021) %>%
  mutate(pl_date_cor = case_when(YEAR == 2018 ~ "2018-05-09", YEAR == 2019 ~ "2019-05-02", YEAR == 2020 ~ "2020-05-19"))

# Make a new column with the heading month and date in calendar format (yyyy-mm-dd).
phenodata_hd <- phenodata_pl %>%
  mutate(HEADING_DATE = case_when(grepl("[0-9][0-9]", phenodata_pl$`HEADING.DATE`) == TRUE ~ as.character(`HEADING.DATE`), grepl("[0-9][0-9]", phenodata_pl$`HEADING.DATE`) == FALSE ~ paste0("0", phenodata_pl$`HEADING.DATE`))) %>%
  mutate(heading_ymd = paste0(YEAR, "-0", HEADING_MONTH, "-", HEADING_DATE))

# Calculate the days to heading. Need difference between pl_date_cor and heading_ymd.
phenodata_dth <- phenodata_hd %>%
  mutate(days_to_heading = as.numeric(as.Date(heading_ymd) - as.Date(pl_date_cor)))

# Clean up the columns and save
phenodata_save <- phenodata_dth[, c(1:11,19)]

phenodata_save <- phenodata_save %>%
  rename(yield = `YIELD_KGHA.1`, testweight = `TESTWEIGHT_hLKG.1`)

write.csv(phenodata_save, "Data/pheno/Phenodata_prelim_hd.csv", row.names = FALSE)

# Check the distributions
hd_dist <- ggplot(data = phenodata_save, mapping = aes(x = factor(YEAR), y = days_to_heading, fill = factor(YEAR))) + geom_boxplot() + geom_jitter() + scale_fill_manual(values = c("#8064A2", "#9BBB59", "#F79646")) + xlab("Cohort (Year of Preliminary Trial)") + ylab("Days to Heading") + theme(legend.position = "none")
hd_dist
