## Format SNP files
# 04/04/2023

# Objectives: Change formatting of chromosome positions so all are consistent with each other 
library(tidyverse)

# Read in data.
prelim_2018_2020_snps_orig <- read.delim("AgCan_Barley_Data/Raw_Data_Copy/Prelim2018_2020_rajabarley.hmp.txt", header = TRUE, sep = "\t")
prelim_2019_snps_orig <- read.delim("AgCan_Barley_Data/Raw_Data_Copy/Prelim2019_rajabarey.hmp.txt", header = TRUE, sep = "\t")

# Fix chromosome positions in Prelim 2019
prelim_2019_snps_fixed <- prelim_2019_snps_orig %>%
  mutate(rs. = gsub("S", "chr", rs.)) %>%
  mutate(rs. = gsub("_", ":", rs.))
write.table(prelim_2019_snps_fixed, "AgCan_Barley_Data/Data_Formatted/prelim_2019_fixed.hmp.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# All SNP files have been converted to vcf using TASSEL. Now can merge, filter SNPs, etc using vcftools. 
# Later, will need to reformat entry names to be consistent.

# Again, note that R recodes headers (entry names), replacing "-" with "." (Phenodata uses "-", although some have extra quotation marks, so these will need to be reformatted). Original Prelim trial files use "-", but unrep file uses ".". Entry names in pretty much all files therefore will need to be reformatted when merging trait and genotype data.