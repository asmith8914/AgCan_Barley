## Assemble genomic prediction data
# 12/14/2023

# Objective: assemble the marker data (first 012 format, convert to -101), combine with trait data

library(tidyverse)

# Assemble marker data ----
# Read in marker data
calls <- read.delim("Data/snps/prelim_merged_5899_012.012", header = FALSE)
lines <- read.delim("Data/snps/prelim_merged_5899_012.012.indv", header = FALSE)
loci <- read.delim("Data/snps/prelim_merged_5899_012.012.pos", header = FALSE)

# Merge into one table (original V1 in calls is an index (rownames), not data, replace with entry. Then rename column names with site positions.)
geno_table <- calls %>%
  mutate(V1 = lines$V1)

colnames(geno_table) <- c("entry", paste0(loci$V1, ":", loci$V2))

# Fix entry names in geno_table
geno_table <- geno_table %>%
  mutate(entry = gsub("\\.", "-", entry))

# Save
write.csv(geno_table, "Data/snps/prelim_5899_012_assembled.csv", row.names = FALSE)

# Subtract 1 to get {-1, 0, 1} from {0, 1, 2}
marker_matx <- as.matrix(geno_table[,-c(1)])
marker_matx <- marker_matx-1
rownames(marker_matx) <- geno_table$entry

# Save
write.csv(marker_matx, "Data/snps/prelim_5899_-101_assembled.csv")

# Combine with categorical and trait data ----
# Read in trait data
pyt_trait <- read.csv("Data/pheno/pyt_emmeans.csv")
pyt_unrep_trait <- read.csv("Data/pheno/pyt_unrep_emmeans.csv")
ayt_trait <- read.csv("Data/pheno/adv_emmeans.csv")

# Read in categorical, fix names
cat <- read.csv("Data/snps/prelim_cat.csv")
cat <- cat %>%
  rename(NAME = entry, Cohort = prelim_cohort)

# Left join categorical, trait data
cat_trait <- left_join(cat, pyt_trait) %>%
  left_join(., pyt_unrep_trait) %>%
  left_join(., ayt_trait)

# Left join categorical and trait with markers
marker_matx2 <- as.data.frame(marker_matx) %>%
  rownames_to_column(var = "NAME")

gp_data <- left_join(cat_trait, marker_matx2)

# Check that this new genomic prediction file has the same trait data as previously.
old_gp_data <- read.csv("C:/Users/alexa/OneDrive - University of Guelph/AgCan_Barley/AgCan_Barley_Data/genomic_prediction/GS_data_PYT_AYT_alltraits_unrep.csv")
old_gp_data <- old_gp_data %>%
  select(colnames(old_gp_data)[c(1:4, 6, 7, 5, 8:21)]) %>%
  rename(mean_PYT_sw = mean_PYT_sw_noblock, mean_PYT_tw = mean_PYT_tw_noblock)

all.equal(old_gp_data[colnames(old_gp_data)], gp_data[colnames(gp_data)[c(1:2, 5:23)]])
# Trait data is the same!

# Save the new genomic prediction data:
write.csv(gp_data, "Data/gp/GS_PYT_PYTu_AYT_alltraits_5899.csv", row.names = FALSE)
