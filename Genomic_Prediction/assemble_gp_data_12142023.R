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

# Combine with trait data ----
