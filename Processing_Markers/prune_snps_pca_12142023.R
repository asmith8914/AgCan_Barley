## Prune SNPs and do PCA
# 12/14/2023

library(SNPRelate)
library(tidyverse)
library(readxl)
library(cowplot)

# Read in vcf and convert to correct format 
vcf <- "Data/snps/prelim_merged_commonsites.vcf"

snpgdsVCF2GDS(vcf, "Data/snps/prelim_merged_commonsites.gds", method = "biallelic.only")

snpgdsSummary("Data/snps/prelim_merged_commonsites.gds")

snps_gds <- snpgdsOpen("Data/snps/prelim_merged_commonsites.gds")

# Prune snps based on LD with threshold 0.99 (only want to remove the entirely colinear ones).
set.seed(1000)
snpset <- snpgdsLDpruning(snps_gds, ld.threshold = 0.99)
# get 5,899 markers.

# First, get all selected snp ids.
snpset.id <- unlist(unname(snpset))

# Do PCAs
pca_set <- snpgdsPCA(snps_gds, snp.id = snpset.id)

# Fix the entry names
pca_set$sample.id <- gsub("\\.", "-", pca_set$sample.id)

# Make tables with PCs
pc_table <- data.frame(entry = pca_set$sample.id, pca_set$eigenvect)
colnames(pc_table)[c(2:33)] <- gsub("X", "PC", colnames(pc_table)[c(2:33)])

# Read in/make the categorical ----
Phenodata_prelim <- read_xlsx("C:/Users/alexa/OneDrive - University of Guelph/AgCan_Barley/AgCan_Barley_Data/May_update/Phenodata.xlsx", sheet = "PYT_2018-19-20-21")
# Fix genotype names
Phenodata_prelim2 <- Phenodata_prelim %>%
  mutate(NAME = gsub("'", "", NAME),
         NAME = gsub("AAC Ling", "AAC-Ling", NAME),
         NAME = gsub("AAC_Ling", "AAC-Ling", NAME))

cat <- Phenodata_prelim2 %>%
  select(c(NAME, YEAR)) %>%
  unique() %>%
  filter(NAME %in% pca_set$sample.id) %>%
  mutate(prelim_cohort = case_when(NAME == "AAC-Ling" ~ "AAC-Ling", NAME == "Leader" ~ "Leader", NAME != c("Leader", "AAC-Ling") ~ as.character(YEAR)),
         cross = case_when(NAME == "AAC-Ling" ~ "AAC-Ling", NAME == "Leader" ~ "Leader", NAME != c("Leader", "AAC-Ling") ~ gsub("-.*", "", NAME))) %>%
  rename(entry = NAME) %>%
  select(c(entry, prelim_cohort, cross)) %>%
  unique()

# Need to add AYT or not to categorical
gs_data <- read.csv("C:/Users/alexa/OneDrive - University of Guelph/AgCan_Barley/AgCan_Barley_Data/genomic_prediction/GS_data_PYT_AYT_alltraits.csv")

ayt_or_no <- !is.na(gs_data$mean_AYT_yield)
names(ayt_or_no) <- gs_data$NAME
ayt_or_no

cat$ayt <- ayt_or_no[cat$entry]

rm(Phenodata_prelim, Phenodata_prelim2, gs_data)

write.csv(cat, "Data/snps/prelim_cat.csv", row.names = FALSE)

# Merge the tables with the categorical
pc_cat <- inner_join(cat, pc_table)

# Plot the PCAs ----
# Coloured by cohort (year of PYT)
pc1pc2y <- ggplot(data = pc_cat, aes(x = PC1, y = PC2, colour = prelim_cohort, shape = prelim_cohort)) + geom_point() + labs(colour = "Cohort", shape = "Cohort") + scale_shape_manual(values = c(15,17,16,3,4)) + scale_colour_manual(values = c("#8064A2", "#9BBB59", "#F79646", "#C0504D", "#4BACC6")) + theme_light()
pc1pc2y
ggsave("Figures/pca/pc1pc2y_5899.pdf", height = 8.5, width = 11)
ggsave("Figures/pca/pc1pc2y_5899.svg", height = 8.5, width = 11)

pc2pc3y <- ggplot(data = pc_cat, aes(x = PC2, y = PC3, colour = prelim_cohort, shape = prelim_cohort)) + geom_point() + labs(colour = "Cohort", shape = "Cohort") + scale_shape_manual(values = c(15,17,16,3,4)) + scale_colour_manual(values = c("#8064A2", "#9BBB59", "#F79646", "#C0504D", "#4BACC6")) + theme_light()
pc2pc3y
ggsave("Figures/pca/pc2pc3y_5899.pdf", height = 8.5, width = 11)
ggsave("Figures/pca/pc2pc3y_5899.svg", height = 8.5, width = 11)

# Make a figure with both plots and one legend.
pc1pc2y_nolegend <- ggplot(data = pc_cat, aes(x = PC1, y = PC2, colour = prelim_cohort, shape = prelim_cohort)) + geom_point() + labs(colour = "Cohort", shape = "Cohort") + scale_shape_manual(values = c(15,17,16,3,4)) + scale_colour_manual(values = c("#8064A2", "#9BBB59", "#F79646", "#C0504D", "#4BACC6")) + theme_light() + theme(legend.position = "none")
pc2pc3y_nolegend <- ggplot(data = pc_cat, aes(x = PC2, y = PC3, colour = prelim_cohort, shape = prelim_cohort)) + geom_point() + labs(colour = "Cohort", shape = "Cohort") + scale_shape_manual(values = c(15,17,16,3,4)) + scale_colour_manual(values = c("#8064A2", "#9BBB59", "#F79646", "#C0504D", "#4BACC6")) + theme_light() + theme(legend.position = "none")
pcs_legend <- get_legend(pc1pc2y)
plot_grid(pc1pc2y_nolegend, pc2pc3y_nolegend, pcs_legend, labels = c("A", "B"), label_size = 12, nrow = 1, rel_widths = c(1,1,0.2))
ggsave("Figures/pca/composite_y_5899.pdf", height = 8.5, width = 14)
ggsave("Figures/pca/composite_y_5899.svg", height = 8.5, width = 14)

# Coloured by cross
pc1pc2c <- ggplot(data = pc_cat, mapping = aes(x = PC1, y = PC2, colour = cross, shape = cross)) + geom_point() + labs(colour = "Cross", shape = "Cross") + scale_shape_manual(values = rep(c(0:8, 15:19), len = 91)) + theme_light()
pc1pc2c
cross_legend <- get_legend(pc1pc2c)
pc1pc2c_noleg <- ggplot(data = pc_cat, mapping = aes(x = PC1, y = PC2, colour = cross, shape = cross)) + geom_point() + labs(colour = "Cross", shape = "Cross") + scale_shape_manual(values = rep(c(0:8, 15:19), len = 91)) + theme_light() + theme(legend.position = "none")
pc2pc3c_noleg <- ggplot(data = pc_cat, mapping = aes(x = PC2, y = PC3, colour = cross, shape = cross)) + geom_point() + labs(colour = "Cross", shape = "Cross") + scale_shape_manual(values = rep(c(0:8, 15:19), len = 91)) + theme_light() + theme(legend.position = "none")
# Composite
plot_grid(pc1pc2c_noleg, pc2pc3c_noleg, cross_legend, labels = c("A", "B"), label_size = 12, nrow = 1, scale = 0.9)
ggsave("Figures/pca/composite_c_5899.pdf", height = 8.5, width = 14)
ggsave("Figures/pca/composite_c_5899.svg", height = 8.5, width = 14)

# Coloured by selection in AYT or not
pc1pc2a <- ggplot(data = pc_cat, aes(x = PC1, y = PC2, colour = ayt, shape = ayt)) + geom_point() + labs(colour = "Advanced to AYT", shape = "Advanced to AYT") + scale_colour_manual(values = c("#4BACC6", "#1F497D"), labels = c("No", "Yes")) + scale_shape_manual(values = c(16,15), labels = c("No", "Yes")) + theme_light()
pc1pc2a
adv_legend <- get_legend(pc1pc2a)
pc1pc2_noleg <- ggplot(data = pc_cat, aes(x = PC1, y = PC2, colour = ayt, shape = ayt)) + geom_point() + labs(colour = "Advanced to AYT", shape = "Advanced to AYT") + scale_colour_manual(values = c("#4BACC6", "#1F497D"), labels = c("No", "Yes")) + scale_shape_manual(values = c(16,15), labels = c("No", "Yes")) + theme_light() + theme(legend.position = "none")
pc2pc3_noleg <- ggplot(data = pc_cat, aes(x = PC2, y = PC3, colour = ayt, shape = ayt)) + geom_point() + labs(colour = "Advanced to AYT", shape = "Advanced to AYT") + scale_colour_manual(values = c("#4BACC6", "#1F497D"), labels = c("No", "Yes")) + scale_shape_manual(values = c(16,15), labels = c("No", "Yes")) + theme_light() + theme(legend.position = "none")
plot_grid(pc1pc2_noleg, pc2pc3_noleg, adv_legend, labels = c("A", "B"), label_size = 12, nrow = 1, scale = 0.9,rel_widths = c(1,1,0.2))
ggsave("Figures/pca/composite_a_5899.pdf", height = 8.5, width = 14)
ggsave("Figures/pca/composite_a_5899.svg", height = 8.5, width = 14)

# Save the % variance explained by each PC ----
pca_set$varprop*100
write.csv(pca_set$varprop*100, "Data/snps/pca_per_var_5899.csv")

# Obtaining the snp set ----
# Going to convert to 012 format for genomic prediction using vcftools. Therefore, make a text file with the list of positions in the reduced snp set.
# Position filtering options in vcftools use either chromosome and position or site id.
library(vcfR)
vcf_all <- read.vcfR("Data/snps/prelim_merged_commonsites.vcf")

vcf_fix_tab <- as.data.frame(vcf_all@fix)
sub_vcf_fix_tab <- vcf_fix_tab %>%
  filter(row.names(vcf_fix_tab) %in% snpset.id)

# Make file with snp ids
pruned_snpids <- sub_vcf_fix_tab %>%
  select(ID)
write_delim(pruned_snpids, "Data/snps/snpids_5899.txt", delim = "\t", col_names = FALSE)
# Make file with snp positions and chromosomes
pruned_snppos <- sub_vcf_fix_tab %>%
  select(CHROM, POS)
write_delim(pruned_snppos, "Data/snps/snppos_5899.txt", delim = "\t", col_names = FALSE)
