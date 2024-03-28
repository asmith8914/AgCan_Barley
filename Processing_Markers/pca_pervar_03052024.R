# Remake PCAs
# 03/05/2024

# Objective: remake the PCAs with percentage of variation on axis labels

library(SNPRelate)
library(tidyverse)
library(readxl)
library(cowplot)

# Read in vcf
vcf <- "Data/snps/prelim_merged_5899.recode.vcf"

snpgdsVCF2GDS(vcf, "Data/snps/prelim_merged_5899.gds", method = "biallelic.only")

snpgdsSummary("Data/snps/prelim_merged_5899.gds")

snps_gds <- snpgdsOpen("Data/snps/prelim_merged_5899.gds")

# Do PCAs
pca<- snpgdsPCA(snps_gds)

# Fix the entry names
pca$sample.id <- gsub("\\.", "-", pca$sample.id)

# Make tables with PCs
pc_table <- data.frame(entry = pca$sample.id, pca$eigenvect)
colnames(pc_table)[c(2:33)] <- gsub("X", "PC", colnames(pc_table)[c(2:33)])

# Read in categorical
cat <- read.csv("Data/snps/prelim_cat.csv")

# Merge the tables with the categorical
pc_cat <- inner_join(cat, pc_table)

# Get the percentages of variation and round
per_var <- sprintf("%.1f", pca$varprop[1:3]*100)

# Plot the PCAs
# Coloured by cohort (year of PYT)
pc1pc2y <- ggplot(data = pc_cat, aes(x = PC1, y = PC2, colour = prelim_cohort, shape = prelim_cohort)) + geom_point() + labs(colour = "Cohort", shape = "Cohort") + xlab(paste0("PC1 (", per_var[1], "%)")) +  ylab(paste0("PC2 (", per_var[2], "%)")) + scale_shape_manual(values = c(15,17,16,3,4)) + scale_colour_manual(values = c("#8064A2", "#9BBB59", "#F79646", "#C0504D", "#4BACC6")) + theme_light()
pc1pc2y
ggsave("Figures/pca/pc1pc2y_5899_pervar.pdf", height = 8.5, width = 11)
ggsave("Figures/pca/pc1pc2y_5899_pervar.svg", height = 8.5, width = 11)

pc2pc3y <- ggplot(data = pc_cat, aes(x = PC2, y = PC3, colour = prelim_cohort, shape = prelim_cohort)) + geom_point() + labs(colour = "Cohort", shape = "Cohort") + xlab(paste0("PC2 (", per_var[2], "%)")) +  ylab(paste0("PC3 (", per_var[3], "%)")) + scale_shape_manual(values = c(15,17,16,3,4)) + scale_colour_manual(values = c("#8064A2", "#9BBB59", "#F79646", "#C0504D", "#4BACC6")) + theme_light()
pc2pc3y
ggsave("Figures/pca/pc2pc3y_5899_pervar.pdf", height = 8.5, width = 11)
ggsave("Figures/pca/pc2pc3y_5899_pervar.svg", height = 8.5, width = 11)

# Make a figure with both plots and one legend.
pcs_legend <- get_legend(pc1pc2y)
plot_grid(pc1pc2y + theme(legend.position = "none"), 
          pc2pc3y + theme(legend.position = "none"), 
          pcs_legend, labels = c("A", "B"), label_size = 12, nrow = 1, rel_widths = c(1,1,0.2))
ggsave("Figures/pca/composite_y_5899_pervar.pdf", height = 8.5, width = 14)
ggsave("Figures/pca/composite_y_5899_pervar.svg", height = 8.5, width = 14)

# Coloured by cross
pc1pc2c <- ggplot(data = pc_cat, mapping = aes(x = PC1, y = PC2, colour = cross, shape = cross)) + geom_point() + labs(colour = "Cross", shape = "Cross") + xlab(paste0("PC1 (", per_var[1], "%)")) +  ylab(paste0("PC2 (", per_var[2], "%)")) + scale_shape_manual(values = rep(c(0:8, 15:19), len = 91)) + theme_light()
pc1pc2c
cross_legend <- get_legend(pc1pc2c)
pc2pc3c <- ggplot(data = pc_cat, mapping = aes(x = PC2, y = PC3, colour = cross, shape = cross)) + geom_point() + labs(colour = "Cross", shape = "Cross") + xlab(paste0("PC2 (", per_var[2], "%)")) +  ylab(paste0("PC3 (", per_var[3], "%)")) + scale_shape_manual(values = rep(c(0:8, 15:19), len = 91)) + theme_light()
# Composite
plot_grid(pc1pc2c + theme(legend.position = "none"), 
          pc2pc3c + theme(legend.position = "none"), 
          cross_legend, labels = c("A", "B"), label_size = 12, nrow = 1, scale = 0.9)
ggsave("Figures/pca/composite_c_5899_pervar.pdf", height = 8.5, width = 14)
ggsave("Figures/pca/composite_c_5899_pervar.svg", height = 8.5, width = 14)

# Coloured by selection in AYT or not
pc1pc2a <- ggplot(data = pc_cat, aes(x = PC1, y = PC2, colour = ayt, shape = ayt)) + geom_point() + labs(colour = "Advanced to AYT", shape = "Advanced to AYT") + xlab(paste0("PC1 (", per_var[1], "%)")) +  ylab(paste0("PC2 (", per_var[2], "%)")) + scale_colour_manual(values = c("#4BACC6", "#1F497D"), labels = c("No", "Yes")) + scale_shape_manual(values = c(16,15), labels = c("No", "Yes")) + theme_light()
pc1pc2a
adv_legend <- get_legend(pc1pc2a)
pc2pc3a <- ggplot(data = pc_cat, aes(x = PC2, y = PC3, colour = ayt, shape = ayt)) + geom_point() + labs(colour = "Advanced to AYT", shape = "Advanced to AYT") + xlab(paste0("PC2 (", per_var[2], "%)")) +  ylab(paste0("PC3 (", per_var[3], "%)")) + scale_colour_manual(values = c("#4BACC6", "#1F497D"), labels = c("No", "Yes")) + scale_shape_manual(values = c(16,15), labels = c("No", "Yes")) + theme_light() + theme(legend.position = "none")
plot_grid(pc1pc2a + theme(legend.position = "none"),
          pc2pc3a + theme(legend.position = "none"), 
          adv_legend, labels = c("A", "B"), label_size = 12, nrow = 1, scale = 0.9,rel_widths = c(1,1,0.2))
ggsave("Figures/pca/composite_a_5899_pervar.pdf", height = 8.5, width = 14)
ggsave("Figures/pca/composite_a_5899_pervar.svg", height = 8.5, width = 14)
