## Trait plots and correlations
# 02/15/2024

# Objectives: 1) xy scatter plots showing correlation between PYTs and AYTs, 2) boxplots showing trait distributions, and 3) correlations of traits with other traits within both PYTs and AYTs, and correlations of traits with themselves between PYTs and AYTs. 
# Using entries' adjusted means. For 3) this will help explanation of trends in response to selection (yield was selected for, but response to selection for all traits was calculated).
# trait_plots_09142023.R and trait_cor_10182023.R

library(tidyverse)
library(cowplot)

# Read in data
gp_data <- read.csv("Data/gp/GS_PYT_PYTu_AYT_alltraits_5899.csv")

# 1) XY Scatters showing correlations between PYTs and AYTs ----
means_pyts_ayts <- gp_data %>%
  filter(!is.na(mean_AYT_yield))

pyt_ayt_y <- ggplot(data = means_pyts_ayts, mapping = aes(x = mean_PYT_yield, y = mean_AYT_yield, colour = Cohort, shape = Cohort)) + geom_point() + xlab("Mean PYT yield (kg/ha)") + ylab("Mean AYT yield (kg/ha)") + scale_shape_manual(values = c(15,17,16,3,4)) + scale_colour_manual(values = c("#8064A2", "#9BBB59", "#F79646", "#C0504D", "#4BACC6")) + theme_light()
pyt_ayt_y

pyt_ayt_ph <- ggplot(data = means_pyts_ayts, mapping = aes(x = mean_PYT_height, y = mean_AYT_height, colour = Cohort, shape = Cohort)) + geom_point() + xlab("Mean PYT height (cm)") + ylab("Mean AYT height (cm)") + scale_shape_manual(values = c(15,17,16,3,4)) + scale_colour_manual(values = c("#8064A2", "#9BBB59", "#F79646", "#C0504D", "#4BACC6")) + theme_light()
pyt_ayt_ph

pyt_ayt_tw <- ggplot(data = means_pyts_ayts, mapping = aes(x = mean_PYT_tw, y = mean_AYT_tw, colour = Cohort, shape = Cohort)) + geom_point() + xlab("Mean PYT test weight (kg/hL)") + ylab("Mean AYT test weight (kg/hL)") + scale_shape_manual(values = c(15,17,16,3,4)) + scale_colour_manual(values = c("#8064A2", "#9BBB59", "#F79646", "#C0504D", "#4BACC6")) + theme_light()
pyt_ayt_tw

pyt_ayt_sw <- ggplot(data = means_pyts_ayts, mapping = aes(x = mean_PYT_sw, y = mean_AYT_sw, colour = Cohort, shape = Cohort)) + geom_point() + xlab("Mean PYT seed weight (g)") + ylab("Mean AYT seed weight (g)") + scale_shape_manual(values = c(15,17,16,3,4)) + scale_colour_manual(values = c("#8064A2", "#9BBB59", "#F79646", "#C0504D", "#4BACC6")) + theme_light()
pyt_ayt_sw

pyt_ayt_hd <- ggplot(data = means_pyts_ayts, mapping = aes(x = mean_PYT_hd, y = mean_AYT_hd, colour = Cohort, shape = Cohort)) + geom_point() + xlab("Mean PYT days to heading") + ylab("Mean AYT days to heading") + scale_shape_manual(values = c(15,17,16,3,4)) + scale_colour_manual(values = c("#8064A2", "#9BBB59", "#F79646", "#C0504D", "#4BACC6")) + theme_light()
pyt_ayt_hd

# Make composite image with all traits
pyt_ayt_leg <- get_legend(pyt_ayt_y)

# Add correlations
trait_cor <- means_pyts_ayts %>%
  select(1:23) %>%
  summarise(cor_y = cor(mean_PYT_yield, mean_AYT_yield),
            cor_ph = cor(mean_PYT_height, mean_AYT_height),
            cor_sw = cor(mean_PYT_sw, mean_AYT_sw),
            cor_tw = cor(mean_PYT_tw, mean_AYT_tw),
            cor_hd = cor(mean_PYT_hd, mean_AYT_hd))

plot_grid(pyt_ayt_y + theme(legend.position = "none") + annotate("text", label = paste0("r = ", sprintf("%.3f", trait_cor$cor_y)), x = 4300, y = 4000),
          pyt_ayt_ph + theme(legend.position = "none") + annotate("text", label = paste0("r = ", sprintf("%.3f", trait_cor$cor_ph)), x = 67, y = 85),
          pyt_ayt_sw + theme(legend.position = "none") + annotate("text", label = paste0("r = ", sprintf("%.3f", trait_cor$cor_sw)), x = 55, y = 40),
          pyt_ayt_tw + theme(legend.position = "none") + annotate("text", label = paste0("r = ", sprintf("%.3f", trait_cor$cor_tw)), x = 64.5, y = 69),
          pyt_ayt_hd + theme(legend.position = "none") + annotate("text", label = paste0("r = ", sprintf("%.3f", trait_cor$cor_hd)), x = 72.5, y = 58),
          pyt_ayt_leg, labels = c("A", "B", "C", "D", "E"))
ggsave("Figures/trait/trait_cor_r.pdf", height = 8.5, width = 11)
ggsave("Figures/trait/trait_cor_r.svg", height = 8.5, width = 11)

# 2) Boxplots with trait distributions ----
# Pivot data longer, remove NAs
for_boxplots <- gp_data %>%
  select(c(1,2,5:9,19:23)) %>%
  pivot_longer(cols = !c(NAME, Cohort), names_to = c("trial", "trait"), names_prefix = "mean_", names_sep = "_", values_to = "mean") %>%
  filter(!Cohort %in% c("AAC-Ling", "Leader")) %>%
  filter(!is.na(mean))

# Make a copy, filter to PYT lines advanced to AYT, pivot data longer, filter to PYT observations only, change trial name
for_boxplots2 <- gp_data %>%
  select(c(1,2,5:9,19:23)) %>%
  filter(!is.na(mean_AYT_yield)) %>%
  pivot_longer(cols = !c(NAME, Cohort), names_to = c("trial", "trait"), names_prefix = "mean_", names_sep = "_", values_to = "mean") %>%
  filter(!Cohort %in% c("AAC-Ling", "Leader")) %>%
  filter(trial == "PYT") %>%
  mutate(trial = "PYT (AYT lines only)")

# Bind them together
for_boxplots3 <- rbind(for_boxplots, for_boxplots2)

# Make a copy to create "Overall" cohort
for_boxplots4 <- for_boxplots3 %>%
  mutate(Cohort = "Overall")

# Bind them together
for_boxplots5 <- rbind(for_boxplots3, for_boxplots4)

# Make a vector with the trial types
trial_types <- c("PYT", "PYT (AYT lines only)", "AYT")

# PYT, PYT selected, AYT boxplots (all traits, each cohort and overall)
box_y <- ggplot(data = for_boxplots5 %>% filter(trait == "yield"), mapping = aes(x = Cohort, y = mean, fill = factor(trial, trial_types))) + labs(fill = "Trial Type") + ylab("Mean Yield (kg/ha)") + geom_boxplot(outlier.shape = NA) + theme_light() + scale_fill_manual(values = c("#9BBB59", "#F79646", "#8064A2")) + geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 1)
box_y

box_ph <- ggplot(data = for_boxplots5 %>% filter(trait == "height"), mapping = aes(x = Cohort, y = mean, fill = factor(trial, trial_types))) + labs(fill = "Trial Type") + ylab("Mean Height (cm)") + geom_boxplot(outlier.shape = NA) + theme_light() + scale_fill_manual(values = c("#9BBB59", "#F79646", "#8064A2")) + geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 1)
box_ph

box_sw <- ggplot(data = for_boxplots5 %>% filter(trait == "sw"), mapping = aes(x = Cohort, y = mean, fill = factor(trial, trial_types))) + labs(fill = "Trial Type") + ylab("Mean Seed Weight (g)") + geom_boxplot(outlier.shape = NA) + theme_light() + scale_fill_manual(values = c("#9BBB59", "#F79646", "#8064A2")) + geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 1)
box_sw

box_tw <- ggplot(data = for_boxplots5 %>% filter(trait == "tw"), mapping = aes(x = Cohort, y = mean, fill = factor(trial, trial_types))) + labs(fill = "Trial Type") + ylab("Mean Test Weight (kg/hL)") + geom_boxplot(outlier.shape = NA) + theme_light() + scale_fill_manual(values = c("#9BBB59", "#F79646", "#8064A2")) + geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 1)
box_tw

box_hd <- ggplot(data = for_boxplots5 %>% filter(trait == "hd"), mapping = aes(x = Cohort, y = mean, fill = factor(trial, trial_types))) + labs(fill = "Trial Type") + ylab("Mean Days to Heading") + geom_boxplot(outlier.shape = NA) + theme_light() + scale_fill_manual(values = c("#9BBB59", "#F79646", "#8064A2")) + geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 1)
box_hd

# Make a composite image!
box_legend <- get_legend(box_y)

plot_grid(box_y + theme(legend.position = "none"),
          box_ph + theme(legend.position = "none"),
          box_sw + theme(legend.position = "none"),
          box_tw + theme(legend.position = "none"),
          box_hd + theme(legend.position = "none"),
          box_legend, labels = c("A", "B", "C", "D", "E"))
ggsave("Figures/trait/trait_boxplots.pdf", height = 8.5, width = 14)
ggsave("Figures/trait/trait_boxplots.svg", height = 8.5, width = 14)

# 3) Trait correlations ----
# Make trait vector
trait_vec <- c("yield", "height", "sw", "tw", "hd")

# Calculate correlations of traits with themselves across PYTs and AYTs
pyt_ayt_cors <- lapply(c(1:5), function(x){
  cor(get("means_pyts_ayts")[paste0("mean_PYT_", trait_vec[x])], get("means_pyts_ayts")[paste0("mean_AYT_", trait_vec[x])])
})

# Calculate correlations of traits with other traits within PYTs and within AYTs.
# "Loop" through each trait
pyt_cors <- lapply(c(1:5), function(x){
  cor(get("gp_data")[paste0("mean_PYT_", trait_vec[x])], get("gp_data")[paste0("mean_PYT_", trait_vec[-x])])
})

ayt_cors <- lapply(c(1:5), function(x){
  cor(get("means_pyts_ayts")[paste0("mean_AYT_", trait_vec[x])], get("means_pyts_ayts")[paste0("mean_AYT_", trait_vec[-x])])
})

# Save all the correlations
sink("Data/pheno/trait_cor_emmeans.txt")
print("Trait correlations using entry adjusted means")
print("Correlations of PYT traits with the same trait in AYTs")
print(pyt_ayt_cors)
print("Correlations of PYT traits with other PYT traits")
print(pyt_cors)
print("Correlations of AYT traits with other AYT traits")
print(ayt_cors)
sink()
