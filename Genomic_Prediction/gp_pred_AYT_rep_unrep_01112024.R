## Genomic predictions and phenotypic correlations between PYTs and AYTs
# 01/11/2024

# Objective: Correlate AYT performance with PYT genomic predictions and phenotypes, for unreplicated and replicated observations. 

library(tidyverse)
library(lme4)
library(emmeans)
library(rrBLUP)
library(cowplot)
library(svglite)

# Read in data
gp_data <- read.csv("Data/gp/GS_PYT_PYTu_AYT_alltraits_5899.csv")

# Separate the data into cohorts
data_2018 <- gp_data %>%
  filter(Cohort == "2018")
data_2019 <- gp_data %>%
  filter(Cohort == "2019")
data_2020 <- gp_data %>%
  filter(Cohort == "2020")

# Marker matrices for cohorts (remove trait data)
markers_2018 <- as.matrix(data_2018[, -c(1:23)])
markers_2019 <- as.matrix(data_2019[, -c(1:23)])
markers_2020 <- as.matrix(data_2020[, -c(1:23)])

# Remove markers from cohort data
data_2018 <- data_2018 %>%
  select(colnames(gp_data)[1:23])
data_2019 <- data_2019 %>%
  select(colnames(gp_data)[1:23])
data_2020 <- data_2020 %>%
  select(colnames(gp_data)[1:23])

# Get indices of AYT lines
vs_ind_2018 <- which(!is.na(data_2018$mean_AYT_yield))
vs_ind_2019 <- which(!is.na(data_2019$mean_AYT_yield))
vs_ind_2020 <- which(!is.na(data_2020$mean_AYT_yield))

# Genomic prediction function
genomic_prediction <- function(data, markers, pyt_var, ayt_var, vs_indices){
  #set parameters
  y <- noquote(pyt_var)
  data[[y]] <- data[[pyt_var]]
  ayt <- noquote(ayt_var)
  data[[ayt]] <- data[[ayt_var]]
  #fit marker effects and make predictions for the fold
  #NOTE: fitting the model on the whole set! ie. not excluding AYT lines.
  marker_eff <- mixed.solve(data[[y]], Z = markers, SE = FALSE)
  #predicted values obtained by multiplying genotype matrix with marker effects and adding mean.
  #NOTE: only need to obtain predicted values for AYT lines
  pred <- markers[vs_indices,] %*% marker_eff$u + c(marker_eff$beta)
  #make dataframe to return predicted and observed values, as well as fold_id, to use for grouping when taking mean.
  res <- data.frame(NAME = data$NAME[vs_indices],
                    Cohort = data$Cohort[vs_indices], 
                    predicted_val = pred,
                    obs_val = data[[y]][vs_indices],
                    ayt_val = data[[ayt]][vs_indices]) %>%
    rownames_to_column("obs_id")
  return(res)}

# Apply genomic prediction to each cohort for each trait:
# Predictions for original traits:
# Cohort vector
cohort_vec <- c("2018", "2019", "2020")
# Trait vector
orig_trait_vec <- c("yield", "height", "sw", "tw", "hd")

orig_trait_res <- lapply(c(1:3), function(x){
  lapply(c(1:5), function(z){
    genomic_prediction(data = get(paste0("data_", cohort_vec[x])),
                       markers = get(paste0("markers_", cohort_vec[x])),
                       pyt_var = paste0("mean_PYT_", orig_trait_vec[z]),
                       ayt_var = paste0("mean_AYT_", orig_trait_vec[z]),
                       vs_indices = get(paste0("vs_ind_", cohort_vec[x])))
  })
})

# Predictions for unrep yield, height, days to heading:
unrep_trait_vec <- paste0("unrep", rep(c("1", "2", "3"), times = 3), "_PYT_", rep(c("yield", "height", "hd"), each = 3))

unrep_trait_res <- lapply(c(1:3), function(x){
  lapply(c(1:9), function(z){
    genomic_prediction(data = get(paste0("data_", cohort_vec[x])),
                       markers = get(paste0("markers_", cohort_vec[x])),
                       pyt_var = unrep_trait_vec[z],
                       ayt_var = paste0("mean_AYT_", gsub("unrep[1-9]_PYT_", "", unrep_trait_vec[z])),
                       vs_indices = get(paste0("vs_ind_", cohort_vec[x])))
  })
})

# Obtain dataframes with correlations
orig_trait_df <- map_dfr(c(1:3), function(x){
  cohort_df <- map_dfr(c(1:5), function(z){
    pluck(orig_trait_res, x) %>%
      pluck(z) %>%
      summarise(gp_cor = cor(predicted_val, ayt_val),
                pheno_cor = cor(obs_val, ayt_val))
  })
  cohort_df$trait <- orig_trait_vec
  cohort_df$cohort <- cohort_vec[x]
  return(cohort_df)
}) %>%
  mutate(rep_level = case_when(trait %in% c("yield", "height", "hd") ~ "rep", trait %in% c("sw", "tw") ~ "unrep"))

unrep_trait_df <- map_dfr(c(1:3), function(x){
  cohort_df <- map_dfr(c(1:9), function(z){
    pluck(unrep_trait_res, x) %>%
      pluck(z) %>%
      summarise(gp_cor = cor(predicted_val, ayt_val),
                pheno_cor = cor(obs_val, ayt_val)) %>%
      mutate(trait = gsub("unrep[1-9]_PYT_", "", unrep_trait_vec[z]))
  })
  cohort_df$cohort <- cohort_vec[x]
  return(cohort_df)
})

# For unrep predictions, want to get mean correlations.
mean_unrep_df <- unrep_trait_df %>%
  group_by(cohort, trait) %>%
  summarise(gp_cor = mean(gp_cor),
            pheno_cor = mean(pheno_cor)) %>%
  mutate(rep_level = "unrep")
# Note, keeping names of correlations the same, even though they are means, because I want them on the same graphs.

# Join together and set up for figures
cor_df <- full_join(orig_trait_df, mean_unrep_df) %>%
  pivot_longer(cols = !c(trait, cohort, rep_level), names_to = "cor_type", values_to = "cor") %>%
  mutate(cor_rep_type = factor(paste0(gsub("_cor", "", cor_type), "_", rep_level), levels = c("pheno_unrep", "gp_unrep", "pheno_rep", "gp_rep")))

write.csv(cor_df, "Data/gp/gp_pa_AYT_rep_unrep.csv", row.names = FALSE)

# 3) Make figures ----
figure_yield <- ggplot(data = cor_df %>% filter(trait == "yield"), aes(x = cohort, y = cor, fill = cor_rep_type)) + geom_bar(stat = "identity", position = "dodge") + labs(fill = "PYT Value", caption = "Yield (kg/ha)") + ylab("Correlation (r) with AYT Entry Adjusted Means") + xlab("Cohort") + scale_fill_manual(values = c("#F79646", "#4BACC6", "#C0504D", "#1F497D"), labels = c(paste0("Unreplicated","\n", "Phenotypic"), paste0("Unreplicated", "\n", "Genomic"), paste0("Replicated","\n", "Phenotypic"), paste0("Replicated", "\n", "Genomic"))) + theme_light()
figure_yield

figure_height <- ggplot(data = cor_df %>% filter(trait == "height"), aes(x = cohort, y = cor, fill = cor_rep_type)) + geom_bar(stat = "identity", position = "dodge") + labs(fill = "PYT Value", caption = "Plant Height (cm)") + ylab("Correlation (r) with AYT Entry Adjusted Means") + xlab("Cohort") + scale_fill_manual(values = c("#F79646", "#4BACC6", "#C0504D", "#1F497D"), labels = c(paste0("Unreplicated","\n", "Phenotypic"), paste0("Unreplicated", "\n", "Genomic"), paste0("Replicated","\n", "Phenotypic"), paste0("Replicated", "\n", "Genomic"))) + ylim(c(0, 0.9)) + theme_light()
figure_height

figure_hd <- ggplot(data = cor_df %>% filter(trait == "hd"), aes(x = cohort, y = cor, fill = cor_rep_type)) + geom_bar(stat = "identity", position = "dodge") + labs(fill = "PYT Value", caption = "Days to Heading") + ylab("Correlation (r) with AYT Entry Adjusted Means") + xlab("Cohort") + scale_fill_manual(values = c("#F79646", "#4BACC6", "#C0504D", "#1F497D"), labels = c(paste0("Unreplicated","\n", "Phenotypic"), paste0("Unreplicated", "\n", "Genomic"), paste0("Replicated","\n", "Phenotypic"), paste0("Replicated", "\n", "Genomic"))) + ylim(c(0, 0.9)) +  theme_light()
figure_hd

figure_sw <- ggplot(data = cor_df %>% filter(trait == "sw"), aes(x = cohort, y = cor, fill = cor_rep_type)) + geom_bar(stat = "identity", position = "dodge") + labs(fill = "PYT Value", caption = "Seed Weight (g)") + ylab("Correlation (r) with AYT Entry Adjusted Means") + xlab("Cohort") + scale_fill_manual(values = c("#F79646", "#4BACC6"), labels = c(paste0("Unreplicated","\n", "Phenotypic"), paste0("Unreplicated", "\n", "Genomic"))) + ylim(c(0, 0.9)) + theme_light()
figure_sw

figure_tw <- ggplot(data = cor_df %>% filter(trait == "tw"), aes(x = cohort, y = cor, fill = cor_rep_type)) + geom_bar(stat = "identity", position = "dodge") + labs(fill = "PYT Value", caption = "Test Weight (kg/hL)") + ylab("Correlation (r) with AYT Entry Adjusted Means") + xlab("Cohort") + scale_fill_manual(values = c("#F79646", "#4BACC6"), labels = c(paste0("Unreplicated","\n", "Phenotypic"), paste0("Unreplicated", "\n", "Genomic"))) + ylim(c(0, 0.9)) + theme_light()
figure_tw

figure_legend <- get_legend(figure_yield)

plot_grid(figure_yield + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          figure_height + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          figure_sw + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          figure_tw + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          figure_hd + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          figure_legend, labels = c("A", "B", "C", "D", "E"))
ggsave("Figures/gp/ayt_pred_rep_unrep.pdf", height = 8.5, width = 11)
ggsave("Figures/gp/ayt_pred_rep_unrep.svg", height = 8.5, width = 11)