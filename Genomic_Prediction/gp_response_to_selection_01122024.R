## Response to Selection with Phenotypic and Genomic Predictions
# 01/12/2024

# Objectives: Fit a model on entry adjusted means with a trial effect, to get emmeans of genotypes in trials (control for differences in means between trial environments). Then calculate the response to selection. 

library(tidyverse)
library(emmeans)
library(rrBLUP)
library(cowplot)

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

# Remove extra columns from cohort data (unrep PYTs and AYT traits, markers)
data_2018 <- data_2018 %>%
  select(colnames(gp_data)[1:9])
data_2019 <- data_2019 %>%
  select(colnames(gp_data)[1:9])
data_2020 <- data_2020 %>%
  select(colnames(gp_data)[1:9])

## Fit model ----
# Make vector of traits to apply functions with
trait_vec <- c("yield", "height", "sw", "tw", "hd")

# Need to reshape data: entry name, if observation is from PYT or AYT, and trait data.
# First exclude markers and unrep values
# Then pivot, separating trial and trait into two variables + values
# Remove NAs from PYT lines with no AYT values
model_data <- gp_data %>%
  select(colnames(gp_data)[1:9], colnames(gp_data)[19:23]) %>%
  pivot_longer(!c(NAME, Cohort, cross, ayt), names_prefix = "mean_", names_to = c("trial_type", "trait"), names_sep = "_", values_to = "value") %>%
  filter(!is.na(value))

# Formula to apply is y = u + trial + geno + error
# Note: also may want to try with a trial x geno interaction
model_lm <- formula(value ~ trial_type + NAME)

# Apply models for all traits
lapply(c(1:5), function(x){
  assign(paste0("lm_", trait_vec[x]),
         lm(formula = model_lm, data = model_data[model_data$trait == trait_vec[x],]),
         envir = .GlobalEnv)
})

# Apply summary and anova to models
lapply(c(1:5), function(x){
  anova(get(paste0("lm_", trait_vec[x])))
})

lapply(c(1:5), function(x){
  summary(get(paste0("lm_", trait_vec[x])))
})

# Get emmeans. Need the value for each line.
lapply(c(1:5), function(x){
  assign(paste0("emmeans_", trait_vec[x]),
         as.data.frame(emmeans(get(paste0("lm_", trait_vec[x])), "NAME")),
         envir = .GlobalEnv)
})

## Response to selection ----
# 1) Make GS and PS "selections" - Use gp_data
# 2) Calculate GS and PS means for selected lines for each cohort (AYT means) - Use emmeans from trial model
# 3) Calculate the means for all lines in each cohort (PYT means) - Use emmeans from trial model
# 4) Calculate difference in means for each method and trait, then plot.

# Remove checks from emmeans
lapply(c(1:5), function(x){
  assign(paste0("nochecks_", trait_vec[x]),
         get(paste0("emmeans_", trait_vec[x])) %>%
           filter(!NAME %in% c("Leader", "AAC-Ling")),
         envir = .GlobalEnv)
})

# Vector of datasets
data_vec <- c("data_2018", "data_2019", "data_2020")

# 1A) Add phenotypic entry rankings (based on yield)
lapply(c(1:length(data_vec)), function(x){
  assign(data_vec[x], 
         get(data_vec[x]) %>% 
           mutate(ps_rank = length(get(data_vec[x])[["NAME"]]) - rank(get(data_vec[x])["mean_PYT_yield"]) + 1) %>%
           select(colnames(get(data_vec[x]))[1:2], "ps_rank", colnames(get(data_vec[x]))[3:length(colnames(get(data_vec[x])))]), 
         envir = .GlobalEnv)
})

# 1B) Do genomic predictions
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

# Apply genomic prediction function
lapply(c(1:length(data_vec)), function(x){
  assign(paste0("gp_", data_vec[x]), 
         genomic_prediction(data = get(data_vec[x]), 
                            markers = get(gsub("data", "markers", data_vec[x])),
                            pyt_var = "mean_PYT_yield",
                            ayt_var = "mean_AYT_yield",
                            vs_indices = which(get(data_vec[x])["Cohort"] == gsub("data_", "", data_vec[x]))),
         envir = .GlobalEnv)
})
# Note: vs_indices is all lines from a cohort (ie. fitting markers and estimating on same lines).

# 1C) Add genomic prediction entry rankings (again, just yield)
# Use the predicted_val to rank entries, then join the rankings with the data.
# Add entry rankings to gp result dataframes
lapply(c(1:length(data_vec)), function(x){
  assign(paste0("gp_", data_vec[x]),
         get(paste0("gp_", data_vec[x])) %>%
           mutate(gp_rank = length(get(paste0("gp_", data_vec[x]))[["NAME"]]) - rank(get(paste0("gp_", data_vec[x]))["predicted_val"]) + 1),
         envir = .GlobalEnv)
})
# Now make named vectors with entry rankings
lapply(c(1:length(data_vec)), function(x){
  assign(paste0("rank_gp_", data_vec[x]),
         as.vector(get(paste0("gp_", data_vec[x]))[["gp_rank"]]),
         envir = .GlobalEnv)
})
lapply(c(1:length(data_vec)), function(x){
  assign(paste0("rank_gp_", data_vec[x]),
         set_names(get(paste0("rank_gp_", data_vec[x])), nm = get(paste0("gp_", data_vec[x]))[["NAME"]]),
         envir = .GlobalEnv)
})

# Add the rankings to the original dataframes
lapply(c(1:length(data_vec)), function(x){
  assign(data_vec[x],
         get(data_vec[x]) %>%
           mutate(gp_rank = get(paste0("rank_gp_", data_vec[x]))[get(data_vec[x])[["NAME"]]]) %>%
           select(colnames(get(data_vec[x]))[1:3], "gp_rank", colnames(get(data_vec[x]))[4:length(colnames(get(data_vec[x])))]),
         envir = .GlobalEnv)
})

# 1D) Make "selections"
# Filter out lines that aren't advanced.
lapply(c(1:length(data_vec)), function(x){
  assign(paste0("ayt_", data_vec[x]),
         get(data_vec[x]) %>%
           filter(!is.na(mean_AYT_yield)),
         envir = .GlobalEnv)
})

# Objective is to get rankings that occur for both GP and PS. Then obtain means of entries with those rankings for GP and PS.
# Get rankings that are common to both GP and PS.
lapply(c(1:length(data_vec)), function(x){
  assign(paste0("comrank_", data_vec[x]),
         unname(get(paste0("ayt_", data_vec[x]))[["ps_rank"]][get(paste0("ayt_", data_vec[x]))[["ps_rank"]] %in% get(paste0("ayt_", data_vec[x]))[["gp_rank"]]]),
         envir = .GlobalEnv)
})

# Filter out the rankings above 30.
lapply(c(1:length(data_vec)), function(x){
  assign(paste0("rankfilt_", data_vec[x]),
         get(paste0("comrank_", data_vec[x]))[get(paste0("comrank_", data_vec[x])) < 31],
         envir = .GlobalEnv)
})

# 1E) Extract line names for PS and GS
lapply(c(1:length(data_vec)), function(x){
  assign(paste0("ps_names_", data_vec[x]),
         get(data_vec[x])[get(data_vec[x])[["ps_rank"]] %in% get(paste0("comrank_", data_vec[x])), 1],
         envir = .GlobalEnv)
})


lapply(c(1:length(data_vec)), function(x){
  assign(paste0("gp_names_", data_vec[x]),
         get(data_vec[x])[get(data_vec[x])[["gp_rank"]] %in% get(paste0("comrank_", data_vec[x])), 1],
         envir = .GlobalEnv)
})

# 2) Calculate "AYT" means for PS and GS (use the trial model emmeans)
ayt_means_ps <- as.data.frame(do.call(rbind, lapply(c(1:5), function(x){
  do.call(rbind, lapply(c(1:3), function(z){
    data.frame(cohort = gsub("data_", "", data_vec[z]),
               trait = trait_vec[x],
               ayt_mean = mean(get(paste0("nochecks_", trait_vec[x]))[get(paste0("nochecks_", trait_vec[x]))[["NAME"]] %in% get(paste0("ps_names_", data_vec[z])),][["emmean"]]))
  }))
})))

ayt_means_gs <- as.data.frame(do.call(rbind, lapply(c(1:5), function(x){
  do.call(rbind, lapply(c(1:3), function(z){
    data.frame(cohort = gsub("data_", "", data_vec[z]),
               trait = trait_vec[x],
               ayt_mean = mean(get(paste0("nochecks_", trait_vec[x]))[get(paste0("nochecks_", trait_vec[x]))[["NAME"]] %in% get(paste0("gp_names_", data_vec[z])),][["emmean"]]))
  }))
})))

# 3) Calculate PYT means (again using trial model emmeans)
pyt_means <- as.data.frame(do.call(rbind, lapply(c(1:5), function(x){
  do.call(rbind, lapply(c(1:3), function(z){
    data.frame(cohort = gsub("data_", "", data_vec[z]),
               trait = trait_vec[x],
               pyt_mean = mean(get(paste0("nochecks_", trait_vec[x]))[get(paste0("nochecks_", trait_vec[x]))[["NAME"]] %in% get(data_vec[z])[["NAME"]],][["emmean"]]))
  }))
})))

# 4A) Calculate change in means between PYT and AYT
ps_results <- left_join(pyt_means, ayt_means_ps) %>%
  mutate(change_mean = ayt_mean - pyt_mean,
         method = "PS")

gs_results <- left_join(pyt_means, ayt_means_gs) %>%
  mutate(change_mean = ayt_mean - pyt_mean,
         method = "GS")

# Merge together
results <- full_join(ps_results, gs_results) %>%
  mutate(method = factor(method, levels = c("PS", "GS")))

# 4B) Plots
fig_rts_y <- ggplot(data = results[results$trait == "yield",], mapping = aes(x = cohort, y = change_mean, fill = method)) + geom_bar(stat = "identity", position = "dodge") + xlab("Cohort") + ylab("Change in Mean Yield (kg/ha)") + labs(fill = "Selection Method") + theme_light() + scale_fill_manual(values = c("#9BBB59", "#4BACC6"))
fig_rts_y

fig_rts_ph <- ggplot(data = results[results$trait == "height",], mapping = aes(x = cohort, y = change_mean, fill = method)) + geom_bar(stat = "identity", position = "dodge") + xlab("Cohort") + ylab("Change in Mean Plant Height (cm)") + labs(fill = "Selection Method") + theme_light() + scale_fill_manual(values = c("#9BBB59", "#4BACC6"))
fig_rts_ph

fig_rts_sw <- ggplot(data = results[results$trait == "sw",], mapping = aes(x = cohort, y = change_mean, fill = method)) + geom_bar(stat = "identity", position = "dodge") + xlab("Cohort") + ylab("Change in Mean Seed Weight (g)") + labs(fill = "Selection Method") + theme_light() + scale_fill_manual(values = c("#9BBB59", "#4BACC6"))
fig_rts_sw

fig_rts_tw <- ggplot(data = results[results$trait == "tw",], mapping = aes(x = cohort, y = change_mean, fill = method)) + geom_bar(stat = "identity", position = "dodge") + xlab("Cohort") + ylab("Change in Mean Test Weight (kg/hL)") + labs(fill = "Selection Method") + theme_light() + scale_fill_manual(values = c("#9BBB59", "#4BACC6"))
fig_rts_tw

fig_rts_hd <- ggplot(data = results[results$trait == "hd",], mapping = aes(x = cohort, y = change_mean, fill = method)) + geom_bar(stat = "identity", position = "dodge") + xlab("Cohort") + ylab("Change in Mean Days to Heading") + labs(fill = "Selection Method") + theme_light() + scale_fill_manual(values = c("#9BBB59", "#4BACC6"))
fig_rts_hd

# Make composite plot
fig_legend <- get_legend(fig_rts_y)

plot_grid(fig_rts_y + theme(legend.position = "none"),
          fig_rts_ph + theme(legend.position = "none"),
          fig_rts_sw + theme(legend.position = "none"),
          fig_rts_tw + theme(legend.position = "none"),
          fig_rts_hd + theme(legend.position = "none"),
          fig_legend, labels = c("A", "B", "C", "D", "E"))
ggsave("Figures/gp/resptosel_wtrialmodel.pdf", height = 8.5, width = 11)
ggsave("Figures/gp/resptosel_wtrialmodel.svg", height = 8.5, width = 11)

# Save the rankings and number of rankings used for each cohort
sink("Data/gp/resptosel_ranks_trialmodel.txt")
print("Ranks and n for each cohort (2018, 2019, 2020) for response to selection.")
print(rankfilt_data_2018)
print(rankfilt_data_2019)
print(rankfilt_data_2020)
print(length(rankfilt_data_2018))
print(length(rankfilt_data_2019))
print(length(rankfilt_data_2020))
print("Ranks for appendix table")
paste(noquote(sort(rankfilt_data_2018)), collapse = ", ")
paste(noquote(sort(rankfilt_data_2019)), collapse = ", ")
paste(noquote(sort(rankfilt_data_2020)), collapse = ", ")
sink()
