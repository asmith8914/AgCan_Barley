## Genomic prediction with random CV vs. stratified by cycle
# 01/08/2024

# Objectives: Predict a PYT using either the other 2 PYTs as the training set (stratified by breeding cycle) or a randomly sampled training set of the same size which includes some lines from validation set PYT, with sub-sampling of the cohort for the validation set.

library(tidyverse)
library(rrBLUP)
library(furrr)
library(cowplot)

# Check RNGkind
RNGkind()
# "Mersenne-Twister" "Inversion"        "Rejection"   
RNGkind("L'Ecuyer-CMRG")
# Set as "L'Ecuyer-CMRG" because genomic prediction functions will use functions that work in parallel, so this is required for reproducibility. 

# Read in data
gp_data <- read.csv("Data/gp/GS_PYT_PYTu_AYT_alltraits_5899.csv")

# Remove checks
gp_data2 <- gp_data %>%
  filter(Cohort %in% c(2018, 2019, 2020))

# Marker matrix for all three cohorts (no checks)
gp_markers <- as.matrix(gp_data2[, -c(1:23)])

## Functions: ----
# Function to assign validation sets
# Inputs: data, validation set cohort (character in quotes), number of folds, number of reps, and seed (integer of format "[int]L").
# Outputs: tibble with observation ids, ids of validation cohort lines, fold ids, and iteration ids.
assign_vs <- function(cv_data, vs_cohort, num_folds = 5, num_reps = 10, set_seed){
  N = nrow(cv_data[cv_data$Cohort == vs_cohort,])
  # vector with folds to sample from
  folds_to_sample <- rep(1:num_folds, length.out = N)
  future_map_dfr(c(1:num_reps), .options = furrr_options(seed = set_seed), function(i){
    print(paste("Iteration", i))
    # randomize the folds to create sets
    sets <- sample(folds_to_sample)
    # output: table with index, iteration
    folds_table <- data.frame(data_id = which(cv_data$Cohort == vs_cohort),
                              fold_id = sets)  %>%
      rownames_to_column("obs_id")
    folds_table$iter <- i
    return(as_tibble(folds_table))}, .progress = TRUE)
}

# Note: a function to assign the training set is not needed for the stratified by breeding cycle scenario, as there is only one possible set.

# Function to perform genomic prediction with training set as two cohorts and provided validation set
# Inputs: data, markers, training set cohorts (as c(character in quotes)), validation sets (output of assign_vs), and phenotype variable (in quotes).
# Outputs: a tibble with entry name, cohort, predicted value, observed value, and fold and iteration ids.
gp_strat_ts_vs <- function(data, markers, ts_cohorts, v_sets, phenotype_var){
  # number of folds and iterations
  num_iter = length(unique(v_sets$iter))
  num_folds = length(unique(v_sets$fold_id))
  # Phenotype data
  y <- noquote(phenotype_var)
  data[[y]] <- data[[phenotype_var]]
  # Get the training set indices
  ts_ind <- which(data$Cohort %in% ts_cohorts)
  # Fit the model on the TS indices
  marker_eff <- mixed.solve(data[[y]][ts_ind], Z = markers[ts_ind,], SE = FALSE)
  # Make predictions for VS indices for each iteration and each fold:
  # First "loop" through iterations:
  c(1:num_iter) %>%
    map_dfr(function(iter){
      print(paste("Iteration", iter))
      #subset the validation sets for each iteration
      v_sets_iter <- v_sets[v_sets$iter == iter, ]
      #then "loop" through folds
      cv_res <- c(1:num_folds) %>% 
        map(function(fold){
          print(paste("Fold", fold))
          #make vector of indices for validation set
          vs_ind <- v_sets_iter$data_id[which(v_sets_iter$fold_id == fold)]
          #predict the validation set values
          pred <- markers[vs_ind,] %*% marker_eff$u + c(marker_eff$beta)
          #make dataframe to return predicted and observed values, as well as fold_id, to use for grouping when taking mean.
          res <- data.frame(NAME = data$NAME[vs_ind],
                            Cohort = data$Cohort[vs_ind],
                            predicted_val = pred,
                            obs_val = data[[y]][vs_ind],
                            fold_id = fold) %>%
            rownames_to_column("obs_id")
          return(res)}) %>%
        bind_rows()
      cv_res$iter <- iter
      return(as_tibble(cv_res))
    })
}
# Note: not the most efficient method, since predicted values are being calculated multiple times for the same entries with the same markers (could just calculate once and then add to the validation set folds, ids, etc.)

# Function to assign randomly sampled training sets (the same size as the stratified by cycle training sets).
# Inputs: data, training set cohorts (as c(character in quotes)), validation sets (output of assign_vs), and seed (integer of format "[int]L").
# Outputs: a dataframe with observation ids, ids of training set lines, fold ids, and iteration ids.
assign_rand_ts <- function(data, ts_cohorts, v_sets, set_seed){
  # number of folds and iterations
  num_iter = length(unique(v_sets$iter))
  num_folds = length(unique(v_sets$fold_id))
  # Determine size of training set used for stratified predictions
  N <- nrow(data[data$Cohort %in% ts_cohorts,])
  # For each iteration and fold:
  future_map_dfr(c(1:num_iter), .options = furrr_options(seed = set_seed), function(iter){
    # subset the training and validation sets for each iteration
    v_sets_iter <- v_sets[v_sets$iter == iter, ]
    obj <- c(1:num_folds) %>%
      map(function(fold){
        # obtain the validation set indices corresponding to the data
        vs_ind <- v_sets_iter$data_id[v_sets_iter$fold_id == fold]
        # get indices of lines that can be sampled for training set
        pos_ts <- which(!(rownames(data) %in% vs_ind))
        #randomly sample a training set which excludes the validation set lines
        ts_ind <- sample(x = pos_ts, size = N, replace = FALSE)
        # dataframe to return
        t_sets <- data.frame(data_id = ts_ind,
                             fold_id = fold) %>%
          rownames_to_column("obs_id")
        t_sets$iter <- iter
        return(t_sets)
      })
  })
}

# Function to perform genomic prediction with assigned randomly sampled training and validation sets.
# Inputs: data, markers, training sets (output of assign_rand_ts), validation sets (output of assign_vs), and phenotype variable (in quotes).
# Outputs: a tibble with entry name, predicted value, observed value, and fold and iteration ids.
gp_rand_ts_vs <- function(data, markers, t_sets, v_sets, phenotype_var){
  # number of folds and iterations
  num_iter = length(unique(v_sets$iter))
  num_folds = length(unique(v_sets$fold_id))
  # Phenotype data
  y <- noquote(phenotype_var)
  data[[y]] <- data[[phenotype_var]]
  # do genomic prediction for each combination of iterations and folds: there is a unique training set and a unique validation set for each iter_fold.
  # first loop through the iterations
  c(1:num_iter) %>%
    map_dfr(function(iter){
      print(paste("Iteration", iter))
      #subset the training and validation sets for each iteration
      v_sets_iter <- v_sets[v_sets$iter == iter, ]
      t_sets_iter <- t_sets[t_sets$iter == iter, ]
      #then loop through folds
      cv_res <- c(1:num_folds) %>% 
        map(function(fold){
          print(paste("Cross-Validating Fold", fold))
          #get indices for training and validation sets
          ts_ind <- t_sets_iter$data_id[t_sets_iter$fold_id == fold]
          vs_ind <- v_sets_iter$data_id[v_sets_iter$fold_id == fold]
          #fit the model on the training set data, subset to the training set indices.
          marker_eff <- mixed.solve(data[[y]][ts_ind], Z = markers[ts_ind,], SE = FALSE)
          #predict the validation set values
          pred <- markers[vs_ind,] %*% marker_eff$u + c(marker_eff$beta)
          #make dataframe to return predicted and observed values, as well as fold_id, to use for grouping when taking mean.
          res <- data.frame(NAME = data$NAME[vs_ind],
                            predicted_val = pred,
                            obs_val = data[[y]][vs_ind],
                            fold_id = fold) %>%
            rownames_to_column("obs_id")
          return(res)}) %>%
        bind_rows()
      cv_res$iter <- iter
      return(as_tibble(cv_res))
    })
}
# Quite slow... 

## Predictions: ----
# Make a vector of traits
trait_vec <- c("yield", "height", "sw", "tw", "hd")

# Assign the validation sets
vs_2018 <- assign_vs(cv_data = gp_data2, vs_cohort = "2018", set_seed = 1L)
vs_2019 <- assign_vs(cv_data = gp_data2, vs_cohort = "2019", set_seed = 2L)
vs_2020 <- assign_vs(cv_data = gp_data2, vs_cohort = "2020", set_seed = 3L)

# Assign the randomly sampled training sets
ts_2018 <- assign_rand_ts(data = gp_data2, ts_cohorts = c("2019", "2020"), v_sets = vs_2018, set_seed = 1L)
ts_2019 <- assign_rand_ts(data = gp_data2, ts_cohorts = c("2018", "2020"), v_sets = vs_2019, set_seed = 2L)
ts_2020 <- assign_rand_ts(data = gp_data2, ts_cohorts = c("2018", "2019"), v_sets = vs_2020, set_seed = 3L)

# Vector of years
year_vec <- c("2018", "2019", "2020")

# Apply the stratified prediction function to all traits and cohorts
lapply(c(1:length(year_vec)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0("strat_", year_vec[x], "_", trait_vec[z]),
           gp_strat_ts_vs(data = gp_data2,
                          markers = gp_markers,
                          v_sets = get(paste0("vs_", year_vec[x])),
                          ts_cohorts = year_vec[-x],
                          phenotype_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
  })
})

# Apply the random sampling prediction function to all traits and cohorts.
# Since function is slow, only going to "loop" through traits, do cohorts separately.
# 2018
lapply(c(1:length(trait_vec)), function(z){
    assign(paste0("rand_", year_vec[1], "_", trait_vec[z]),
           gp_rand_ts_vs(data = gp_data2,
                         markers = gp_markers,
                         v_sets = get(paste0("vs_", year_vec[1])),
                         t_sets = get(paste0("ts_", year_vec[1])),
                         phenotype_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
})
# 2019
lapply(c(1:length(trait_vec)), function(z){
  assign(paste0("rand_", year_vec[2], "_", trait_vec[z]),
         gp_rand_ts_vs(data = gp_data2,
                       markers = gp_markers,
                       v_sets = get(paste0("vs_", year_vec[2])),
                       t_sets = get(paste0("ts_", year_vec[2])),
                       phenotype_var = paste0("mean_PYT_", trait_vec[z])),
         envir = .GlobalEnv)
})
# 2020
lapply(c(1:length(trait_vec)), function(z){
  assign(paste0("rand_", year_vec[3], "_", trait_vec[z]),
         gp_rand_ts_vs(data = gp_data2,
                       markers = gp_markers,
                       v_sets = get(paste0("vs_", year_vec[3])),
                       t_sets = get(paste0("ts_", year_vec[3])),
                       phenotype_var = paste0("mean_PYT_", trait_vec[z])),
         envir = .GlobalEnv)
})

## Obtain correlations and make figures ----
# Function to calculate correlations
gp_cor <- function(gp_result){
  gp_result %>%
    group_by(iter, fold_id) %>%
    summarise(fold_pa = cor(predicted_val, obs_val)) %>%
    group_by(iter) %>%
    summarise(iter_pa = mean(fold_pa)) %>%
    summarise(overall_pa = mean(iter_pa), stdev = sqrt(var(iter_pa)))
}

# Apply function to results of prediction with stratified and random training sets
strat_res <- do.call(rbind, lapply(c(1:length(year_vec)), function(x){
  do.call(rbind, lapply(c(1:length(trait_vec)), function(z){
    gp_cor(get(paste0("strat_", year_vec[x], "_", trait_vec[z])))
  })) %>%
    mutate(trait = trait_vec, vs = year_vec[x], method = "strat") %>%
    select(vs, trait, overall_pa, stdev, method)
}))

rand_res <- do.call(rbind, lapply(c(1:length(year_vec)), function(x){
  do.call(rbind, lapply(c(1:length(trait_vec)), function(z){
    gp_cor(get(paste0("rand_", year_vec[x], "_", trait_vec[z])))
  })) %>%
    mutate(trait = trait_vec, vs = year_vec[x], method = "rand") %>%
    select(vs, trait, overall_pa, stdev, method)
}))

# Combine
results <- rbind(strat_res, rand_res)

# Make figures
fig_strat_rand_y <- ggplot(data = results %>% filter(trait == "yield"), mapping = aes(x = vs, y = overall_pa, fill = method)) + geom_bar(stat = "identity", position = "dodge") + xlab("Validation Set Size/Validation Set") + ylab("Mean Correlation") + labs(fill = "TS/VS Composition") + scale_fill_manual(values = c("#C0504D", "#1F497D"), labels = c(paste0("Randomly Sampled", "\n", "(Size Controlled)"), paste0("Stratified by", "\n", "Breeding Cycle"))) + geom_errorbar(aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + theme_light()
fig_strat_rand_y

fig_strat_rand_ph <- ggplot(data = results %>% filter(trait == "height"), mapping = aes(x = vs, y = overall_pa, fill = method)) + geom_bar(stat = "identity", position = "dodge") + xlab("Validation Set Size/Validation Set") + ylab("Mean Correlation") + labs(fill = "TS/VS Composition") + scale_fill_manual(values = c("#C0504D", "#1F497D"), labels = c(paste0("Randomly Sampled", "\n", "(Size Controlled)"), paste0("Stratified by", "\n", "Breeding Cycle"))) + geom_errorbar(aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + theme_light()
fig_strat_rand_ph

fig_strat_rand_sw <- ggplot(data = results %>% filter(trait == "sw"), mapping = aes(x = vs, y = overall_pa, fill = method)) + geom_bar(stat = "identity", position = "dodge") + xlab("Validation Set Size/Validation Set") + ylab("Mean Correlation") + labs(fill = "TS/VS Composition") + scale_fill_manual(values = c("#C0504D", "#1F497D"), labels = c(paste0("Randomly Sampled", "\n", "(Size Controlled)"), paste0("Stratified by", "\n", "Breeding Cycle"))) + geom_errorbar(aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + theme_light()
fig_strat_rand_sw

fig_strat_rand_tw <- ggplot(data = results %>% filter(trait == "tw"), mapping = aes(x = vs, y = overall_pa, fill = method)) + geom_bar(stat = "identity", position = "dodge") + xlab("Validation Set Size/Validation Set") + ylab("Mean Correlation") + labs(fill = "TS/VS Composition") + scale_fill_manual(values = c("#C0504D", "#1F497D"), labels = c(paste0("Randomly Sampled", "\n", "(Size Controlled)"), paste0("Stratified by", "\n", "Breeding Cycle"))) + geom_errorbar(aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + theme_light()
fig_strat_rand_tw

fig_strat_rand_hd <- ggplot(data = results %>% filter(trait == "hd"), mapping = aes(x = vs, y = overall_pa, fill = method)) + geom_bar(stat = "identity", position = "dodge") + xlab("Validation Set Size/Validation Set") + ylab("Mean Correlation") + labs(fill = "TS/VS Composition") + scale_fill_manual(values = c("#C0504D", "#1F497D"), labels = c(paste0("Randomly Sampled", "\n", "(Size Controlled)"), paste0("Stratified by", "\n", "Breeding Cycle"))) + geom_errorbar(aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + theme_light()
fig_strat_rand_hd

strat_rand_legend <- get_legend(fig_strat_rand_hd)

plot_grid(fig_strat_rand_y + ylim(0,0.75) + labs(caption = "Yield (kg/ha)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)), 
          fig_strat_rand_ph + ylim(0,0.75) + labs(caption = "Plant Height (cm)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)), 
          fig_strat_rand_sw + ylim(0,0.75) + labs(caption = "Seed Weight (g)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)), 
          fig_strat_rand_tw + ylim(0,0.75) + labs(caption = "Test Weight (kg/hL)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)), 
          fig_strat_rand_hd + ylim(0,0.75) + labs(caption = "Days to Heading") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          strat_rand_legend,
          labels = c("A", "B", "C", "D", "E"))
ggsave("Figures/gp/strat_vs_rand.pdf", height = 8.5, width = 11)
ggsave("Figures/gp/strat_vs_rand.svg", height = 8.5, width = 11)
