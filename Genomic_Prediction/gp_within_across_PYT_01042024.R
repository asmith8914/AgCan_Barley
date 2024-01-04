## Genomic prediction within and across PYTs
# 01/04/2024

# Objective: Perform genomic predictions within and across preliminary trials using 5899 marker set.

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

## Genomic prediction functions ----
# 2 tasks: assign folds and perform predictions

# This function assigns folds for multiple runs of cross validation within a dataset.
# Inputs: the data, number of folds (default 5), number of reps (default 10), the seed.
# Outputs: a tibble with the observation id (index to use for matching with data), the fold id (used to separate data into folds), and the iteration id (used to separate iterations).
# Note that results will be the same if seed and dataset size are the same (should change seed between datasets).
assign_sets_within <- function(cv_data, num_folds = 5, num_reps = 10, set_seed){
  N = nrow(cv_data)
  # vector with folds to sample from
  folds_to_sample <- rep(1:num_folds, length.out = N)
  future_map_dfr(c(1:num_reps), .options = furrr_options(seed = set_seed), function(i){
    print(paste("Iteration", i))
    # randomize the folds to create sets
    sets <- sample(folds_to_sample)
    # output: table with index, iteration
    folds_table <- data.frame(fold_id = sets)  %>%
      rownames_to_column("obs_id")
    folds_table$iter <- i
    return(as_tibble(folds_table))}, .progress = TRUE)
}

# This function creates unique training sets of the appropriate size for each fold of each iteration given a dataset and the validation sets output from assign_sets_within.
# Inputs: training set data, validation set output from assign_sets_within, seed (integer of format "[int]L")
# Outputs: a dataframe with __
# Note that results will be the same if seed and dataset size are the same (should change seed between datasets).
assign_sets_across <- function(ts_data, v_sets, set_seed){
  # determine the sizes that the training sets should be
  N = length(unique(v_sets$obs_id))
  num_folds = length(unique(v_sets$fold_id))
  folds <- rep(1:num_folds, length.out = N)
  num_iter = length(unique(v_sets$iter))
  # for each fold, sample the correct number of lines from the ts_data
  future_map_dfr(c(1:num_folds), .options = furrr_options(seed = set_seed), function(fold){
    #determine the size of the training set for each fold as validation set
    ts_ind <- which(folds != fold)
    ts_size <- length(ts_ind)
    #then randomly sample from ts_data to use as the training set.
    #for each fold, make a unique training set for each iteration.
    obj <- c(1:num_iter) %>%
      map(function(iter){
        t_sets <- data.frame(ts = sample(x = 1:nrow(ts_data), size = ts_size),
                             fold_id = fold) %>%
          rownames_to_column("obs_id")
        t_sets$iter <- iter
        return(t_sets)})
  })
}

# This function performs genomic prediction within datasets
# Inputs: data, validation set output from assign_sets_within, marker matrix, name of phenotype variable
# Outputs: a dataframe with
gp_within <- function(cv_data, v_sets, markers, phenotype_var){
  # parameters
  num_iter = length(unique(v_sets$iter))
  num_folds = length(unique(v_sets$fold_id))
  y <- noquote(phenotype_var)
  cv_data[[y]] <- cv_data[[phenotype_var]]
  # need to group by iteration
  obj <- c(1:num_iter) %>%
    map_dfr(function(iter){
      print(paste("Iteration", iter))
      v_sets_iter <- v_sets[v_sets$iter == iter, ]
      # then do cross-validation as usual
      cv_res <- c(1:num_folds) %>% 
        map(function(fold){
          print(paste("Cross-Validating Fold", fold))
          #make vector of indices for training set
          trainIndex <- which(v_sets_iter$fold_id != fold)
          #make a copy of the data
          data1 <- cv_data
          #mask validation set values
          data1[[y]][-trainIndex] <- NA
          #fit marker effects and make predictions for the fold
          marker_eff <- mixed.solve(data1[[y]][trainIndex], Z = markers[trainIndex,], SE = FALSE)
          #predicted values obtained by multiplying genotype matrix with marker effects and adding mean.
          pred <- markers[-trainIndex,] %*% marker_eff$u + c(marker_eff$beta)
          #make dataframe to return predicted and observed values, as well as fold_id, to use for grouping when taking mean.
          res <- data.frame(NAME = cv_data$NAME[-trainIndex],
                            predicted_val = pred,
                            obs_val = cv_data[[y]][-trainIndex],
                            fold_id = fold) %>%
            rownames_to_column("obs_id")
          return(res)}) %>%
        bind_rows()
      cv_res$iter <- iter
      return(as_tibble(cv_res))})
}

# This function performs genomic prediction across datasets
# Inputs: training set data, validation set data, training set markers, validation set markers, training set output from assign_sets_across, validation set output from assign_sets_within, name of phenotype variable.
# Outputs: a dataframe with
gp_across <- function(ts_data, vs_data, ts_markers, vs_markers, t_sets, v_sets, phenotype_var){
  # parameters:
  # number of folds and iterations
  num_iter = length(unique(v_sets$iter))
  num_folds = length(unique(v_sets$fold_id))
  # phenotype data
  y <- noquote(phenotype_var)
  ts_data[[y]] <- ts_data[[phenotype_var]]
  vs_data[[y]] <- vs_data[[phenotype_var]]
  
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
          #make vector of indices for training set
          ts_ind <- t_sets_iter$ts[t_sets_iter$fold_id == fold]
          #make vector of indices for validation set
          vs_ind <- which(v_sets_iter$fold_id == fold)
          #fit the model on the training set data, subset to the training set indices.
          marker_eff <- mixed.solve(ts_data[[y]][ts_ind], Z = ts_markers[ts_ind,], SE = FALSE)
          #predict the validation set values
          pred <- vs_markers[vs_ind,] %*% marker_eff$u + c(marker_eff$beta)
          #make dataframe to return predicted and observed values, as well as fold_id, to use for grouping when taking mean.
          res <- data.frame(NAME = vs_data$NAME[vs_ind],
                            predicted_val = pred,
                            obs_val = vs_data[[y]][vs_ind],
                            fold_id = fold) %>%
            rownames_to_column("obs_id")
          return(res)}) %>%
        bind_rows()
      cv_res$iter <- iter
      return(as_tibble(cv_res))})
}

## Within PYT predictions ----
# Assign the sets 
within_2018_sets <- assign_sets_within(cv_data = data_2018, set_seed = 1L)
within_2019_sets <- assign_sets_within(cv_data = data_2019, set_seed = 2L)
within_2020_sets <- assign_sets_within(cv_data = data_2020, set_seed = 3L)

# Make trait and data vectors
trait_vec <- c("yield", "height", "sw", "tw", "hd")
data_vec <- c("data_2018", "data_2019", "data_2020")

# Apply the prediction function for each phenotype variable and each cohort
lapply(c(1:length(data_vec)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0("within", gsub("data", "", data_vec[x]), "_", trait_vec[z]),
           gp_within(cv_data = get(data_vec[x]),
                     markers = get(gsub("data", "markers", data_vec[x])),
                     v_sets = get(paste0("within", gsub("data", "", data_vec[x]), "_sets")),
                     phenotype_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
  })
})
# This step is a bit slow and could likely be optimized.

## Across PYT predictions ----
# Assign the sets
across_ts20_vs18_sets <- assign_sets_across(v_sets = within_2018_sets, ts_data = data_2020, set_seed = 1L)
across_ts19_vs18_sets <- assign_sets_across(v_sets = within_2018_sets, ts_data = data_2019, set_seed = 2L)
across_ts20_vs19_sets <- assign_sets_across(v_sets = within_2019_sets, ts_data = data_2020, set_seed = 3L)
across_ts18_vs19_sets <- assign_sets_across(v_sets = within_2019_sets, ts_data = data_2018, set_seed = 4L)
across_ts18_vs20_sets <- assign_sets_across(v_sets = within_2020_sets, ts_data = data_2018, set_seed = 5L)
across_ts19_vs20_sets <- assign_sets_across(v_sets = within_2020_sets, ts_data = data_2019, set_seed = 6L)

# A dataset vector for predictions across trials
year_comb <- t(combn(c("2018", "2019", "2020"), 2))
across_vec <- c(paste0(year_comb[,1], "_", year_comb[,2]), paste0(year_comb[,2], "_", year_comb[,1]))
rm(year_comb)

# Apply the prediction function for each phenotype variable and all combinations of cohorts
# Resulting objects have ts first then vs given in name
lapply(c(1:length(across_vec)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0("across_", across_vec[x], "_", trait_vec[z]),
           gp_across(ts_data = get(paste0("data_", gsub("_[0-9]{4}", "", across_vec[x]))),
                     ts_markers = get(paste0("markers_", gsub("_[0-9]{4}", "", across_vec[x]))),
                     t_sets = get(paste0("across_ts", gsub("^20", "", gsub("_[0-9]{4}", "", across_vec[x])), "_vs", gsub("[0-9]{4}_[0-9]{2}", "", across_vec[x]), "_sets")),
                     vs_data = get(paste0("data_", gsub("[0-9]{4}_", "", across_vec[x]))),
                     vs_markers = get(paste0("markers_", gsub("[0-9]{4}_", "", across_vec[x]))),
                     v_sets = get(paste0("within_", gsub("[0-9]{4}_", "", across_vec[x]), "_sets")),
                     phenotype_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
  })
})
# Also slow and likely can be optimized.

## Calculate prediction accuracies and assemble figure ----
# Function to calculate prediction accuracies for given results
gp_pa <- function(gp_result){
  gp_result %>%
    group_by(iter, fold_id) %>%
    summarise(fold_pa = cor(predicted_val, obs_val)) %>%
    group_by(iter) %>%
    summarise(iter_pa = mean(fold_pa)) %>%
    summarise(overall_pa = mean(iter_pa), stdev = sqrt(var(iter_pa)))
}

# Apply prediction accuracy function to within trial results
res_within <- do.call(rbind, lapply(c(1:length(data_vec)), function(x){
  do.call(rbind, lapply(c(1:length(trait_vec)), function(z){
    gp_pa(get(paste0("within", gsub("data", "", data_vec[x]), "_", trait_vec[z])))
  })) %>%
    mutate(trait = trait_vec, ts = gsub("data_", "", data_vec[x]), vs = gsub("data_", "", data_vec[x])) %>%
    select(ts, vs, trait, overall_pa, stdev)
}))

# Apply prediction accuracy function to across trial results
res_across <- do.call(rbind, lapply(c(1:length(across_vec)), function(x){
  do.call(rbind, lapply(c(1:length(trait_vec)), function(z){
    gp_pa(get(paste0("across_", across_vec[x], "_", trait_vec[z])))
  })) %>%
    mutate(trait = trait_vec, ts = gsub("_[0-9]{4}", "", across_vec[x]), vs = gsub("[0-9]{4}_", "", across_vec[x])) %>%
    select(ts, vs, trait, overall_pa, stdev)
}))

# Combine results
results <- rbind(res_within, res_across)

# Make figures
lapply(c(1:length(trait_vec)), function(x){
  assign(paste0("figure_", trait_vec[x]),
         ggplot(data = results %>% filter(trait == trait_vec[x]), mapping = aes(x = vs, y = reorder(ts, desc(ts)), fill = overall_pa)) + geom_tile(colour = "black") + geom_text(data = results %>% filter(trait == trait_vec[x]), aes(label = paste0(sprintf("%.3f", overall_pa), "\n", "± ", sprintf("%.3f", stdev))), colour = "white") + labs(fill = "Mean PA") + xlab("Validation Set") + ylab("Training Set") + coord_fixed() + theme(panel.background = element_blank(), axis.text = element_text(size = 11, colour = "black"), axis.ticks = element_blank()) + scale_x_discrete(position = "top") + scale_fill_gradient(low = "#4BACC6", high = "#1F497D", limits=c(-0.01,0.75)),
         envir = .GlobalEnv)
})

# Combine into one plot
plot_grid(figure_yield + labs(caption = "Yield (kg/ha)") + theme(plot.caption = element_text(hjust = 0.5, size = 14)),
          figure_height + labs(caption = "Plant Height (cm)") + theme(plot.caption = element_text(hjust = 0.5, size = 14)), 
          figure_sw + labs(caption = "Seed Weight (g)") + theme(plot.caption = element_text(hjust = 0.5, size = 14)), 
          figure_tw + labs(caption = "Test Weight (kg/hL)") + theme(plot.caption = element_text(hjust = 0.5, size = 14)), 
          figure_hd + labs(caption = "Days to Heading") + theme(plot.caption = element_text(hjust = 0.5, size = 14)), 
          labels = c("A", "B", "C", "D", "E"))
ggsave("Figures/gp/within_across_pyt_alltraits.pdf", height = 8.5, width = 14)
ggsave("Figures/gp/within_across_pyt_alltraits.svg", height = 8.5, width = 14)
