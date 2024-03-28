## Genomic prediction reliabilities
# 03/12/2024

# Objective: Estimate reliabilities of genomic predictions

library(tidyverse)
library(furrr)
library(rrBLUP)

# Read in data
gp_data <- read.csv("Data/gp/GS_PYT_PYTu_AYT_alltraits_5899.csv")

# Markers
gp_markers <- as.matrix(gp_data[, -c(1:23)])

# Genomic prediction functions:
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

gp_pa <- function(gp_result){
  gp_result %>%
    group_by(iter, fold_id) %>%
    summarise(fold_pa = cor(predicted_val, obs_val)) %>%
    group_by(iter) %>%
    summarise(iter_pa = mean(fold_pa)) %>%
    summarise(overall_pa = mean(iter_pa), stdev = sqrt(var(iter_pa)))
}

# Trait vector
trait_vec <- c("yield", "height", "sw", "tw", "hd")

# For each trait, apply cross validation and genomic prediction to the whole dataset
sets <- assign_sets_within(cv_data = gp_data, set_seed = 2024L)

lapply(c(1:length(trait_vec)), function(z){
  assign(paste0("gp_", trait_vec[z]),
         gp_within(cv_data = gp_data,
                            v_sets = sets,
                            markers = gp_markers,
                            phenotype_var = paste0("mean_PYT_", trait_vec[z])),
         envir = .GlobalEnv)
})

# Obtain correlations of predicted and observed values
gp_cors <- sapply(c(1:length(trait_vec)), function(z){
  gp_pa(get(paste0("gp_", trait_vec[z])))
})

gp_cor <- unlist(gp_cors[1,])
names(gp_cor) <- trait_vec

# Read in the heritabilities
brd_h2 <- read.csv("Data/pheno/broad_sense_heritabilities.csv")

# Reliabilities = cor(GEBV,Y)^2/h2
gp_rel <- (gp_cor^2)/(brd_h2[brd_h2$h2 == "pyt", -1])

write.csv(gp_rel, "Data/gp/gp_reliabilities2.csv", row.names = FALSE)
