## Genomic Prediction Trained with 1, 2, or 3 Cohorts (Size Controlled)
# 01/15/2024

# Objectives: Compare predictions using 1, 2, or 3 cohorts in training set to predict the same validation sets (size controlled).

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

# Function to randomly assign training sets from specified cohorts
# Inputs:
# Outputs:
assign_ts_cohorts <- function(data, ts_cohorts, v_sets, set_seed){
  # number of folds and iterations
  num_iter = length(unique(v_sets$iter))
  num_folds = length(unique(v_sets$fold_id))
  # For each iteration and fold:
  future_map_dfr(c(1:num_iter), .options = furrr_options(seed = set_seed), function(iter){
    # subset the training and validation sets for each iteration
    v_sets_iter <- v_sets[v_sets$iter == iter, ]
    obj <- c(1:num_folds) %>%
      map(function(fold){
        # Set size of training set (size of validation set cohort minus current validation set)
        N <- length(v_sets_iter$fold_id) - length(v_sets_iter$fold_id[v_sets_iter$fold_id == fold])
        # obtain the validation set indices corresponding to the data
        vs_ind <- v_sets_iter$data_id[v_sets_iter$fold_id == fold]
        # get indices of lines that can be sampled for training set (from the cohorts to be included but not in validation set)
        pos_ts <- which(!(rownames(data) %in% vs_ind) & (data$Cohort %in% ts_cohorts))
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

# Function to apply genomic prediction with specified training and validation sets
# Inputs:
# Outputs:
gp_ts_vs <- function(data, markers, t_sets, v_sets, phenotype_var){
  # number of folds and iterations
  num_iter = length(unique(v_sets$iter))
  num_folds = length(unique(v_sets$fold_id))
  # Phenotype data
  y <- noquote(phenotype_var)
  data[[y]] <- data[[phenotype_var]]
  # do genomic prediction for each combination of iterations and folds.
  # first "loop" through the iterations
  c(1:num_iter) %>%
    map_dfr(function(iter){
      print(paste("Iteration", iter))
      #subset the training and validation sets for each iteration
      v_sets_iter <- v_sets[v_sets$iter == iter, ]
      t_sets_iter <- t_sets[t_sets$iter == iter, ]
      #then "loop" through folds
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

## Predictions ----
# Make a vector of traits
trait_vec <- c("yield", "height", "sw", "tw", "hd")

# Assign the validation sets
vs_2018 <- assign_vs(cv_data = gp_data2, vs_cohort = "2018", set_seed = 1L)
vs_2019 <- assign_vs(cv_data = gp_data2, vs_cohort = "2019", set_seed = 2L)
vs_2020 <- assign_vs(cv_data = gp_data2, vs_cohort = "2020", set_seed = 3L)

# Assign the training sets: for each validation set, need to use each other cohort, both other cohorts, all three cohorts. 
ts_2018_pyt_1 <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2019"), v_sets = vs_2018, set_seed = 101L)
ts_2018_pyt_2 <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2020"), v_sets = vs_2018, set_seed = 102L)
ts_2018_2pyt <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2019", "2020"), v_sets = vs_2018, set_seed = 103L)
ts_2018_3pyt <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2018", "2019", "2020"), v_sets = vs_2018, set_seed = 104L)

ts_2019_pyt_1 <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2018"), v_sets = vs_2019, set_seed = 101L)
ts_2019_pyt_2 <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2020"), v_sets = vs_2019, set_seed = 102L)
ts_2019_2pyt <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2018", "2020"), v_sets = vs_2019, set_seed = 103L)
ts_2019_3pyt <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2018", "2019", "2020"), v_sets = vs_2019, set_seed = 104L)

ts_2020_pyt_1 <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2018"), v_sets = vs_2020, set_seed = 101L)
ts_2020_pyt_2 <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2019"), v_sets = vs_2020, set_seed = 102L)
ts_2020_2pyt <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2018", "2019"), v_sets = vs_2020, set_seed = 103L)
ts_2020_3pyt <- assign_ts_cohorts(data = gp_data2, ts_cohorts = c("2018", "2019", "2020"), v_sets = vs_2020, set_seed = 104L)

# Make vectors of training sets to use
ts_2018_vec <- c("ts_2018_pyt_1", "ts_2018_pyt_2", "ts_2018_2pyt", "ts_2018_3pyt")
ts_2019_vec <- gsub("2018", "2019", ts_2018_vec)
ts_2020_vec <- gsub("2018", "2020", ts_2018_vec)

# Apply genomic prediction function to each validation set with each training set
lapply(c(1:length(ts_2018_vec)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0(gsub("ts", "gp", ts_2018_vec[x]), "_", trait_vec[z]),
           gp_ts_vs(data = gp_data2,
                    markers = gp_markers,
                    t_sets = get(ts_2018_vec[x]),
                    v_sets = vs_2018,
                    phenotype_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
  })
})

lapply(c(1:length(ts_2019_vec)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0(gsub("ts", "gp", ts_2019_vec[x]), "_", trait_vec[z]),
           gp_ts_vs(data = gp_data2,
                    markers = gp_markers,
                    t_sets = get(ts_2019_vec[x]),
                    v_sets = vs_2019,
                    phenotype_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
  })
})

lapply(c(1:length(ts_2020_vec)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0(gsub("ts", "gp", ts_2020_vec[x]), "_", trait_vec[z]),
           gp_ts_vs(data = gp_data2,
                    markers = gp_markers,
                    t_sets = get(ts_2020_vec[x]),
                    v_sets = vs_2020,
                    phenotype_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
  })
})

## Obtain correlations and make figure ----
# Function to calculate correlations
gp_cor <- function(gp_result){
  gp_result %>%
    group_by(iter, fold_id) %>%
    summarise(fold_pa = cor(predicted_val, obs_val)) %>%
    group_by(iter) %>%
    summarise(iter_pa = mean(fold_pa)) %>%
    summarise(overall_pa = mean(iter_pa), stdev = sqrt(var(iter_pa)))
}

# Apply function to results
res_2018 <- do.call(rbind, lapply(c(1:length(ts_2018_vec)), function(x){
  do.call(rbind, lapply(c(1:length(trait_vec)), function(z){
    gp_cor(get(paste0(gsub("ts", "gp", ts_2018_vec[x]), "_", trait_vec[z]))) %>%
      mutate(Cohort = "2018", trait = trait_vec[z], ts = gsub("ts_2018_", "", ts_2018_vec[x]))
  }))
}))

res_2019 <- do.call(rbind, lapply(c(1:length(ts_2019_vec)), function(x){
  do.call(rbind, lapply(c(1:length(trait_vec)), function(z){
    gp_cor(get(paste0(gsub("ts", "gp", ts_2019_vec[x]), "_", trait_vec[z]))) %>%
      mutate(Cohort = "2019", trait = trait_vec[z], ts = gsub("ts_2019_", "", ts_2019_vec[x]))
  }))
}))

res_2020 <- do.call(rbind, lapply(c(1:length(ts_2020_vec)), function(x){
  do.call(rbind, lapply(c(1:length(trait_vec)), function(z){
    gp_cor(get(paste0(gsub("ts", "gp", ts_2020_vec[x]), "_", trait_vec[z]))) %>%
      mutate(Cohort = "2020", trait = trait_vec[z], ts = gsub("ts_2020_", "", ts_2020_vec[x]))
  }))
}))

# Combine
results <- rbind(res_2018, res_2019, res_2020) %>%
  mutate(ts = factor(case_when(.$ts == "2pyt" ~ "Two Cohorts", .$ts == "3pyt" ~ "Three Cohorts", .$ts == "pyt_1" ~ "One Cohort", .$ts == "pyt_2" ~ "One Cohort (Other)"), levels = c("One Cohort", "One Cohort (Other)", "Two Cohorts", "Three Cohorts")))

# Make figures
lapply(c(1:length(trait_vec)), function(x){
  assign(paste0("figure_", trait_vec[x]),
         ggplot(data = results %>% filter(trait == trait_vec[x]), mapping = aes(x = Cohort, y = overall_pa, fill = ts)) + geom_bar(stat = "identity", position = "dodge") + xlab("Validation Set") + ylab("Mean Correlation") + labs(fill = "Cohorts in Training set") + scale_fill_manual(values = c("#4BACC6", "#1F497D", "#8064A2", "#9BBB59")) + geom_errorbar(aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + theme_light(),
         envir = .GlobalEnv)
})

figure_legend <- get_legend(figure_yield)

# Combine into one plot
plot_grid(figure_yield + ylim(-0.05,0.75) + labs(caption = "Yield (kg/ha)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          figure_height + ylim(-0.05,0.75) + labs(caption = "Plant Height (cm)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          figure_sw + ylim(-0.05,0.75) + labs(caption = "Seed Weight (g)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          figure_tw + ylim(-0.05,0.75) + labs(caption = "Test Weight (kg/hL)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          figure_hd + ylim(-0.05,0.75) + labs(caption = "Days to Heading") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          figure_legend,
          labels = c("A", "B", "C", "D", "E"))
ggsave("Figures/gp/gp_123cohort_size_control.pdf", height = 8.5, width = 11)
ggsave("Figures/gp/gp_123cohort_size_control.svg", height = 8.5, width = 11)
