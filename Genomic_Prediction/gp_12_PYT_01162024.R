## Genomic Prediction Trained with 1, 2, or 3 Cohorts (Not Size Controlled)
# 01/16/2024

# Objectives: Compare predictions using 1, 2, or 3 cohorts in training set to predict the same validation sets (not size controlled).

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

gp_set_ts_vs <- function(data, markers, ts_cohorts, v_sets, phenotype_var){
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

## Predictions ----
# Make a vector of traits
trait_vec <- c("yield", "height", "sw", "tw", "hd")

# Assign the validation sets
vs_2018 <- assign_vs(cv_data = gp_data2, vs_cohort = "2018", set_seed = 1L)
vs_2019 <- assign_vs(cv_data = gp_data2, vs_cohort = "2019", set_seed = 2L)
vs_2020 <- assign_vs(cv_data = gp_data2, vs_cohort = "2020", set_seed = 3L)

# Apply the genomic prediction function to each validation set, using each other cohort and both other cohorts for training.
ts_2018_list <- list("2019", "2020", c("2019", "2020"))

lapply(c(1:length(ts_2018_list)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0("gp_2018_", x, "_", trait_vec[z]),
           gp_set_ts_vs(data = gp_data2,
                        markers = gp_markers,
                        ts_cohorts = ts_2018_list[[x]],
                        v_sets = vs_2018,
                        phenotype_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
  })
})

ts_2019_list <- list("2018", "2020", c("2018", "2020"))

lapply(c(1:length(ts_2019_list)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0("gp_2019_", x, "_", trait_vec[z]),
           gp_set_ts_vs(data = gp_data2,
                        markers = gp_markers,
                        ts_cohorts = ts_2019_list[[x]],
                        v_sets = vs_2019,
                        phenotype_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
  })
})

ts_2020_list <- list("2018", "2019", c("2018", "2019"))

lapply(c(1:length(ts_2020_list)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0("gp_2020_", x, "_", trait_vec[z]),
           gp_set_ts_vs(data = gp_data2,
                        markers = gp_markers,
                        ts_cohorts = ts_2020_list[[x]],
                        v_sets = vs_2020,
                        phenotype_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
  })
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

# Make vector with result names
res_vec <- paste0("gp_20", rep(18:20, each = 15), "_", rep(1:3, each = 5), "_", c("yield", "height", "sw", "tw", "hd"))


# Apply function to results
results <- do.call(rbind, lapply(c(1:length(res_vec)), function(x){
  gp_cor(get(res_vec[x])) %>%
    mutate(source = res_vec[x])
}))

# Extract info from source
results2 <- results %>%
  mutate(trait = gsub("^gp_[0-9]{4}_[1-3]{1}_", "", .$source),
         Cohort = gsub("_[1-3]{1}_[a-z]*$", "", gsub("^gp_", "", .$source)),
         ts_type = gsub("_[a-z]*$", "", gsub("^gp_[0-9]{4}_", "", .$source)))

results3 <- results2 %>%
  mutate(ts = factor(case_when(ts_type == "1" ~ "One Cohort", ts_type == "2" ~ "One Cohort (Other)", ts_type == "3" ~ "Two Cohorts"), levels = c("One Cohort", "One Cohort (Other)", "Two Cohorts"))) %>%
  select(Cohort, trait, ts, overall_pa, stdev)

# Make figures
lapply(c(1:length(trait_vec)), function(x){
  assign(paste0("figure_", trait_vec[x]),
         ggplot(data = results3 %>% filter(trait == trait_vec[x]), mapping = aes(x = Cohort, y = overall_pa, fill = ts)) + geom_bar(stat = "identity", position = "dodge") + xlab("Validation Set") + ylab("Mean Correlation") + labs(fill = "Cohorts in Training set") + scale_fill_manual(values = c("#4BACC6", "#1F497D", "#8064A2")) + geom_errorbar(aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + theme_light(),
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
ggsave("Figures/gp/gp_12cohort_not_sz_cont.pdf", height = 8.5, width = 11)
ggsave("Figures/gp/gp_12cohort_not_sz_cont.svg", height = 8.5, width = 11)
