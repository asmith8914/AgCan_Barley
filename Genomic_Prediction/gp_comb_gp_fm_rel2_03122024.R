## Combined Family Means and Genomic Predictions
# 03/12/2024

# Objective: Compare genomic predictions, family means, and combined family means and genomic predictions using reliabilities as weights.

library(tidyverse)
library(rrBLUP)
library(furrr)
library(cowplot)
library(ggpubr)
library(rstatix)

# Check RNGkind
RNGkind()
# "Mersenne-Twister" "Inversion"        "Rejection"   
RNGkind("L'Ecuyer-CMRG")
# Set as "L'Ecuyer-CMRG" because genomic prediction functions will use functions that work in parallel, so this is required for reproducibility.

# Read in data
gp_data <- read.csv("Data/gp/GS_PYT_PYTu_AYT_alltraits_5899.csv")

# Remove lines from crosses with only 1 line
cross_include <- table(gp_data$cross)
cross_include <- cross_include[cross_include > 1]
gp_cross <- gp_data %>%
  filter(cross %in% names(cross_include))

# Separate the data into cohorts
data_2018 <- gp_cross %>%
  filter(Cohort == "2018")
data_2019 <- gp_cross %>%
  filter(Cohort == "2019")
data_2020 <- gp_cross %>%
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

## Functions ----
# Function to assign training and validation sets
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

# Genomic prediction function to match family means training and validation sets. Note that this function keeps training set lines whose cross isn't represented in the validation set (whereas information from those lines is not used for family means - training set and validation sets are still the same though, gp just uses more information from training set than family means).
gp_fam_inc <- function(cv_data, v_sets, markers, phenotype_var){
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
      # then do cross-validation
      cv_res <- c(1:num_folds) %>% 
        map(function(fold){
          print(paste("Cross-Validating Fold", fold))
          #make vector of indices for training and validation sets
          vs_ind <- which(v_sets_iter$fold_id == fold)
          ts_ind <- which(v_sets_iter$fold_id != fold)
          #get entry names for training and validation sets
          vs_entry <- cv_data[vs_ind, 1]
          ts_entry <- cv_data[ts_ind, 1]
          #get training set crosses
          ts_cross <- gsub("-.*", "", ts_entry)
          #make a copy of the data, removing lines whose crosses aren't in the training set
          data1 <- cv_data[cv_data$cross %in% ts_cross,]
          #update training and validation set indices: need indices that refer to both copied and original datasets
          ts_ind_copy <- which(data1$NAME %in% ts_entry)
          vs_ind_copy <- which(data1$NAME %in% vs_entry)
          ts_ind_orig <- which(cv_data$NAME %in% data1$NAME[data1$NAME %in% ts_entry])
          vs_ind_orig <- which(cv_data$NAME %in% data1$NAME[data1$NAME %in% vs_entry])
          #mask validation set values
          data1[[y]][vs_ind_copy] <- NA
          #fit marker effects and make predictions for the fold
          marker_eff <- mixed.solve(data1[[y]][ts_ind_copy], Z = markers[ts_ind_orig,], SE = FALSE)
          #predicted values obtained by multiplying genotype matrix with marker effects and adding mean.
          pred <- markers[vs_ind_orig,] %*% marker_eff$u + c(marker_eff$beta)
          #make dataframe to return predicted and observed values, as well as fold_id, to use for grouping when taking mean.
          res <- data.frame(NAME = cv_data$NAME[vs_ind_orig],
                            gp_pred_val = pred,
                            obs_val = cv_data[[y]][vs_ind_orig],
                            fold_id = fold) %>%
            rownames_to_column("obs_id")
          return(res)}) %>%
        bind_rows()
      cv_res$iter <- iter
      return(as_tibble(cv_res))})
}

# Function to calculate family means. Note that validation set entries whose cross isn't represented in training set are removed.
family_means <- function(data, v_sets, pheno_var){
  # parameters
  num_iter = length(unique(v_sets$iter))
  num_folds = length(unique(v_sets$fold_id))
  y <- noquote(pheno_var)
  data[[y]] <- data[[pheno_var]]
  # need to group by iteration
  obj <- c(1:num_iter) %>%
    map_dfr(function(iter){
      print(paste("Iteration", iter))
      v_sets_iter <- v_sets[v_sets$iter == iter, ]
      # now get means in folds in iter
      cv_res <- c(1:num_folds) %>% 
        map(function(fold){
          #make vector of indices for training and validation sets
          vs_ind <- which(v_sets_iter$fold_id == fold)
          ts_ind <- which(v_sets_iter$fold_id != fold)
          #get entry names for training and validation sets
          vs_entry <- data[vs_ind, 1]
          ts_entry <- data[ts_ind, 1]
          #subset data to training set entries
          ts_data <- data[data$NAME %in% ts_entry,]
          #calculate cross means for filtered training set
          ts_cross_means <- ts_data %>%
            group_by(cross) %>%
            summarise(cross_mean = mean(.data[[y]])) %>%
            mutate(fold_id = fold)
          #make dataframe to return entry names, cross means, observed values, and fold and iter ids.
          vs_data <- data.frame(NAME = data$NAME[vs_ind],
                                cross = data$cross[vs_ind],
                                obs_val = data[[y]][vs_ind],
                                fold_id = fold) 
          vs_data2 <- left_join(vs_data, ts_cross_means) %>%
            filter(!is.na(cross_mean)) %>%
            rownames_to_column("obs_id") %>%
            #return results
            return(vs_data2)}) %>%
        bind_rows()
      cv_res$iter <- iter
      return(as_tibble(cv_res))
    })
}

# Function to calculate prediction accuracy for specified prediction type (column with prediction)
iter_pa_var <- function(input, pred_col){
  y <- noquote(pred_col)
  input$y <- input[[pred_col]]
  input %>%
    group_by(iter, fold_id) %>%
    summarise(fold_pa = cor(y, obs_val)) %>%
    group_by(iter) %>%
    summarise(iter_pa = mean(fold_pa))
} 

## Predictions ----
# Assign the sets 
within_2018_sets <- assign_sets_within(cv_data = data_2018, set_seed = 1L)
within_2019_sets <- assign_sets_within(cv_data = data_2019, set_seed = 2L)
within_2020_sets <- assign_sets_within(cv_data = data_2020, set_seed = 3L)

# Make trait and data vectors
trait_vec <- c("yield", "height", "sw", "tw", "hd")
data_vec <- c("data_2018", "data_2019", "data_2020")

# Apply the genomic prediction function for each phenotype variable and each cohort
lapply(c(1:length(data_vec)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0("gp", gsub("data", "", data_vec[x]), "_", trait_vec[z]),
           gp_fam_inc(cv_data = get(data_vec[x]),
                      markers = get(gsub("data", "markers", data_vec[x])),
                      v_sets = get(paste0("within", gsub("data", "", data_vec[x]), "_sets")),
                      phenotype_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
  })
})

# Apply the family means prediction function for each phenotype variable and each cohort
lapply(c(1:length(data_vec)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0("fm", gsub("data", "", data_vec[x]), "_", trait_vec[z]),
           family_means(data = get(data_vec[x]),
                        v_sets = get(paste0("within", gsub("data", "", data_vec[x]), "_sets")),
                        pheno_var = paste0("mean_PYT_", trait_vec[z])),
           envir = .GlobalEnv)
  })
})

# Weights (reliabilities) of combined predictions are: h2 of family means, and cor(GEBV,y)^2/h2 for genomic prediction
# Read in the heritabilities
brd_h2 <- read.csv("Data/pheno/broad_sense_heritabilities.csv")

# Read in the reliabilities of genomic predictions
gp_rel <- read.csv("Data/gp/gp_reliabilities2.csv")

# Extract family means heritability (remove rownames to allow column names to be set later)
fam_rel <- brd_h2[brd_h2$h2 == "pyt_fam", -1]
row.names(fam_rel) <- NULL

# Get sum of reliabilties
w1_w2 <- gp_rel + fam_rel

# Weights:
weights <- data.frame(w_gp = t(gp_rel/w1_w2),
                      w_fm = t(fam_rel/w1_w2)) %>%
  rownames_to_column(var = "trait")

# Make a data frame with both types of results, and add a column with combined prediction
# Combined prediction = (w1 * GEBV) + (w2 * family mean), where w1 = w_gp (rel_gp / (rel_gp + fam_rel))and w2 = w_fam (fam_rel / (rel_gp + fam_rel)).
lapply(c(1:length(data_vec)), function(x){
  lapply(c(1:length(trait_vec)), function(z){
    assign(paste0("comb", gsub("data", "", data_vec[x]), "_", trait_vec[z]),
           full_join(get(paste0("gp", gsub("data", "", data_vec[x]), "_", trait_vec[z])),
                     get(paste0("fm", gsub("data", "", data_vec[x]), "_", trait_vec[z]))) %>%
             mutate(comb_pred = (weights$w_gp[weights$trait == trait_vec[z]]*gp_pred_val)+(weights$w_fm[weights$trait == trait_vec[z]]*cross_mean)),
           envir = .GlobalEnv)
  })
})

# Calculate correlations:
# First make list of results dataframes and vector of types of predicted values
res_list <- mget(paste0(rep(paste0("comb", gsub("data", "", data_vec)), each = 5), "_", trait_vec))

# Apply to each result and type of predicted values (adding method column to make joining easier)
gp_iter_res <- lapply(res_list, function(x){
  iter_pa_var(input = x, pred_col = "gp_pred_val") %>%
    mutate(method = as.factor("gp"))})

fm_iter_res <- lapply(res_list, function(x){
  iter_pa_var(input = x, pred_col = "cross_mean")  %>%
    mutate(method = as.factor("fm"))})

comb_iter_res <- lapply(res_list, function(x){
  iter_pa_var(input = x, pred_col = "comb_pred") %>%
    mutate(method = as.factor("comb"))})

# Calculate the overall correlations (method was dropped, added back in)
gp_over_res <- sapply(c(1:15), function(x){
  gp_iter_res[[x]] %>%
    summarise(overall_pa = mean(iter_pa), stdev = sqrt(var(iter_pa))) %>%
    mutate(source = names(gp_iter_res)[x], 
           method = "gp")
})

fm_over_res <- sapply(c(1:15), function(x){
  fm_iter_res[[x]] %>%
    summarise(overall_pa = mean(iter_pa), stdev = sqrt(var(iter_pa))) %>%
    mutate(source = names(fm_iter_res)[x],
           method = "fm")
})

comb_over_res <- sapply(c(1:15), function(x){
  comb_iter_res[[x]] %>%
    summarise(overall_pa = mean(iter_pa), stdev = sqrt(var(iter_pa))) %>%
    mutate(source = names(comb_iter_res)[x],
           method = "comb")
})

# Transpose results and bind together (fix classes from transposed list)
overall_res <- data.frame(rbind(t(gp_over_res), t(fm_over_res), t(comb_over_res))) %>%
  mutate(Cohort = gsub("_[a-z]*$", "", gsub("^[a-z]{2,}_", "", .$source)),
         trait = gsub("^[a-z]{2,}_[0-9]{4}_", "", .$source),
         method = case_when(.$method == "gp" ~ "Genomic Prediction", .$method == "fm" ~ "Family Means", .$method == "comb" ~ "Combined Prediction"),
         stdev = as.numeric(.$stdev),
         overall_pa = as.numeric(.$overall_pa)) %>%
  select(method, Cohort, trait, overall_pa, stdev)

write.csv(overall_res, "Data/gp/gp_pa_fm_comb_rel2.csv", row.names = FALSE)

## Add statistical tests ----
# Since there are multiple groups but data is still paired, need to do pairwise Bonferroni corrected paired t-tests.
# Make list of dataframes (cohort + trait) with iteration PAs and types of predictions
all_iter_res <- lapply(c(1:15), function(x){
  rbind(gp_iter_res[[x]], fm_iter_res[[x]], comb_iter_res[[x]])
})

# Add names back to list
all_iter_res2 <- all_iter_res
names(all_iter_res2) <- gsub("comb_", "", names(comb_iter_res))  

# Turn list into dataframe (specify levels so order is correct on plot)
all_iter_df <- c(1:15) %>% map(function(x){
  as.data.frame(all_iter_res2[[x]]) %>%
    mutate(Cohort = gsub("_[a-z]*", "", names(all_iter_res2)[x]),
           trait = gsub("[0-9]{4}_", "", names(all_iter_res2)[x]))
}) %>%
  bind_rows() %>%
  mutate(method = factor(.$method, levels = c("gp", "fm", "comb")))

# Apply pairwise paired t-tests with Bonferroni correction
pair_ttest <- all_iter_df %>%
  group_by(Cohort, trait) %>%
  pairwise_t_test(iter_pa ~ method, p.adjust.method = "bonferroni", paired = TRUE) %>%
  add_xy_position(x = "Cohort", dodge = 0.9, step.increase = 0.1)

write.csv(pair_ttest[,1:12], "Data/gp/fam_gp_comb_ttest_res2.csv", row.names = FALSE)

# Remake the overall result table (need method to match groups in t-test results, will make labels more informative on figure later)
overall_res2 <- data.frame(rbind(t(gp_over_res), t(fm_over_res), t(comb_over_res))) %>%
  mutate(Cohort = gsub("_[a-z]*$", "", gsub("^[a-z]{2,}_", "", .$source)),
         trait = gsub("^[a-z]{2,}_[0-9]{4}_", "", .$source),
         method = factor(as.character(.$method), levels = c("gp", "fm", "comb")),
         stdev = as.numeric(.$stdev),
         overall_pa = as.numeric(.$overall_pa)) %>%
  select(method, Cohort, trait, overall_pa, stdev)

# Figures:
fig_ttest_yield <- ggplot(data = overall_res2 %>% filter(trait == trait_vec[1]), mapping = aes(x = Cohort, y = overall_pa, fill = method)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylab("Prediction Accurary (r)") + 
  labs(fill = "Prediction Method") + 
  scale_fill_manual(values = c("#4BACC6", "#C0504D", "#8064A2"), labels = c("Genomic Prediction", "Family Means", "Combined Prediction")) + 
  geom_errorbar(data = overall_res2 %>% filter(trait == trait_vec[1]), mapping = aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + 
  stat_pvalue_manual(pair_ttest %>% filter(trait == trait_vec[1]), label = "p.adj.signif", inherit.aes = FALSE, hide.ns = TRUE) + 
  ylim(c(0,1)) + 
  theme_light()
fig_ttest_yield

fig_ttest_height <- ggplot(data = overall_res2 %>% filter(trait == trait_vec[2]), mapping = aes(x = Cohort, y = overall_pa, fill = method)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylab("Prediction Accurary (r)") + 
  labs(fill = "Prediction Method") + 
  scale_fill_manual(values = c("#4BACC6", "#C0504D", "#8064A2"), labels = c("Genomic Prediction", "Family Means", "Combined Prediction")) + 
  geom_errorbar(data = overall_res2 %>% filter(trait == trait_vec[2]), mapping = aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + 
  stat_pvalue_manual(pair_ttest %>% filter(trait == trait_vec[2]), label = "p.adj.signif", inherit.aes = FALSE, hide.ns = TRUE) + 
  ylim(c(0,1)) + 
  theme_light()
fig_ttest_height

fig_ttest_sw <- ggplot(data = overall_res2 %>% filter(trait == trait_vec[3]), mapping = aes(x = Cohort, y = overall_pa, fill = method)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylab("Prediction Accurary (r)") + 
  labs(fill = "Prediction Method") + 
  scale_fill_manual(values = c("#4BACC6", "#C0504D", "#8064A2"), labels = c("Genomic Prediction", "Family Means", "Combined Prediction")) + 
  geom_errorbar(data = overall_res2 %>% filter(trait == trait_vec[3]), mapping = aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + 
  stat_pvalue_manual(pair_ttest %>% filter(trait == trait_vec[3]), label = "p.adj.signif", inherit.aes = FALSE, hide.ns = TRUE) + 
  ylim(c(0,1)) + 
  theme_light()
fig_ttest_sw

fig_ttest_tw <- ggplot(data = overall_res2 %>% filter(trait == trait_vec[4]), mapping = aes(x = Cohort, y = overall_pa, fill = method)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylab("Prediction Accurary (r)") + 
  labs(fill = "Prediction Method") + 
  scale_fill_manual(values = c("#4BACC6", "#C0504D", "#8064A2"), labels = c("Genomic Prediction", "Family Means", "Combined Prediction")) + 
  geom_errorbar(data = overall_res2 %>% filter(trait == trait_vec[4]), mapping = aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + 
  stat_pvalue_manual(pair_ttest %>% filter(trait == trait_vec[4]), label = "p.adj.signif", inherit.aes = FALSE, hide.ns = TRUE) + 
  ylim(c(0,1)) + 
  theme_light()
fig_ttest_tw

fig_ttest_hd <- ggplot(data = overall_res2 %>% filter(trait == trait_vec[5]), mapping = aes(x = Cohort, y = overall_pa, fill = method)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylab("Prediction Accurary (r)") + 
  labs(fill = "Prediction Method") + 
  scale_fill_manual(values = c("#4BACC6", "#C0504D", "#8064A2"), labels = c("Genomic Prediction", "Family Means", "Combined Prediction")) + 
  geom_errorbar(data = overall_res2 %>% filter(trait == trait_vec[5]), mapping = aes(ymin = overall_pa-stdev, ymax = overall_pa+stdev), width = 0.2, position = position_dodge(0.9)) + 
  stat_pvalue_manual(pair_ttest %>% filter(trait == trait_vec[5]), label = "p.adj.signif", inherit.aes = FALSE, hide.ns = TRUE) + 
  ylim(c(0,1)) + 
  theme_light()
fig_ttest_hd

fig_ttest_legend <- cowplot::get_legend(fig_ttest_hd)

plot_grid(fig_ttest_yield + labs(caption = "Yield (kg/ha)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)),
          fig_ttest_height + labs(caption = "Plant Height (cm)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)), 
          fig_ttest_sw + labs(caption = "Seed Weight (g)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)), 
          fig_ttest_tw + labs(caption = "Test Weight (kg/hL)") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)), 
          fig_ttest_hd + labs(caption = "Days to Heading") + theme(legend.position = "none", plot.caption = element_text(hjust = 0.5, size = 14)), 
          fig_ttest_legend,
          labels = c("A", "B", "C", "D", "E"))
ggsave("Figures/gp/fam_gp_comb_rel2_ttest.svg", height = 8.5, width = 14)
ggsave("Figures/gp/fam_gp_comb_rel2_ttest.pdf", height = 8.5, width = 14)
