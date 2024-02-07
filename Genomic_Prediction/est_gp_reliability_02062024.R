## Genomic prediction reliabilities
# 02/06/2024

# Objective: Estimate reliabilities of genomic predictions

library(tidyverse)
library(rrBLUP)

# Read in data
gp_data <- read.csv("Data/gp/GS_PYT_PYTu_AYT_alltraits_5899.csv")

# Markers
gp_markers <- as.matrix(gp_data[, -c(1:23)])

# Function to apply mixed.solve to calculate marker effects and calculate GEBVs.
ms_GEBVs <- function(data, phenotype_var, markers){
  #set parameters
  y <- noquote(phenotype_var)
  data[[y]] <- data[[phenotype_var]]
  #fit marker effects
  marker_eff <- mixed.solve(data[[y]], Z = markers, SE = FALSE)
  #predicted values obtained by multiplying genotype matrix with marker effects and adding mean.
  pred <- markers %*% marker_eff$u + c(marker_eff$beta)
  #make dataframe to return predicted and observed values
  res <- data.frame(NAME = data$NAME,
                    Cohort = data$Cohort, 
                    predicted_val = pred,
                    obs_val = data[[y]]) %>%
    rownames_to_column("obs_id")
  return(res)}

# Trait vector
trait_vec <- c("yield", "height", "sw", "tw", "hd")

# For each trait, apply mixed.solve to the whole dataset
lapply(c(1:length(trait_vec)), function(z){
  assign(paste0("gp_", trait_vec[z]),
         ms_GEBVs(data = gp_data,
                  phenotype_var = paste0("mean_PYT_", trait_vec[z]),
                  markers = gp_markers),
         envir = .GlobalEnv)
})

# Obtain correlations of predicted and observed values
gp_cors <- sapply(c(1:length(trait_vec)), function(z){
  cor(get(paste0("gp_", trait_vec[z]))["predicted_val"],
      get(paste0("gp_", trait_vec[z]))["obs_val"])
})
names(gp_cors) <- trait_vec

# Read in the heritabilities
brd_h2 <- read.csv("Data/pheno/broad_sense_heritabilities.csv")

# Reliabilities = cor(GEBV,Y)^2/h2
gp_rel <- (gp_cors^2)/(brd_h2[brd_h2$h2 == "pyt", -1])
# Note: reliabilities are quite high, height is above 1.

write.csv(gp_rel, "Data/gp/gp_reliabilities.csv", row.names = FALSE)