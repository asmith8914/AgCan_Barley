## Obtain AYT entry means
# 12/15/2023

# Objective: fit models to the AYT traits and obtain the genotype means.

library(tidyverse)
library(emmeans)
library(lme4)

# Read in data, change column names, set numeric categories as factors
adv_pheno_fixed <- read.csv("Data/pheno/Phenodata_adv_hd_cor.csv")
pheno_adv <- adv_pheno_fixed %>%
  rename(tw = testweight, sw = SEEDWEIGHT_G, height = HT_CM, hd = days_to_heading) %>%
  mutate(YEAR = as.factor(YEAR), REP = as.factor(REP))

# Make vector of traits to apply functions with
trait_vec <- c("yield", "height", "sw", "tw", "hd")

# Pivot the data so there is a trait column and a value column
data_adv <- pheno_adv %>%
  pivot_longer(!colnames(pheno_adv)[1:6], names_to = "trait", values_to = "value")

# Model to fit is: y = mean + geno + year + location:year + block:location:year + geno:location:year + e
model_adv <- formula(value ~ NAME + YEAR + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))

# Apply models for all traits
lapply(c(1:5), function(x){
  assign(paste0("lm_", trait_vec[x]),
         lmer(formula = model_adv, data = data_adv[data_adv$trait == trait_vec[x],]),
         envir = .GlobalEnv)
})

# Apply anova and summary to models, show variance of random effects
lapply(c(1:5), function(x){
  anova(get(paste0("lm_", trait_vec[x])))
})

lapply(c(1:5), function(x){
  summary(get(paste0("lm_", trait_vec[x])))
})

lapply(c(1:5), function(x){
  data.frame(VarCorr(get(paste0("lm_", trait_vec[x]))))
})

# Get emmeans for entries
lapply(c(1:5), function(x){
  assign(paste0("emmeans_", trait_vec[x]),
         as.data.frame(emmeans(get(paste0("lm_", trait_vec[x])), "NAME")),
         envir = .GlobalEnv)
})

# Save the emmeans
adv_emmeans <- emmeans_yield %>%
  select(NAME, emmean) %>%
  rename(mean_AYT_yield = emmean) %>%
  mutate(mean_AYT_height = emmeans_height$emmean,
         mean_AYT_sw = emmeans_sw$emmean,
         mean_AYT_tw = emmeans_tw$emmean,
         mean_AYT_hd = emmeans_hd$emmean)

write.csv(adv_emmeans, "Data/pheno/adv_emmeans.csv", row.names = FALSE)
