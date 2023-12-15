## Obtain PYT entry means
# 12/15/2023

# Objective: fit models to the PYT traits and obtain the genotype means.

library(tidyverse)
library(emmeans)
library(lme4)

# Read in data
pheno_pyt <- read.csv("Data/pheno/Phenodata_prelim_hd.csv")

# Clean up names and set numerical categories as factors
pheno_pyt <- pheno_pyt %>%
  rename(tw = testweight, sw = SEEDWEIGHT_G, height = HT_CM, hd = days_to_heading) %>%
  mutate(YEAR = as.factor(YEAR), BLOCK = as.factor(BLOCK), REP = as.factor(REP))

# Trait vectors
swtw <- c("sw", "tw")
yphhd <- c("yield", "height", "hd")

# Pivot the data so there is a trait column and a value column
data_pyt <- pheno_pyt %>%
  pivot_longer(!colnames(pheno_pyt)[1:7], names_to = "trait", values_to = "value")

# Model for yield, height, and hd is y = mean + geno + year + block:year + error
model_yphhd <- formula(value ~ NAME + YEAR + (1|BLOCK:YEAR))

# Model for seed weight and test weight is y = mean + geno + year + error
model_swtw <- formula(value ~ NAME + YEAR)
# Note this is also model for treating replicated traits as unreplicated.

# Apply models
lapply(c(1:3), function(x){
  assign(paste0("lm_", yphhd[x]),
         lmer(formula = model_yphhd, data = data_pyt[data_pyt$trait == yphhd[x],]),
         envir = .GlobalEnv)
})

lapply(c(1:2), function(x){
  assign(paste0("lm_", swtw[x]),
         lm(formula = model_swtw, data = data_pyt[data_pyt$trait == swtw[x],]),
         envir = .GlobalEnv)
})

# Apply anova and summary to models, show variance of random effects for yield, ph, hd.
lapply(c(1:2), function(x){
  anova(get(paste0("lm_", swtw[x])))
})

lapply(c(1:2), function(x){
  summary(get(paste0("lm_", swtw[x])))
})

lapply(c(1:3), function(x){
  anova(get(paste0("lm_", yphhd[x])))
})

lapply(c(1:3), function(x){
  summary(get(paste0("lm_", yphhd[x])))
})

lapply(c(1:3), function(x){
  data.frame(VarCorr(get(paste0("lm_", yphhd[x]))))
})

# Get emmeans for entries
lapply(c(1:3), function(x){
  assign(paste0("emmeans_", yphhd[x]),
         as.data.frame(emmeans(get(paste0("lm_", yphhd[x])), "NAME")),
         envir = .GlobalEnv)
})

lapply(c(1:2), function(x){
  assign(paste0("emmeans_", swtw[x]),
         as.data.frame(emmeans(get(paste0("lm_", swtw[x])), "NAME")),
         envir = .GlobalEnv)
})

# Save the emmeans
pyt_emmeans <- emmeans_yield %>%
  select(NAME, emmean) %>%
  rename(mean_PYT_yield = emmean) %>%
  mutate(mean_PYT_height = emmeans_height$emmean,
         mean_PYT_sw = emmeans_sw$emmean,
         mean_PYT_tw = emmeans_tw$emmean,
         mean_PYT_hd = emmeans_hd$emmean)

write.csv(pyt_emmeans, "Data/pheno/pyt_emmeans.csv", row.names = FALSE)

# Treat yield, plant height and days to heading as unreplicated ----
# Apply the unrep model (lapply for each rep within lapply for each trait)
lapply(c(1:3), function(x){
  assign(paste0("lm_unrep_", yphhd[x]),
         lapply(c(1:3), function(z){
           lm(formula = model_swtw, data = data_pyt %>% filter(trait == yphhd[x], REP == z))
           }),
         envir = .GlobalEnv)
})

# Apply anova and summary
lapply(lm_unrep_yield, anova)
lapply(lm_unrep_height, anova)
lapply(lm_unrep_hd, anova)

lapply(lm_unrep_yield, summary)
lapply(lm_unrep_height, summary)
lapply(lm_unrep_hd, summary)

# Obtain the unrep emmeans
lapply(c(1:3), function(x){
  assign(paste0("emmeans_unrep_", yphhd[x]),
         lapply(c(1:3), function(z){
           as.data.frame(emmeans(get(paste0("lm_unrep_", yphhd[x]))[[z]], "NAME"))
         }),
         envir = .GlobalEnv)
})

# Make dataframe
pyt_unrep_emmeans <- emmeans_unrep_yield[[1]] %>%
  select(NAME, emmean) %>%
  rename(unrep1_PYT_yield = emmean) %>%
  mutate(unrep2_PYT_yield = emmeans_unrep_yield[[2]]$emmean,
         unrep3_PYT_yield = emmeans_unrep_yield[[3]]$emmean,
         unrep1_PYT_height = emmeans_unrep_height[[1]]$emmean,
         unrep2_PYT_height = emmeans_unrep_height[[2]]$emmean,
         unrep3_PYT_height = emmeans_unrep_height[[3]]$emmean,
         unrep1_PYT_hd = emmeans_unrep_hd[[1]]$emmean,
         unrep2_PYT_hd = emmeans_unrep_hd[[2]]$emmean,
         unrep3_PYT_hd = emmeans_unrep_hd[[3]]$emmean)

# Save
write.csv(pyt_unrep_emmeans, "Data/pheno/pyt_unrep_emmeans.csv", row.names = FALSE)
