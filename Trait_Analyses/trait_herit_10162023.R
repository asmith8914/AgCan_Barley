## Trait heritabilities
# Orig. 10/16/2023, tidied 03/28/2024
# Note: paths in this script are from original directory, entry adjusted means did not change so re-running this script was not necessary.

# Objective: Obtain heritabilities for traits.

library(tidyverse)
library(lme4)
library(emmeans)

# Piepho and Mohring 2007 ad hoc measure of broad sense heritability (equation 19): H = genetic variance / (genetic variance + (mean variance of differences between BLUEs/2))
# Schmidt et al 2019 (https://doi.org/10.2135/cropsci2018.06.0376) includes a github link which shows how to calculate heritability this way: https://github.com/PaulSchmidtGit/Heritability
# Also note that this method of calculating heritability reduces to the standard method if data is balanced.

# Read in data ----
pyt_data <- read.csv("AgCan_Barley_Data/phenodata/Phenodata_prelim_hd.csv")
ayt_data <- read.csv("AgCan_Barley_Data/phenodata/Phenodata_adv_hd.csv")

# Fix names + column names in ayt data
ayt_data <- ayt_data %>%
  mutate(NAME = gsub("'", "", NAME),
         NAME = gsub("AAC Ling", "AAC-Ling", NAME),
         NAME = gsub("AAC_Ling", "AAC-Ling", NAME),
         NAME = gsub("AAC Synergy", "AAC-Synergy", NAME)) %>%
  rename(yield = `YIELD_KGHA.1`, testweight = `TESTWEIGHT_hLKG.1`) %>%
  mutate(YEAR = as.factor(YEAR), REP = as.factor(REP))

# Fix column names in pyt data
pyt_data <- pyt_data %>%
  rename(yield = `YIELD_KGHA.1`, testweight = `TESTWEIGHT_hLKG.1`) %>%
  mutate(YEAR = as.factor(YEAR), REP = as.factor(REP))

# Fit fully random models (to obtain genetic variance) ----
## Models
# PYTs
m_r_pyt_y <- formula(yield ~ (1|NAME) + (1|YEAR) + (1|BLOCK:YEAR))
m_r_pyt_ph <- formula(HT_CM ~ (1|NAME) + (1|YEAR) + (1|BLOCK:YEAR))
m_r_pyt_hd <- formula(days_to_heading ~ (1|NAME) + (1|YEAR) + (1|BLOCK:YEAR))
m_r_pyt_sw <- formula(SEEDWEIGHT_G ~ (1|NAME) + (1|YEAR))
m_r_pyt_tw <- formula(testweight ~ (1|NAME) + (1|YEAR))

# AYTs
m_r_ayt_y <- formula(yield ~ (1|NAME) + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
m_r_ayt_ph <- formula(HT_CM ~ (1|NAME) + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
m_r_ayt_sw <- formula(SEEDWEIGHT_G ~ (1|NAME) + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
m_r_ayt_tw <- formula(testweight ~ (1|NAME) + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
m_r_ayt_hd <- formula(days_to_heading ~ (1|NAME) + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))

## Apply models
# PYTs
lmr_pyt_y <- lmer(data = pyt_data, formula = m_r_pyt_y)
lmr_pyt_ph <- lmer(data = pyt_data, formula = m_r_pyt_ph)
lmr_pyt_hd <- lmer(data = pyt_data, formula = m_r_pyt_hd)
lmr_pyt_sw <- lmer(data = pyt_data, formula = m_r_pyt_sw)
lmr_pyt_tw <- lmer(data = pyt_data, formula = m_r_pyt_tw)

# AYTs
lmr_ayt_y <- lmer(data = ayt_data, formula = m_r_ayt_y)
lmr_ayt_ph <- lmer(data = ayt_data, formula = m_r_ayt_ph)
lmr_ayt_hd <- lmer(data = ayt_data, formula = m_r_ayt_hd)
lmr_ayt_sw <- lmer(data = ayt_data, formula = m_r_ayt_sw)
lmr_ayt_tw <- lmer(data = ayt_data, formula = m_r_ayt_tw)

## Extract genetic variance
# Trait vector
trait_vec <- c("y", "ph", "hd", "sw", "tw")

# PYTs
pyt_gvar <- as.data.frame(do.call(cbind, lapply(c(1:5), function(x){
  get(paste0("lmr_pyt_", trait_vec[x])) %>% 
    VarCorr %>% as_tibble %>% filter(grp == "NAME") %>% pull(vcov)
})))
colnames(pyt_gvar) <- paste0("pyt_", trait_vec)

# AYTs
ayt_gvar <- as.data.frame(do.call(cbind, lapply(c(1:5), function(x){
  get(paste0("lmr_ayt_", trait_vec[x])) %>% 
    VarCorr %>% as_tibble %>% filter(grp == "NAME") %>% pull(vcov)
})))
colnames(ayt_gvar) <- paste0("ayt_", trait_vec)

# Fit model with fixed genotypes (to obtain mean variance of differences between BLUEs) ----

# Year was treated as fixed when obtaining BLUEs for genomic predictions on the basis that there are only 3 years, so not enough levels to treat as random.
# Should year be treated as fixed or random for calculation of mean variance of a difference of two genotypic BLUEs for heritability??
# Using random year for now, since comparing with genetic variance where year was random.

## Models
# PYTs
m_f_pyt_y <- formula(yield ~ NAME + (1|YEAR) + (1|BLOCK:YEAR))
m_f_pyt_ph <- formula(HT_CM ~ NAME + (1|YEAR) + (1|BLOCK:YEAR))
m_f_pyt_hd <- formula(days_to_heading ~ NAME + (1|YEAR) + (1|BLOCK:YEAR))
m_f_pyt_sw <- formula(SEEDWEIGHT_G ~ NAME + (1|YEAR))
m_f_pyt_tw <- formula(testweight ~ NAME + (1|YEAR))

# AYTs
m_f_ayt_y <- formula(yield ~ NAME + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
m_f_ayt_ph <- formula(HT_CM ~ NAME + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
m_f_ayt_sw <- formula(SEEDWEIGHT_G ~ NAME + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
m_f_ayt_tw <- formula(testweight ~ NAME + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
m_f_ayt_hd <- formula(days_to_heading ~ NAME + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))

## Apply models
# PYTs
lmf_pyt_y <- lmer(data = pyt_data, formula = m_f_pyt_y)
lmf_pyt_ph <- lmer(data = pyt_data, formula = m_f_pyt_ph)
lmf_pyt_hd <- lmer(data = pyt_data, formula = m_f_pyt_hd)
lmf_pyt_sw <- lmer(data = pyt_data, formula = m_f_pyt_sw)
lmf_pyt_tw <- lmer(data = pyt_data, formula = m_f_pyt_tw)

# AYTs
lmf_ayt_y <- lmer(data = ayt_data, formula = m_f_ayt_y)
lmf_ayt_ph <- lmer(data = ayt_data, formula = m_f_ayt_ph)
lmf_ayt_hd <- lmer(data = ayt_data, formula = m_f_ayt_hd)
lmf_ayt_sw <- lmer(data = ayt_data, formula = m_f_ayt_sw)
lmf_ayt_tw <- lmer(data = ayt_data, formula = m_f_ayt_tw)

## Obtain the mean variances of a difference of two genotypic BLUEs
# AYTs
ayt_mvd <- as.data.frame(do.call(cbind, lapply(c(1:5), function(x){
  get(paste0("lmf_ayt_", trait_vec[x])) %>% 
    emmeans(pairwise ~ NAME) %>% pluck("contrasts") %>% as_tibble %>% mutate(Var = SE^2) %>% pull(Var) %>% mean
})))
colnames(ayt_mvd) <- paste0("ayt_", trait_vec)

# PYTs
pyt_mvd <- as.data.frame(do.call(cbind, lapply(c(1:5), function(x){
  get(paste0("lmf_pyt_", trait_vec[x])) %>% 
    emmeans(pairwise ~ NAME) %>% pluck("contrasts") %>% as_tibble %>% mutate(Var = SE^2) %>% pull(Var) %>% mean
})))
# NOTE: VERY SLOW (ran it overnight)
colnames(pyt_mvd) <- paste0("pyt_", trait_vec)

# Calculate broad sense heritability ----
ayt_herit <- as.data.frame(do.call(cbind, lapply(c(1:5), function(x){
  get("ayt_gvar")[paste0("ayt_", trait_vec[x])] / (get("ayt_gvar")[paste0("ayt_", trait_vec[x])] + (get("ayt_mvd")[paste0("ayt_", trait_vec[x])]/2))
})))

pyt_herit <- as.data.frame(do.call(cbind, lapply(c(1:5), function(x){
  get("pyt_gvar")[paste0("pyt_", trait_vec[x])] / (get("pyt_gvar")[paste0("pyt_", trait_vec[x])] + (get("pyt_mvd")[paste0("pyt_", trait_vec[x])]/2))
})))

# Save the genetic variances, MVDs, and heritabilities:
sink("AgCan_Barley_Data/phenodata/gvar_mvd_herit.txt")
print("PYT genetic variances")
print(pyt_gvar)
print("AYT genetic variances")
print(ayt_gvar)
print("PYT mean variances of difference of genotypic BLUEs")
print(pyt_mvd)
print("AYT mean variances of difference of genotypic BLUEs")
print(ayt_mvd)
print("PYT broad sense heritability")
print(pyt_herit)
print("AYT broad sense heritability")
print(ayt_herit)
sink()
