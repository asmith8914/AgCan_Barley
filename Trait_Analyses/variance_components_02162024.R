## Variance components
# 10/12/2023

# Objectives: Re-calculate variance components and round for tables

library(tidyverse)
library(lme4)

# 1) Calculate variance components ----
# Read in data
pyt_data <- read.csv("Data/pheno/Phenodata_prelim_hd.csv")
ayt_data <- read.csv("Data/pheno/Phenodata_adv_hd_cor.csv")

# Set factors
ayt_data <- ayt_data %>%
  mutate(YEAR = as.factor(YEAR), REP = as.factor(REP))

pyt_data <- pyt_data %>%
  mutate(YEAR = as.factor(YEAR), REP = as.factor(REP), BLOCK = as.factor(BLOCK))

# Models used for adjusted entry means across PYTs were: trait ~ NAME + YEAR + (1|BLOCK:YEAR) and trait ~ NAME + YEAR
# Make genotype random but otherwise the same models
model_pyt_y <- formula(yield ~ (1|NAME) + (1|YEAR) + (1|BLOCK:YEAR))
model_pyt_ph <- formula(HT_CM ~ (1|NAME) + (1|YEAR) + (1|BLOCK:YEAR))
model_pyt_hd <- formula(days_to_heading ~ (1|NAME) + (1|YEAR) + (1|BLOCK:YEAR))
model_pyt_sw <- formula(SEEDWEIGHT_G ~ (1|NAME) + (1|YEAR))
model_pyt_tw <- formula(testweight ~ (1|NAME) + (1|YEAR))

# Model used for adjusted entry means across AYTs was: trait ~ NAME + YEAR + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR)
# Make genotype random but otherwise the same model
model_adv_y <- formula(yield ~ (1|NAME) + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
model_adv_ph <- formula(HT_CM ~ (1|NAME) + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
model_adv_sw <- formula(SEEDWEIGHT_G ~ (1|NAME) + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
model_adv_tw <- formula(testweight ~ (1|NAME) + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))
model_adv_hd <- formula(days_to_heading ~ (1|NAME) + (1|YEAR) + (1|LOCATION:YEAR) + (1|REP:LOCATION:YEAR) + (1|NAME:LOCATION:YEAR))

# Fit the models
lm_pyt_y <- lmer(data = pyt_data, formula = model_pyt_y)
lm_pyt_ph <- lmer(data = pyt_data, formula = model_pyt_ph)
lm_pyt_hd <- lmer(data = pyt_data, formula = model_pyt_hd)
lm_pyt_sw <- lmer(data = pyt_data, formula = model_pyt_sw)
lm_pyt_tw <- lmer(data = pyt_data, formula = model_pyt_tw)

lm_ayt_y <- lmer(data = ayt_data, formula = model_adv_y)
lm_ayt_ph <- lmer(data = ayt_data, formula = model_adv_ph)
lm_ayt_hd <- lmer(data = ayt_data, formula = model_adv_hd)
lm_ayt_sw <- lmer(data = ayt_data, formula = model_adv_sw)
lm_ayt_tw <- lmer(data = ayt_data, formula = model_adv_tw)

summary(lm_pyt_y)
summary(lm_pyt_ph)
summary(lm_pyt_hd)
summary(lm_pyt_tw)
summary(lm_pyt_sw)

summary(lm_ayt_y)
summary(lm_ayt_ph)
summary(lm_ayt_hd)
summary(lm_ayt_tw)
summary(lm_ayt_sw)

# Obtain the variance components, make rounding consistent.

sink("Data/pheno/var_comp_rounded.txt")
print("PYT yield", quote = FALSE)
print(model_pyt_y)
data.frame(VarCorr(lm_pyt_y)) %>%
  mutate(total_var = sum(vcov),
         prop_var = vcov/total_var,
         perc_var = sprintf(prop_var*100, fmt = "%.3f"),
         total_perc = sum(as.numeric(perc_var)))

print("PYT height", quote = FALSE)
print(model_pyt_ph)
data.frame(VarCorr(lm_pyt_ph)) %>%
  mutate(total_var = sum(vcov),
         prop_var = vcov/total_var,
         perc_var = sprintf(prop_var*100, fmt = "%.3f"),
         total_perc = sum(as.numeric(perc_var)))

print("PYT days to heading", quote = FALSE)
print(model_pyt_hd)
data.frame(VarCorr(lm_pyt_hd)) %>%
  mutate(total_var = sum(vcov),
         prop_var = vcov/total_var,
         perc_var = sprintf(prop_var*100, fmt = "%.3f"),
         total_perc = sum(as.numeric(perc_var)))

print("PYT seed weight", quote = FALSE)
print(model_pyt_sw)
data.frame(VarCorr(lm_pyt_sw)) %>%
  mutate(total_var = sum(vcov),
         prop_var = vcov/total_var,
         perc_var = sprintf(prop_var*100, fmt = "%.3f"),
         total_perc = sum(as.numeric(perc_var)))

print("PYT test weight", quote = FALSE)
print(model_pyt_tw)
data.frame(VarCorr(lm_pyt_tw)) %>%
  mutate(total_var = sum(vcov),
         prop_var = vcov/total_var,
         perc_var = sprintf(prop_var*100, fmt = "%.3f"),
         total_perc = sum(as.numeric(perc_var)))

print("AYT yield", quote = FALSE)
print(model_adv_y)
data.frame(VarCorr(lm_ayt_y)) %>%
  mutate(total_var = sum(vcov),
         prop_var = vcov/total_var,
         perc_var = sprintf(prop_var*100, fmt = "%.3f"),
         total_perc = sum(as.numeric(perc_var)))

print("AYT height", quote = FALSE)
print(model_adv_ph)
data.frame(VarCorr(lm_ayt_ph)) %>%
  mutate(total_var = sum(vcov),
         prop_var = vcov/total_var,
         perc_var = sprintf(prop_var*100, fmt = "%.3f"),
         total_perc = sum(as.numeric(perc_var)))

print("AYT days to heading", quote = FALSE)
print(model_adv_hd)
data.frame(VarCorr(lm_ayt_hd)) %>%
  mutate(total_var = sum(vcov),
         prop_var = vcov/total_var,
         perc_var = sprintf(prop_var*100, fmt = "%.3f"),
         total_perc = sum(as.numeric(perc_var)))

print("AYT seed weight", quote = FALSE)
print(model_adv_sw)
data.frame(VarCorr(lm_ayt_sw)) %>%
  mutate(total_var = sum(vcov),
         prop_var = vcov/total_var,
         perc_var = sprintf(prop_var*100, fmt = "%.3f"),
         total_perc = sum(as.numeric(perc_var)))

print("AYT test weight", quote = FALSE)
print(model_adv_tw)
data.frame(VarCorr(lm_ayt_tw)) %>%
  mutate(total_var = sum(vcov),
         prop_var = vcov/total_var,
         perc_var = sprintf(prop_var*100, fmt = "%.3f"),
         total_perc = sum(as.numeric(perc_var)))
sink()
