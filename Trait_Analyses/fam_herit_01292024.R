## Obtain heritabilities of family means
# 01/29/2024

# Objective: Obtain the heritability of family means for each trait, to use as reliability (weight for combined predictions).

library(tidyverse)
library(lme4)
library(emmeans)

# Calculate broad sense heritability of family means
# Use the same equation for unbalanced designs that was used for individual genotypes: H(fam) = genetic variance(fam) / (genetic variance(fam) + (mean variance of differences between BLUEs(fam)/2)) (Piepho and Mohring 2007, Holland et al. 2003, pg. 64).
# Schmidt et al 2019 (https://doi.org/10.2135/cropsci2018.06.0376) includes a github link which shows how to calculate heritability this way: https://github.com/PaulSchmidtGit/Heritability

# Read in data
pyt_data <- read.csv("Data/pheno/Phenodata_prelim_hd.csv")

# Make rep, block, year factors in pyt data
pyt_data <- pyt_data %>%
  mutate(YEAR = as.factor(YEAR), REP = as.factor(REP), BLOCK = as.factor(BLOCK))

# Add column for cross
pyt_data <- pyt_data %>%
  mutate(cross = case_when(NAME == "AAC-Ling" ~ "AAC-Ling", NAME == "Leader" ~ "Leader", NAME != c("Leader", "AAC-Ling") ~ gsub("-.*", "", NAME)))

# Fit fully random models (to obtain genetic variance of family)
# Models are: y = mean + family + line:family + year + block:year (if applicable)
m_r_pyt_y <- formula(yield ~ (1|cross) + (1|NAME:cross) + (1|YEAR) + (1|BLOCK:YEAR))
m_r_pyt_ph <- formula(HT_CM ~ (1|cross) + (1|NAME:cross) + (1|YEAR) + (1|BLOCK:YEAR))
m_r_pyt_hd <- formula(days_to_heading ~ (1|cross) + (1|NAME:cross) + (1|YEAR) + (1|BLOCK:YEAR))
m_r_pyt_sw <- formula(SEEDWEIGHT_G ~ (1|cross) + (1|NAME:cross) + (1|YEAR))
m_r_pyt_tw <- formula(testweight ~ (1|cross) + (1|NAME:cross) + (1|YEAR))

# Apply models
lmr_pyt_y <- lmer(data = pyt_data, formula = m_r_pyt_y)
lmr_pyt_ph <- lmer(data = pyt_data, formula = m_r_pyt_ph)
lmr_pyt_hd <- lmer(data = pyt_data, formula = m_r_pyt_hd)
lmr_pyt_sw <- lmer(data = pyt_data, formula = m_r_pyt_sw)
lmr_pyt_tw <- lmer(data = pyt_data, formula = m_r_pyt_tw)

# Trait vector
trait_vec <- c("y", "ph", "hd", "sw", "tw")

# Model summaries:
lapply(c(1:5), function(x){
  summary(get(paste0("lmr_pyt_", trait_vec[x])))
})
# Note: yield is still singular... Year variance is zero.

# Extract genetic variance
pyt_gvar <- as.data.frame(do.call(cbind, lapply(c(1:5), function(x){
  get(paste0("lmr_pyt_", trait_vec[x])) %>% 
    VarCorr %>% as_tibble %>% filter(grp == "cross") %>% pull(vcov)
})))
colnames(pyt_gvar) <- paste0("pyt_", trait_vec)

# Fit model with fixed genotypes (to obtain mean variance of differences between BLUEs)
# Models
m_f_pyt_y <- formula(yield ~ cross + (1|NAME:cross) + (1|YEAR) + (1|BLOCK:YEAR))
m_f_pyt_ph <- formula(HT_CM ~ cross + (1|NAME:cross) + (1|YEAR) + (1|BLOCK:YEAR))
m_f_pyt_hd <- formula(days_to_heading ~ cross + (1|NAME:cross) + (1|YEAR) + (1|BLOCK:YEAR))
m_f_pyt_sw <- formula(SEEDWEIGHT_G ~ cross + (1|NAME:cross) + (1|YEAR))
m_f_pyt_tw <- formula(testweight ~ cross + (1|NAME:cross) + (1|YEAR))

# Apply models
lmf_pyt_y <- lmer(data = pyt_data, formula = m_f_pyt_y)
lmf_pyt_ph <- lmer(data = pyt_data, formula = m_f_pyt_ph)
lmf_pyt_hd <- lmer(data = pyt_data, formula = m_f_pyt_hd)
lmf_pyt_sw <- lmer(data = pyt_data, formula = m_f_pyt_sw)
lmf_pyt_tw <- lmer(data = pyt_data, formula = m_f_pyt_tw)

# Obtain the mean variances of a difference of two genotypic BLUEs
pyt_mvd <- as.data.frame(do.call(cbind, lapply(c(1:5), function(x){
  get(paste0("lmf_pyt_", trait_vec[x])) %>% 
    emmeans(pairwise ~ cross) %>% pluck("contrasts") %>% as_tibble %>% mutate(Var = SE^2) %>% pull(Var) %>% mean
})))
# Note: quite slow (>15 min)
colnames(pyt_mvd) <- paste0("pyt_", trait_vec)

# Calculate broad sense heritability
pyt_herit <- as.data.frame(do.call(cbind, lapply(c(1:5), function(x){
  get("pyt_gvar")[paste0("pyt_", trait_vec[x])] / (get("pyt_gvar")[paste0("pyt_", trait_vec[x])] + (get("pyt_mvd")[paste0("pyt_", trait_vec[x])]/2))
})))

# Save the genetic variances, MVDs, and heritabilities:
sink("Data/pheno/fam_mean_var_mvd_herit.txt")
print("PYT genetic variances of family means")
print(pyt_gvar)
print("PYT mean variances of difference of family BLUEs")
print(pyt_mvd)
print("PYT broad sense heritability of family means")
print(pyt_herit)
sink()

# Should also save as tables:
write.csv(pyt_herit, "Data/pheno/fam_mean_herit.csv", row.names = FALSE)
write.csv(pyt_mvd, "Data/pheno/fam_mean_mvd.csv", row.names = FALSE)
write.csv(pyt_gvar, "Data/pheno/fam_mean_var.csv", row.names = FALSE)
