
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for OECD using 2011 wave
##      3) run intensive margin regressions: electricity expenditure on climate + covariates
##         using Dubin and McFadden (1984) approach

rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Load packages
library(data.table)
library(plyr)
library(dplyr)
library(FSA)
library(haven)
library(readstata13)
library(stringr)
library(tidyverse)
library(sandwich)
library(lmtest)
library(ResourceSelection)
library(multiwayvcov)
library(msm) # https://stats.oarc.ucla.edu/r/faq/how-can-i-estimate-the-standard-error-of-transformed-regression-parameters-in-r-using-the-delta-method/
library(margins) 
library(texreg)
library(xtable)
library(survey)
library(fixest)
library(tibble)
library(marginaleffects)


# Set directory
user <- 'fp'
#user <- 'gf'

if (user=='fp') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

if (user=='gf') {
  stub <- "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/"
}

house <- paste(stub,'6-Projections/data/household/', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')
script <- paste(stub,'6-Projections/rscripts/dmcf/regressions/country/with_continuous_urbanisation/', sep='')

# Load global data
global <- readRDS(paste(house,'global.rds', sep=''))

# Interaction prices
global <- global %>% mutate(ln_ely_p_cdd = ln_ely_p*mean_CDD18_db,
                            ln_ely_p_cdd2 = ln_ely_p*(mean_CDD18_db^2),
                            ln_ely_p_own = ln_ely_p*ownership_d,
                            ln_ely_p_nme = ln_ely_p*n_members,
                            mean_CDD18_db2 = mean_CDD18_db^2,
                            mean_CDD18_db_exp = ln_total_exp_usd_2011*mean_CDD18_db,
                            mean_CDD18_db2_exp = ln_total_exp_usd_2011*(mean_CDD18_db^2),
                            curr_CDD18_db2 = curr_CDD18_db^2,
                            edu_head_2 = as.factor(edu_head_2))

# Check
global <- global[complete.cases(global$ln_ely_q), ]
global <- global[complete.cases(global$ac), ]
global <- global[complete.cases(global$ln_total_exp_usd_2011), ]
global <- global[complete.cases(global$mean_CDD18_db), ]
global <- global[complete.cases(global$ownership_d), ]
global <- global[complete.cases(global$n_members), ]
global <- global[complete.cases(global$age_head), ]
global <- global[complete.cases(global$country), ]
global <- global[complete.cases(global$weight), ]
global <- global[complete.cases(global$sex_head), ]
global <- global[complete.cases(global$urban_sh), ]
global <- global[complete.cases(global$ln_ely_p), ]
global <- global[complete.cases(global$curr_CDD18_db), ]
global <- global[complete.cases(global$curr_HDD18_db), ]
global <- global[complete.cases(global$adm1), ]
global <- global %>% filter(ln_ely_q > 0)
global <- global %>% filter(weight > 0)

# Choosing cluster
HH_Europe <- global %>% filter(country == "Sweden" | country == "Spain" | country == "Netherlands" | country == "Switzerland " |
                                 country == "France")
HH_NonEurope <- global %>% filter(country == "Canada" | country == "Australia" | country == "Japan")

# Ref category for education
HH_Europe$edu_head_2 <- relevel(HH_Europe$edu_head_2,"1") # we do so since there are no edu = 0
HH_NonEurope$edu_head_2 <- relevel(HH_NonEurope$edu_head_2,"1") # we do so since there are no edu = 0

# Scale variable
HH_Europe <- HH_Europe %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD18_db)),
                              std_CDD = as.numeric(scale(curr_CDD18_db)),
                              std_elyp = as.numeric(scale(ln_ely_p)),
                              std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                              std_HDD = as.numeric(scale(curr_HDD18_db)),
                              std_urban_sh = as.numeric(scale(urban_sh)),
                              std_n_members = as.numeric(scale(n_members)),
                              std_age_head = as.numeric(scale(age_head)))

HH_NonEurope <- HH_NonEurope %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD18_db)),
                                  std_CDD = as.numeric(scale(curr_CDD18_db)),
                                  std_elyp = as.numeric(scale(ln_ely_p)),
                                  std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                  std_HDD = as.numeric(scale(curr_HDD18_db)),
                                  std_urban_sh = as.numeric(scale(urban_sh)),
                                  std_n_members = as.numeric(scale(n_members)),
                                  std_age_head = as.numeric(scale(age_head)))


# AC formula for OECD
ac_formula_oecd <- ac ~ std_CDD_mean + I(std_CDD_mean^2) + std_CDD_mean*std_texp + I(std_CDD_mean^2)*std_texp + std_texp + std_CDD + I(std_CDD^2) +
  std_elyp + std_elyp*std_CDD_mean + std_elyp*I(std_CDD_mean^2) + std_elyp*ownership_d + std_elyp*std_n_members + 
  std_urban_sh + std_n_members + ownership_d + edu_head_2 + std_age_head + sex_head | country

## Logistic regression of AC on covariates
# Europe
reg_ac_eu <- feglm(ac_formula_oecd, family = binomial(link = "logit"), 
                   data = HH_Europe, weights = ~weight, cluster = c("adm1"))

# Average marginal effects (AMEs)
ac_margins_eu <- summary(avg_slopes(reg_ac_eu, wts = HH_Europe$weight))
gc()

# Predicted probabilities
HH_Europe$phat0_obs <- as.numeric(predict(reg_ac_eu, type="response"))

# Non Europe
reg_ac_noneu <- feglm(ac_formula_oecd, family = binomial(link = "logit"), 
                      data = HH_NonEurope, weights = ~weight, cluster = c("adm1"))

# Average marginal effects (AMEs)
ac_margins_noneu <- summary(avg_slopes(reg_ac_noneu, wts = HH_NonEurope$weight))
gc()

# Predicted probabilities
HH_NonEurope$phat0_obs <- as.numeric(predict(reg_ac_noneu, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Europe$ac_obs <- ifelse(HH_Europe$phat0_obs>0.5 & !is.na(HH_Europe$phat0_obs), 1 , 0)
HH_NonEurope$ac_obs <- ifelse(HH_NonEurope$phat0_obs>0.5 & !is.na(HH_NonEurope$phat0_obs), 1 , 0)

# Selection term
HH_Europe$xb_noac = 1-HH_Europe$phat0_obs               
HH_Europe$selection = ifelse(HH_Europe$ac==1, 
                             (HH_Europe$xb_noac*log(HH_Europe$xb_noac)/HH_Europe$phat0_obs) + log(HH_Europe$phat0_obs), 
                             (HH_Europe$phat0_obs*log(HH_Europe$phat0_obs)/HH_Europe$xb_noac) + log(HH_Europe$xb_noac))

HH_NonEurope$xb_noac = 1-HH_NonEurope$phat0_obs               
HH_NonEurope$selection = ifelse(HH_NonEurope$ac==1, 
                                (HH_NonEurope$xb_noac*log(HH_NonEurope$xb_noac)/HH_NonEurope$phat0_obs) + log(HH_NonEurope$phat0_obs), 
                                (HH_NonEurope$phat0_obs*log(HH_NonEurope$phat0_obs)/HH_NonEurope$xb_noac) + log(HH_NonEurope$xb_noac))


# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + ac*std_CDD + ac*I(std_CDD^2) + std_CDD + I(std_CDD^2) + 
  std_texp + std_HDD + I(std_HDD^2) + std_elyp +
  std_urban_sh + std_n_members + 
  ownership_d + edu_head_2 + std_age_head + sex_head + selection | country

## Europe
# With selection
model <- feols(ely_formula_oecd, data = HH_Europe, weights = ~weight, cluster = c("adm1")); summary(model)

# Marginal effects
ely_margins_eu <- summary(avg_slopes(model, slope = "dydx", wts = HH_Europe$weight))

## Non-Europe
# With selection
model <- feols(ely_formula_oecd, data = HH_NonEurope, weights = ~weight, cluster = c("adm1")); summary(model)

# Marginal effects
ely_margins_noneu <- summary(avg_slopes(model, slope = "dydx", wts = HH_NonEurope$weight))

# Export
ely_margins <- ely_margins_eu
ac_margins <- ac_margins_eu

save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/oecdeu_dmcf.RData', sep=''))

ely_margins <- ely_margins_noneu
ac_margins <- ac_margins_noneu

save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/oecdnoneu_dmcf.RData', sep=''))
