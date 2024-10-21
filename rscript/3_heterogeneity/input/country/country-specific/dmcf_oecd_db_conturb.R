
# .rs.restartR()
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Load packages
library(data.table)
library(plyr)
library(dplyr)
library(FSA)
library(haven)
library(stringr)
library(tidyverse)
library(sandwich)
library(lmtest)
library(ResourceSelection)
library(multiwayvcov)
library(texreg)
library(xtable)
library(stargazer)
library(effects)
library(fixest)
library(marginaleffects)

# Set users
user <- 'user'

if (user=='user') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/'
}

house <- paste(stub,'data/household/', sep='')
interm <- paste(stub,'results/regressions/for_graphs/subsamples/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/'
output <- paste(stub,'output/figures/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/figures/'
script <- 'C:/Users/Standard/Documents/Github/acglobal/rscript/3_heterogeneity/input/country/country-specific/'


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


### Instead of pooling the countries: same regression but using clusters of country
# AC formula for OECD
ac_formula_oecd <- as.numeric(as.character(ac)) ~ mean_CDD18_db + mean_CDD18_db2 + curr_CDD18_db + curr_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 +
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +  
  curr_HDD18_db + I(curr_HDD18_db^2) +
  n_members + edu_head_2 + age_head + sex_head | country

## Logistic regression of AC on covariates
# Europe
reg_ac_eu <- feglm(ac_formula_oecd, family = binomial(link = "logit"), 
                   data = HH_Europe, weights = ~weight, cluster = c("adm1"))

# Predicted probabilities
HH_Europe$phat0_obs <- as.numeric(predict(reg_ac_eu, type="response"))

# Non Europe
reg_ac_noneu <- feglm(ac_formula_oecd, family = binomial(link = "logit"), 
                      data = HH_NonEurope, weights = ~weight, cluster = c("adm1"))

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

## Europe
# Formula without selection
ely_formula_oecd <- ln_ely_q ~ ac +   
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) | country

# Without selection
model0 <- feols(ely_formula_oecd, data = HH_Europe, weights = ~weight, cluster = c("adm1")); summary(model0)

# Formula with selection
ely_formula_oecd <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head | country

# With selection
model1 <- feols(ely_formula_oecd, data = HH_Europe, weights = ~weight, cluster = c("adm1")); summary(model1)

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model2 <- feols(ely_formula_oecd, data = HH_Europe, weights = ~weight, cluster = c("adm1")); summary(model2)

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model3 <- feols(ely_formula_oecd, data = HH_Europe, weights = ~weight, cluster = c("adm1")); summary(model3)

#  Marginal effect of AC
ac_eff_eu <- avg_slopes(model3, variables = "ac", slope = "dydx", wts = HH_Europe$weight)

# Save coefficients in data frame
dydx_ac_eu <- ac_eff_eu


## Non-Europe
# Formula without selection
ely_formula_oecd <- ln_ely_q ~ ac +   
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) | country

# Without selection
model4 <- feols(ely_formula_oecd, data = HH_NonEurope, weights = ~weight, cluster = c("adm1")); summary(model4)

# Formula with selection
ely_formula_oecd <- ln_ely_q ~ ac +   
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head  | country

# With selection
model5 <- feols(ely_formula_oecd, data = HH_NonEurope, weights = ~weight, cluster = c("adm1")); summary(model5)

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head  + selection | country

# With selection
model6 <- feols(ely_formula_oecd, data = HH_NonEurope, weights = ~weight, cluster = c("adm1")); summary(model6)

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) +
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head  + selection | country

# With selection
model7 <- feols(ely_formula_oecd, data = HH_NonEurope, weights = ~weight, cluster = c("adm1")); summary(model7)

#  Marginal effect of AC
ac_eff_neu <- avg_slopes(model7, variables = "ac", slope = "dydx", wts = HH_NonEurope$weight)

# Save coefficients in data frame
dydx_ac_neu <- ac_eff_neu


# Save the R Environment will be used for the projections
dydx_ac <- dydx_ac_eu
save(list = c("reg_ac_eu", "HH_Europe", "model3", "dydx_ac"), 
     file = paste(interm,'oecdeu_dmcf.RData', sep=''))
dydx_ac <- dydx_ac_neu
save(list = c("reg_ac_noneu", "HH_NonEurope", "model7", "dydx_ac"), 
     file = paste(interm,'oecdnoneu_dmcf.RData', sep=''))
