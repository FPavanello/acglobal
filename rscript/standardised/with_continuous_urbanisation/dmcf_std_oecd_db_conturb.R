
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
library(stargazer)
library(lm.beta)
library(reghelper)
library(survey)


# Set users
user <- 'fp'
#user <- 'gf'

if (user=='fp') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

if (user=='gf') {
  stub <- 'H:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

house <- paste(stub,'6-Projections/data/household/', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')
script <- paste(stub,'6-Projections/rscripts/dmcf/regressions/standardised/with_continuous_urbanisation/', sep='')

# Load Household data
HH_OECD <- readRDS(paste(house,'OECD/EPIC/Household_OECD_mod.rds', sep=''))

# Add urbanisation share
source(paste0(stub, "6-Projections/rscripts/process_raw_data/add_urban/add_urban_oecd.R"))
HH_OECD$geometry <- NULL
HH_OECD$adm1 <- as.character(HH_OECD$state)

# Interaction prices
HH_OECD$mean_CDD18_db <- HH_OECD$meanpy_CDD18_db
HH_OECD$mean_hDD18_db <- HH_OECD$meanpy_hDD18_db
HH_OECD <- HH_OECD %>% mutate(ln_ely_p = log(ely_p_usd_2011),
                              ln_ely_p_cdd = ln_ely_p*mean_CDD18_db,
                              ln_ely_p_cdd2 = ln_ely_p*(mean_CDD18_db^2),
                              ln_ely_p_own = ln_ely_p*as.numeric(as.character(ownership_d)),
                              ln_ely_p_nme = ln_ely_p*n_members,
                              mean_CDD18_db2 = mean_CDD18_db^2,
                              mean_CDD18_db_exp = ln_total_exp_usd_2011*mean_CDD18_db,
                              mean_CDD18_db2_exp = ln_total_exp_usd_2011*(mean_CDD18_db^2),
                              curr_CDD18_db2 = curr_CDD18_db^2)

# Only those with not missing values 
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ac), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ln_ely_q), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$mean_CDD18_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$mean_HDD18_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$curr_CDD18_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$curr_HDD18_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ln_total_exp_usd_2011), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$urban_sh), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$n_members), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$share_under18), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$edu_head_2), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ownership_d), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ely_q_impute_5_95), ] # 3648 observations
HH_OECD <- HH_OECD[complete.cases(HH_OECD$age_head), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$sex_head), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ln_ely_p), ]
HH_OECD <- HH_OECD %>% filter(ln_ely_q > 0)
HH_OECD <- HH_OECD %>% filter(weight > 0)

# Choosing cluster
HH_Europe <- HH_OECD %>% filter(country == "Sweden" | country == "Spain" | country == "Netherlands" | country == "Switzerland " |
                                  country == "France")
HH_NonEurope <- HH_OECD %>% filter(country == "Canada" | country == "Australia" | country == "Japan")

HH_Europe$country2 <- as.factor(HH_Europe$country)
HH_NonEurope$country2 <- as.factor(HH_NonEurope$country)

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

# Survey
HH_Europe_svy <- svydesign(data = HH_Europe, ids = ~adm1, weights = ~ weight)
HH_NonEurope_svy <- svydesign(data = HH_NonEurope, ids = ~adm1, weights = ~ weight)


##################################

#        Extensive margin        #

##################################

### Instead of pooling the countries: same regression but using clusters of country
# AC formula for OECD
ac_formula_oecd <- ac ~ std_CDD_mean + I(std_CDD_mean^2) + std_CDD_mean*std_texp + I(std_CDD_mean^2)*std_texp + std_texp + std_CDD + I(std_CDD^2) +
  std_elyp + std_elyp*std_CDD_mean + std_elyp*I(std_CDD_mean^2) + std_elyp*ownership_d + std_elyp*std_n_members + 
  std_urban_sh + std_n_members + ownership_d + edu_head_2 + std_age_head + sex_head + country2

## Logistic regression of AC on covariates
# Europe
reg_ac_eu <- svyglm(ac_formula_oecd, design = HH_Europe_svy, family = binomial(logit), na.action=na.omit); summary(reg_ac_eu)

# Save AME results
margins <- margins(reg_ac_eu, design = HH_Europe_svy)
ac_margins_eu <- summary(margins)

# Predicted probabilities
HH_Europe$phat0_obs <- as.numeric(predict(reg_ac_eu, type="response"))

# Non Europe
reg_ac_noneu <- svyglm(ac_formula_oecd, design = HH_NonEurope_svy, family = binomial(logit), na.action=na.omit); summary(reg_ac_noneu)

# Save AME results
margins <- margins(reg_ac_noneu, design = HH_NonEurope_svy)
ac_margins_noneu <- summary(margins)

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

# Survey - re-run to add new variable
HH_Europe_svy <- svydesign(data = HH_Europe, ids = ~adm1, weights = ~ weight)
HH_NonEurope_svy <- svydesign(data = HH_NonEurope, ids = ~adm1, weights = ~ weight)


#################################################################

#     Intensive margin projections at 2040 + the role of AC

#     1) We interact AC with a set of variables to understand
#     how AC affects the adoption based on different charact.

#     2) We compute coefficients not only at the average, but
#     also based on specific values of our variables.
#     For instance, we compute the coefficients by decile, and
#     not only for the average household

#     Somehow point 1) is similar to a CDA, but without the
#     other appliances. For simplicity, I am going to interact
#     AC only with climate

#################################################################

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + ac*std_CDD + ac*I(std_CDD^2) + std_CDD + I(std_CDD^2) + 
  std_texp + std_HDD + I(std_HDD^2) + std_elyp +
  std_urban_sh + std_n_members + 
  ownership_d + edu_head_2 + std_age_head + sex_head + 
  country2 + selection

## Europe
# With selection
model <- svyglm(ely_formula_oecd, design = HH_Europe_svy, na.action=na.omit); summary(model)

# Marginal effects
ely_margins_eu <- summary(margins(model, design = HH_Europe_svy))

## Non-Europe
# With selection
model <- svyglm(ely_formula_oecd, design = HH_NonEurope_svy, na.action=na.omit); summary(model)

# Marginal effects
ely_margins_noneu <- summary(margins(model, design = HH_NonEurope_svy))

# Export
ely_margins <- ely_margins_eu
ac_margins <- ac_margins_eu

save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/oecdeu_dmcf.RData', sep=''))

ely_margins <- ely_margins_noneu
ac_margins <- ac_margins_noneu

save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/oecdnoneu_dmcf.RData', sep=''))
