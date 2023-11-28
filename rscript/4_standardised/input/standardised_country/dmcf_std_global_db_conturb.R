
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for the global data set
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
library(tibble)

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

# Load global data
global <- readRDS(paste(house,'global.rds', sep=''))

# Interaction prices
global$mean_CDD18_db <- global$meanpy_CDD18_db
global$mean_hDD18_db <- global$meanpy_hDD18_db
global <- global %>% mutate(ln_ely_p_cdd = ln_ely_p*mean_CDD18_db,
                            ln_ely_p_cdd2 = ln_ely_p*(mean_CDD18_db^2),
                            ln_ely_p_own = ln_ely_p*as.numeric(as.character(ownership_d)),
                            ln_ely_p_nme = ln_ely_p*n_members,
                            mean_CDD18_db2 = mean_CDD18_db^2,
                            mean_CDD18_db_exp = ln_total_exp_usd_2011*mean_CDD18_db,
                            mean_CDD18_db2_exp = ln_total_exp_usd_2011*(mean_CDD18_db^2),
                            curr_CDD18_db2 = curr_CDD18_db^2)

# Check
global <- global[complete.cases(global$ln_ely_q), ]
global <- global[complete.cases(global$ac), ]
global <- global[complete.cases(global$ln_total_exp_usd_2011), ]
global <- global[complete.cases(global$mean_CDD18_db), ]
global <- global[complete.cases(global$urban), ]
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

# Scale variable
global <- global %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD18_db)),
                            std_CDD = as.numeric(scale(curr_CDD18_db)),
                            std_elyp = as.numeric(scale(ln_ely_p)),
                            std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                            std_HDD = as.numeric(scale(curr_HDD18_db)),
                            std_urban_sh = as.numeric(scale(urban_sh)),
                            std_n_members = as.numeric(scale(n_members)),
                            std_age_head = as.numeric(scale(age_head)))

# Survey
global_svy <- svydesign(data = global, ids = ~ adm1, weights = ~ weight)


##################################

#        Extensive margin        #

##################################

# AC formula for global
ac_formula <- ac ~ std_CDD_mean + I(std_CDD_mean^2) + std_CDD_mean*std_texp + I(std_CDD_mean^2)*std_texp + std_texp +
  std_elyp + std_elyp*std_CDD_mean + std_elyp*I(std_CDD_mean^2) + std_CDD + I(std_CDD^2) + std_HDD + I(std_HDD^2) +
  std_elyp*ownership_d + std_elyp*std_n_members + std_urban_sh + ownership_d + 
  std_n_members + edu_head_2 + std_age_head + sex_head + country

# Logistic regression of AC on covariates
reg_ac <- svyglm(ac_formula, design = global_svy, 
                 family = binomial(logit), na.action=na.omit); summary(reg_ac)

# Save AME results
margins <- margins(reg_ac, design = global_svy)
ac_margins <- summary(margins)

# Predicted probabilities
global$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
global$ac_obs <- ifelse(global$phat0_obs>0.5 & !is.na(global$phat0_obs), 1 , 0)

# Selection term
global$xb_noac = 1-global$phat0_obs               
global$selection = ifelse(global$ac==1, 
                          (global$xb_noac*log(global$xb_noac)/global$phat0_obs) + log(global$phat0_obs), 
                          (global$phat0_obs*log(global$phat0_obs)/global$xb_noac) + log(global$xb_noac))

# Survey - re-run to add new variable
global_svy <- svydesign(data = global, ids = ~ adm1, weights = ~ weight)


#################################################################

#     Intensive margin + the role of AC

#     1) We interact AC with a set of variables to understand
#     how AC affects the adoption based on different charact.

#     2) We compute coefficients not only at the averages, but
#     also based on specific values of our variables.
#     For instance, we compute the coefficients by decile, and
#     not only for the average household

#     Somehow point 1) is similar to a CDA, but without the
#     other appliances. For simplicity, I am going to interact
#     AC only with climate

#################################################################

# Formula electricity expenditure without selection
ely_formula  <- ln_ely_q ~ ac + ac*std_CDD + ac*I(std_CDD^2) + std_CDD + I(std_CDD^2) + 
  std_texp + std_HDD + I(std_HDD^2) + std_elyp + std_urban_sh + ownership_d + 
  std_n_members + edu_head_2 + std_age_head + sex_head + selection + country

# Without selection
model <- svyglm(ely_formula, design = global_svy, na.action=na.omit); summary(model)

# Marginal effect of AC
ely_margins <- summary(margins(model, design = global_svy))

# Export
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/global_dmcf.RData', sep=''))
