
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for China using 2016 wave
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
script <- paste(stub,'6-Projections/rscripts/dmcf/regressions/standardised/with_continuous_urbanisation/', sep='')

# Load Household data
HH_China <- readRDS(paste(house,'China/china.rds', sep=''))

# Add urbanisation share
source(paste0(stub, "6-Projections/rscripts/process_raw_data/add_urban/add_urban_chn.R"))
HH_China$adm1 <- as.character(HH_China$state)

# Interaction prices
HH_China$mean_CDD18_db <- HH_China$meanpy_CDD18_db
HH_China$mean_hDD18_db <- HH_China$meanpy_hDD18_db
HH_China <- HH_China %>% mutate(ln_ely_p = log(ely_p_usd_2011),
                                ln_ely_p_cdd = ln_ely_p*mean_CDD18_db,
                                ln_ely_p_cdd2 = ln_ely_p*(mean_CDD18_db^2),
                                ln_ely_p_own = ln_ely_p*as.numeric(as.character(ownership_d)),
                                ln_ely_p_nme = ln_ely_p*n_members,
                                mean_CDD18_db2 = mean_CDD18_db^2,
                                mean_CDD18_db_exp = ln_total_exp_usd_2011*mean_CDD18_db,
                                mean_CDD18_db2_exp = ln_total_exp_usd_2011*(mean_CDD18_db^2),
                                curr_CDD18_db2 = curr_CDD18_db^2)

# Only those with not missing values 
HH_China <- HH_China[complete.cases(HH_China$ac), ]
HH_China <- HH_China[complete.cases(HH_China$mean_CDD18_db), ]
HH_China <- HH_China[complete.cases(HH_China$mean_HDD18_db), ]
HH_China <- HH_China[complete.cases(HH_China$curr_CDD18_db), ]
HH_China <- HH_China[complete.cases(HH_China$curr_HDD18_db), ]
HH_China <- HH_China[complete.cases(HH_China$ln_total_exp_usd_2011), ]
HH_China <- HH_China[complete.cases(HH_China$urban_sh), ]
HH_China <- HH_China[complete.cases(HH_China$n_members), ]
HH_China <- HH_China[complete.cases(HH_China$sh_under16), ]
HH_China <- HH_China[complete.cases(HH_China$housing_index_lab), ]
HH_China <- HH_China[complete.cases(HH_China$ownership_d), ]
HH_China <- HH_China[complete.cases(HH_China$edu_head_2), ]
HH_China <- HH_China[complete.cases(HH_China$age_head), ]
HH_China <- HH_China[complete.cases(HH_China$sex_head), ]
HH_China <- HH_China[complete.cases(HH_China$weight), ]
HH_China <- HH_China[complete.cases(HH_China$ln_ely_p), ]

# Scale variable
HH_China <- HH_China %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD18_db)),
                                std_CDD = as.numeric(scale(curr_CDD18_db)),
                                std_elyp = as.numeric(scale(ln_ely_p)),
                                std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                std_HDD = as.numeric(scale(curr_HDD18_db)),
                                std_urban_sh = as.numeric(scale(urban_sh)),
                                std_n_members = as.numeric(scale(n_members)),
                                std_age_head = as.numeric(scale(age_head)),
                                std_sh_under16 = as.numeric(scale(sh_under16)))

# Survey
HH_China_svy <- svydesign(data = HH_China, ids = ~adm1, weights = ~ weight)


##################################

#        Extensive margin        #

##################################

# AC formula for China
ac_formula_chn <- ac ~ std_CDD_mean + I(std_CDD_mean^2) + std_CDD_mean*std_texp + I(std_CDD_mean^2)*std_texp + std_texp + std_CDD + I(std_CDD^2) +
  std_elyp + std_elyp*std_CDD_mean + std_elyp*I(std_CDD_mean^2) + std_elyp*ownership_d + std_elyp*std_n_members + 
  std_urban_sh + std_n_members + ownership_d + edu_head_2 + std_age_head + sex_head + housing_index_lab + macroarea

# Logistic regression of AC on covariates
reg_ac <- svyglm(ac_formula_chn, design = HH_China_svy,
                 family = binomial(logit), na.action=na.omit); summary(reg_ac)

# Save AME results
margins <- margins(reg_ac, design = HH_China_svy)
ac_margins <- summary(margins)

# Predicted probabilities
HH_China$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_China$ac_obs <- ifelse(HH_China$phat0_obs>0.5 & !is.na(HH_China$phat0_obs), 1 , 0)

# Selection term for intensive margin part
HH_China$xb_noac = 1-HH_China$phat0_obs               
HH_China$selection = ifelse(HH_China$ac==1, 
                            (HH_China$xb_noac*log(HH_China$xb_noac)/HH_China$phat0_obs) + log(HH_China$phat0_obs), 
                            (HH_China$phat0_obs*log(HH_China$phat0_obs)/HH_China$xb_noac) + log(HH_China$xb_noac))

# Survey - re-run to add new variable
HH_China_svy <- svydesign(data = HH_China, ids = ~adm1, weights = ~ weight)


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
ely_formula_chn <- ln_ely_q ~ ac + ac*std_CDD + ac*I(std_CDD^2) + std_CDD + I(std_CDD^2) + 
  std_texp + std_HDD + I(std_HDD^2) + std_elyp + 
  std_urban_sh + std_n_members + ownership_d + edu_head_2 + 
  housing_index_lab + std_age_head + sex_head + selection + macroarea

# With selection
model <- svyglm(ely_formula_chn, design = HH_China_svy, na.action=na.omit); summary(model)

# Marginal effect of AC
ely_margins <- summary(margins(model, design = HH_China_svy))

# Export
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/chn_dmcf.RData', sep=''))
