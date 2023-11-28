
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for Africa using 2016 wave
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
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}


house <- paste(stub,'6-Projections/data/household/', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')
script <- paste(stub,'6-Projections/rscripts/dmcf/regressions/standardised/with_continuous_urbanisation/', sep='')

## Load Household data - Africa
HH_Nigeria <- readRDS(paste(house,'Nigeria/nga/Nigeria_ghs.rds', sep=''))
HH_Tanzania <- readRDS(paste(house,'Tanzania/HBS/Tanzania_hbs.rds', sep=''))
HH_Kenya <- readRDS(paste(house,'Kenya/IHBS/kenya_ihbs.rds', sep=''))
HH_Ghana <- readRDS(paste(house,'Ghana/ghana.rds', sep=''))
HH_Niger <- readRDS(paste(house,'Niger/HLCA/niger_hlca.rds', sep=''))
HH_Malawi <- readRDS(paste(house,'Malawi/IHS/malawi_ihs.rds', sep=''))
HH_Burkina <- readRDS(paste(house,'Burkina Faso/EMC/burkina_emc.rds', sep=''))

# Adjust gender head for Ghana
HH_Ghana$sex_head <- ifelse(HH_Ghana$sex_head==1, 0, 1)

# Country
HH_Nigeria <- HH_Nigeria %>% mutate(country = "Nigeria")
HH_Tanzania <- HH_Tanzania %>% mutate(country ="Tanzania")
HH_Kenya <- HH_Kenya %>% mutate(country = "Kenya")
HH_Ghana <- HH_Ghana %>% mutate(country = "Ghana")
HH_Burkina <- HH_Burkina %>% mutate(country = "Burkina Faso")
HH_Malawi <- HH_Malawi %>% mutate(country = "Malawi")
HH_Niger <- HH_Niger %>% mutate(country = "Niger")

# Append
HH_Africa <- rbind.fill(HH_Nigeria, HH_Tanzania)
HH_Africa <- rbind.fill(HH_Africa, HH_Kenya)
HH_Africa <- rbind.fill(HH_Africa, HH_Ghana)
HH_Africa <- rbind.fill(HH_Africa, HH_Burkina)
HH_Africa <- rbind.fill(HH_Africa, HH_Malawi)
HH_Africa <- rbind.fill(HH_Africa, HH_Niger)

# Arrange some vars
HH_Africa$mean_CDD_db <- as.numeric(HH_Africa$mean_CDD_db)
HH_Africa$mean_HDD_db <- as.numeric(HH_Africa$mean_HDD_db)
HH_Africa$curr_CDD_db <- as.numeric(HH_Africa$curr_CDD_db)
HH_Africa$curr_HDD_db <- as.numeric(HH_Africa$curr_HDD_db)
HH_Africa$curr_CDD18_db <- as.numeric(HH_Africa$curr_CDD18_db)
HH_Africa$mean_CDD18_db <- as.numeric(HH_Africa$mean_CDD18_db)
HH_Africa$curr_HDD18_db <- as.numeric(HH_Africa$curr_HDD18_db)
HH_Africa$mean_HDD18_db <- as.numeric(HH_Africa$mean_HDD18_db)
HH_Africa$state <- as.factor(HH_Africa$state)
HH_Africa$urban <- as.factor(as.character(HH_Africa$urban))
HH_Africa$ac <- as.factor(as.character(HH_Africa$ac))
HH_Africa$ref <- as.factor(as.character(HH_Africa$ref))
HH_Africa$pc <- as.factor(as.character(HH_Africa$pc))
HH_Africa$tv <- as.factor(as.character(HH_Africa$tv))
HH_Africa$wshm <- as.factor(as.character(HH_Africa$wshm))
HH_Africa$fan <- as.factor(as.character(HH_Africa$fan))
HH_Africa$ownership <- as.factor(as.character(HH_Africa$ownership))
HH_Africa$edu_head_2 <- as.factor(as.character(HH_Africa$edu_head_2))
HH_Africa$sex_head <- as.factor(as.character(HH_Africa$sex_head))
HH_Africa$housing_index_lab <- as.factor(as.character(HH_Africa$housing_index_lab))

# Adjust some variables
HH_Africa$ln_total_exp_usd_2011 <- ifelse(HH_Africa$ln_total_exp_usd_2011>13, NA, HH_Africa$ln_total_exp_usd_2011)
HH_Africa$age_head <- ifelse(HH_Africa$age_head ==999, NA, HH_Africa$age_head)
HH_Africa$ln_total_exp_usd_2011[HH_Africa$total_exp_usd_2011 == 0] <- NA
HH_Africa$ln_ely_q[HH_Africa$ln_ely_q < 0] <- NA
HH_Africa$housing_index_lab <- as.factor(HH_Africa$housing_index_lab)

# Minimum sub-national level
HH_Africa <- HH_Africa %>% mutate(subnat = as.character(state)) 
HH_Africa$subnat[HH_Africa$country == "Malawi"] <- as.character(HH_Africa$district[HH_Africa$country == "Malawi"])
HH_Africa$subnat[HH_Africa$country == "Burkina Faso"] <- as.character(HH_Africa$district[HH_Africa$country == "Burkina Faso"])
HH_Africa$subnat[HH_Africa$country == "Nigeria"] <- as.character(HH_Africa$lga[HH_Africa$country == "Nigeria"])
HH_Africa$adm1 <- as.character(HH_Africa$state)

# Add urbanisation share
source(paste0(stub, "6-Projections/rscripts/process_raw_data/add_urban/add_urban_AFR.R"))

# Only those with not missing values 
HH_Africa <- HH_Africa[complete.cases(HH_Africa$ac), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$mean_CDD_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$mean_HDD_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$curr_CDD_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$curr_HDD_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$mean_CDD18_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$mean_HDD18_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$curr_CDD18_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$curr_HDD18_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$ln_total_exp_usd_2011), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$urban_sh), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$n_members), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$housing_index_lab), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$ownership_d), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$edu_head_2), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$age_head), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$sex_head), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$country), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$adm1), ]

# CDD in 100s
HH_Africa <- HH_Africa %>% mutate(mean_CDD_db = mean_CDD_db/100,
                                  mean_HDD_db = mean_HDD_db/100,
                                  curr_CDD_db = curr_CDD_db/100,
                                  curr_HDD_db = curr_HDD_db/100,
                                  mean_CDD18_db = mean_CDD18_db/100,
                                  mean_HDD18_db = mean_HDD18_db/100,
                                  curr_CDD18_db = curr_CDD18_db/100,
                                  curr_HDD18_db = curr_HDD18_db/100,
                                  meanpy_CDD_db = meanpy_CDD_db/100,
                                  meanpy_HDD_db = meanpy_HDD_db/100,
                                  meanpy_CDD18_db = meanpy_CDD18_db/100,
                                  meanpy_HDD18_db = meanpy_HDD18_db/100)

# Check
HH_Africa <- HH_Africa[(is.finite(HH_Africa$ln_ely_q) & !is.na(HH_Africa$ln_ely_q)),] # only those who have access to electricity
unique.default(as.factor(HH_Africa$country))
summary(as.factor(HH_Africa$country))
rm(HH_Burkina, HH_Ghana, HH_Niger, HH_Kenya, HH_Nigeria, HH_Malawi, HH_Tanzania)

# Interaction prices
HH_Africa$mean_CDD18_db <- HH_Africa$meanpy_CDD18_db
HH_Africa$mean_hDD18_db <- HH_Africa$meanpy_hDD18_db
HH_Africa <- HH_Africa %>% mutate(ln_ely_p = log(ely_p_usd_2011),
                                  ln_ely_p_cdd = ln_ely_p*mean_CDD18_db,
                                  ln_ely_p_cdd2 = ln_ely_p*(mean_CDD18_db^2),
                                  ln_ely_p_own = ln_ely_p*as.numeric(as.character(ownership_d)),
                                  ln_ely_p_nme = ln_ely_p*n_members,
                                  mean_CDD18_db2 = mean_CDD18_db^2,
                                  mean_CDD18_db_exp = ln_total_exp_usd_2011*mean_CDD18_db,
                                  mean_CDD18_db2_exp = ln_total_exp_usd_2011*(mean_CDD18_db^2),
                                  curr_CDD18_db2 = curr_CDD18_db^2)

# Select variables
HH_Africa <- dplyr::select(HH_Africa, ac, ln_ely_q, mean_CDD18_db, mean_HDD18_db, ln_total_exp_usd_2011, urban,
                           n_members, sh_under16, ownership_d, edu_head_2, age_head, sex_head, housing_index_lab, 
                           country, curr_CDD18_db, curr_CDD18_db2, curr_HDD18_db, hhid, state, weight, mean_CDD18_db, mean_HDD18_db, curr_CDD18_db, 
                           curr_HDD18_db, fan, wshm, ref, pc, tv, ely_p, ely_p_usd_2011, ely_exp, ely_exp_usd_2011,
                           ln_ely_p_cdd, ln_ely_p_cdd2, ln_ely_p_own, ln_ely_p_nme, mean_CDD18_db2, mean_CDD18_db_exp,
                           mean_CDD18_db2_exp, ln_ely_p, urban_sh, adm1)

# Scale variable
HH_Africa <- HH_Africa %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD18_db)),
                                  std_elyp = as.numeric(scale(ln_ely_p)),
                                  std_CDD = as.numeric(scale(curr_CDD18_db)),
                                  std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                  std_HDD = as.numeric(scale(curr_HDD18_db)),
                                  std_urban_sh = as.numeric(scale(urban_sh)),
                                  std_n_members = as.numeric(scale(n_members)),
                                  std_age_head = as.numeric(scale(age_head)),
                                  std_sh_under16 = as.numeric(scale(sh_under16)))

# Survey
HH_Africa_svy <- svydesign(data = HH_Africa, ids = ~adm1, weights = ~ weight)


##################################

#        Extensive margin        #

##################################

# AC formula for Africa
ac_formula_afr <- ac ~ std_CDD_mean + I(std_CDD_mean^2) + std_CDD_mean*std_texp + I(std_CDD_mean^2)*std_texp + std_texp +
  std_elyp + std_elyp*std_CDD_mean + std_elyp*I(std_CDD_mean^2) + std_CDD + I(std_CDD^2) + std_elyp*ownership_d + std_elyp*std_n_members + 
  std_urban_sh + std_n_members + ownership_d + edu_head_2 + std_age_head + sex_head + 
  housing_index_lab + country

# Logistic regression of AC on covariates
reg_ac <- svyglm(ac_formula_afr, design = HH_Africa_svy,
                 family = binomial(logit), na.action=na.omit); summary(reg_ac)

# Save AME results
margins <- margins(reg_ac, design = HH_Africa_svy)
ac_margins <- summary(margins)

# Predicted probabilities
HH_Africa$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Africa$ac_obs <- ifelse(HH_Africa$phat0_obs>0.5 & !is.na(HH_Africa$phat0_obs), 1 , 0)

# Selection term
HH_Africa$xb_noac = 1-HH_Africa$phat0_obs               
HH_Africa$selection = ifelse(HH_Africa$ac==1, 
                             (HH_Africa$xb_noac*log(HH_Africa$xb_noac)/HH_Africa$phat0_obs) + log(HH_Africa$phat0_obs), 
                             (HH_Africa$phat0_obs*log(HH_Africa$phat0_obs)/HH_Africa$xb_noac) + log(HH_Africa$xb_noac))

# Survey - re-run to add new variable
HH_Africa_svy <- svydesign(data = HH_Africa, ids = ~adm1, weights = ~ weight)


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

# Formula electricity expenditure for electricity expenditure
ely_formula_afr  <- ln_ely_q ~ ac + ac*std_CDD + ac*I(std_CDD^2) + std_CDD + I(std_CDD^2) + 
  std_texp + std_HDD + I(std_HDD^2) + std_elyp +
  std_urban_sh + ownership_d + std_n_members + edu_head_2 + std_age_head + sex_head + 
  housing_index_lab + country + selection

# With selection
model <- svyglm(ely_formula_afr, design = HH_Africa_svy, na.action=na.omit); summary(model)

# Marginal effect of AC
ely_margins <- summary(margins(model, design = HH_Africa_svy))

# Save the R Environment will be used for the projections
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/afr_dmcf.RData', sep=''))
