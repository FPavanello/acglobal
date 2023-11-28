
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts logit regressions for the global data set
##      3) run intensive margin regressions: electricity expenditure on climate + covariates
##         using Dubin and McFadden (1984) approach
##      4) ROBUSTNESS CHECKS

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
library(effects)
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

# Load global data
global <- readRDS(paste(house,'global.rds', sep=''))
global_housing <- readRDS(paste(house,'further_data/global_housing.rds', sep=''))
global_italy <- readRDS(paste(house,'further_data/global_italy.rds', sep=''))

# Check
global <- global[complete.cases(global$ln_ely_q), ]
global <- global[complete.cases(global$ac), ]
global <- global[complete.cases(global$ln_total_exp_usd_2011), ]
global <- global[complete.cases(global$mean_CDD_db), ]
global <- global[complete.cases(global$urban_sh), ]
global <- global[complete.cases(global$ownership_d), ]
global <- global[complete.cases(global$n_members), ]
global <- global[complete.cases(global$age_head), ]
global <- global[complete.cases(global$country), ]
global <- global[complete.cases(global$weight), ]

# Check
global_housing <- global_housing[complete.cases(global_housing$ln_ely_q), ]
global_housing <- global_housing[complete.cases(global_housing$ac), ]
global_housing <- global_housing[complete.cases(global_housing$ln_total_exp_usd_2011), ]
global_housing <- global_housing[complete.cases(global_housing$mean_CDD_db), ]
global_housing <- global_housing[complete.cases(global_housing$urban_sh), ]
global_housing <- global_housing[complete.cases(global_housing$ownership_d), ]
global_housing <- global_housing[complete.cases(global_housing$n_members), ]
global_housing <- global_housing[complete.cases(global_housing$age_head), ]
global_housing <- global_housing[complete.cases(global_housing$country), ]
global_housing <- global_housing[complete.cases(global_housing$weight), ]

# Check
global_italy <- global_italy[complete.cases(global_italy$ln_ely_q), ]
global_italy <- global_italy[complete.cases(global_italy$ac), ]
global_italy <- global_italy[complete.cases(global_italy$ln_total_exp_usd_2011), ]
global_italy <- global_italy[complete.cases(global_italy$mean_CDD_db), ]
global_italy <- global_italy[complete.cases(global_italy$n_members), ]
global_italy <- global_italy[complete.cases(global_italy$age_head), ]
global_italy <- global_italy[complete.cases(global_italy$country), ]
global_italy <- global_italy[complete.cases(global_italy$urban_sh), ]


####################################################################

#                       ROBUSTNESS CHECKS        

#       (1) Add Italy
#       (2) Drop Germany
#       (3) Include only countries with housing information
#       (4) Only climate and income variables
#       (5) Unweighted regressions

####################################################################

###################

#  (1) Add Italy  #

###################

# AC formula for global
ac_formula <- ac ~ mean_CDD_db + I(mean_CDD_db^2) + mean_CDD_db*ln_total_exp_usd_2011 + 
  I(mean_CDD_db^2)*ln_total_exp_usd_2011 + n_members + edu_head_2 + age_head + country

# Logistic regression of AC on covariates
reg_ac <- glm(ac_formula, data = global_italy, family = binomial(logit), na.action=na.omit); summary(reg_ac)
vcov <- vcovCL(reg_ac, cluster = global_italy$hhid) # clusterise SEs
coeftest(reg_ac, vcov = vcovCL(reg_ac, cluster = global_italy$hhid)) # clustered standard errors

# Predicted probabilities
global_italy$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
global_italy$ac_obs <- ifelse(global_italy$phat0_obs>0.5 & !is.na(global_italy$phat0_obs), 1 , 0)

table(global_italy$ac, global_italy$ac_obs)

# Selection term
global_italy$xb_noac = 1-global_italy$phat0_obs               
global_italy$selection = ifelse(global_italy$ac==1, 
                             (global_italy$xb_noac*log(global_italy$xb_noac)/global_italy$phat0_obs) + log(global_italy$phat0_obs), 
                             (global_italy$phat0_obs*log(global_italy$phat0_obs)/global_italy$xb_noac) + log(global_italy$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD_db + I(curr_CDD_db^2) + curr_CDD_db*ln_total_exp_usd_2011 + I(curr_CDD_db^2)*ln_total_exp_usd_2011 + 
  curr_HDD_db + I(curr_HDD_db^2) + curr_HDD_db*ln_total_exp_usd_2011 + I(curr_HDD_db^2)*ln_total_exp_usd_2011 + 
  urban_sh + n_members + edu_head_2 + age_head + country + selection

# With selection
model01 <- lm(ely_formula, data = global_italy, na.action=na.omit); summary(model01)
vcov01 <- cluster.vcov(model01, global_italy$hhid) # clusterise SEs
coeftest(model01, vcov=vcov01)
model01_cl <- coeftest(model01, vcov=vcov01) # with clustered SEs

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD_db + ac*I(curr_CDD_db^2) + 
                           curr_CDD_db + I(curr_CDD_db^2) + curr_CDD_db*ln_total_exp_usd_2011 + I(curr_CDD_db^2)*ln_total_exp_usd_2011 + 
                           curr_HDD_db + I(curr_HDD_db^2) + curr_HDD_db*ln_total_exp_usd_2011 + I(curr_HDD_db^2)*ln_total_exp_usd_2011 + 
                           urban_sh + n_members + edu_head_2 + age_head + country + selection

# With selection
model1 <- lm(ely_formula, data = global_italy, na.action=na.omit); summary(model1)
vcov1 <- cluster.vcov(model1, global_italy$hhid) # clusterise SEs
coeftest(model1, vcov=vcov1)
model1_cl <- coeftest(model1, vcov=vcov1) # with clustered SEs

# Marginal effect of AC
#ac_eff <- margins(model1, variables = "ac", vcov = vcov1)
#summary(ac_eff)


######################

#  (2) Drop Germany  #

######################

# Drop Germany
global_nodeu <- global %>% filter(country != "Germany")

# Survey
global_nodeu_svy <- svydesign(data = global_nodeu, ids = ~1, weights = ~ weight)

# AC formula for global
ac_formula <- ac ~ mean_CDD_db + I(mean_CDD_db^2) + mean_CDD_db*ln_total_exp_usd_2011 + 
  I(mean_CDD_db^2)*ln_total_exp_usd_2011 + urban_sh + ownership_d + n_members + edu_head_2 + age_head + country

# Logistic regression of AC on covariates
reg_ac <- svyglm(ac_formula, design = global_nodeu_svy, family = binomial(logit), na.action=na.omit); summary(reg_ac)

# Predicted probabilities
global_nodeu$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
global_nodeu$ac_obs <- ifelse(global_nodeu$phat0_obs>0.5 & !is.na(global_nodeu$phat0_obs), 1 , 0)

# Selection term
global_nodeu$xb_noac = 1-global_nodeu$phat0_obs               
global_nodeu$selection = ifelse(global_nodeu$ac==1, 
                          (global_nodeu$xb_noac*log(global_nodeu$xb_noac)/global_nodeu$phat0_obs) + log(global_nodeu$phat0_obs), 
                          (global_nodeu$phat0_obs*log(global_nodeu$phat0_obs)/global_nodeu$xb_noac) + log(global_nodeu$xb_noac))

# Survey - re-run to add new variable
global_nodeu_svy <- svydesign(data = global_nodeu, ids = ~1, weights = ~ weight)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + 
  curr_CDD_db + I(curr_CDD_db^2) + curr_CDD_db*ln_total_exp_usd_2011 + I(curr_CDD_db^2)*ln_total_exp_usd_2011 + 
  curr_HDD_db + I(curr_HDD_db^2) + curr_HDD_db*ln_total_exp_usd_2011 + I(curr_HDD_db^2)*ln_total_exp_usd_2011 + 
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + country + selection

# With selection
model02_cl <- svyglm(ely_formula, design = global_nodeu_svy, na.action=na.omit); summary(model02_cl)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD_db + ac*I(curr_CDD_db^2) + 
  curr_CDD_db + I(curr_CDD_db^2) + curr_CDD_db*ln_total_exp_usd_2011 + I(curr_CDD_db^2)*ln_total_exp_usd_2011 + 
  curr_HDD_db + I(curr_HDD_db^2) + curr_HDD_db*ln_total_exp_usd_2011 + I(curr_HDD_db^2)*ln_total_exp_usd_2011 + 
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + country + selection

# With selection
model2_cl <- svyglm(ely_formula, design = global_nodeu_svy, na.action=na.omit); summary(model2_cl)

# Marginal effect of AC
#ac_eff <- margins(model2_cl, variables = "ac", design = global_nodeu_svy)
#summary(ac_eff)


##########################################

#  (3) Only countries with housing info  #

##########################################

# Survey
global_housing_svy <- svydesign(data = global_housing, ids = ~1, weights = ~ weight)

# AC formula for global
ac_formula <- ac ~ mean_CDD_db + I(mean_CDD_db^2) + mean_CDD_db*ln_total_exp_usd_2011 + 
                   I(mean_CDD_db^2)*ln_total_exp_usd_2011 + urban_sh + ownership_d + n_members + 
                   edu_head_2 + age_head + housing_index_lab + country

# Logistic regression of AC on covariates
reg_ac <- svyglm(ac_formula, design = global_housing_svy, family = binomial(logit), na.action=na.omit); summary(reg_ac)

# Predicted probabilities
global_housing$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
global_housing$ac_obs <- ifelse(global_housing$phat0_obs>0.5 & !is.na(global_housing$phat0_obs), 1 , 0)

# Selection term
global_housing$xb_noac = 1-global_housing$phat0_obs               
global_housing$selection = ifelse(global_housing$ac==1, 
                                (global_housing$xb_noac*log(global_housing$xb_noac)/global_housing$phat0_obs) + log(global_housing$phat0_obs), 
                                (global_housing$phat0_obs*log(global_housing$phat0_obs)/global_housing$xb_noac) + log(global_housing$xb_noac))

# Survey
global_housing_svy <- svydesign(data = global_housing, ids = ~1, weights = ~ weight)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + 
  curr_CDD_db + I(curr_CDD_db^2) + curr_CDD_db*ln_total_exp_usd_2011 + I(curr_CDD_db^2)*ln_total_exp_usd_2011 + 
  curr_HDD_db + I(curr_HDD_db^2) + curr_HDD_db*ln_total_exp_usd_2011 + I(curr_HDD_db^2)*ln_total_exp_usd_2011 + 
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + housing_index_lab + country + selection

# With selection
model03_cl <- svyglm(ely_formula, design = global_housing_svy, na.action=na.omit); summary(model03_cl)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD_db + ac*I(curr_CDD_db^2) + 
  curr_CDD_db + I(curr_CDD_db^2) + curr_CDD_db*ln_total_exp_usd_2011 + I(curr_CDD_db^2)*ln_total_exp_usd_2011 + 
  curr_HDD_db + I(curr_HDD_db^2) + curr_HDD_db*ln_total_exp_usd_2011 + I(curr_HDD_db^2)*ln_total_exp_usd_2011 + 
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + housing_index_lab + country + selection

# With selection
model3_cl <- svyglm(ely_formula, design = global_housing_svy, na.action=na.omit); summary(model3_cl)

# Marginal effect of AC
#ac_eff <- margins(model3_cl, variables = "ac", design = global_housing_svy)
#summary(ac_eff)


#################################

#  (4) Only income and climate  #

#################################

# Dataframe
global_ci <- global

# Survey
global_ci_svy <- svydesign(data = global_ci, ids = ~1, weights = ~ weight)

# AC formula for global
ac_formula <- ac ~ mean_CDD_db + I(mean_CDD_db^2) + mean_CDD_db*ln_total_exp_usd_2011 + 
                   I(mean_CDD_db^2)*ln_total_exp_usd_2011 + country

# Logistic regression of AC on covariates
reg_ac <- svyglm(ac_formula, design = global_ci_svy, family = binomial(logit), na.action=na.omit); summary(reg_ac)

# Predicted probabilities
global_ci$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
global_ci$ac_obs <- ifelse(global_ci$phat0_obs>0.5 & !is.na(global_ci$phat0_obs), 1 , 0)

# Selection term
global_ci$xb_noac = 1-global_ci$phat0_obs               
global_ci$selection = ifelse(global_ci$ac==1, 
                                (global_ci$xb_noac*log(global_ci$xb_noac)/global_ci$phat0_obs) + log(global_ci$phat0_obs), 
                                (global_ci$phat0_obs*log(global_ci$phat0_obs)/global_ci$xb_noac) + log(global_ci$xb_noac))

# Survey
global_ci_svy <- svydesign(data = global_ci, ids = ~1, weights = ~ weight)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD_db + I(curr_CDD_db^2) + curr_CDD_db*ln_total_exp_usd_2011 + I(curr_CDD_db^2)*ln_total_exp_usd_2011 + 
  curr_HDD_db + I(curr_HDD_db^2) + curr_HDD_db*ln_total_exp_usd_2011 + I(curr_HDD_db^2)*ln_total_exp_usd_2011 + 
  country + selection

# With selection
model04_cl <- svyglm(ely_formula, design = global_ci_svy, na.action=na.omit); summary(model04_cl)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD_db + ac*I(curr_CDD_db^2) + 
  curr_CDD_db + I(curr_CDD_db^2) + curr_CDD_db*ln_total_exp_usd_2011 + I(curr_CDD_db^2)*ln_total_exp_usd_2011 + 
  curr_HDD_db + I(curr_HDD_db^2) + curr_HDD_db*ln_total_exp_usd_2011 + I(curr_HDD_db^2)*ln_total_exp_usd_2011 + 
  country + selection

# With selection
model4_cl <- svyglm(ely_formula, design = global_ci_svy, na.action=na.omit); summary(model4_cl)

# Marginal effect of AC
#ac_eff <- margins(model4_cl, variables = "ac", design = global_ci_svy)
#summary(ac_eff)


#################################

#  (5) Unweighted regressions   #

#################################

# AC formula for global
ac_formula <- ac ~ mean_CDD_db + I(mean_CDD_db^2) + mean_CDD_db*ln_total_exp_usd_2011 + 
  I(mean_CDD_db^2)*ln_total_exp_usd_2011 + urban_sh + ownership_d + n_members + edu_head_2 + age_head + country

# Logistic regression of AC on covariates
reg_ac <- glm(ac_formula, data = global, family = binomial(logit), na.action=na.omit); summary(reg_ac)

# Predicted probabilities
global$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
global$ac_obs <- ifelse(global$phat0_obs>0.5 & !is.na(global$phat0_obs), 1 , 0)

# Selection term
global$xb_noac = 1-global$phat0_obs               
global$selection = ifelse(global$ac==1, 
                                (global$xb_noac*log(global$xb_noac)/global$phat0_obs) + log(global$phat0_obs), 
                                (global$phat0_obs*log(global$phat0_obs)/global$xb_noac) + log(global$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD_db + I(curr_CDD_db^2) + curr_CDD_db*ln_total_exp_usd_2011 + I(curr_CDD_db^2)*ln_total_exp_usd_2011 + 
  curr_HDD_db + I(curr_HDD_db^2) + curr_HDD_db*ln_total_exp_usd_2011 + I(curr_HDD_db^2)*ln_total_exp_usd_2011 +  
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + country + selection

# With selection
model05 <- lm(ely_formula, data = global, na.action=na.omit); summary(model05)
vcov05 <- cluster.vcov(model05, global$hhid) # clusterise SEs
coeftest(model05, vcov=vcov05)
model05_cl <- coeftest(model05, vcov=vcov05) # with clustered SEs

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD_db + ac*I(curr_CDD_db^2) + 
  curr_CDD_db + I(curr_CDD_db^2) + curr_CDD_db*ln_total_exp_usd_2011 + I(curr_CDD_db^2)*ln_total_exp_usd_2011 + 
  curr_HDD_db + I(curr_HDD_db^2) + curr_HDD_db*ln_total_exp_usd_2011 + I(curr_HDD_db^2)*ln_total_exp_usd_2011 +  
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + country + selection

# With selection
model5 <- lm(ely_formula, data = global, na.action=na.omit); summary(model5)
vcov5 <- cluster.vcov(model5, global$hhid) # clusterise SEs
coeftest(model5, vcov=vcov5)
model5_cl <- coeftest(model5, vcov=vcov5) # with clustered SEs

# Marginal effect of AC
#ac_eff <- margins(model5, variables = "ac", vcov = vcov5)
#summary(ac_eff)


## Export
# Compare the models
screenreg(list(model01_cl, model1_cl, model02_cl, model2_cl, model03_cl, model3_cl, 
               model04_cl, model4_cl, model05_cl, model5_cl), digits = 3, 
          caption = "Robustness Checks",
          stars = c(0.1, 0.05, 0.01), 
          custom.model.names = c("Include Italy", "Include Italy", 
                                 "Drop Germany", "Drop Germany", "Housing", "Housing", 
                                 "Expenditure-Climate", "Expenditure-Climate", 
                                 "Unweighted", "Unweighted"),
          custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD_db" = "AC $\\times$ CDD", 
                                 "ac1:I(curr_CDD_db^2)" = "AC $\\times$ CDD$^2$"),
          custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "Country FE" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "Countries" = c("22", "22", "20", "20", "10", "10", "21", "21", "21", "21")))

# Export
texreg(list(model01_cl, model1_cl, model02_cl, model2_cl, model03_cl, model3_cl, 
            model04_cl, model4_cl, model05_cl, model5_cl), digits = 3, 
       caption = "Robustness Checks",
       stars = c(0.1, 0.05, 0.01), 
       custom.model.names = c("Include Italy", "Include Italy", 
                              "Drop Germany", "Drop Germany", "Housing", "Housing", 
                              "Expenditure-Climate", "Expenditure-Climate", 
                              "Unweighted", "Unweighted"),
       custom.note = "Robust std. errors at the district level in parentheses. 
       $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$. Regressions are conducted using survey weights.", 
       file = paste(output,'electricity/robustness/Global_robustness.tex', sep=''), append=F,  float.pos = "htbp", label = "rob: ely_global",
       custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD_db" = "AC $\\times$ CDD", 
                              "ac1:I(curr_CDD_db^2)" = "AC $\\times$ CDD$^2$"),
       custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Country FE" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Countries" = c("22", "22", "20", "20", "10", "10", "21", "21", "21", "21")), caption.above = TRUE)

