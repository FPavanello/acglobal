
## This R-script:
##      1) exploits CDD-dry bulb 18 deg and HDD-dry bulb 18 deg
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
library(fixest)
library(marginaleffects)

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


####################################################################

#       ROBUSTNESS CHECKS (Appendix)       

#       (1) Country-level Fixed effects
#       (2) Sub-national Fixed effects
#       (3) CDD 24 and HDD 15
#       (4) Drop Electricity Prices
#       (5) Electricity Prices interact with Income Deciles
#       (6) Correction term squared
#       (7) Correction term x interaction
#       (8) Winsorizing the sample
#       (9) Trimming the sample
#       (10) Unweighted Regressions


#       ROBUSTNESS CHECKS (Supplementary information)

#       (11) In levels
#       (12) Controlling for fuel-sepcific shares and 
#            potential solar output

####################################################################



###########################################################################

#  (1) Subnational Fixed Effects (smallest ADM available in each survey)  #

###########################################################################

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | subnat

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),] # > 50k observations lost

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | subnat

# With selection
model0_sfe <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0_sfe)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | subnat

# With selection
model_sfe <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_sfe)

# Mean electricity quantity
mean_sfe <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean_sfe

# Number of countries
csfe <- length(unique.default(sec$country))

# Clean
rm(fs, sec)
gc()


###############################

#  (2) Country Fixed Effects  #

###############################

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | country

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), 
            data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),]

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model0_cfe <- feols(ely_formula, data = sec, weights = ~weight, 
                    cluster = c("adm1")); summary(model0_cfe)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_cfe <- feols(ely_formula, data = sec, weights = ~weight, 
                   cluster = c("adm1")); summary(model_cfe)

# Mean electricity quantity
mean_cfe <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean_cfe

# Number of countries
ccfe <- length(unique.default(sec$country))

# Clean
rm(fs, sec)
gc()

###################################

#  (3) CDD 24 and HDD 15 degress  #

###################################

# AC formula for global
ac_formula <- ac ~ mean_CDD_db + I(mean_CDD_db^2) +
  ln_total_exp_usd_2011*mean_CDD_db + ln_total_exp_usd_2011*I(mean_CDD_db^2) + ln_total_exp_usd_2011 + curr_CDD_db + I(curr_CDD_db^2) +  
  ln_ely_p + ln_ely_p*mean_CDD_db + ln_ely_p*I(mean_CDD_db^2) + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), 
            data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),]

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD_db + I(curr_CDD_db^2) + curr_HDD_db + I(curr_HDD_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model0_24 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0_24)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD_db + ac*I(curr_CDD_db^2) + 
  curr_CDD_db + I(curr_CDD_db^2) + curr_HDD_db + I(curr_HDD_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_24 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_24)

# Mean electricity quantity
mean_24 <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean_24

# Number of countries
c24 <- length(unique.default(sec$country))

# Clean
rm(fs, sec)
gc()

##########################################

#      (4) Drop electricity prices       #

##########################################

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p*mean_CDD_db + ln_ely_p*I(mean_CDD_db^2) + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), 
            data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),]

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model0_np <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0_np)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_np <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_np)

# Mean electricity quantity
mean_np <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean_np

# Number of countries
cnp <- length(unique.default(sec$country))

# Clean
rm(fs, sec)
gc()


##################################################

# (5) Electricity prices interacted with deciles #

##################################################

# Decile of total expenditure
setDT(global)[,dec_inc := cut(total_exp_usd_2011, breaks = quantile(total_exp_usd_2011, probs = seq(0, 1, 0.1)),
                              labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), include.lowest = TRUE)]

# AC formula for global
ac_formula <- as.numeric(as.character(ac)) ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + ln_ely_p*dec_inc + dec_inc +
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), 
            data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),]

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p*dec_inc +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model0_dec <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0_dec)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p*dec_inc +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_dec <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_dec)

# Mean electricity quantity
mean_dec <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean_dec

# Number of countries
cdec <- length(unique.default(sec$country))

# Clean
rm(fs, sec)
gc()


##################################

#  (6) Correction term squared  #

#################################

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),] # 52k observations lost

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Squared
sec$selection_sq <- sec$selection*sec$selection

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection + selection_sq | adm1

# With selection
model0_ct2 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0_ct2)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection + selection_sq | adm1

# With selection
model_ct2 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_ct2)

# Mean electricity quantity
mean_ct2 <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean_ct2

# Number of countries
cct2 <- length(unique.default(sec$country))

# Clean
rm(fs, sec)
gc()


#####################################

#  (7) Correction term interaction #

####################################

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),] # 52k observations lost

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Linear Interaction
sec$selection_int <- sec$curr_CDD18_db*sec$selection
sec$selection_intsq <- sec$selection*sec$curr_CDD18_db*sec$curr_CDD18_db

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection + selection_int + selection_intsq | adm1

# With selection
model0_ctintsq <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0_ctintsq)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection + selection_int + selection_intsq | adm1

# With selection
model_ctintsq <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_ctintsq)

# Mean electricity quantity
mean_ctintsq <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean_ctintsq

# Number of countries
cctintsq <- length(unique.default(sec$country))

# Clean
rm(fs, sec)
gc()


#####################################

#  (8) Winsorizing the sample 5-95  #

#####################################

# Load
library(DescTools)

# Winsorize based on electricity consumption, income and CDD
global$ely_qw <- DescTools::Winsorize(global$ely_q, probs = c(0.05, 0.95), na.rm = TRUE)
global$ln_ely_qw <- DescTools::Winsorize(global$ln_ely_q, probs = c(0.05, 0.95), na.rm = TRUE)
global$ln_total_exp_usd_2011w <- Winsorize(global$total_exp_usd_2011, probs = c(0.05, 0.95), na.rm = TRUE)
global$curr_CDD18_dbw <- Winsorize(global$curr_CDD18_db, probs = c(0.05, 0.95), na.rm = TRUE)
global$mean_CDD18_dbw <- Winsorize(global$mean_CDD18_db, probs = c(0.05, 0.95), na.rm = TRUE)
gc()


# AC formula for global
ac_formula <- ac ~ mean_CDD18_dbw + I(mean_CDD18_dbw^2) +
  ln_total_exp_usd_2011w*mean_CDD18_dbw + ln_total_exp_usd_2011w*I(mean_CDD18_dbw^2) + ln_total_exp_usd_2011w + curr_CDD18_dbw + I(curr_CDD18_dbw^2) +  
  ln_ely_p + ln_ely_p*mean_CDD18_dbw + ln_ely_p*I(mean_CDD18_dbw^2) + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head  | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),] # 52k observations lost

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))
gc()

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_qw ~ ac +
  curr_CDD18_dbw + I(curr_CDD18_dbw^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011w + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model0_win <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0_win)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_qw ~ ac + ac*curr_CDD18_dbw + ac*I(curr_CDD18_dbw^2) + 
  curr_CDD18_dbw + I(curr_CDD18_dbw^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011w + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_win <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_win)

# Mean electricity quantity
mean_win <- weighted.mean(exp(sec$ln_ely_qw), sec$weight)
mean_win

# Number of countries
cwin <- length(unique.default(sec$country))

# Clean
rm(fs, sec)
gc()


####################################

#  (9) Trimmering the sample 5-95  #

####################################

# Winsorize based on electricity consumption
cut_point_top <- quantile(global$ely_q, 0.95)
cut_point_bottom <- quantile(global$ely_q, 0.05)

# Filter
global_trim <- filter(global, ely_q >= cut_point_bottom & ely_q <= cut_point_top)

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global_trim, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global_trim[obs(fs),] # 52k observations lost

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))
gc()

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model0_trim <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0_trim)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_trim <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_trim)

# Mean electricity quantity
mean_trim <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean_trim

# Number of countries
ctrim <- length(unique.default(sec$country))

# Clean
rm(fs, sec)
gc()


#################################

#  (10) Unweighted regressions   #

#################################

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression of AC on covariates
fs <- feglm(ac_formula, family = binomial(link = "logit"), 
            data = global, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),] # 52k observations lost

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model0_unw <- feols(ely_formula, data = sec, cluster = c("adm1")); summary(model0_unw)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_unw <- feols(ely_formula, data = sec, cluster = c("adm1")); summary(model_unw)

# Mean electricity quantity
mean_unw <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean_unw

# Number of countries
cunw <- length(unique.default(sec$country))

# Clean
rm(fs, sec)
gc()



#################################

#  (11) Electricity in Levels  #

################################

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),] # 52k observations lost

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev0 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_lev0)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev1 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_lev1)


# Winsorize based on electricity consumption, income and CDD
global$ely_qw <- DescTools::Winsorize(global$ely_q, probs = c(0.05, 0.95), na.rm = TRUE)
global$ln_total_exp_usd_2011w <- Winsorize(global$total_exp_usd_2011, probs = c(0.05, 0.95), na.rm = TRUE)
global$curr_CDD18_dbw <- Winsorize(global$curr_CDD18_db, probs = c(0.05, 0.95), na.rm = TRUE)
global$mean_CDD18_dbw <- Winsorize(global$mean_CDD18_db, probs = c(0.05, 0.95), na.rm = TRUE)
gc()

# AC formula for global
ac_formula <- ac ~ mean_CDD18_dbw + I(mean_CDD18_dbw^2) +
  ln_total_exp_usd_2011w*mean_CDD18_dbw + ln_total_exp_usd_2011w*I(mean_CDD18_dbw^2) + ln_total_exp_usd_2011w + curr_CDD18_dbw + I(curr_CDD18_dbw^2) +  
  ln_ely_p + ln_ely_p*mean_CDD18_dbw + ln_ely_p*I(mean_CDD18_dbw^2) + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head  | adm1

# Logistic regression
fsw <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fsw)

# Save data set for which there are obs both in first and second stage
secw <- global[obs(fsw),] # 52k observations lost

# Predicted probabilities
secw$phat0_obs <- as.numeric(predict(fsw, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
secw$ac_obs <- ifelse(secw$phat0_obs>0.5 & !is.na(secw$phat0_obs), 1 , 0)

# Selection term
secw$xb_noac = 1-secw$phat0_obs               
secw$selection = ifelse(secw$ac==1, 
                       (secw$xb_noac*log(secw$xb_noac)/secw$phat0_obs) + log(secw$phat0_obs), 
                       (secw$phat0_obs*log(secw$phat0_obs)/secw$xb_noac) + log(secw$xb_noac))
gc()

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_qw ~ ac +
  curr_CDD18_dbw + I(curr_CDD18_dbw^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011w + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev2 <- feols(ely_formula, data = secw, weights = ~weight, cluster = c("adm1")); summary(model_lev2)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_qw ~ ac + ac*curr_CDD18_dbw + ac*I(curr_CDD18_dbw^2) + 
  curr_CDD18_dbw + I(curr_CDD18_dbw^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011w + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev3 <- feols(ely_formula, data = secw, weights = ~weight, cluster = c("adm1")); summary(model_lev3)


# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fst <- feglm(ac_formula, family = binomial(link = "logit"), data = global_trim, weights = ~weight, cluster = c("adm1")); summary(fst)

# Save data set for which there are obs both in first and second stage
sect <- global_trim[obs(fst),] # 52k observations lost

# Predicted probabilities
sect$phat0_obs <- as.numeric(predict(fst, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sect$ac_obs <- ifelse(sect$phat0_obs>0.5 & !is.na(sect$phat0_obs), 1 , 0)

# Selection term
sect$xb_noac = 1-sect$phat0_obs               
sect$selection = ifelse(sect$ac==1, 
                       (sect$xb_noac*log(sect$xb_noac)/sect$phat0_obs) + log(sect$phat0_obs), 
                       (sect$phat0_obs*log(sect$phat0_obs)/sect$xb_noac) + log(sect$xb_noac))
gc()

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev4 <- feols(ely_formula, data = sect, weights = ~weight, cluster = c("adm1")); summary(model_lev4)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev5 <- feols(ely_formula, data = sect, weights = ~weight, cluster = c("adm1")); summary(model_lev5)


# Mean electricity quantity
mean_lev1 <- weighted.mean(sec$ely_q, sec$weight)
mean_lev1
mean_lev2 <- weighted.mean(secw$ely_qw, secw$weight)
mean_lev2
mean_lev3 <- weighted.mean(sect$ely_qw, sect$weight)
mean_lev3

# Number of countries
cwin1 <- length(unique.default(sec$country))
cwin2 <- length(unique.default(secw$country))
cwin3 <- length(unique.default(sect$country))

# Clean
rm(fs, sec, fsw, secw, global_trim, sect, fst)
gc()


##########################################################################

#  (12) Controlling for fuel-sepcific shares and potential solar output  #

##########################################################################

# Shares
ivshare <- global %>% dplyr::group_by(country, subnat) %>% dplyr::summarise(cap_solar = mean(cap_solar, na.rm = TRUE), cap_wind = mean(cap_wind, na.rm = TRUE), 
                                                                            cap_hydro = mean(cap_hydro, na.rm = TRUE), cap_nuclear  = mean(cap_nuclear, na.rm = TRUE), 
                                                                            cap_otherres = mean(cap_otherres, na.rm = TRUE), cap_coal = mean(cap_coal, na.rm = TRUE),
                                                                            cap_gas = mean(cap_gas, na.rm = TRUE), cap_oil = mean(cap_oil, na.rm = TRUE), cap_other = mean(cap_other, na.rm = TRUE))
ivshare <- ivshare %>% mutate(cap_tot = cap_solar + cap_wind + cap_hydro + cap_nuclear + cap_otherres + cap_coal + cap_gas + cap_oil + cap_other)
ivshare <- ivshare %>% mutate(sh_solar = ifelse(cap_solar == 0, 0, cap_solar/cap_tot),
                              sh_wind = ifelse(cap_wind == 0, 0, cap_wind/cap_tot),
                              sh_hydro = ifelse(cap_hydro == 0, 0, cap_hydro/cap_tot),
                              sh_nuclear = ifelse(cap_nuclear == 0, 0, cap_nuclear/cap_tot),
                              sh_otherres = ifelse(cap_otherres == 0, 0, cap_otherres/cap_tot),
                              sh_coal = ifelse(cap_coal == 0, 0, cap_coal/cap_tot),
                              sh_gas = ifelse(cap_gas == 0, 0, cap_gas/cap_tot),
                              sh_oil = ifelse(cap_oil == 0, 0, cap_oil/cap_tot),
                              sh_other = ifelse(cap_other == 0, 0, cap_other/cap_tot),
                              sh_ren = sh_solar + sh_wind + sh_otherres)
ivshare <- ivshare %>% dplyr::select(-c(cap_solar, cap_wind, cap_nuclear, cap_otherres, cap_hydro, cap_coal, cap_gas, cap_oil, cap_other))
gc()

# Merge
global <- merge(global, ivshare, by = c("country", "subnat"))

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head + sh_coal + sh_hydro + sh_gas + sh_ren | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),] # 52k observations lost

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection + sh_coal + sh_hydro + sh_gas + sh_ren | adm1

# With selection
model_cont0 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_cont0)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection + sh_coal + sh_hydro + sh_gas + sh_ren | adm1

# With selection
model_cont1 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_cont1)


## Export
# Compare the models
screenreg(list(model0_sfe, model_sfe, model0_cfe, model_cfe, model0_24, model_24, model0_np, model_np, 
               model0_dec, model_dec), digits = 3, 
          caption = "Robustness Checks",
          stars = c(0.1, 0.05, 0.01), 
          custom.model.names = c("Subnational FE", "Subnational FE", "Country FE", "Country FE", 
                                 "CDD 24 - HDD 15", "CDD 24 - HDD 15", "No Elec. Price", "No Elec. Price", 
                                 "Price Interactions", "Price Interactions"),
          custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                                 "ac:curr_CDD_db" = "AC $\\times$ CDD", 
                                 "ac:I(curr_CDD_db^2)" = "AC $\\times$ CDD$^2$"),
          custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "Country FE" = c("NO", "NO", "NO", "YES", "YES", "NO", "NO", "NO", "NO", "NO"), 
                                 "Subnational FE" = c("YES", "YES", "NO", "NO", "NO", "NO", "NO", "NO", "NO", "NO"), 
                                 "ADM-1 FE" = c("NO", "NO", "NO", "NO", "YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "Mean Outcome (kWh)" = c(mean_sfe, mean_sfe, mean_cfe, mean_cfe, mean_24, mean_24, mean_np, mean_np, mean_dec, mean_dec), 
                                 "Countries" = c(csfe, csfe, ccfe, ccfe, c24, c24, cnp, cnp, cdec, cdec)))

screenreg(list(model0_ct2, model_ct2, model0_ctintsq, model_ctintsq, model0_win, model_win, model0_trim, model_trim, 
               model0_unw, model_unw), digits = 3, 
          caption = "Robustness Checks",
          stars = c(0.1, 0.05, 0.01), 
          custom.model.names = c("Squared CT", "Squared CT", "Interaction CT with CDD", "Interaction CT with CDD", 
                                 "Winsorized", "Winsorized", "Trimmed", "Trimmed", 
                                 "Unweighted", "Unweighted"),
          custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", 
                                 "ac:I(curr_CDD18_dbw^2)" = "AC $\\times$ CDD$^2$"),
          custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "ADM-1 FE" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "Correction Term$^2$" = c("YES", "YES", "NO", "NO", "NO", "NO", "NO", "NO", "NO", "NO"), 
                                 "Correction Term $\\times$ f(CDD)" = c("NO", "NO", "NO", "NO", "YES", "YES", "NO", "NO", "NO", "NO"), 
                                 "Mean Outcome (kWh)" = c(mean_ct2, mean_ct2, mean_ctintsq, mean_ctintsq, mean_win, mean_win, mean_trim, mean_trim, mean_unw, mean_unw), 
                                 "Countries" = c(cct2, cct2, cctintsq, cctintsq, cwin, cwin, ctrim, ctrim, cunw, cunw)))

screenreg(list(model_lev0, model_lev1, model_lev2, model_lev3, model_lev4, model_lev5), digits = 3, 
          caption = "Robustness Checks - Electricity in Levels",
          stars = c(0.1, 0.05, 0.01), 
          custom.model.names = c("Full", "Full", 
                                 "Winsorized", "Winsorized", "Trimmed", "Trimmed"),
          custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac:curr_CDD18_dbw" = "AC $\\times$ CDD", 
                                 "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", 
                                 "ac:I(curr_CDD18_dbw^2)" = "AC $\\times$ CDD$^2$"),
          custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "ADM-1 FE" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "Mean Outcome (kWh)" = c(mean_lev1, mean_lev1, mean_lev2, mean_lev2, mean_lev3, mean_lev3), 
                                 "Countries" = c(cwin1, cwin1, cwin2, cwin2, cwin3, cwin3)))


# Export
texreg(list(model0_sfe, model_sfe, model0_cfe, model_cfe, model0_24, model_24, model0_np, model_np, 
            model0_dec, model_dec), digits = 3, 
       caption = "Robustness Checks",
       stars = c(0.1, 0.05, 0.01), 
       custom.model.names = c("Subnational FE", "Subnational FE", "Country FE", "Country FE", 
                              "CDD 24 - HDD 15", "CDD 24 - HDD 15", "No Elec. Price", "No Elec. Price", 
                              "Price Interactions", "Price Interactions"),
       custom.note = "'Subnational' means at the most disaggregated geographical information for each country. Regressions are conducted using survey weights.
       Standard errors are clustered at the ADM1 level. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$.", 
       file = paste(output,'electricity/robustness/Global_robustness_part1.tex', sep=''), append=F,  float.pos = "H", 
       label = "rob: ely_global",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "ac:curr_CDD_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD_db^2)" = "AC $\\times$ CDD$^2$"),
       custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Country FE" = c("NO", "NO", "NO", "YES", "YES", "NO", "NO", "NO", "NO", "NO"), 
                              "Subnational FE" = c("YES", "YES", "NO", "NO", "NO", "NO", "NO", "NO", "NO", "NO"), 
                              "ADM-1 FE" = c("NO", "NO", "NO", "NO", "YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Mean Outcome (kWh)" = c(mean_sfe, mean_sfe, mean_cfe, mean_cfe, mean_24, mean_24, mean_np, mean_np, mean_dec, mean_dec), 
                              "Countries" = c(csfe, csfe, ccfe, ccfe, c24, c24, cnp, cnp, cdec, cdec)), 
       caption.above = TRUE)

texreg(list(model0_ct2, model_ct2, model0_ctintsq, model_ctintsq, model0_win, model_win, model0_trim, model_trim, 
            model0_unw, model_unw), digits = 3, 
       caption = "Robustness Checks",
       stars = c(0.1, 0.05, 0.01), 
       custom.model.names = c("Squared CT", "Squared CT", "Interaction CT with CDD", "Interaction CT with CDD", 
                              "Winsorized", "Winsorized", "Trimmed", "Trimmed", 
                              "Unweighted", "Unweighted"),
       custom.note = "'Subnational' means at the most disaggregated geographical information for each country. Regressions are conducted using survey weights.
       Standard errors are clustered at the ADM1 level. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$.", 
       file = paste(output,'electricity/robustness/Global_robustness_part2.tex', sep=''), append=F,  float.pos = "H", label = "rob: ely_global",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "ac:curr_CDD18_dbw" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_dbw^2)" = "AC $\\times$ CDD$^2$"),
       custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                              "ADM-1 FE" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Correction Term$^2$" = c("YES", "YES", "NO", "NO", "NO", "NO", "NO", "NO", "NO", "NO"), 
                              "Correction Term $\\times$ f(CDD)" = c("NO", "NO", "NO", "NO", "YES", "YES", "NO", "NO", "NO", "NO"), 
                              "Mean Outcome (kWh)" = c(mean_ct2, mean_ct2, mean_ctintsq, mean_ctintsq, mean_win, mean_win, mean_trim, mean_trim, mean_unw, mean_unw), 
                              "Countries" = c(cct2, cct2, cctintsq, cctintsq, cwin, cwin, ctrim, ctrim, cunw, cunw)), 
       caption.above = TRUE)

texreg(list(model_lev0, model_lev1, model_lev2, model_lev3, model_lev4, model_lev5), digits = 3, 
          caption = "Robustness Checks - Electricity in Levels",
          stars = c(0.1, 0.05, 0.01), 
          custom.model.names = c("Full", "Full", 
                                 "Winsorized", "Winsorized", "Trimmed", "Trimmed"),
          custom.note = "'Subnational' means at the most disaggregated geographical information for each country. Regressions are conducted using survey weights.
          Standard errors are clustered at the ADM1 level. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$.", 
          file = paste(output,'electricity/robustness/Global_robustness_level.tex', sep=''), append=F,  float.pos = "H", label = "rob: ely_global_level",
          custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac:curr_CDD18_dbw" = "AC $\\times$ CDD",
                                 "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", 
                                 "ac:I(curr_CDD18_dbw^2)" = "AC $\\times$ CDD$^2$"),
          custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "ADM-1 FE" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "Mean Outcome (kWh)" = c(mean_lev1, mean_lev1, mean_lev2, mean_lev2, mean_lev3, mean_lev3), 
                                 "Countries" = c(cwin1, cwin1, cwin2, cwin2, cwin3, cwin3)), 
          caption.above = TRUE)
