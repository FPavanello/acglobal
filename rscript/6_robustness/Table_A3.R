
#############################################################

#                          Table A3

#############################################################

# Free memory
.rs.restartR()
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
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'data/household/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/'
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/tables/'


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


###########################################################################

#  (1) Subnational Fixed Effects (smallest ADM available in each survey)  #

###########################################################################

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  curr_HDD18_db + I(curr_HDD18_db^2) +
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
  curr_HDD18_db + I(curr_HDD18_db^2) +
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
  curr_HDD_db + I(curr_HDD_db^2) +
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
  curr_HDD18_db + I(curr_HDD18_db^2) +
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
  curr_HDD18_db + I(curr_HDD18_db^2) +
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
  curr_HDD18_db + I(curr_HDD18_db^2) +
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
  curr_HDD18_db + I(curr_HDD18_db^2) +
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
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_ely_p + ln_ely_p*mean_CDD18_dbw + ln_ely_p*I(mean_CDD18_dbw^2) + ln_ely_p_nme + ln_ely_p_own + 
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
  curr_HDD18_db + I(curr_HDD18_db^2) +
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
  curr_HDD18_db + I(curr_HDD18_db^2) +
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
       file = paste(output,'TableA3_1.tex', sep=''), append=F,  float.pos = "H", 
       label = "main: tableA3",
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
       file = paste(output,'TableA3_2.tex', sep=''), append=F,  float.pos = "H", label = "main: tableA3",
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
