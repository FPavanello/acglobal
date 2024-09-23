
##########################################

#                 Table S16

##########################################

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
interm <- "C:/Users/Standard/Documents/Github/acglobal/interm/"
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/supplementary/'

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

# List of countries
countries <- unique.default(global$country)


## Weighted
# Prepare empty data frames to store the coefficients
coeff1_df <- data.frame()
coeff2_df <- data.frame()

# Loop over all countries
for (country_exclude in countries) {

# Drop country
global_filtered <- global %>% filter(country != country_exclude)

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global_filtered, weights = ~weight, cluster = c("adm1")); summary(fs)
gc()

# Save data set for which there are obs both in first and second stage
sec <- global_filtered[obs(fs),]

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

# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model1 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model1)

# AC effect
coeff1 <- data.frame(ac = model1[["coefficients"]][1], country = country_exclude)

# Append
coeff1_df <- rbind(coeff1_df, coeff1)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model2 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model2)

# AC effect
coeff2 <- data.frame(ac = model2[["coefficients"]][1], ac_cdd = model2[["coefficients"]][17], ac_cdd2 = model2[["coefficients"]][18], country = country_exclude)

# Append
coeff2_df <- rbind(coeff2_df, coeff2)

# Clean
gc()

}

# Summary of the coefficients
summary(coeff1_df$ac)

summary(coeff2_df$ac)
summary(coeff2_df$ac_cdd)
summary(coeff2_df$ac_cdd2)

# Results
print(xtable(t(summary(coeff1_df$ac)), digits=c(4,4,4,4,4,4,4)), include.rownames=FALSE, 
      file = paste(output,'TableS16_1.tex', sep=''))
print(xtable(t(summary(coeff2_df$ac)), digits=c(4,4,4,4,4,4,4)), include.rownames=FALSE,
      file = paste(output,'TableS16_2.tex', sep=''))
print(xtable(t(summary(coeff2_df$ac_cdd)), digits=c(4,4,4,4,4,4,4)),  include.rownames=FALSE,
      file = paste(output,'TableS16_3.tex', sep=''))
print(xtable(t(summary(coeff2_df$ac_cdd2)), digits=c(4,4,4,4,4,4,4)), include.rownames=FALSE,
      file = paste(output,'TableS16_4.tex', sep=''))


## Unweighted
# Prepare empty data frames to store the coefficients
coeff1_df <- data.frame()
coeff2_df <- data.frame()

# Loop over all countries
for (country_exclude in countries) {
  
  # Drop country
  global_filtered <- global %>% filter(country != country_exclude)
  
  # AC formula for global
  ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
    mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
    ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
    urban_sh + ownership_d + 
    n_members + edu_head_2 + age_head + sex_head | adm1
  
  # Logistic regression
  fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global_filtered, cluster = c("adm1")); summary(fs)
  gc()
  
  # Save data set for which there are obs both in first and second stage
  sec <- global_filtered[obs(fs),]
  
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
  
  # Formula electricity expenditure without interaction
  ely_formula <- ln_ely_q ~ ac + 
    curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
    urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1
  
  # With selection
  model1 <- feols(ely_formula, data = sec, cluster = c("adm1")); summary(model1)
  
  # AC effect
  coeff1 <- data.frame(ac = model1[["coefficients"]][1], country = country_exclude)
  
  # Append
  coeff1_df <- rbind(coeff1_df, coeff1)
  
  # Formula electricity expenditure for electricity expenditure
  ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
    curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
    urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1
  
  # With selection
  model2 <- feols(ely_formula, data = sec, cluster = c("adm1")); summary(model2)
  
  # AC effect
  coeff2 <- data.frame(ac = model2[["coefficients"]][1], ac_cdd = model2[["coefficients"]][17], ac_cdd2 = model2[["coefficients"]][18], country = country_exclude)
  
  # Append
  coeff2_df <- rbind(coeff2_df, coeff2)
  
  # Clean
  gc()
  
}

# Summary of the coefficients
summary(coeff1_df$ac)

summary(coeff2_df$ac)
summary(coeff2_df$ac_cdd)
summary(coeff2_df$ac_cdd2)

# Results
print(xtable(t(summary(coeff1_df$ac)), digits=c(4,4,4,4,4,4,4)), include.rownames=FALSE,
      file = paste(output,'TableS16_5.tex', sep=''))
print(xtable(t(summary(coeff2_df$ac)), digits=c(4,4,4,4,4,4,4)), include.rownames=FALSE,
      file = paste(output,'TableS16_6.tex', sep=''))
print(xtable(t(summary(coeff2_df$ac_cdd)), digits=c(4,4,4,4,4,4,4)),  include.rownames=FALSE,
      file = paste(output,'TableS16_7.tex', sep=''))
print(xtable(t(summary(coeff2_df$ac_cdd2)), digits=c(4,4,4,4,4,4,4)), include.rownames=FALSE,
      file = paste(output,'TableS16_8.tex', sep=''))
