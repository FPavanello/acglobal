
################################################

#    Table A6, Table A7, Table A8, Table A9

################################################

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
interm <- paste(stub,'results/regressions/for_graphs/subsamples/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/'
output <- paste(stub,'output/tables/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/tables/'

# Load global data
global <- readRDS(paste(house,'global.rds', sep=''))

# Interaction prices
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

# Appliances
appl <- global %>% group_by(country) %>% summarise(ac = mean(ac, na.rm = TRUE),
                                                   ref = mean(ref, na.rm = TRUE),
                                                   wshm = mean(wshm, na.rm = TRUE),
                                                   tv = mean(tv, na.rm = TRUE),
                                                   pc = mean(pc, na.rm = TRUE))

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | adm1

# REF formula for global
ref_formula <- ref ~ mean_CDD18_db + mean_CDD18_db2 + mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | adm1

# TV formula for global
tv_formula <- tv ~ mean_CDD18_db + mean_CDD18_db2 + mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | adm1

# PC formula for global
pc_formula <- pc ~ mean_CDD18_db + mean_CDD18_db2 + mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | adm1

# Washing machine formula for global
wshm_formula <- wshm ~ mean_CDD18_db + mean_CDD18_db2 + mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | adm1


########################################

#        Including refrigerators       #

########################################

# Logistic regression
fs_ac <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs_ac)

# Save data set for which there are obs both in first and second stage
sec_ref <- global[obs(fs_ac),]

# Predicted probabilities
sec_ref$phat0_ac <- as.numeric(predict(fs_ac, type="response"))

# Logistic regression
fs_ref <- feglm(ref_formula, family = binomial(link = "logit"), data = sec_ref, weights = ~weight, cluster = c("adm1")); summary(fs_ref)

# Save data set for which there are obs both in first and second stage
sec_ref <- sec_ref[obs(fs_ref),]

# Predicted probabilities
sec_ref$phat0_ref <- as.numeric(predict(fs_ref, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec_ref$ac_obs <- ifelse(sec_ref$phat0_ac>0.5 & !is.na(sec_ref$phat0_ac), 1 , 0)
sec_ref$ref_obs <- ifelse(sec_ref$phat0_ref>0.5 & !is.na(sec_ref$phat0_ref), 1 , 0)

# Selection term
sec_ref$xb_noac = 1-sec_ref$phat0_ac               
sec_ref$selectionac = ifelse(sec_ref$ac==1, 
                       (sec_ref$xb_noac*log(sec_ref$xb_noac)/sec_ref$phat0_ac) + log(sec_ref$phat0_ac), 
                       (sec_ref$phat0_ac*log(sec_ref$phat0_ac)/sec_ref$xb_noac) + log(sec_ref$xb_noac))

sec_ref$xb_noref = 1-sec_ref$phat0_ref              
sec_ref$selectionref = ifelse(sec_ref$ref==1, 
                           (sec_ref$xb_noref*log(sec_ref$xb_noref)/sec_ref$phat0_ref) + log(sec_ref$phat0_ref), 
                           (sec_ref$phat0_ref*log(sec_ref$phat0_ref)/sec_ref$xb_noref) + log(sec_ref$xb_noref))


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac | adm1

# With selection
model_ref1 <- feols(ely_formula, data = sec_ref, weights = ~weight, cluster = c("adm1")); summary(model_ref1)


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + ref +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionref | adm1

# With selection
model_ref2 <- feols(ely_formula, data = sec_ref, weights = ~weight, cluster = c("adm1")); summary(model_ref2)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + ref + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionref | adm1

# With selection
model_ref3 <- feols(ely_formula, data = sec_ref, weights = ~weight, cluster = c("adm1")); summary(model_ref3)


# Testing whether the impact of ref changes with temperature
# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ref + ref*curr_CDD18_db + ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionref | adm1

# With selection
model_ref4 <- feols(ely_formula, data = sec_ref, weights = ~weight, cluster = c("adm1")); summary(model_ref4)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ref + ref*curr_CDD18_db + ref*I(curr_CDD18_db^2) + ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionref | adm1

# With selection
model_ref5 <- feols(ely_formula, data = sec_ref, weights = ~weight, cluster = c("adm1")); summary(model_ref5)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ref + ref*curr_CDD18_db + ref*I(curr_CDD18_db^2) + ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionref | adm1

# With selection
model_ref6 <- feols(ely_formula, data = sec_ref, weights = ~weight, cluster = c("adm1")); summary(model_ref6)

# Clean
rm(fs_ac, fs_ref)
gc()



########################################

#        Including televisions       #

########################################

# Logistic regression
fs_ac <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs_ac)

# Save data set for which there are obs both in first and second stage
sec_tv <- global[obs(fs_ac),]

# Predicted probabilities
sec_tv$phat0_ac <- as.numeric(predict(fs_ac, type="response"))

# Logistic regression
fs_tv <- feglm(tv_formula, family = binomial(link = "logit"), data = sec_tv, weights = ~weight, cluster = c("adm1")); summary(fs_tv)

# Save data set for which there are obs both in first and second stage
sec_tv <- sec_tv[obs(fs_tv),]

# Predicted probabilities
sec_tv$phat0_tv <- as.numeric(predict(fs_tv, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec_tv$ac_obs <- ifelse(sec_tv$phat0_ac>0.5 & !is.na(sec_tv$phat0_ac), 1 , 0)
sec_tv$tv_obs <- ifelse(sec_tv$phat0_tv>0.5 & !is.na(sec_tv$phat0_tv), 1 , 0)

# Selection term
sec_tv$xb_noac = 1-sec_tv$phat0_ac               
sec_tv$selectionac = ifelse(sec_tv$ac==1, 
                             (sec_tv$xb_noac*log(sec_tv$xb_noac)/sec_tv$phat0_ac) + log(sec_tv$phat0_ac), 
                             (sec_tv$phat0_ac*log(sec_tv$phat0_ac)/sec_tv$xb_noac) + log(sec_tv$xb_noac))

sec_tv$xb_notv = 1-sec_tv$phat0_tv              
sec_tv$selectiontv = ifelse(sec_tv$tv==1, 
                              (sec_tv$xb_notv*log(sec_tv$xb_notv)/sec_tv$phat0_tv) + log(sec_tv$phat0_tv), 
                              (sec_tv$phat0_tv*log(sec_tv$phat0_tv)/sec_tv$xb_notv) + log(sec_tv$xb_notv))


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac | adm1

# With selection
model_tv1 <- feols(ely_formula, data = sec_tv, weights = ~weight, cluster = c("adm1")); summary(model_tv1)


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + tv +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectiontv | adm1

# With selection
model_tv2 <- feols(ely_formula, data = sec_tv, weights = ~weight, cluster = c("adm1")); summary(model_tv2)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + tv + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectiontv | adm1

# With selection
model_tv3 <- feols(ely_formula, data = sec_tv, weights = ~weight, cluster = c("adm1")); summary(model_tv3)


# Testing whether the impact of tv changes with temperature
# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ tv + tv*curr_CDD18_db + ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectiontv | adm1

# With selection
model_tv4 <- feols(ely_formula, data = sec_tv, weights = ~weight, cluster = c("adm1")); summary(model_tv4)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ tv + tv*curr_CDD18_db + tv*I(curr_CDD18_db^2) + ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectiontv | adm1

# With selection
model_tv5 <- feols(ely_formula, data = sec_tv, weights = ~weight, cluster = c("adm1")); summary(model_tv5)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ tv + tv*curr_CDD18_db + tv*I(curr_CDD18_db^2) + ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectiontv | adm1

# With selection
model_tv6 <- feols(ely_formula, data = sec_tv, weights = ~weight, cluster = c("adm1")); summary(model_tv6)

# Clean
rm(fs_ac, fs_tv)
gc()


#####################################

#        Including computers       #

####################################

# Logistic regression
fs_ac <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs_ac)

# Save data set for which there are obs both in first and second stage
sec_pc <- global[obs(fs_ac),]

# Predicted probabilities
sec_pc$phat0_ac <- as.numeric(predict(fs_ac, type="response"))

# Logistic regression
fs_pc <- feglm(pc_formula, family = binomial(link = "logit"), data = sec_pc, weights = ~weight, cluster = c("adm1")); summary(fs_pc)

# Save data set for which there are obs both in first and second stage
sec_pc <- sec_pc[obs(fs_pc),]

# Predicted probabilities
sec_pc$phat0_pc <- as.numeric(predict(fs_pc, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec_pc$ac_obs <- ifelse(sec_pc$phat0_ac>0.5 & !is.na(sec_pc$phat0_ac), 1 , 0)
sec_pc$pc_obs <- ifelse(sec_pc$phat0_pc>0.5 & !is.na(sec_pc$phat0_pc), 1 , 0)

# Selection term
sec_pc$xb_noac = 1-sec_pc$phat0_ac               
sec_pc$selectionac = ifelse(sec_pc$ac==1, 
                            (sec_pc$xb_noac*log(sec_pc$xb_noac)/sec_pc$phat0_ac) + log(sec_pc$phat0_ac), 
                            (sec_pc$phat0_ac*log(sec_pc$phat0_ac)/sec_pc$xb_noac) + log(sec_pc$xb_noac))

sec_pc$xb_nopc = 1-sec_pc$phat0_pc              
sec_pc$selectionpc = ifelse(sec_pc$pc==1, 
                            (sec_pc$xb_nopc*log(sec_pc$xb_nopc)/sec_pc$phat0_pc) + log(sec_pc$phat0_pc), 
                            (sec_pc$phat0_pc*log(sec_pc$phat0_pc)/sec_pc$xb_nopc) + log(sec_pc$xb_nopc))


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac | adm1

# With selection
model_pc1 <- feols(ely_formula, data = sec_pc, weights = ~weight, cluster = c("adm1")); summary(model_pc1)


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + pc +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionpc | adm1

# With selection
model_pc2 <- feols(ely_formula, data = sec_pc, weights = ~weight, cluster = c("adm1")); summary(model_pc2)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + pc + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionpc | adm1

# With selection
model_pc3 <- feols(ely_formula, data = sec_pc, weights = ~weight, cluster = c("adm1")); summary(model_pc3)


# Testing whether the impact of pc changes with temperature
# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ pc + pc*curr_CDD18_db + ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionpc | adm1

# With selection
model_pc4 <- feols(ely_formula, data = sec_pc, weights = ~weight, cluster = c("adm1")); summary(model_pc4)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ pc + pc*curr_CDD18_db + pc*I(curr_CDD18_db^2) + ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionpc | adm1

# With selection
model_pc5 <- feols(ely_formula, data = sec_pc, weights = ~weight, cluster = c("adm1")); summary(model_pc5)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ pc + pc*curr_CDD18_db + pc*I(curr_CDD18_db^2) + ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionpc | adm1

# With selection
model_pc6 <- feols(ely_formula, data = sec_pc, weights = ~weight, cluster = c("adm1")); summary(model_pc6)

# Clean
rm(fs_ac, fs_pc)
gc()


###########################################

#        Including washing machines       #

###########################################

# Logistic regression
fs_ac <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs_ac)

# Save data set for which there are obs both in first and second stage
sec_wshm <- global[obs(fs_ac),]

# Predicted probabilities
sec_wshm$phat0_ac <- as.numeric(predict(fs_ac, type="response"))

# Logistic regression
fs_wshm <- feglm(wshm_formula, family = binomial(link = "logit"), data = sec_wshm, weights = ~weight, cluster = c("adm1")); summary(fs_wshm)

# Save data set for which there are obs both in first and second stage
sec_wshm <- sec_wshm[obs(fs_wshm),]

# Predicted probabilities
sec_wshm$phat0_wshm <- as.numeric(predict(fs_wshm, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec_wshm$ac_obs <- ifelse(sec_wshm$phat0_ac>0.5 & !is.na(sec_wshm$phat0_ac), 1 , 0)
sec_wshm$wshm_obs <- ifelse(sec_wshm$phat0_wshm>0.5 & !is.na(sec_wshm$phat0_wshm), 1 , 0)

# Selection term
sec_wshm$xb_noac = 1-sec_wshm$phat0_ac               
sec_wshm$selectionac = ifelse(sec_wshm$ac==1, 
                            (sec_wshm$xb_noac*log(sec_wshm$xb_noac)/sec_wshm$phat0_ac) + log(sec_wshm$phat0_ac), 
                            (sec_wshm$phat0_ac*log(sec_wshm$phat0_ac)/sec_wshm$xb_noac) + log(sec_wshm$xb_noac))

sec_wshm$xb_nowshm = 1-sec_wshm$phat0_wshm              
sec_wshm$selectionwshm = ifelse(sec_wshm$wshm==1, 
                            (sec_wshm$xb_nowshm*log(sec_wshm$xb_nowshm)/sec_wshm$phat0_wshm) + log(sec_wshm$phat0_wshm), 
                            (sec_wshm$phat0_wshm*log(sec_wshm$phat0_wshm)/sec_wshm$xb_nowshm) + log(sec_wshm$xb_nowshm))


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac | adm1

# With selection
model_wshm1 <- feols(ely_formula, data = sec_wshm, weights = ~weight, cluster = c("adm1")); summary(model_wshm1)


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + wshm +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionwshm | adm1

# With selection
model_wshm2 <- feols(ely_formula, data = sec_wshm, weights = ~weight, cluster = c("adm1")); summary(model_wshm2)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + wshm + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionwshm | adm1

# With selection
model_wshm3 <- feols(ely_formula, data = sec_wshm, weights = ~weight, cluster = c("adm1")); summary(model_wshm3)


# Testing whether the impact of wshm changes with temperature
# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ wshm + wshm*curr_CDD18_db + ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionwshm | adm1

# With selection
model_wshm4 <- feols(ely_formula, data = sec_wshm, weights = ~weight, cluster = c("adm1")); summary(model_wshm4)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ wshm + wshm*curr_CDD18_db + wshm*I(curr_CDD18_db^2) + ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionwshm | adm1

# With selection
model_wshm5 <- feols(ely_formula, data = sec_wshm, weights = ~weight, cluster = c("adm1")); summary(model_wshm5)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ wshm + wshm*curr_CDD18_db + wshm*I(curr_CDD18_db^2) + ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selectionac + selectionwshm | adm1

# With selection
model_wshm6 <- feols(ely_formula, data = sec_wshm, weights = ~weight, cluster = c("adm1")); summary(model_wshm6)

# Clean
rm(fs_ac, fs_wshm)
gc()


## Exporting
# Mean electricity quantity
mean_ref <- weighted.mean(exp(sec_ref$ln_ely_q), sec_ref$weight)
mean_ref
mean_tv <- weighted.mean(exp(sec_tv$ln_ely_q), sec_tv$weight)
mean_tv
mean_pc <- weighted.mean(exp(sec_pc$ln_ely_q), sec_pc$weight)
mean_pc
mean_wshm <- weighted.mean(exp(sec_wshm$ln_ely_q), sec_wshm$weight)
mean_wshm

# Countries
cref <- length(unique.default(sec_ref$country))
cref
ctv <- length(unique.default(sec_tv$country))
ctv
cpc <- length(unique.default(sec_pc$country))
cpc
cwshm <- length(unique.default(sec_wshm$country))
cwshm


# Refrigerator
texreg(list(model_ref1, model_ref2, model_ref3, model_ref4, model_ref5, model_ref6), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - Refrigerators",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF", "DMF", "DMF", "DMF"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'TableA6.tex', sep=''), append=F,  
       float.pos = "htbp", label = "main: tableA6",
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "curr_CDD18_db:ac" = "AC $\\times$ CDD", 
                              "I(curr_CDD18_db^2):ac" = "AC $\\times$ CDD$^2$",
                              "ref"= "Refrigerator", "ref:curr_CDD18_db" = "Refrigerator $\\times$ CDD", 
                              "ref:I(curr_CDD18_db^2)" = "Refrigerator $\\times$ CDD$^2$"),
       custom.gof.rows = list("Correction Term (AC)" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Correction Term (Refrigerator)" = c("NO", "YES", "YES", "YES", "YES", "YEST"), 
                              "ADM-1 FE" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Mean Outcome" = c(mean_ref, mean_ref, mean_ref, mean_ref, mean_ref, mean_ref), 
                              "Countries" = c(cref, cref, cref, cref, cref, cref)), 
       caption.above = TRUE)

# Television
texreg(list(model_tv1, model_tv2, model_tv3, model_tv4, model_tv5, model_tv6), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - Television",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF", "DMF", "DMF", "DMF"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'TableA7.tex', sep=''), append=F,  
       float.pos = "htbp", label = "main: tableA7",
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "curr_CDD18_db:ac" = "AC $\\times$ CDD", 
                              "I(curr_CDD18_db^2):ac" = "AC $\\times$ CDD$^2$",
                              "tv"= "TV", "tv:curr_CDD18_db" = "TV $\\times$ CDD", 
                              "tv:I(curr_CDD18_db^2)" = "TV $\\times$ CDD$^2$"),
       custom.gof.rows = list("Correction Term (AC)" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Correction Term (TV)" = c("NO", "YES", "YES", "YES", "YES", "YES"), 
                              "ADM-1 FE" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Mean Outcome" = c(mean_tv, mean_tv, mean_tv, mean_tv, mean_tv, mean_tv), 
                              "Countries" = c(ctv, ctv, ctv, ctv, ctv, ctv)), 
       caption.above = TRUE)

# PC
texreg(list(model_pc1, model_pc2, model_pc3, model_pc4, model_pc5, model_pc6), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - PC",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF", "DMF", "DMF", "DMF"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'TableA8.tex', sep=''), append=F,  
       float.pos = "htbp", label = "main: tableA8",
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "curr_CDD18_db:ac" = "AC $\\times$ CDD", 
                              "I(curr_CDD18_db^2):ac" = "AC $\\times$ CDD$^2$",
                              "pc"= "PC", "pc:curr_CDD18_db" = "PC $\\times$ CDD", 
                              "pc:I(curr_CDD18_db^2)" = "PC $\\times$ CDD$^2$"),
       custom.gof.rows = list("Correction Term (AC)" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Correction Term (pc)" = c("NO", "YES", "YES", "YES", "YES", "YES"), 
                              "ADM-1 FE" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Mean Outcome" = c(mean_pc, mean_pc, mean_pc, mean_pc, mean_pc, mean_pc), 
                              "Countries" = c(cpc, cpc, cpc, cpc, cpc, cpc)), 
       caption.above = TRUE)

# Washing Machine
texreg(list(model_wshm1, model_wshm2, model_wshm3, model_wshm4, model_wshm5, model_wshm6), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - Washing Machine",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF", "DMF", "DMF", "DMF"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'TableA9.tex', sep=''), append=F,  
       float.pos = "htbp", label = "main: tableA9",
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "curr_CDD18_db:ac" = "AC $\\times$ CDD", 
                              "I(curr_CDD18_db^2):ac" = "AC $\\times$ CDD$^2$",
                              "wshm"= "Washing Machine", "wshm:curr_CDD18_db" = "Washing Machine $\\times$ CDD", 
                              "wshm:I(curr_CDD18_db^2)" = "Washing Machine $\\times$ CDD$^2$"),
       custom.gof.rows = list("Correction Term (AC)" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Correction Term (Washing Machine)" = c("NO", "YES", "YES", "YES", "YES", "YES"), 
                              "ADM-1 FE" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Mean Outcome" = c(mean_wshm, mean_wshm, mean_wshm, mean_wshm, mean_wshm, mean_wshm), 
                              "Countries" = c(cwshm, cwshm, cwshm, cwshm, cwshm, cwshm)), 
       caption.above = TRUE)

