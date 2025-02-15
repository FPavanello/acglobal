
##########################################

#                 Table S18

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
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'data/household/', sep='')
interm <- paste(stub,'results/regressions/for_graphs/subsamples/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/'
output <- paste(stub,'output/supplementary/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/supplementary/'


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
                            curr_CDD18_db2 = curr_CDD18_db^2,
                            edu_head_2 = as.factor(edu_head_2))

# Median 
global  <- global %>% mutate(pvgen = pvcap*pvout, # pvgen in MW
                             pvgen_kWh = pvcap*1000*pvout, # kW/(kWh/kW) -> kWh
                             dpvint = ifelse(pvgen > median(pvgen, na.rm = TRUE), 1, 0),
                             pvgen = asinh(pvgen),
                             pvcap = asinh(pvcap))

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


## 1) Continuous variable: PVGEN
# Formula first-stage
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  curr_HDD18_db + I(curr_HDD18_db^2) +
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + pvgen | country

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

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

# Mean AC
mean_ac <- weighted.mean(sec$ac, sec$weight)
mean_ac


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*pvgen +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model1 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model1)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ln_ely_p*pvgen +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model2 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model2)

# Mean electricity quantity
meang <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
meang

# Country
cntryg <- length(unique.default(sec$country))
cntryg

# Clean
gc()


## 2) Dummy variable: D(PVOUT X PVCAP) 
# Formula first-stage
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  curr_HDD18_db + I(curr_HDD18_db^2) +
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + dpvint | country

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

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

# Average
meanpv <- global %>% group_by(country) %>% summarise(pvcap = mean(pvcap, na.rm = TRUE), pvout = mean(pvout, na.rm = TRUE), dpvint = mean(dpvint, na.rm = TRUE)) 
print(meanpv, n = 25)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*dpvint + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model3 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model3)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ln_ely_p*dpvint +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model4 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model4)

# Mean electricity quantity
meand <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
meand

# Country
cntryd <- length(unique.default(sec$country))
cntryd

# Clean
gc()


## 3) PV Capacity
# Formula first-stage
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  curr_HDD18_db + I(curr_HDD18_db^2) +
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + pvcap | country

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

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

# Average
meanpv <- global %>% group_by(country) %>% summarise(pvcap = mean(pvcap, na.rm = TRUE), pvout = mean(pvout, na.rm = TRUE), dpvint = mean(dpvint, na.rm = TRUE)) 
print(meanpv, n = 25)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*pvcap + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model5 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model5)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ln_ely_p*pvcap +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model6 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model6)

# Mean electricity quantity
meanc <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
meanc

# Country
cntryc <- length(unique.default(sec$country))
cntryc

# Clean
gc()


## 4) PV Output
# Formula first-stage
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  curr_HDD18_db + I(curr_HDD18_db^2) +
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + pvout | country

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

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

# Average
meanpv <- global %>% group_by(country) %>% summarise(pvcap = mean(pvcap, na.rm = TRUE), pvout = mean(pvout, na.rm = TRUE), dpvint = mean(dpvint, na.rm = TRUE)) 
print(meanpv, n = 25)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*pvout + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model7 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model7)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ln_ely_p*pvout +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model8 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model8)

# Mean electricity quantity
meano <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
meano

# Country
cntryo <- length(unique.default(sec$country))
cntryo

# Clean
gc()


## Export
# Full sample
texreg(list(model1, model2, model3, model4, model5, model6, model7, model8), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - PV (Country FE)",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF", "DMF", "DMF", "DMF", "DMF", "DMF"),
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'TableS18.tex', sep=''), append=F,  
       float.pos = "htbp", label = "si: tableS18",
       custom.coef.map = list("ac"= "AC",
                              "pvgen" = "asinh(PV Generation)",
                              "dpvint" = "$\\mathds{1}$(PV Gen. $\\times$ PVOUT > Median)",
                              "pvcap" = "asinh(PV Capacity)",
                              "pvout" = "PVOUT",
                              "ac:pvgen" = "AC $\\times$ asinh(PV Generation)",
                              "ac:dpvint" = "AC $\\times$ $\\mathds{1}$(PV Gen. $\\times$ PVOUT > Median)",
                              "ac:pvcap" = "AC $\\times$ asinh(PV Capacity)",
                              "ac:pvout" = "AC $\\times$ PVOUT",
                              "ln_ely_p" = "Log(P)",
                              "ln_ely_p:pvgen" = "Log(P) $\\times$ asinh(PV Generation)",
                              "ln_ely_p:dpvint" = "Log(P) $\\times$ $\\mathds{1}$(PV Gen > Median)",
                              "ln_ely_p:pvcap" = "Log(P) $\\times$ asinh(PV Capacity)",
                              "ln_ely_p:pvout" = "Log(P) $\\times$ PVOUT"),
       custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Country FE" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Controls" = c("YES", "YES", "YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Mean Outcome" = c(meang, meang, meand, meand, meanc, meanc, meano, meano), 
                              "Countries" = c(cntryg, cntryg, cntryd, cntryd, cntryc, cntryc, cntryo, cntryo)), 
       caption.above = TRUE)
