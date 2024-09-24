
#################################################################

# Additional analysis:

#   - Controlling for potential solar output x solar capacity
#   - Interacting PVCAP X PVOUT with air-conditioning
#   - Interacting PVCAP X PVOUT with AC X f(CDD)
#   - Interacting PVCAP X PVOUT with electricity prices

# We use both continuous and binary/categorical specifications

#################################################################

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
  stub <- 'F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

house <- paste(stub,'6-Projections/data/household/', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')

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
global  <- global %>% mutate(pvcap_kw = pvcap*1000, # kW
                             pvint = pvcap*pvout,
                             pvint_kWh = pvcap_kw*pvout, # kW/(kWh/kW) -> kWh
                             pvint_MWh = pvint/1000, # MWh
                             dpvint = ifelse(pvint > median(pvint, na.rm = TRUE), 1, 0),
                             log_pvint_kWh = log(pvint_kWh+1),
                             as_pvint_kWh = asinh(pvint_kWh))

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



################################################

#    1) Continuous variable: PVOUT X PVCAP     #

################################################

# Formula first-stage
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + pvcap*pvout | adm1

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
meanpv <- global %>% group_by(country) %>% summarise(pvcap = mean(pvcap, na.rm = TRUE), pvout = mean(pvout, na.rm = TRUE)) 
print(meanpv, n = 25)

# Mean AC
mean_ac <- weighted.mean(sec$ac, sec$weight)
mean_ac


# Formula electricity q without selection
ely_formula  <- ln_ely_q ~ ac + pvcap*pvout +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head | adm1

# Without selection
model0 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0)


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + pvcap*pvout +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model1 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model1)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*pvcap*pvout + pvcap*pvout +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model2 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model2)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ln_ely_p*pvcap*pvout + pvcap*pvout +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model3 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model3)

# Mean electricity quantity
mean <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean


## Export
# Compare the models
screenreg(list(model1, model2, model3), digits = 3, 
          caption = "The Effect of Air-conditioning on Residential Electricity Quantity - PV Potential Output and PV Capacity",
          stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF"),
          omit.coef = "(country)|(Intercept)|(selection)",
          custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                                 "pvout" = "PVOUT",
                                 "pvcap" = "PV Capacity",
                                 "pvcap:pvout" = "PV Capacity $\\times$ PVOUT",
                                 "ac:pvout" = "AC $\\times$ PVOUT",
                                 "ac:pvcap" = "AC $\\times$ PV Capacity",
                                 "ac:pvcap:pvout" = "AC $\\times$ PV Capacity $\\times$ PVOUT",
                                 "ln_ely_p" = "Log(P)",
                                 "ln_ely_p:pvout" = "Log(P) $\\times$ PVOUT",
                                 "ln_ely_p:pvcap" = "Log(P) $\\times$ PV Capacity",
                                 "ln_ely_p:pvcap:pvout" = "Log(P) $\\times$ PV Capacity $\\times$ PVOUT"),
          custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES"), 
                                 "ADM-1 FE" = c("YES", "YES", "YES"), 
                                 "Mean Outcome" = c(mean, mean, mean), 
                                 "Countries" = c("25", "25", "25")))

# Full sample
texreg(list(model1, model2, model3), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - PV Potential Output and PV Capacity (continuous)",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF"),
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/additional/Global_wgt_pvoutxcap.tex', sep=''), append=F,  
       float.pos = "htbp", label = "app: ely_global_pvoutxcap",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "pvout" = "PVOUT",
                              "pvcap" = "PV Capacity",
                              "pvcap:pvout" = "PV Capacity $\\times$ PVOUT",
                              "ac:pvout" = "AC $\\times$ PVOUT",
                              "ac:pvcap" = "AC $\\times$ PV Capacity",
                              "ac:pvcap:pvout" = "AC $\\times$ PV Capacity $\\times$ PVOUT",
                              "ln_ely_p" = "Log(P)",
                              "ln_ely_p:pvout" = "Log(P) $\\times$ PVOUT",
                              "ln_ely_p:pvcap" = "Log(P) $\\times$ PV Capacity",
                              "ln_ely_p:pvcap:pvout" = "Log(P) $\\times$ PV Capacity $\\times$ PVOUT"),
       custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES"), 
                              "ADM-1 FE" = c("YES", "YES", "YES"), 
                              "Mean Outcome" = c(mean, mean, mean), 
                              "Countries" = c("25", "25", "25")), 
       caption.above = TRUE)

# Air-conditioning
texreg(list(fs), digits = 3, caption = "Logit Regression for Air-conditioning Ownership - PV Potential Output and PV Capacity (continuous)",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("Logit"),
       custom.note = "\\textbf{Notes}: Dependent variable is air-conditioning (0,1). Clustered std. errors at the ADM-1 level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$. Regressions are conducted using survey weights.", 
       file = paste(output,'airconditioning/additional/Global_wgt_coeff_pvoutxpvcap.tex', sep=''), append=F,  
       float.pos = "H", label = "app: ac_global_pvoutxpvcap",
       omit.coef = "(country)|(Intercept)", 
       custom.coef.map = list("pvout" = "PVOUT",
                              "pvcap" = "PV Capacity",
                              "pvcap:pvout" = "PV Capacity $\\times$ PVOUT"),
       custom.gof.rows = list("ADM-1 FE" = c("YES"), 
                              "Mean Outcome" = c(mean_ac), 
                              "Countries" = c("25")), caption.above = TRUE)

# Clean
rm(model0, model1, model2, model3, fs, sec)
gc()



################################################

#     2) Dummy variable: D(PVOUT X PVCAP)       

# We define a dummy variable for above the 
# median value of such interaction

################################################

# Formula first-stage
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + dpvint | adm1

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


# Formula electricity q without selection
ely_formula  <- ln_ely_q ~ ac + dpvint +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head | adm1

# Without selection
model4 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model4)


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + dpvint +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model5 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model5)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*dpvint + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model6 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model6)

# Plot
plot_slopes(model6, variables = "ac", slope = "dydx", condition = c("dpvint"))


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*dpvint*curr_CDD18_db + ac*dpvint*I(curr_CDD18_db^2) +  
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model7 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model7)

# Plot
plot_slopes(model7, variables = "ac", condition = list("curr_CDD18_db" = 0:30, "dpvint"))


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ln_ely_p*dpvint +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model8 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model8)

# Mean electricity quantity
mean <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean

# Clean
rm(sec)
gc()


## Export
# Compare the models
screenreg(list(model5, model6, model8), digits = 3, 
          caption = "The Effect of Air-conditioning on Residential Electricity Quantity - PV Potential Output and PV Capacity (dummy)",
          stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF"),
          omit.coef = "(country)|(Intercept)|(selection)",
          custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                                 "dpvint" = "$\\mathds{1}$(PV Gen. $\\times$ PVOUT > Median)",
                                 "ac:dpvint" = "AC $\\times$ $\\mathds{1}$(PV Gen. $\\times$ PVOUT > Median)",
                                 "ln_ely_p" = "Log(P)",
                                 "ln_ely_p:dpvint" = "Log(P) $\\times$ $\\mathds{1}$(PV Gen > Median)"),
          custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES"), 
                                 "ADM-1 FE" = c("YES", "YES", "YES"), 
                                 "Mean Outcome" = c(mean, mean, mean), 
                                 "Countries" = c("25", "25", "25")))

# Full sample
texreg(list(model5, model6, model8), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - PV Potential Output and PV Capacity (dummy)",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF"),
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/additional/Global_wgt_dpvint.tex', sep=''), append=F,  
       float.pos = "htbp", label = "app: ely_global_dpvint",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "dpvint" = "$\\mathds{1}$(PV Gen. $\\times$ PVOUT > Median)",
                              "ac:dpvint" = "AC $\\times$ $\\mathds{1}$(PV Gen. $\\times$ PVOUT > Median)",
                              "ln_ely_p" = "Log(P)",
                              "ln_ely_p:dpvint" = "Log(P) $\\times$ $\\mathds{1}$(PV Gen > Median)"),
       custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES"), 
                              "ADM-1 FE" = c("YES", "YES", "YES"), 
                              "Mean Outcome" = c(mean, mean, mean), 
                              "Countries" = c("25", "25", "25")), 
       caption.above = TRUE)

# Air-conditioning
texreg(list(fs), digits = 3, caption = "Logit Regression for Air-conditioning Ownership - PV Potential Output and PV Capacity (dummy)",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("Logit"),
       custom.note = "\\textbf{Notes}: Dependent variable is air-conditioning (0,1). Clustered std. errors at the ADM-1 level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$. Regressions are conducted using survey weights.", 
       file = paste(output,'airconditioning/additional/Global_wgt_coeff_dpvint.tex', sep=''), append=F,  
       float.pos = "H", label = "app: ac_global_dpvint",
       omit.coef = "(country)|(Intercept)", 
       custom.coef.map = list("dpvint" = "$\\mathds{1}$(PV Gen. $\\times$ PVOUT > Median)"),
       custom.gof.rows = list("ADM-1 FE" = c("YES"), 
                              "Mean Outcome" = c(mean_ac), 
                              "Countries" = c("25")), caption.above = TRUE)



################################################

#               3) PVGEN: in MWh               #   

################################################

# Formula first-stage
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + pvint_MWh | adm1

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


# Formula electricity q without selection
ely_formula  <- ln_ely_q ~ ac + pvint_MWh +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head | adm1

# Without selection
model9 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model9)


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + pvint_MWh +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model10 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model10)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*pvint_MWh + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model11 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model11)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ln_ely_p*pvint_MWh + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model12 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model12)

# Mean electricity quantity
mean <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean

# Clean
rm(sec)
gc()


## Export
# Air-conditioning
texreg(list(fs), digits = 3, caption = "Logit Regression for Air-conditioning Ownership - PV Gen.",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("Logit"),
       custom.note = "\\textbf{Notes}: Dependent variable is air-conditioning (0,1). Clustered std. errors at the ADM-1 level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$. Regressions are conducted using survey weights.", 
       file = paste(output,'airconditioning/additional/Global_wgt_coeff_pvgen.tex', sep=''), append=F,  
       float.pos = "H", label = "app: ac_global_pvgen",
       omit.coef = "(country)|(Intercept)", 
       custom.coef.map = list("pvint_MWh" = "PV Gen. (MWh)"),
       custom.gof.rows = list("ADM-1 FE" = c("YES"), 
                              "Mean Outcome" = c(mean_ac), 
                              "Countries" = c("25")), caption.above = TRUE)



################################################

#         4) PVGEN: in asinh/log+1             #   

################################################

# Formula first-stage
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + as_pvint_kWh | adm1

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


# Formula electricity q without selection
ely_formula  <- ln_ely_q ~ ac + as_pvint_kWh +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head | adm1

# Without selection
model13 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model13)


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + as_pvint_kWh +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model14 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model14)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*as_pvint_kWh + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model15 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model15)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ln_ely_p*as_pvint_kWh + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model16 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model16)

# Mean electricity quantity
mean <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean

# Clean
rm(sec)
gc()


## Export
# Compare the models
screenreg(list(model10, model11, model12, model14, model15, model16), digits = 3, 
          caption = "The Effect of Air-conditioning on Residential Electricity Quantity - PV Gen",
          stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF", "DMF", "DMF", "DMF"),
          omit.coef = "(country)|(Intercept)|(selection)",
          custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                                 "pvint_MWh" = "PV Gen.",
                                 "ac:pvint_MWh" = "AC $\\times$ PV Gen.",
                                 "ln_ely_p" = "Log(P)",
                                 "ln_ely_p:pvint_MWh" = "Log(P) $\\times$ PV Gen.",
                                 "as_pvint_kWh" = "asinh(PV Gen.)",
                                 "ac:as_pvint_kWh" = "AC $\\times$ asinh(PV Gen.)",
                                 "ln_ely_p:as_pvint_kWh" = "Log(P) $\\times$ asinh(PV Gen.)"),
          custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "ADM-1 FE" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                                 "Mean Outcome" = c(mean, mean, mean, mean, mean, mean), 
                                 "Countries" = c("25", "25", "25", "25", "25", "25")))

# Full sample
texreg(list(model10, model11, model12, model14, model15, model16), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - PV Gen.",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF", "DMF", "DMF", "DMF"),
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/additional/Global_wgt_pvgen.tex', sep=''), append=F,  
       float.pos = "htbp", label = "app: ely_global_pvgen",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "pvint_MWh" = "PV Gen. (MWh)",
                              "ac:pvint_MWh" = "AC $\\times$ PV Gen.",
                              "ln_ely_p" = "Log(P)",
                              "ln_ely_p:pvint_MWh" = "Log(P) $\\times$ PV Gen.",
                              "as_pvint_kWh" = "asinh(PV Gen.)",
                              "ac:as_pvint_kWh" = "AC $\\times$ asinh(PV Gen.)",
                              "ln_ely_p:as_pvint_kWh" = "Log(P) $\\times$ asinh(PV Gen.)"),
       custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "ADM-1 FE" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Mean Outcome" = c(mean, mean, mean, mean, mean, mean), 
                              "Countries" = c("25", "25", "25", "25", "25", "25")), 
       caption.above = TRUE)

# Air-conditioning
texreg(list(fs), digits = 3, caption = "Logit Regression for Air-conditioning Ownership - PV Gen. (asinh)",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("Logit"),
       custom.note = "\\textbf{Notes}: Dependent variable is air-conditioning (0,1). Clustered std. errors at the ADM-1 level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$. Regressions are conducted using survey weights.", 
       file = paste(output,'airconditioning/additional/Global_wgt_coeff_aspvgen.tex', sep=''), append=F,  
       float.pos = "H", label = "app: ac_global_aspvgen",
       omit.coef = "(country)|(Intercept)", 
       custom.coef.map = list("as_pvint_kWh" = "asinh(PV Gen.)"),
       custom.gof.rows = list("ADM-1 FE" = c("YES"), 
                              "Mean Outcome" = c(mean_ac), 
                              "Countries" = c("25")), caption.above = TRUE)


