
##########################################

#      Table 5, Table A13, Table A14

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
                            curr_CDD18_db2 = curr_CDD18_db^2,
                            edu_head_2 = as.factor(edu_head_2))

# Median 
global  <- global %>% mutate(pvgen = pvcap*pvout, # pvgen in MW
                             pvgen_kWh = pvcap*1000*pvout, # kW/(kWh/kW) -> kWh
                             dpvint = ifelse(pvgen > median(pvgen, na.rm = TRUE), 1, 0),
                             pvgen = asinh(pvgen))

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
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + pvgen | adm1

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


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + pvgen +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model1 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model1)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*pvgen +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model2 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model2)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ln_ely_p*pvgen +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model3 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model3)

# Mean electricity quantity
mean <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean

# Country
cntry <- length(unique.default(sec$country))
cntry


## Export
# Full sample
texreg(list(model1, model2, model3), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - PV Generation",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF"),
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'TableA8.tex', sep=''), append=F,  
       float.pos = "htbp", label = "main: tableA8",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "pvgen" = "asinh(PV Generation)",
                              "ac:pvgen" = "AC $\\times$ PV Capacity $\\times$ asinh(PV Generation)",
                              "ln_ely_p" = "Log(P)",
                              "ln_ely_p:pvgen" = "Log(P) $\\times$ asinh(PV Generation)"),
       custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES"), 
                              "ADM-1 FE" = c("YES", "YES", "YES"), 
                              "Mean Outcome" = c(mean, mean, mean), 
                              "Countries" = c(cntry, cntry, cntry)), 
       caption.above = TRUE)

# Air-conditioning
texreg(list(fs), digits = 3, caption = "Logit Regression for Air-conditioning Ownership - PV Potential Output and PV Capacity (continuous)",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("Logit"),
       custom.note = "\\textbf{Notes}: Dependent variable is air-conditioning (0,1). Clustered std. errors at the ADM-1 level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$. Regressions are conducted using survey weights.", 
       file = paste(output,'TableA9_1.tex', sep=''), append=F,  
       float.pos = "H", label = "main: tableA9",
       omit.coef = "(country)|(Intercept)", 
       custom.coef.map = list("pvout" = "PVOUT",
                              "pvcap" = "PV Capacity",
                              "pvgen" = "asinh(PV Generation)"),
       custom.gof.rows = list("ADM-1 FE" = c("YES"), 
                              "Mean Outcome" = c(mean_ac), 
                              "Countries" = c(cntry)), caption.above = TRUE)

# Clean
rm(model1, model2, model3, fs, sec)
gc()



## 2) Dummy variable: D(PVOUT X PVCAP) 
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


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + dpvint +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model4 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model4)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*dpvint + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model5 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model5)

# Plot
plot_slopes(model5, variables = "ac", slope = "dydx", condition = c("dpvint"))


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ln_ely_p*dpvint +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model6 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model6)

# Mean electricity quantity
mean <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean

# Country
cntry <- length(unique.default(sec$country))
cntry

# Clean
gc()


## Export
# Full sample
texreg(list(model4, model5, model6), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - PV Potential Output and PV Capacity (dummy)",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF", "DMF", "DMF"),
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'Table4.tex', sep=''), append=F,  
       float.pos = "htbp", label = "main: table4",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "dpvint" = "$\\mathds{1}$(PV Gen. $\\times$ PVOUT > Median)",
                              "ac:dpvint" = "AC $\\times$ $\\mathds{1}$(PV Gen. $\\times$ PVOUT > Median)",
                              "ln_ely_p" = "Log(P)",
                              "ln_ely_p:dpvint" = "Log(P) $\\times$ $\\mathds{1}$(PV Gen > Median)"),
       custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES"), 
                              "ADM-1 FE" = c("YES", "YES", "YES"), 
                              "Mean Outcome" = c(mean, mean, mean), 
                              "Countries" = c(cntry, cntry, cntry)), 
       caption.above = TRUE)

# Air-conditioning
texreg(list(fs), digits = 3, caption = "Logit Regression for Air-conditioning Ownership - PV Potential Output and PV Capacity (dummy)",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("Logit"),
       custom.note = "\\textbf{Notes}: Dependent variable is air-conditioning (0,1). Clustered std. errors at the ADM-1 level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$. Regressions are conducted using survey weights.", 
       file = paste(output,'tableA9_2.tex', sep=''), append=F,  
       float.pos = "H", label = "main: tableA9_2",
       omit.coef = "(country)|(Intercept)", 
       custom.coef.map = list("dpvint" = "$\\mathds{1}$(PV Gen. $\\times$ PVOUT > Median)"),
       custom.gof.rows = list("ADM-1 FE" = c("YES"), 
                              "Mean Outcome" = c(mean_ac), 
                              "Countries" = c(cntry)), caption.above = TRUE)


