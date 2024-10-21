
#############################################################

#               Table 3, Table A1, Table A2

#############################################################

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
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

house <- paste(stub,'6-Projections/data/household/', sep='')
interm <- paste(stub,'results/regressions/for_projections/', sep='')
interm <- "C:/Users/Standard/Documents/Github/acglobal/interm/"
output <- paste(stub,'6-Projections/results/regressions/', sep='')
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
                            curr_HDD18_db2 = curr_HDD18_db^2,
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




##################################

#     Extensive margin choice

##################################

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  curr_HDD18_db + curr_HDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | adm1

# Linear probability model
reg_ac00 <- feols(ac_formula, data = global, weights = ~weight, cluster = c("adm1")); summary(reg_ac00)

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  curr_HDD18_db + curr_HDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | adm1

# Linear probability model
reg_ac0 <- feols(ac_formula, data = global, weights = ~weight, cluster = c("adm1")); summary(reg_ac0)

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)
gc()

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),] # dropped obs for which logit predicts perfectly 0 or 1

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), # [P_0 * lnP_0 / (1 - P_0)] + lnP_1 = [P_0 * lnP_0 / P_1] + lnP_1
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac)) # [P_1 * lnP_1 / (1 - P_1)] + lnP_0 = [P_1 * lnP_1 / P_0] + lnP_0

# Mean AC
mean_ac <- weighted.mean(sec$ac, sec$weight)
mean_ac

# Number of countries
ncoun <- length(unique.default(sec$country))
ncoun

# Table B2
texreg(list(reg_ac00, reg_ac0, fs), digits = 3, caption = "Logit Regression for Air-conditioning Ownership --- Global",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("LPM", "LPM", "Logit"),
       custom.note = "\\textbf{Notes}: Dependent variable is air-conditioning (0,1). Clustered std. errors at the ADM-1 level in parentheses. Column (4) shows the average marginal effects (AMEs) from the logit regression. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$. Regressions are conducted using survey weights.", 
       file = paste(output,'TableA2_1to3.tex', sep=''), append=F,  
       float.pos = "H", label = "main: ac_global",
       omit.coef = "(country)|(Intercept)", 
       custom.coef.map = list("mean_CDD18_db" = "$\\overline{CDD}$", "mean_CDD18_db2" = "$\\overline{CDD}^2$", 
                              "mean_CDD18_db_exp" = "$\\overline{CDD}$ $\\times$ Log(Exp)", 
                              "mean_CDD18_db2_exp" = "$\\overline{CDD}^2$ $\\times$ Log(Exp)",
                              "curr_CDD18_db" = "$CDD$", "curr_CDD18_db2" = "$CDD^2$",
                              "curr_HDD18_db" = "$HDD$", "curr_HDD18_db2" = "$HDD^2$",
                              "ln_ely_p_cdd" = "$\\overline{CDD}$ $\\times$ Log(P)", 
                              "ln_ely_p_cdd2" = "$\\overline{CDD}^2$ $\\times$ Log(P)",
                              "ln_total_exp_usd_2011" = "Log(Exp)",
                              "ln_ely_p" = "Log(P)",
                              "ln_ely_p_nme" = "Log(P) $\\times$ Household Size",
                              "ln_ely_p_own" = "Log(P) $\\times$ House Ownership",
                              "urban_sh" = "Urbanisation (\\%)", 
                              "ownership_d1" = "House Ownership (Yes = 1)", 
                              "n_members" = "Household Size", 
                              "edu_head_21" = "Primary Edu.", 
                              "edu_head_22" = "Secondary Edu.", 
                              "edu_head_23" = "Post Edu.", 
                              "age_head" = "Age (Head)", 
                              "sex_head" = "Female (Yes = 1)"),
       custom.gof.rows = list("ADM-1 FE" = c("YES", "YES", "YES"), 
                              "Mean Outcome" = c(mean_ac, mean_ac, mean_ac), 
                              "Countries" = c(ncoun, ncoun, ncoun)), caption.above = TRUE)

# Average marginal effects (AMEs)
rm(reg_ac0, reg_ac00)
margins <- avg_slopes(fs, wts = sec$weight, newdata = sec)
summary(margins)
gc()

# Table B2
xtable(summary(margins), display=rep('g', 9), caption = "Logit Regression for Air-conditioning Ownership - Global", digits = 3)
print(xtable(summary(margins), display=rep('g', 9), caption = "Logit Regression for Air-conditioning Ownership - Global"), 
      file= paste(output,'TableA2_4.tex', sep=''),append=F, table.placement = "htbp",
      caption.placement="top")

# Clean
rm(margins)
gc()


# Get full marginal effects
# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + I(mean_CDD18_db^2) + 
  mean_CDD18_db*ln_total_exp_usd_2011 + I(mean_CDD18_db^2)*ln_total_exp_usd_2011 + ln_total_exp_usd_2011 + 
  curr_CDD18_db + I(curr_CDD18_db^2) +
  curr_HDD18_db + I(curr_HDD18_db^2) +
  ln_ely_p + ln_ely_p*mean_CDD18_db + ln_ely_p*I(mean_CDD18_db^2) + ln_ely_p*n_members + ln_ely_p*ownership_d + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | adm1

# LPM
flpm <- feols(ac_formula, data = global, weights = ~weight, cluster = c("adm1")); summary(flpm)
flpm_marg <- avg_slopes(flpm, wts = sec$weight, newdata = sec)
summary(flpm_marg) # 100 CDD -> 4.1, 1% Exp -> 0.08
gc()

# Logit regression
ffs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(ffs)
ffs_marg <- avg_slopes(ffs, wts = sec$weight, newdata = sec)
summary(ffs_marg) # 100 CDD -> 5.7, 1% Exp -> 0.08
gc()


########################################

#      Intensive margin decision

########################################

# Formula electricity expenditure without selection
ely_formula  <- ln_ely_q ~ ac | 
  adm1

# Without selection
model00 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model00)


# Formula electricity expenditure without selection
ely_formula  <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head | adm1

# Without selection
model0 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0)


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model1 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model1)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model2 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model2)

# Mean electricity quantity
mean <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean

# Number of countries
ncoun <- length(unique.default(sec$country))
ncoun


## Export
# Compare the models
screenreg(list(model00, model0, model1, model2), digits = 3, 
          caption = "The Effect of Air-conditioning on Residential Electricity Quantity",
          stars = c(0.1, 0.05, 0.01), custom.model.names = c("OLS", "OLS", "DMF", "DMF"),
          omit.coef = "(country)|(Intercept)|(selection)",
          custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", "curr_CDD18_db" = "CDD", 
                                 "I(curr_CDD18_db^2)" = "CDD$^2$", "curr_HDD18_db" = "HDD", "I(curr_HDD18_db^2)" = "HDD$^2$",
                                 "ln_total_exp_usd_2011" = "Log(Exp)",
                                 "ln_ely_p" = "Log(P)",
                                 "urban_sh" = "Urbanisation (\\%)", 
                                 "ownership_d" = "House Ownership (Yes = 1)", "n_members" = "Household Size", 
                                 "edu_head_21" = "Primary Edu.", "edu_head_22" = "Secondary Edu.", "edu_head_23" = "Post Edu.", 
                                 "age_head" = "Age (Head)", "sex_head" = "Female (Yes = 1)",
                                 "selection" = "$\\hat{\\zeta}$"),
          custom.gof.rows = list("Correction Term" = c("NO", "NO", "YES", "YES"), 
                                 "ADM-1 FE" = c("YES", "YES", "YES", "YES"), 
                                 "Mean Outcome (kWh)" = c(mean, mean, mean, mean), 
                                 "Countries" = c("25", "25", "25", "25")))

# Table A1
texreg(list(model00, model0, model1, model2), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("OLS", "OLS", "DMF", "DMF"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'TableA1.tex', sep=''), append=F,  
       float.pos = "htbp", label = "main: ely_global",
       omit.coef = "(country)|(Intercept)",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", "curr_CDD18_db" = "CDD", 
                              "I(curr_CDD18_db^2)" = "CDD$^2$", "curr_HDD18_db" = "HDD", "I(curr_HDD18_db^2)" = "HDD$^2$",
                              "ln_total_exp_usd_2011" = "Log(Exp)",
                              "ln_ely_p" = "Log(P)",
                              "urban_sh" = "Urbanisation (\\%)", 
                              "ownership_d" = "House Ownership (Yes = 1)", "n_members" = "Household Size", 
                              "edu_head_21" = "Primary Edu.", "edu_head_22" = "Secondary Edu.", "edu_head_23" = "Post Edu.", 
                              "age_head" = "Age (Head)", "sex_head" = "Female (Yes = 1)",
                              "selection" = "$\\hat{\\zeta}$"),
       custom.gof.rows = list("ADM-1 FE" = c("YES", "YES", "YES", "YES"), 
                              "Mean Outcome (kWh)" = c(mean, mean, mean, mean), 
                              "Countries" = c("25", "25", "25", "25")), 
       caption.above = TRUE)

# Table 3
texreg(list(model00, model0, model1, model2), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("OLS", "OLS", "DMF", "DMF"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'Table3.tex', sep=''), append=F,  
       float.pos = "htbp", label = "main: ely_global",
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$"),
       custom.gof.rows = list("Controls" = c("NO", "YES", "YES", "YES"),
                              "Correction Term" = c("NO", "NO", "YES", "YES"), 
                              "ADM-1 FE" = c("YES", "YES", "YES", "YES"), 
                              "Mean Outcome (kWh)" = c(mean, mean, mean, mean), 
                              "Countries" = c("25", "25", "25", "25")), 
       caption.above = TRUE)


# Save the R Environment will be used for the projections
reg_ac <- fs
reg_ely <- model2
save(list = c("reg_ac", "reg_ely", "sec", "global"), 
     file = paste(interm,'global_wgt_dmcf.RData', sep=''))

