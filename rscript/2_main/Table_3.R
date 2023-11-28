
## This R-script:
##      1) exploits CDD-dry bulb 18 deg and HDD-dry bulb 18 deg
##      2) conducts logit regressions for the global data set
##      3) run intensive margin regressions: electricity expenditure on climate + covariates
##         using Dubin and McFadden (1984) approach
##      4) apply survey weights to run regressions 

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
global$mean_CDD18_db <- global$meanpy_CDD18_db
global$mean_hDD18_db <- global$meanpy_hDD18_db
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
global <- global[complete.cases(global$urban), ]
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

# Survey
global_svy <- svydesign(data = global, ids = ~ adm1, weights = ~ weight)


##################################

#        Extensive margin        #

##################################

# AC formula for global
ac_formula <- as.numeric(as.character(ac)) ~ mean_CDD18_db + mean_CDD18_db2 +
  ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | country

# Linear probability model
reg_ac00 <- feols(ac_formula, data = global, weights = ~weight, cluster = c("adm1")); summary(reg_ac00)

# AC formula for global
ac_formula <- as.numeric(as.character(ac)) ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | country

# Linear probability model
reg_ac0 <- feols(ac_formula, data = global, weights = ~weight, cluster = c("adm1")); summary(reg_ac0)

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save coefficient results
texreg(list(reg_ac00, reg_ac0, fs), digits = 3, caption = "Logit Regression for Air-conditioning Ownership --- Global",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("LPM", "LPM", "Logit"),
       custom.note = "\\textbf{Notes}: Dependent variable is air-conditioning (0,1). Clustered std. errors at the ADM-1 level in parentheses. Column (4) shows the average marginal effects (AMEs) from the logit regression. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$. Regressions are conducted using survey weights.", 
       file = paste(output,'airconditioning/main/Global_wgt_coeff.tex', sep=''), append=F,  
       float.pos = "H", label = "main: ac_global",
       omit.coef = "(country)|(Intercept)", 
       custom.coef.map = list("mean_CDD18_db" = "$\\overline{CDD}$", "mean_CDD18_db2" = "$\\overline{CDD}^2$", 
                              "mean_CDD18_db_exp" = "$\\overline{CDD}$ $\\times$ Log(Exp)", 
                              "mean_CDD18_db2_exp" = "$\\overline{CDD}^2$ $\\times$ Log(Exp)",
                              "curr_CDD18_db" = "$CDD$", "curr_CDD18_db2" = "$CDD^2$",
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
                              "sex_head1" = "Female (Yes = 1)"),
       custom.gof.rows = list("Country FE" = c("YES", "YES", "YES"), 
                              "Countries" = c("25", "25", "25")), caption.above = TRUE)

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + country

# Logistic regression
reg_ac <- svyglm(ac_formula, design = global_svy, 
                 family = binomial(logit), na.action=na.omit); summary(reg_ac)

# Save AME results
margins <- margins(reg_ac, design = global_svy)
summary(margins)
xtable(summary(margins), display=rep('g', 8), caption = "Logit Regression for Air-conditioning Ownership - Global")
print(xtable(summary(margins), display=rep('g', 8), caption = "Logit Regression for Air-conditioning Ownership - Global"), 
      file= paste(output,'airconditioning/main/Global_wgt_ame.tex', sep=''),append=F, table.placement = "htbp",
      caption.placement="top")

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

# Survey
global_svy <- svydesign(data = sec, ids = ~ adm1, weights = ~ weight)


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

# Formula electricity expenditure without selection
ely_formula  <- ln_ely_q ~ ac | 
  country

# Without selection
model00 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model00)


# Formula electricity expenditure without selection
ely_formula  <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head | country

# Without selection
model0 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0)


# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model1 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model1)


# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

# With selection
model2 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model2)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection + country

# With selection
model3 <- svyglm(ely_formula, design = global_svy, na.action=na.omit); summary(model3)
reg_ely <- model3


# Compare the models
screenreg(list(model00, model0, model1, model2), digits = 3, 
          caption = "The Effect of Air-conditioning on Residential Electricity Quantity",
          stars = c(0.1, 0.05, 0.01), custom.model.names = c("OLS", "OLS", "DMF", "DMF"),
          omit.coef = "(country)|(Intercept)|(selection)",
          custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac1:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", "curr_CDD18_db" = "CDD", 
                                 "I(curr_CDD18_db^2)" = "CDD$^2$", "curr_HDD18_db" = "HDD", "I(curr_HDD18_db^2)" = "HDD$^2$",
                                 "ln_total_exp_usd_2011" = "Log(Exp)",
                                 "ln_ely_p" = "Log(P)",
                                 "urban_sh" = "Urbanisation (\\%)", 
                                 "ownership_d1" = "House Ownership (Yes = 1)", "n_members" = "Household Size", 
                                 "edu_head_21" = "Primary Edu.", "edu_head_22" = "Secondary Edu.", "edu_head_23" = "Post Edu.", 
                                 "age_head" = "Age (Head)", "sex_head1" = "Female (Yes = 1)"),
          custom.gof.rows = list("Correction Term" = c("NO", "NO", "YES", "YES"), 
                                 "Country FE" = c("YES", "YES", "YES", "YES"), 
                                 "Countries" = c("25", "25", "25", "25")))

## Export
# Full table
texreg(list(model00, model0, model1, model2), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("OLS", "OLS", "DMF", "DMF"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/main/Global_wgt_ft.tex', sep=''), append=F,  
       float.pos = "htbp", label = "main: ely_global",
       omit.coef = "(country)|(Intercept)",
       custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac1:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", "curr_CDD18_db" = "CDD", 
                              "I(curr_CDD18_db^2)" = "CDD$^2$", "curr_HDD18_db" = "HDD", "I(curr_HDD18_db^2)" = "HDD$^2$",
                              "ln_total_exp_usd_2011" = "Log(Exp)",
                              "ln_ely_p" = "Log(P)",
                              "urban_sh" = "Urbanisation (\\%)", 
                              "ownership_d1" = "House Ownership (Yes = 1)", "n_members" = "Household Size", 
                              "edu_head_21" = "Primary Edu.", "edu_head_22" = "Secondary Edu.", "edu_head_23" = "Post Edu.", 
                              "age_head" = "Age (Head)", "sex_head1" = "Female (Yes = 1)",
                              "selection" = "$\\hat{\\zeta}$"),
       custom.gof.rows = list("Country FE" = c("YES", "YES", "YES", "YES"), 
                              "Countries" = c("25", "25", "25", "25")), 
       caption.above = TRUE)

# Only AC
texreg(list(model00, model0, model1, model2), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("OLS", "OLS", "DMF", "DMF"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/main/Global_wgt.tex', sep=''), append=F,  
       float.pos = "htbp", label = "main: ely_global",
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac1:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$"),
       custom.gof.rows = list("Controls" = c("NO", "YES", "YES", "YES"),
                              "Correction Term" = c("NO", "NO", "YES", "YES"), 
                              "Country FE" = c("YES", "YES", "YES", "YES"), 
                              "Countries" = c("25", "25", "25", "25")), 
       caption.above = TRUE)


# Save the R Environment will be used for the projections
save(list = c("reg_ac", "reg_ely", "sec"), 
     file = paste(output,'/for_projections/global_wgt_dmcf.RData', sep=''))

