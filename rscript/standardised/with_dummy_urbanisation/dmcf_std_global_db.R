
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for the global data set
##      3) run intensive margin regressions: electricity expenditure on climate + covariates
##         using Dubin and McFadden (1984) approach

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
library(lm.beta)
library(reghelper)


# Set users
user <- 'fp'
#user <- 'gf'

if (user=='fp') {
  stub <- 'G:/Il mio Drive/'
}

if (user=='gf') {
  stub <- 'H:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

house <- paste(stub,'6-Projections/data/household/', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')

# Load global data
global <- readRDS(paste(house,'global.rds', sep=''))

# Scale variable
global <- global %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD_db)),
                            std_CDD = as.numeric(scale(curr_CDD_db)),
                            std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                            std_HDD = as.numeric(scale(curr_HDD_db)),
                            std_n_members = as.numeric(scale(n_members)),
                            std_age_head = as.numeric(scale(age_head)))


##################################

#        Extensive margin        #

##################################

# AC formula for global
ac_formula <- ac ~ std_CDD_mean*std_texp + std_n_members + edu_head_2 + std_age_head + country

# Logistic regression of AC on covariates
reg_ac <- glm(ac_formula, data = global, family = binomial(logit), na.action=na.omit); summary(reg_ac)
vcov <- vcovCL(reg_ac, cluster = global$hhid) # clusterise SEs
coeftest(reg_ac, vcov = vcovCL(reg_ac, cluster = global$hhid)) # clustered standard errors

# Save AME results
margins <- margins(reg_ac, vcov = vcov)
ac_margins <- summary(margins)
xtable(summary(margins), display=rep('g', 8), 
       caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Global")
print(xtable(summary(margins), display=rep('g', 8), 
             caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Global"), 
      file= paste(output,'airconditioning/standardised/Global.tex', sep=''),append=F, table.placement = "htbp",
      caption.placement="top")

# Predicted probabilities
global$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
global$ac_obs <- ifelse(global$phat0_obs>0.5 & !is.na(global$phat0_obs), 1 , 0)

table(global$ac, global$ac_obs)

# Selection term
global$xb_noac = 1-global$phat0_obs               
global$selection = ifelse(global$ac==1, 
                          (global$xb_noac*log(global$xb_noac)/global$phat0_obs) + log(global$phat0_obs), 
                          (global$phat0_obs*log(global$phat0_obs)/global$xb_noac) + log(global$xb_noac))


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
ely_formula  <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp + std_n_members + std_age_head +
                edu_head_2 + country

# Without selection
model0 <- lm(ely_formula, data = global, na.action=na.omit); summary(model0)
vcov0 <- cluster.vcov(model0, global$hhid) # clusterise SEs
coeftest(model0, vcov=vcov0)
model0_cl <- coeftest(model0, vcov=vcov0) # with clustered SEs

# Marginal effect of AC
ac_eff <- margins(model0, variables = "ac", vcov = vcov0)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(global$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(global$ely_q) # Average effect 282.56 kWh

# Formula electricity expenditure without interaction
ely_formula <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp + std_n_members + std_age_head +
               edu_head_2 + country + selection

# With selection
model1 <- lm(ely_formula, data = global, na.action=na.omit); summary(model1)
vcov1 <- cluster.vcov(model1, global$hhid) # clusterise SEs
coeftest(model1, vcov=vcov1)
model1_cl <- coeftest(model1, vcov=vcov1) # with clustered SEs

# Marginal effect of AC
ac_eff <- margins(model1, variables = "ac", vcov = vcov1)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(global$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(global$ely_q) # Average effect 322.95 kWh

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*std_CDD + std_CDD*std_texp + std_HDD*std_texp + std_n_members + std_age_head +
                           edu_head_2 + country + selection

# With selection
model2 <- lm(ely_formula, data = global, na.action=na.omit); summary(model2)
vcov2 <- cluster.vcov(model2, global$hhid) # clusterise SEs
coeftest(model2, vcov=vcov2)
model2_cl <- coeftest(model2, vcov=vcov2) # with clustered SEs

# Marginal effect of AC
ac_eff <- margins(model2, variables = "ac", vcov = vcov2)
ely_margins <- summary(margins(model2, vcov = vcov2))
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(global$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(global$ely_q) # Average effect 12.59 kWh


# Compare the models
screenreg(list(model0_cl, model1_cl, model2_cl), digits = 3, 
          caption = "Air-conditioning Impact on Electricity Quantity using Standardised Covariates - Global",
          stars = c(0.1, 0.05, 0.01), custom.model.names = c("No Correction Term", "Correction Term", "Interaction"),
          omit.coef = "(country)|(Intercept)", reorder.coef = c(1, 13, 2, 4, 3, 10, 11, 5, 6, 7, 8, 9, 12),
          custom.coef.names = c("AC", "CDD", "Log(Exp)", "HDD", "Household Size", "Primary Edu.", "Secondary Edu.", "Post Edu.",
                                "Age (Head)", "CDD $\\times$ Log(Exp)", "HDD $\\times$ Log(Exp)", "Correction Term", "AC $\\times$ CDD"),
          custom.gof.rows = list("Country FE" = c("YES", "YES", "YES")))

texreg(list(model0_cl, model1_cl, model2_cl), digits = 3, 
       caption = "Air-conditioning Impact on Electricity Quantity using Standardised Covariates - Global",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("No Correction Term", "Correction Term", "Interaction"),
       omit.coef = "(country)|(Intercept)", reorder.coef = c(1, 13, 2, 4, 3, 10, 11, 5, 6, 7, 8, 9, 12),
       custom.coef.names = c("AC", "CDD", "Log(Exp)", "HDD", "Household Size", "Primary Edu.", "Secondary Edu.", "Post Edu.",
                             "Age (Head)", "CDD $\\times$ Log(Exp)", "HDD $\\times$ Log(Exp)", "Correction Term", "AC $\\times$ CDD"),
       custom.gof.rows = list("Country FE" = c("YES", "YES", "YES")))

# Export
texreg(list(model0_cl, model1_cl, model2_cl), digits = 3, 
       caption = "Air-conditioning Impact on Electricity Quantity using Standardised Covariates - Global",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("No Correction Term", "Correction Term", "Interaction"),
       custom.note = "Clustered std. errors at the district level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/standardised/Global.tex', sep=''), append=F,  float.pos = "htbp", label = "main: ely_afr",
       omit.coef = "(country)|(Intercept)", reorder.coef = c(1, 13, 2, 4, 3, 10, 11, 5, 6, 7, 8, 9, 12),
       custom.coef.names = c("AC", "CDD", "Log(Exp)", "HDD", "Household Size", "Primary Edu.", "Secondary Edu.", "Post Edu.",
                             "Age (Head)", "CDD $\\times$ Log(Exp)", "HDD $\\times$ Log(Exp)", "Correction Term", "AC $\\times$ CDD"),
       custom.gof.rows = list("Country FE" = c("YES", "YES", "YES")), caption.above = TRUE)

# Export
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/global_dmcf.RData', sep=''))
