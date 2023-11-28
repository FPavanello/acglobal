
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for Italy using 2016 wave
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
user <- 'gf'

if (user=='fp') {
  stub <- 'G:/Il mio Drive/'
}

if (user=='gf') {
  stub <- 'H:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}


house <- paste(stub,'6-Projections/data/household/', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')

# Load Household data
HH_Italy <- readRDS(paste(house,'Italy/Italy_hbs.rds', sep=''))

# Arrange some vars
HH_Italy$state <- as.factor(HH_Italy$state)
HH_Italy$macroarea <- as.factor(HH_Italy$macroarea)
HH_Italy$ac <- as.factor(HH_Italy$ac)
HH_Italy$edu_head_2 <- as.factor(HH_Italy$edu_head_2)
HH_Italy$sex_head <- as.factor(HH_Italy$sex_head)

# Only those with not missing values 
HH_Italy <- HH_Italy[complete.cases(HH_Italy$ac), ]
HH_Italy <- HH_Italy[complete.cases(HH_Italy$mean_CDD_db), ]
HH_Italy <- HH_Italy[complete.cases(HH_Italy$mean_HDD_db), ]
HH_Italy <- HH_Italy[complete.cases(HH_Italy$curr_CDD_db), ]
HH_Italy <- HH_Italy[complete.cases(HH_Italy$curr_HDD_db), ]
HH_Italy <- HH_Italy[complete.cases(HH_Italy$ln_total_exp_usd_2011), ]
HH_Italy <- HH_Italy[complete.cases(HH_Italy$n_members), ]
HH_Italy <- HH_Italy[complete.cases(HH_Italy$edu_head_2), ]
HH_Italy <- HH_Italy[complete.cases(HH_Italy$age_head), ]
HH_Italy <- HH_Italy[complete.cases(HH_Italy$sex_head), ]

# CDD in 100s
HH_Italy <- HH_Italy %>% mutate(mean_CDD_db = mean_CDD_db/100,
                                mean_HDD_db = mean_HDD_db/100,
                                curr_CDD_db = curr_CDD_db/100,
                                curr_HDD_db = curr_HDD_db/100)

# Check
HH_Italy <- HH_Italy[(is.finite(HH_Italy$ln_ely_q) & !is.na(HH_Italy$ln_ely_q)),]

# Scale variable
HH_Italy <- HH_Italy %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD_db)),
                                  std_CDD = as.numeric(scale(curr_CDD_db)),
                                  std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                  std_HDD = as.numeric(scale(curr_HDD_db)),
                                std_n_members = as.numeric(scale(n_members)),
                                std_age_head = as.numeric(scale(age_head)))


##################################

#        Extensive margin        #

##################################

# AC formula for Italy - no ownership_id at the moment
ac_formula_ITA <- ac ~ std_CDD_mean*std_texp + std_n_members + 
  as.factor(edu_head_2) + std_age_head + as.factor(sex_head) +
  as.factor(macroarea)

# Logistic regression of AC on covariates
reg_ac <- glm(ac_formula_ITA, data = HH_Italy, family = binomial(logit), na.action=na.omit); summary(reg_ac)
vcov <- vcovCL(reg_ac, cluster = HH_Italy$hhid) # clusterise SEs
coeftest(reg_ac, vcov = vcovCL(reg_ac, cluster = HH_Italy$hhid)) # clustered standard errors

# Save AME results
margins <- margins(reg_ac, vcov = vcov)
summary(margins)
ac_margins <- summary(margins)
xtable(summary(margins), display=rep('g', 8), 
       caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Italy")
print(xtable(summary(margins), display=rep('g', 8), 
             caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Italy"), 
      file= paste(output,'airconditioning/standardised/ITA.tex', sep=''),append=F, table.placement = "htbp",
      caption.placement="top")

# Predicted probabilities
HH_Italy$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Italy$ac_obs <- ifelse(HH_Italy$phat0_obs>0.5 & !is.na(HH_Italy$phat0_obs), 1 , 0)

table(HH_Italy$ac, HH_Italy$ac_obs)

# Selection term for intensive margin part
HH_Italy$xb_noac = 1-HH_Italy$phat0_obs               
HH_Italy$selection = ifelse(HH_Italy$ac==1, 
                            (HH_Italy$xb_noac*log(HH_Italy$xb_noac)/HH_Italy$phat0_obs) + log(HH_Italy$phat0_obs), 
                            (HH_Italy$phat0_obs*log(HH_Italy$phat0_obs)/HH_Italy$xb_noac) + log(HH_Italy$xb_noac))


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
ely_formula_ITA <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp +
  std_n_members + edu_head_2 + std_age_head + sex_head + macroarea

# With selection
model0 <- lm(ely_formula_ITA, data = HH_Italy, na.action=na.omit); summary(model0)
vcov0 <- cluster.vcov(model0, HH_Italy$hhid) # clusterise SEs
coeftest(model0, vcov=vcov0)
model0_cl <- coeftest(model0, vcov=vcov0) # with clustered SEs

# Save betas
betas <- coef(model0)

# Marginal effect of AC
ac_eff <- margins(model0, variables = c("ac", "curr_CDD_db"), vcov = vcov0)
summary(ac_eff)
mean(ac_eff$dydx_ac1)
mean(HH_Italy$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Italy$ely_q) # Average effect 1221.62 kWh

# Formula electricity expenditure without interactions
ely_formula_ITA <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp +
  std_n_members + edu_head_2 + std_age_head + sex_head + macroarea + 
  selection

# With selection
model1 <- lm(ely_formula_ITA, data = HH_Italy, na.action=na.omit); summary(model1)
vcov1 <- cluster.vcov(model1, HH_Italy$hhid) # clusterise SEs
coeftest(model1, vcov=vcov1)
model1_cl <- coeftest(model1, vcov=vcov1) # with clustered SEs

# Save betas
betas <- coef(model1)

# Marginal effect of AC
ac_eff <- margins(model1, variables = "ac", vcov = vcov1)
summary(ac_eff)
mean(ac_eff$dydx_ac1)
mean(HH_Italy$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Italy$ely_q) # Average effect 1213.49 kWh

# Formula electricity expenditure with interactions
ely_formula_ITA <- ln_ely_q ~ ac + ac*std_CDD + 
  std_CDD*std_texp + std_HDD*std_texp +
  std_n_members + edu_head_2 + std_age_head + sex_head + macroarea + selection

# With selection
model2 <- lm(ely_formula_ITA, data = HH_Italy, na.action=na.omit); summary(model2)
vcov2 <- cluster.vcov(model2, HH_Italy$hhid) # clusterise SEs
coeftest(model2, vcov=vcov2)
model2_cl <- coeftest(model2, vcov=vcov2) # with clustered SEs

# Save betas
betas <- coef(model2)

# Marginal effect of AC
ac_eff <- margins(model2, variables = "ac", vcov = vcov2)
summary(ac_eff)
ely_margins <- summary(margins(model2, vcov = vcov2))
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Italy$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Italy$ely_q) # Average effect 1046.00 kWh

# Export
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/ita_dmcf.RData', sep=''))
