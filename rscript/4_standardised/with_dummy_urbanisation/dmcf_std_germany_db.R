
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for Germany using 2016 wave
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
HH_Germany <- readRDS(paste(house,'Germany/Germany.rds', sep=''))

# Only those with not missing values 
HH_Germany <- HH_Germany[complete.cases(HH_Germany$ac), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$mean_CDD_db), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$mean_HDD_db), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$curr_CDD_db), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$curr_HDD_db), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$ln_total_exp_usd_2011), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$urban), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$n_members), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$sh_under16), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$ownership_d), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$edu_head_2), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$age_head), ]
HH_Germany <- HH_Germany[complete.cases(HH_Germany$sex_head), ]

# CDD in 100s
HH_Germany <- HH_Germany %>% mutate(mean_CDD_db = mean_CDD_db/100,
                                    mean_HDD_db = mean_HDD_db/100,
                                    curr_CDD_db = curr_CDD_db/100,
                                    curr_HDD_db = curr_HDD_db/100)

# Macro-region
HH_Germany$macroarea <- ifelse((HH_Germany$state2==15 | HH_Germany$state2==8 | HH_Germany$state2==6| HH_Germany$state2==5| HH_Germany$state2==9| HH_Germany$state2==4 | HH_Germany$state2==3| HH_Germany$state2==14), 1, 
                               ifelse((HH_Germany$state2==10 | HH_Germany$state2==7 | HH_Germany$state2==16 | HH_Germany$state2==13), 2, 3))

# Arrange some vars
HH_Germany$state <- as.factor(HH_Germany$state)
HH_Germany$macroarea <- as.factor(HH_Germany$macroarea)
HH_Germany$urban <- as.factor(as.character(HH_Germany$urban))
HH_Germany$ac <- as.factor(as.character(HH_Germany$ac))
HH_Germany$ownership_d <- as.factor(as.character(HH_Germany$ownership_d))
HH_Germany$edu_head_2 <- as.factor(as.character(HH_Germany$edu_head_2))
HH_Germany$sex_head <- as.factor(as.character(HH_Germany$sex_head))

# Check
HH_Germany <- HH_Germany[(is.finite(HH_Germany$ln_ely_q) & !is.na(HH_Germany$ln_ely_q)),]

# Scale variable
HH_Germany <- HH_Germany %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD_db)),
                                    std_CDD = as.numeric(scale(curr_CDD_db)),
                                    std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                    std_HDD = as.numeric(scale(curr_HDD_db)),
                                    std_n_members = as.numeric(scale(n_members)),
                                    std_age_head = as.numeric(scale(age_head)),
                                    std_sh_under16 = as.numeric(scale(sh_under16)))


##################################

#        Extensive margin        #

##################################

# AC formula for Germany - no ownership_id at the moment
ac_formula_DEU <- ac ~ std_CDD_mean*std_texp + as.factor(urban) + std_n_members + 
  as.factor(edu_head_2) + std_age_head + as.factor(sex_head) + std_sh_under16 +
  as.factor(macroarea)

# Logistic regression of AC on covariates
reg_ac <- glm(ac_formula_DEU, data = HH_Germany, family = binomial(logit), na.action=na.omit); summary(reg_ac)
vcov <- vcovCL(reg_ac, cluster = HH_Germany$hhid) # clusterise SEs
coeftest(reg_ac, vcov = vcovCL(reg_ac, cluster = HH_Germany$hhid)) # clustered standard errors

# Save AME results
margins <- margins(reg_ac, vcov = vcov)
summary(margins)
ac_margins <- summary(margins)
xtable(summary(margins), display=rep('g', 8), 
       caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Germany")
print(xtable(summary(margins), display=rep('g', 8), 
             caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Germany"), 
      file= paste(output,'airconditioning/standardised/DEU.tex', sep=''),append=F, table.placement = "htbp",
      caption.placement="top")

# Predicted probabilities
HH_Germany$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Germany$ac_obs <- ifelse(HH_Germany$phat0_obs>0.5 & !is.na(HH_Germany$phat0_obs), 1 , 0)

table(HH_Germany$ac, HH_Germany$ac_obs)

# Selection term for intensive margin part
HH_Germany$xb_noac = 1-HH_Germany$phat0_obs               
HH_Germany$selection = ifelse(HH_Germany$ac==1, 
                              (HH_Germany$xb_noac*log(HH_Germany$xb_noac)/HH_Germany$phat0_obs) + log(HH_Germany$phat0_obs), 
                              (HH_Germany$phat0_obs*log(HH_Germany$phat0_obs)/HH_Germany$xb_noac) + log(HH_Germany$xb_noac))


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
ely_formula_DEU <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp + 
  urban + std_n_members + std_sh_under16 + edu_head_2 +
  std_age_head + sex_head + macroarea

# With selection
model0 <- lm(ely_formula_DEU, data = HH_Germany, na.action=na.omit); summary(model0)
vcov0 <- cluster.vcov(model0, HH_Germany$hhid) # clusterise SEs
coeftest(model0, vcov=vcov0)
model0_cl <- coeftest(model0, vcov=vcov0) # with clustered SEs

# Save betas
betas <- coef(model0)

# Marginal effect of AC
ac_eff <- margins(model0, variables = "ac", vcov = vcov0)
summary(ac_eff)
mean(ac_eff$dydx_ac1)
mean(HH_Germany$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Germany$ely_q) # Average effect 1221.62 kWh

# Formula electricity expenditure without interactions
ely_formula_DEU <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp + 
  urban + std_n_members + std_sh_under16 + edu_head_2 +
  std_age_head + sex_head + macroarea + selection

# With selection
model1 <- lm(ely_formula_DEU, data = HH_Germany, na.action=na.omit); summary(model1)
vcov1 <- cluster.vcov(model1, HH_Germany$hhid) # clusterise SEs
coeftest(model1, vcov=vcov1)
model1_cl <- coeftest(model1, vcov=vcov1) # with clustered SEs

# Save betas
betas <- coef(model1)

# Marginal effect of AC
ac_eff <- margins(model1, variables = "ac", vcov = vcov1)
summary(ac_eff)
mean(ac_eff$dydx_ac1)
mean(HH_Germany$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Germany$ely_q) # Average effect 1213.49 kWh

# Formula electricity expenditure with interactions
ely_formula_DEU <- ln_ely_q ~ ac + ac*std_CDD +
  std_CDD*std_texp + std_HDD*std_texp + 
  urban + std_n_members + std_sh_under16 + edu_head_2 +
  std_age_head + sex_head + macroarea + selection

# With selection
model2 <- lm(ely_formula_DEU, data = HH_Germany, na.action=na.omit); summary(model2)
vcov2 <- cluster.vcov(model2, HH_Germany$hhid) # clusterise SEs
coeftest(model2, vcov=vcov2)
model2_cl <- coeftest(model2, vcov=vcov2) # with clustered SEs

# Save betas
betas <- coef(model2)

# Marginal effect of AC
ac_eff <- margins(model2,  variables = c("ac", "std_CDD"), vcov = vcov2)
summary(ac_eff)
ely_margins <- summary(margins(model2, vcov = vcov2))
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Germany$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Germany$ely_q) # Average effect 931.31 kWh

# Export
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/deu_dmcf.RData', sep=''))
