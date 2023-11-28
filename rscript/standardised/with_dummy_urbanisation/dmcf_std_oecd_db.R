
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for OECD using 2011 wave
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
HH_OECD <- readRDS(paste(house,'OECD/EPIC/Household_OECD_mod.rds', sep=''))

# Arrange some vars
HH_OECD$urban <- as.factor(HH_OECD$urban)
HH_OECD$country2 <- as.factor(HH_OECD$country)
HH_OECD$ac <- as.factor(HH_OECD$ac)
HH_OECD$ownership_d <- as.factor(HH_OECD$ownership_d)
HH_OECD$sex_head <- as.factor(HH_OECD$sex_head)

# Log electricity quantity
HH_OECD$ln_ely_q <- log(HH_OECD$ely_q_impute_5_95)
HH_OECD$ln_ely_q[HH_OECD$ely_q_impute_5_95 == 0] <- NA

# Only those with not missing values 
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ac), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$mean_CDD_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$mean_HDD_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$curr_CDD_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$curr_HDD_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ln_total_exp_usd_2011), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$urban), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$n_members), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$share_under18), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$year_postschool), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ownership_d), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ely_q_impute_5_95), ] # 3648 observations
HH_OECD <- HH_OECD[complete.cases(HH_OECD$age_head), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$sex_head), ]

# CDD in 100s
HH_OECD <- HH_OECD %>% mutate(mean_CDD_db = mean_CDD_db/100,
                              mean_HDD_db = mean_HDD_db/100,
                              curr_HDD_db = curr_HDD_db/100,
                              curr_CDD_db = curr_CDD_db/100)
# Scale variable
HH_OECD <- HH_OECD %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD_db)),
                                  std_CDD = as.numeric(scale(curr_CDD_db)),
                                  std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                  std_HDD = as.numeric(scale(curr_HDD_db)),
                                  std_n_members = as.numeric(scale(n_members)),
                                  std_age_head = as.numeric(scale(age_head)),
                                  std_sh_under16 = as.numeric(scale(share_under18)),
                                  std_year_postschool = as.numeric(scale(year_postschool)))


### Instead of pooling the countries: same regression but using clusters of country
# AC formula for OECD
ac_formula_oecd <- ac ~ std_CDD_mean*std_texp + as.factor(urban) + std_n_members + std_sh_under16 + 
  as.factor(ownership_d) + std_year_postschool + std_age_head + as.factor(sex_head) + as.factor(country)

# Choosing cluster
HH_Europe <- HH_OECD %>% filter(country == "Sweden" | country == "Spain" | country == "Netherlands" | country == "France")
HH_NonEurope <- HH_OECD %>% filter(country == "Canada" | country == "Australia" | country == "Japan")

## Logistic regression of AC on covariates
# Europe
reg_ac_eu <- glm(ac_formula_oecd, data = HH_Europe, family = binomial(logit), na.action=na.omit); summary(reg_ac)
coeftest(reg_ac_eu, vcov = vcovCL(reg_ac_eu, cluster = HH_Europe$hhid)) # robust standard errors
vcov <- vcovCL(reg_ac_eu, cluster = HH_Europe$hhid) # clusterise SEs

# Save AME results
margins <- margins(reg_ac_eu, vcov = vcov)
summary(margins)
ac_margins_eu <- summary(margins)
xtable(summary(margins), display=rep('g', 8), 
       caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - OECD Europe")
print(xtable(summary(margins), display=rep('g', 8), 
             caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - OECD Europe"), 
      file= paste(output,'airconditioning/standardised/OECD_Europe.tex', sep=''),append=F, table.placement = "htbp",
      caption.placement="top")

# Predicted probabilities
HH_Europe$phat0_obs <- as.numeric(predict(reg_ac_eu, type="response"))

# Non Europe
reg_ac_noneu <- glm(ac_formula_oecd, data = HH_NonEurope, family = binomial(logit), na.action=na.omit); summary(reg_ac)
coeftest(reg_ac_noneu, vcov = vcovCL(reg_ac_noneu, cluster = HH_NonEurope$hhid)) # robust standard errors
vcov <- vcovCL(reg_ac_noneu, cluster = HH_NonEurope$hhid) # clusterise SEs

# Save AME results
margins <- margins(reg_ac_noneu, vcov = vcov)
summary(margins)
ac_margins_noneu <- summary(margins)
xtable(summary(margins), display=rep('g', 8), 
       caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - OECD Non-Europe")
print(xtable(summary(margins), display=rep('g', 8), 
             caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - OECD Non-Europe"), 
      file= paste(output,'airconditioning/standardised/OECD_NonEurope.tex', sep=''),append=F, table.placement = "htbp",
      caption.placement="top")

# Predicted probabilities
HH_NonEurope$phat0_obs <- as.numeric(predict(reg_ac_noneu, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Europe$ac_obs <- ifelse(HH_Europe$phat0_obs>0.5 & !is.na(HH_Europe$phat0_obs), 1 , 0)
table(HH_Europe$ac, HH_Europe$ac_obs)

HH_NonEurope$ac_obs <- ifelse(HH_NonEurope$phat0_obs>0.5 & !is.na(HH_NonEurope$phat0_obs), 1 , 0)
table(HH_NonEurope$ac, HH_NonEurope$ac_obs)


# Selection term
HH_Europe$xb_noac = 1-HH_Europe$phat0_obs               
HH_Europe$selection = ifelse(HH_Europe$ac==1, 
                             (HH_Europe$xb_noac*log(HH_Europe$xb_noac)/HH_Europe$phat0_obs) + log(HH_Europe$phat0_obs), 
                             (HH_Europe$phat0_obs*log(HH_Europe$phat0_obs)/HH_Europe$xb_noac) + log(HH_Europe$xb_noac))

HH_NonEurope$xb_noac = 1-HH_NonEurope$phat0_obs               
HH_NonEurope$selection = ifelse(HH_NonEurope$ac==1, 
                                (HH_NonEurope$xb_noac*log(HH_NonEurope$xb_noac)/HH_NonEurope$phat0_obs) + log(HH_NonEurope$phat0_obs), 
                                (HH_NonEurope$phat0_obs*log(HH_NonEurope$phat0_obs)/HH_NonEurope$xb_noac) + log(HH_NonEurope$xb_noac))


#################################################################

#     Intensive margin projections at 2040 + the role of AC

#     1) We interact AC with a set of variables to understand
#     how AC affects the adoption based on different charact.

#     2) We compute coefficients not only at the average, but
#     also based on specific values of our variables.
#     For instance, we compute the coefficients by decile, and
#     not only for the average household

#     Somehow point 1) is similar to a CDA, but without the
#     other appliances. For simplicity, I am going to interact
#     AC only with climate

#################################################################

### Europe vs NonEurope

## Europe
# Formula without selection
ely_formula_oecd <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp +
  urban + std_n_members + std_sh_under16 + 
  ownership_d + std_year_postschool + std_age_head + sex_head + 
  country2

# Without selection
model0 <- lm(ely_formula_oecd, data = HH_Europe, na.action=na.omit); summary(model0)
vcov0 <- cluster.vcov(model0, HH_Europe$hhid) # clusterise SEs
coeftest(model0, vcov=vcov0)
model0_cl <- coeftest(model0, vcov=vcov0) # with clustered SEs

# Save betas
betas <- coef(model0)

# Marginal effect of AC
ac_eff_eu <- margins(model0, variables = c("ac", "std_CDD"), vcov = vcov0)
summary(ac_eff_eu)
mean(ac_eff_eu$dydx_acYes) # same as #213
mean(HH_Europe$ely_q_impute_5_95)*exp(mean(ac_eff_eu$dydx_acYes)) - mean(HH_Europe$ely_q_impute_5_95) # average effect 1884.04 kWh


# Formula with selection
ely_formula_oecd <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp +
  urban + std_n_members + std_sh_under16 + 
  ownership_d + std_year_postschool + std_age_head + sex_head + 
  country2 + selection

# With selection
model1 <- lm(ely_formula_oecd, data = HH_Europe, na.action=na.omit); summary(model1)
vcov1 <- cluster.vcov(model1, HH_Europe$hhid) # clusterise SEs
coeftest(model1, vcov=vcov1)
model1_cl <- coeftest(model1, vcov=vcov1) # with clustered SEs

# Save betas
betas <- coef(model1)

# Marginal effect of AC
ac_eff_eu <- margins(model1, variables = "ac", vcov = vcov1)
summary(ac_eff_eu)
mean(ac_eff_eu$dydx_acYes) # same as #213
mean(HH_Europe$ely_q_impute_5_95)*exp(mean(ac_eff_eu$dydx_acYes)) - mean(HH_Europe$ely_q_impute_5_95) # average effect 1924.43 kWh


# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + ac*std_CDD + 
  std_CDD*std_texp + std_HDD*std_texp +
  urban + std_n_members + std_sh_under16 + 
  ownership_d + std_year_postschool + std_age_head + sex_head + 
  country2 + selection

# With selection
model2 <- lm(ely_formula_oecd, data = HH_Europe, na.action=na.omit); summary(model2)
vcov2 <- cluster.vcov(model2, HH_Europe$hhid) # clusterise SEs
coeftest(model2, vcov=vcov2)
model2_cl <- coeftest(model2, vcov=vcov2) # with clustered SEs

# Save betas
betas <- coef(model2)

# Marginal effect of AC
ac_eff_eu <- margins(model2, variables = "ac", vcov = vcov2)
summary(ac_eff_eu)
ely_margins_eu <- summary(margins(model2, variables = c("ac", "std_CDD", "std_texp", "std_HDD",
                                                        "urban", "std_n_members", "std_sh_under16", 
                                                        "ownership_d", "std_year_postschool", "std_age_head", "sex_head", 
                                                        "selection"), vcov = vcov2))
mean(ac_eff_eu$dydx_acYes) # same as #213
mean(HH_Europe$ely_q_impute_5_95)*exp(mean(ac_eff_eu$dydx_acYes)) - mean(HH_Europe$ely_q_impute_5_95) # average effect 1944.21 kWh


## Non-Europe
# Formula without selection
ely_formula_oecd <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp +
  urban + std_n_members + std_sh_under16 + 
  ownership_d + std_year_postschool + std_age_head + sex_head + 
  country2

# Without selection
model3 <- lm(ely_formula_oecd, data = HH_NonEurope, na.action=na.omit); summary(model3)
vcov3 <- cluster.vcov(model3, HH_NonEurope$hhid) # clusterise SEs
coeftest(model3, vcov=vcov3)
model3_cl <- coeftest(model3, vcov=vcov3) # with clustered SEs

# Save betas
betas <- coef(model3)

# Marginal effect of AC
ac_eff_neu <- margins(model3, variables = "ac", vcov = vcov3)
summary(ac_eff_neu)
mean(ac_eff_neu$dydx_acYes) # same as #213
mean(HH_NonEurope$ely_q_impute_5_95)*exp(mean(ac_eff_neu$dydx_acYes)) - mean(HH_NonEurope$ely_q_impute_5_95) # average effect 880.63 kWh


# Formula with selection
ely_formula_oecd <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp +
  urban + std_n_members + std_sh_under16 + 
  ownership_d + std_year_postschool + std_age_head + sex_head + 
  country2 + selection

# With selection
model4 <- lm(ely_formula_oecd, data = HH_NonEurope, na.action=na.omit); summary(model4)
vcov4 <- cluster.vcov(model4, HH_NonEurope$hhid) # clusterise SEs
coeftest(model4, vcov=vcov4) # with clustered SEs
model4_cl <- coeftest(model4, vcov=vcov4) # with clustered SEs

# Save betas
betas <- coef(model4)

# Marginal effect of AC
ac_eff_neu <- margins(model4, variables = "ac", vcov = vcov4)
summary(ac_eff_neu)
mean(ac_eff_neu$dydx_acYes) # same as #213
mean(HH_NonEurope$ely_q_impute_5_95)*exp(mean(ac_eff_neu$dydx_acYes)) - mean(HH_NonEurope$ely_q_impute_5_95) # average effect 877.31 kWh


# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + ac*std_CDD + 
  std_CDD*std_texp + std_HDD*std_texp +
  urban + std_n_members + std_sh_under16 + 
  ownership_d + std_year_postschool + std_age_head + sex_head + 
  country2 + selection

# With selection
model5 <- lm(ely_formula_oecd, data = HH_NonEurope, na.action=na.omit); summary(model5)
vcov5 <- cluster.vcov(model5, HH_NonEurope$hhid) # clusterise SEs
coeftest(model5, vcov=vcov5)
model5_cl <- coeftest(model5, vcov=vcov5) # with clustered SEs

# Save betas
betas <- coef(model5)

# Marginal effect of AC
ac_eff_neu <- margins(model5, variables = "ac", vcov = vcov5)
summary(ac_eff_neu)
ely_margins_noneu <- summary(margins(model5, 
                                     variables = c("ac", "std_CDD", "std_texp", "std_HDD",
                                                   "urban", "std_n_members", "std_sh_under16", 
                                                   "ownership_d", "std_year_postschool", "std_age_head", "sex_head", 
                                                   "selection"), vcov = vcov5))
mean(ac_eff_neu$dydx_acYes) # same as #213
mean(HH_NonEurope$ely_q_impute_5_95)*exp(mean(ac_eff_neu$dydx_acYes)) - mean(HH_NonEurope$ely_q_impute_5_95) # average effect 1090.13 kWh

# Export
ely_margins <- ely_margins_eu
ac_margins <- ac_margins_eu

save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/oecdeu_dmcf.RData', sep=''))

ely_margins <- ely_margins_noneu
ac_margins <- ac_margins_noneu

save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/oecdnoneu_dmcf.RData', sep=''))
