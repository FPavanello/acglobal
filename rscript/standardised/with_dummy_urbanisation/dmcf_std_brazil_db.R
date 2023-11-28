
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for Brazil using 2016 wave
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


house <- paste(stub,'6-Projections/data/household/Fourcountries', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')

# Load Household data
HH_Brazil <- readRDS(paste(house,'/bra_pof.rds', sep=''))

# Arrange some vars
HH_Brazil$state <- as.factor(HH_Brazil$state)
HH_Brazil$urban <- as.factor(HH_Brazil$urban)
HH_Brazil$ac <- as.factor(HH_Brazil$ac)
HH_Brazil$fan <- as.factor(HH_Brazil$fan)
HH_Brazil$refrigerator <- as.factor(HH_Brazil$refrigerator)
HH_Brazil$ownership <- as.factor(HH_Brazil$ownership)
HH_Brazil$edu_head_2 <- as.factor(HH_Brazil$edu_head_2)
HH_Brazil$sex_head <- as.factor(HH_Brazil$sex_head)
HH_Brazil$housing_index_lab <- as.factor(HH_Brazil$housing_index_lab)
HH_Brazil$occupation_head <- as.factor(HH_Brazil$occupation_head)
HH_Brazil$region3 <- as.factor(HH_Brazil$region3)

# log total exp
HH_Brazil$ln_total_exp_usd_2011 <- log(HH_Brazil$total_exp_usd_2011)
HH_Brazil$ln_total_exp_usd_2011[HH_Brazil$total_exp_usd_2011 == 0] <- NA

# Log electricity quantity
HH_Brazil$ln_ely_q <- log(HH_Brazil$ely_q)

# Only those with not missing values 
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$ac), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$mean_CDD_db), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$mean_HDD_db), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$curr_CDD_db), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$curr_HDD_db), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$ln_total_exp_usd_2011), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$urban), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$n_members), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$sh_under16), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$housing_index_lab), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$ownership_d), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$edu_head_2), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$age_head), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$sex_head), ]

# Only who has electricity access
HH_Brazil <- HH_Brazil %>% filter(ely_access == 1)

# Keep only urban stratums
#HH_Brazil <- HH_Brazil %>% filter(stratum_groups=="capital" | stratum_groups=="other_urban") 

# CDD in 100s
HH_Brazil <- HH_Brazil %>% mutate(mean_CDD_db = mean_CDD_db/100)
HH_Brazil <- HH_Brazil %>% mutate(curr_CDD_db = curr_CDD_db/100)
HH_Brazil <- HH_Brazil %>% mutate(mean_HDD_db = mean_HDD_db/100)
HH_Brazil <- HH_Brazil %>% mutate(curr_HDD_db = curr_HDD_db/100)

# Check
HH_Brazil <- HH_Brazil[(is.finite(HH_Brazil$ln_ely_q) & !is.na(HH_Brazil$ln_ely_q)),]

# Filter variables and rows
HH_Brazil <- dplyr::select(HH_Brazil, hhid, ln_ely_q, ely_q, ac, curr_CDD_db, ln_total_exp_usd_2011, curr_HDD_db, ln_total_exp_usd_2011, 
                           urban, n_members, sh_under16, ownership_d, edu_head_2,
                           housing_index_lab, housing_index, age_head , sex_head, state, mean_CDD_db, mean_HDD_db, 
                           state_district, state3, region3, stratum_groups, state_code, weight)

# Scale variable
HH_Brazil <- HH_Brazil %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD_db)),
                                  std_CDD = as.numeric(scale(curr_CDD_db)),
                                  std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                  std_HDD = as.numeric(scale(curr_HDD_db)),
                                  std_n_members = as.numeric(scale(n_members)),
                                  std_age_head = as.numeric(scale(age_head)),
                                  std_sh_under16 = as.numeric(scale(sh_under16)))


##################################

#        Extensive margin        #

##################################

# AC formula for Brazil
ac_formula_bra <- ac ~ std_CDD_mean*std_texp + as.factor(urban) + std_n_members + std_sh_under16 + as.factor(ownership_d) + 
  as.factor(edu_head_2) + std_age_head + as.factor(sex_head) + as.factor(housing_index_lab) + 
  as.factor(region3)

# Logistic regression of AC on covariates
reg_ac <- glm(ac_formula_bra, data = HH_Brazil, family = binomial(logit), na.action=na.omit); summary(reg_ac)
vcov <- vcovCL(reg_ac, cluster = HH_Brazil$hhid) # clusterise SEs
coeftest(reg_ac, vcov = vcovCL(reg_ac, cluster = HH_Brazil$hhid)) # clustered standard errors

# Save AME results
margins <- margins(reg_ac, vcov = vcov)
summary(margins)
ac_margins <- summary(margins)
xtable(summary(margins), display=rep('g', 8), 
       caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Brazil")
print(xtable(summary(margins), display=rep('g', 8), 
             caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Brazil"), 
      file= paste(output,'airconditioning/standardised/BRA.tex', sep=''),append=F, table.placement = "htbp",
      caption.placement="top")

# Predicted probabilities
HH_Brazil$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Brazil$ac_obs <- ifelse(HH_Brazil$phat0_obs>0.5 & !is.na(HH_Brazil$phat0_obs), 1 , 0)

table(HH_Brazil$ac, HH_Brazil$ac_obs)

# Selection term
HH_Brazil$xb_noac = 1-HH_Brazil$phat0_obs               
HH_Brazil$selection = ifelse(HH_Brazil$ac==1, 
                             (HH_Brazil$xb_noac*log(HH_Brazil$xb_noac)/HH_Brazil$phat0_obs) + log(HH_Brazil$phat0_obs), 
                             (HH_Brazil$phat0_obs*log(HH_Brazil$phat0_obs)/HH_Brazil$xb_noac) + log(HH_Brazil$xb_noac))


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

# Formula electricity quantity without interactions
ely_formula_bra <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp + 
  urban + std_n_members + std_sh_under16 + ownership_d + edu_head_2 + 
  housing_index_lab + std_age_head + sex_head + region3 

# Without selection
model0 <- lm(ely_formula_bra, data = HH_Brazil, na.action=na.omit); summary(model0)
vcov0 <- cluster.vcov(model0, HH_Brazil$hhid) # clusterise SEs
coeftest(model0, vcov=vcov0) # with clustered SEs
model0_cl <- coeftest(model0, vcov=vcov0) # with clustered SEs

# Save betas
betas <- coef(model0)

# Marginal effect of AC
ac_eff <- margins(model0, variables = "ac", vcov = vcov0)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Brazil$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Brazil$ely_q) # Average effect 820.51 kWh

# Formula electricity quantity without interactions
ely_formula_bra <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp + 
  urban + std_n_members + std_sh_under16 + ownership_d + edu_head_2 + 
  housing_index_lab + std_age_head + sex_head + region3 + selection

# With selection
model1 <- lm(ely_formula_bra, data = HH_Brazil, na.action=na.omit); summary(model1)
vcov1 <- cluster.vcov(model1, HH_Brazil$hhid) # clusterise SEs
coeftest(model1, vcov=vcov1) # with clustered SEs
model1_cl <- coeftest(model1, vcov=vcov1) # with clustered SEs

# Save betas
betas <- coef(model1)

# Marginal effect of AC
ac_eff <- margins(model1, variables = "ac", vcov = vcov1)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Brazil$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Brazil$ely_q) # Average effect 933.20 kWh

# Formula electricity quantity with interactions
ely_formula_bra <- ln_ely_q ~ ac + ac*std_CDD +
  std_CDD*std_texp + std_HDD*std_texp + 
  urban + std_n_members + std_sh_under16 + ownership_d + edu_head_2 + 
  housing_index_lab + std_age_head + sex_head + region3 + 
  selection

# With selection
model2 <- lm(ely_formula_bra, data = HH_Brazil, na.action=na.omit); summary(model2)
vcov2 <- cluster.vcov(model2, HH_Brazil$hhid) # clusterise SEs
coeftest(model2, vcov=vcov2) # with clustered SEs
model2_cl <- coeftest(model2, vcov=vcov2) # with clustered SEs

# Save betas
betas <- coef(model2)

# Marginal effect of AC
ac_eff <- margins(model2,  variables = c("ac", "std_CDD"), vcov = vcov2)
summary(ac_eff)
ely_margins <- summary(margins(model2, vcov = vcov2))
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Brazil$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Brazil$ely_q) # Average effect 931.31 kWh

# Export
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/bra_dmcf.RData', sep=''))
