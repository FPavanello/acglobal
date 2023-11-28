
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for Spain (EPIC 2011)
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

# Filter - Spain
HH_Spain <- HH_OECD %>% filter(country == "Spain")

# Log electricity quantity
HH_Spain$ln_ely_q <- log(HH_Spain$ely_q_impute_5_95)
HH_Spain$ln_ely_q[HH_Spain$ely_q_impute_5_95 == 0] <- NA

# Only those with not missing values 
HH_Spain <- HH_Spain[complete.cases(HH_Spain$ac), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$mean_CDD_db), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$mean_HDD_db), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$curr_CDD_db), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$curr_HDD_db), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$ln_total_exp_usd_2011), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$urban), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$n_members), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$share_under18), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$year_postschool), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$ownership_d), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$ely_q_impute_5_95), ] # 3648 observations
HH_Spain <- HH_Spain[complete.cases(HH_Spain$age_head), ]
HH_Spain <- HH_Spain[complete.cases(HH_Spain$sex_head), ]

# CDD in 100s
HH_Spain <- HH_Spain %>% mutate(mean_CDD_db = mean_CDD_db/100,
                                mean_HDD_db = mean_HDD_db/100,
                                curr_HDD_db = curr_HDD_db/100,
                                curr_CDD_db = curr_CDD_db/100)
# Scale variable
HH_Spain <- HH_Spain %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD_db)),
                                std_CDD = as.numeric(scale(curr_CDD_db)),
                                std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                std_HDD = as.numeric(scale(curr_HDD_db)),
                                std_n_members = as.numeric(scale(n_members)),
                                std_age_head = as.numeric(scale(age_head)),
                                std_sh_under16 = as.numeric(scale(share_under18)),
                                std_year_postschool = as.numeric(scale(year_postschool)))

# Arrange some vars
HH_Spain$urban <- as.factor(HH_Spain$urban)
HH_Spain$country2 <- as.factor(HH_Spain$country)
HH_Spain$ac <- as.factor(HH_Spain$ac)
HH_Spain$ownership_d <- as.factor(HH_Spain$ownership_d)
HH_Spain$sex_head <- as.factor(HH_Spain$sex_head)
HH_Spain$region <- as.factor(HH_Spain$region)

# AC formula for OECD
ac_formula_oecd <- ac ~ std_CDD_mean*std_texp + as.factor(urban) + std_n_members + std_sh_under16 + 
  as.factor(ownership_d) + std_year_postschool + std_age_head + as.factor(sex_head) + as.factor(region)

## Logistic regression of AC on covariates
# Spain
reg_ac <- glm(ac_formula_oecd, data = HH_Spain, family = binomial(logit), na.action=na.omit); summary(reg_ac)
coeftest(reg_ac, vcov = vcovCL(reg_ac, cluster = HH_Spain$hhid)) # robust standard errors
vcov <- vcovCL(reg_ac, cluster = HH_Spain$hhid) # clusterise SEs

# Save AME results
margins <- margins(reg_ac, vcov = vcov)
summary(margins)
ac_margins <- summary(margins)

# Predicted probabilities
HH_Spain$phat0_obs <- as.numeric(predict(reg_ac, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Spain$ac_obs <- ifelse(HH_Spain$phat0_obs>0.5 & !is.na(HH_Spain$phat0_obs), 1 , 0)
table(HH_Spain$ac, HH_Spain$ac_obs)

# Selection term
HH_Spain$xb_noac = 1-HH_Spain$phat0_obs               
HH_Spain$selection = ifelse(HH_Spain$ac==1, 
                            (HH_Spain$xb_noac*log(HH_Spain$xb_noac)/HH_Spain$phat0_obs) + log(HH_Spain$phat0_obs), 
                            (HH_Spain$phat0_obs*log(HH_Spain$phat0_obs)/HH_Spain$xb_noac) + log(HH_Spain$xb_noac))

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

# Formula without selection
ely_formula_oecd <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp +
  urban + std_n_members + std_sh_under16 + 
  ownership_d + std_year_postschool + std_age_head + sex_head + 
  region

# Without selection
model0 <- lm(ely_formula_oecd, data = HH_Spain, na.action=na.omit); summary(model0)
vcov0 <- cluster.vcov(model0, HH_Spain$hhid) # clusterise SEs
coeftest(model0, vcov=vcov0)
model0_cl <- coeftest(model0, vcov=vcov0) # with clustered SEs

# Formula with selection
ely_formula_oecd <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp +
  urban + std_n_members + std_sh_under16 + 
  ownership_d + std_year_postschool + std_age_head + sex_head + 
  region + selection

# With selection
model1 <- lm(ely_formula_oecd, data = HH_Spain, na.action=na.omit); summary(model1)
vcov1 <- cluster.vcov(model1, HH_Spain$hhid) # clusterise SEs
coeftest(model1, vcov=vcov1)
model1_cl <- coeftest(model1, vcov=vcov1) # with clustered SEs

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + ac*std_CDD + 
  std_CDD*std_texp + std_HDD*std_texp +
  urban + std_n_members + std_sh_under16 + 
  ownership_d + std_year_postschool + std_age_head + sex_head + 
  region + selection

# With selection
model2 <- lm(ely_formula_oecd, data = HH_Spain, na.action=na.omit); summary(model2)
vcov2 <- cluster.vcov(model2, HH_Spain$hhid) # clusterise SEs
coeftest(model2, vcov=vcov2)
model2_cl <- coeftest(model2, vcov=vcov2) # with clustered SEs

# Marginal effects
ely_margins <- summary(margins(model2, variables = c("ac", "std_CDD", "std_texp", "std_HDD",
                                                     "urban", "std_n_members", "std_sh_under16", 
                                                     "ownership_d", "std_year_postschool", "std_age_head", "sex_head", 
                                                     "selection"), vcov = vcov2))
mean(HH_Spain$ely_q_impute_5_95)*exp(mean(ely_margins$dydx_acYes)) - mean(HH_Spain$ely_q_impute_5_95) # average effect 1944.21 kWh

# Export
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/esp_dmcf.RData', sep=''))
