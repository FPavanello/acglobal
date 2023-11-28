
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for Africa using 2016 wave
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
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}


house <- paste(stub,'6-Projections/data/household/', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')

# Load Household data
HH_Nigeria <- readRDS(paste(house,'Nigeria/nga/Nigeria_ghs.rds', sep=''))
HH_Tanzania <- readRDS(paste(house,'Tanzania/HBS/Tanzania_hbs.rds', sep=''))
HH_Kenya <- readRDS(paste(house,'Kenya/IHBS/kenya_ihbs.rds', sep=''))
HH_Ghana <- readRDS(paste(house,'Ghana/ghana.rds', sep=''))

HH_Ghana$sex_head <- ifelse(HH_Ghana$sex_head==1, 0, 1)

# Macroarea
HH_Kenya$macroarea <- 1 # Coast
HH_Kenya$macroarea[HH_Kenya$state == "Garissa" | HH_Kenya$state == "Wajir" | HH_Kenya$state == "Mandera"] <- 2 # North-Eastern
HH_Kenya$macroarea[HH_Kenya$state == "Marsabit" | HH_Kenya$state == "Isiolo" | HH_Kenya$state == "Meru" |
                     HH_Kenya$state == "Tharaka-Nithi" | HH_Kenya$state == "Embu" | HH_Kenya$state == "Kitui" |  
                     HH_Kenya$state == "Machakos" | HH_Kenya$state == "Makueni"] <- 3 # Eastern
HH_Kenya$macroarea[HH_Kenya$state == "Nyandarua" | HH_Kenya$state == "Nyeri" | HH_Kenya$state == "Kirinyaga" |
                     HH_Kenya$state == "Murang'a" | HH_Kenya$state == "Kiambu"] <- 3 # Central
HH_Kenya$macroarea[HH_Kenya$state == "Turkana" | HH_Kenya$state == "West Pokot" | HH_Kenya$state == "Samburu" |
                     HH_Kenya$state == "Trans Nzoia" | HH_Kenya$state == "Uasin Gishu" | HH_Kenya$state == "Samburu" |
                     HH_Kenya$state == "Elgeyo-Marakwet" | HH_Kenya$state == "Nandi" | HH_Kenya$state == "Baringo" |
                     HH_Kenya$state == "Laikipia" | HH_Kenya$state == "Nakuru" | HH_Kenya$state == "Narok" |
                     HH_Kenya$state == "Kajiado" | HH_Kenya$state == "Kericho" | HH_Kenya$state == "Bomet"] <- 4 # Rift Valley
HH_Kenya$macroarea[HH_Kenya$state == "Kakamega" | HH_Kenya$state == "Vihiga" | HH_Kenya$state == "Bungoma" |
                     HH_Kenya$state == "Busia"] <- 5 # Western
HH_Kenya$macroarea[HH_Kenya$state == "Siaya" | HH_Kenya$state == "Kisumu" | HH_Kenya$state == "Homa Bay" |
                     HH_Kenya$state == "Migori" | HH_Kenya$state == "Kisii" | HH_Kenya$state == "Nyamira"] <- 6 # Nyanza
HH_Kenya$macroarea[HH_Kenya$state == "Nairobi"] <- 7 # Nairobi

HH_Tanzania$macroarea <- 1 # North
HH_Tanzania$macroarea[HH_Tanzania$state == "Taafr" | HH_Tanzania$state == "Manyara" | HH_Tanzania$state == "Singida" |
                        HH_Tanzania$state == "Singida" | HH_Tanzania$state == "Tabora" | HH_Tanzania$state == "Kilimanjaro" |
                        HH_Tanzania$state == "Katavi" | HH_Tanzania$state == "Rukwa" | HH_Tanzania$state == "Songwe" |
                        HH_Tanzania$state == "Mbeya" | HH_Tanzania$state == "Dodoma" | HH_Tanzania$state == "Taafr" |
                        HH_Tanzania$state == "Iriafr"] <- 2
HH_Tanzania$macroarea[HH_Tanzania$state == "Njombe" | HH_Tanzania$state == "Morogoro" | HH_Tanzania$state == "Pwani" |
                        HH_Tanzania$state == "Lindi" | HH_Tanzania$state == "Ruvuma" | HH_Tanzania$state == "Dar es Salaam" |
                        HH_Tanzania$state == "Mtwara"] <- 3

# HH_Ghana$macroarea[HH_Ghana$state == "Western" | HH_Ghana$state == "Central" | HH_Ghana$state == "Greater Accra"] <- 1 
# HH_Ghana$macroarea[HH_Ghana$state == "Volta" | HH_Ghana$state == "Eastern" | HH_Ghana$state == "Ashanti"] <- 2
# HH_Ghana$macroarea[HH_Ghana$state == "Brong Ahafo" | HH_Ghana$state == "Northern" | HH_Ghana$state == "Upper East" | HH_Ghana$state == "Upper West"] <- 4 

# Country
HH_Nigeria <- HH_Nigeria %>% mutate(country = "Nigeria")
HH_Tanzania <- HH_Tanzania %>% mutate(country ="Tanzania")
HH_Kenya <- HH_Kenya %>% mutate(country = "Kenya")
HH_Ghana <- HH_Ghana %>% mutate(country = "Ghana")

# Append
HH_Africa <- rbind.fill(HH_Nigeria, HH_Tanzania)
HH_Africa <- rbind.fill(HH_Africa, HH_Kenya)
HH_Africa <- rbind.fill(HH_Africa, HH_Ghana)

# Arrange some vars
HH_Africa$state <- as.factor(HH_Africa$state)
HH_Africa$urban <- as.factor(as.character(HH_Africa$urban))
HH_Africa$ac <- as.factor(as.character(HH_Africa$ac))
HH_Africa$ownership <- as.factor(as.character(HH_Africa$ownership))
HH_Africa$edu_head_2 <- as.factor(as.character(HH_Africa$edu_head_2))
HH_Africa$sex_head <- as.factor(as.character(HH_Africa$sex_head))
HH_Africa$housing_index_lab <- as.factor(as.character(HH_Africa$housing_index_lab))

# log total exp
HH_Africa$ln_total_exp_usd_2011[HH_Africa$total_exp_usd_2011 == 0] <- NA
HH_Africa$ln_ely_q[HH_Africa$ln_ely_q < 0] <- NA

# Only those with not missing values 
HH_Africa <- HH_Africa[complete.cases(HH_Africa$ac), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$mean_CDD_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$mean_HDD_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$curr_CDD_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$curr_HDD_db), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$urban), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$n_members), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$sh_under16), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$housing_index_lab), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$ownership_d), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$edu_head_2), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$age_head), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$sex_head), ]
HH_Africa <- HH_Africa[complete.cases(HH_Africa$country), ]

HH_Africa$mean_CDD_db <- as.numeric(HH_Africa$mean_CDD_db)

# CDD in 100s
HH_Africa <- HH_Africa %>% mutate(mean_CDD_db = mean_CDD_db/100,
                                  mean_HDD_db = mean_HDD_db/100,
                                  curr_CDD_db = curr_CDD_db/100,
                                  curr_HDD_db = curr_HDD_db/100)

# Only HH with electricity access
#HH_Africa <- HH_Africa %>% filter(ely_access == 1)

# Check
HH_Africa <- HH_Africa[(is.finite(HH_Africa$ln_ely_q) & !is.na(HH_Africa$ln_ely_q)),]

# Scale variable
HH_Africa <- HH_Africa %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD_db)),
                                  std_CDD = as.numeric(scale(curr_CDD_db)),
                                  std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                  std_HDD = as.numeric(scale(curr_HDD_db)),
                                  std_n_members = as.numeric(scale(n_members)),
                                  std_age_head = as.numeric(scale(age_head)),
                                  std_sh_under16 = as.numeric(scale(sh_under16)))


##################################

#        Extensive margin        #

##################################

HH_Africa <- dplyr::select(HH_Africa, ac, ln_ely_q, std_CDD_mean, std_texp, urban, std_n_members, 
                           std_sh_under16, ownership_d, edu_head_2, std_age_head, sex_head, housing_index_lab, country, 
                           std_CDD, std_HDD, hhid)
HH_Africa <- na.omit(HH_Africa)

# AC formula for Africa
ac_formula_afr <- ac ~ std_CDD_mean*std_texp + as.factor(urban) + std_n_members + std_sh_under16 + 
  as.factor(ownership_d) + as.factor(edu_head_2) + std_age_head + as.factor(sex_head) + 
  as.factor(housing_index_lab) + as.factor(country)

# Logistic regression of AC on covariates
reg_ac <- glm(ac_formula_afr, data = HH_Africa, family = binomial(logit), na.action=na.omit); summary(reg_ac)
vcov <- vcovCL(reg_ac, cluster = HH_Africa$hhid) # clusterise SEs
coeftest(reg_ac, vcov = vcovCL(reg_ac, cluster = HH_Africa$hhid)) # clustered standard errors

# Save AME results
margins <- margins(reg_ac, vcov = vcov)
summary(margins)
ac_margins <- summary(margins)
xtable(summary(margins), display=rep('g', 8), 
       caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Africa")
print(xtable(summary(margins), display=rep('g', 8), 
             caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Africa"), 
      file= paste(output,'airconditioning/standardised/AFR.tex', sep=''),append=F, table.placement = "htbp",
      caption.placement="top")


# Predicted probabilities
HH_Africa$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Africa$ac_obs <- ifelse(HH_Africa$phat0_obs>0.5 & !is.na(HH_Africa$phat0_obs), 1 , 0)

table(HH_Africa$ac, HH_Africa$ac_obs)

# Selection term
HH_Africa$xb_noac = 1-HH_Africa$phat0_obs               
HH_Africa$selection = ifelse(HH_Africa$ac==1, 
                             (HH_Africa$xb_noac*log(HH_Africa$xb_noac)/HH_Africa$phat0_obs) + log(HH_Africa$phat0_obs), 
                             (HH_Africa$phat0_obs*log(HH_Africa$phat0_obs)/HH_Africa$xb_noac) + log(HH_Africa$xb_noac))


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
ely_formula_afr  <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp + 
  urban + std_n_members + std_sh_under16 + ownership_d + edu_head_2 + 
  housing_index_lab + std_age_head + sex_head + country

# Without selection
model0 <- lm(ely_formula_afr, data = HH_Africa, na.action=na.omit); summary(model0)
vcov0 <- cluster.vcov(model0, HH_Africa$hhid) # clusterise SEs
coeftest(model0, vcov=vcov0)
model0_cl <- coeftest(model0, vcov=vcov0) # with clustered SEs

# Save betas
betas <- coef(model0)

# Marginal effect of AC
ac_eff <- margins(model0, variables = "ac", vcov = vcov0)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Africa$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Africa$ely_q) # Average effect 282.56 kWh

# Formula electricity expenditure without interaction
ely_formula_afr  <- ln_ely_q ~ ac + std_CDD*std_texp + std_HDD*std_texp + 
  urban + std_n_members + std_sh_under16 + ownership_d + edu_head_2 + 
  housing_index_lab + std_age_head + sex_head + country + selection

# With selection
model1 <- lm(ely_formula_afr, data = HH_Africa, na.action=na.omit); summary(model1)
vcov1 <- cluster.vcov(model1, HH_Africa$hhid) # clusterise SEs
coeftest(model1, vcov=vcov1)
model1_cl <- coeftest(model1, vcov=vcov1) # with clustered SEs

# Save betas
betas <- coef(model1)

# Marginal effect of AC
ac_eff <- margins(model1, variables = "ac", vcov = vcov1)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Africa$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Africa$ely_q) # Average effect 322.95 kWh

# Formula electricity expenditure for electricity expenditure
ely_formula_afr  <- ln_ely_q ~ ac + ac*std_CDD + 
  std_CDD*std_texp + std_HDD*std_texp + 
  urban + std_n_members + std_sh_under16 + ownership_d + edu_head_2 + 
  housing_index_lab + std_age_head + sex_head + country + selection

# With selection
model2 <- lm(ely_formula_afr, data = HH_Africa, na.action=na.omit); summary(model2)
vcov2 <- cluster.vcov(model2, HH_Africa$hhid) # clusterise SEs
coeftest(model2, vcov=vcov2)
model2_cl <- coeftest(model2, vcov=vcov2) # with clustered SEs

# Save betas
betas <- coef(model2)

# Marginal effect of AC
ac_eff <- margins(model2, variables = "ac", vcov = vcov2)
summary(ac_eff)
ely_margins <- summary(margins(model2, vcov = vcov2))
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Africa$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Africa$ely_q) # Average effect 12.59 kWh


# Save the R Environment will be used for the projections
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/afr_dmcf.RData', sep=''))
