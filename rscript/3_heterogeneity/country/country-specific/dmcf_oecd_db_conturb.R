
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts logit regressions for OECD countries using EPIC 2011 wave
##      3) run intensive margin regressions: electricity expenditure on climate + covariates
##         using Dubin and McFadden (1984) approach

rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Load package_heads
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
library(survey)
library(fixest)
library(tibble)

# Set directory
user <- 'fp'
#user <- 'gf'

if (user=='fp') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

if (user=='gf') {
  stub <- "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/"
}

house <- paste(stub,'6-Projections/data/household/', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')
script <- paste(stub,'6-Projections/rscripts/dmcf/regressions/country/with_continuous_urbanisation/', sep='')

# Load Household data
HH_OECD <- readRDS(paste(house,'OECD/EPIC/Household_OECD_mod.rds', sep=''))

# Add urbanisation share
source(paste0(stub, "6-Projections/rscripts/process_raw_data/add_urban/add_urban_oecd.R"))
HH_OECD$geometry <- NULL
HH_OECD$adm1 <- as.character(HH_OECD$state)

# Interaction prices
HH_OECD$mean_CDD18_db <- HH_OECD$meanpy_CDD18_db
HH_OECD$mean_hDD18_db <- HH_OECD$meanpy_hDD18_db
HH_OECD <- HH_OECD %>% mutate(ln_ely_p = log(ely_p_usd_2011),
                                  ln_ely_p_cdd = ln_ely_p*mean_CDD18_db,
                                  ln_ely_p_cdd2 = ln_ely_p*(mean_CDD18_db^2),
                                  ln_ely_p_own = ln_ely_p*as.numeric(as.character(ownership_d)),
                                  ln_ely_p_nme = ln_ely_p*n_members,
                                  mean_CDD18_db2 = mean_CDD18_db^2,
                                  mean_CDD18_db_exp = ln_total_exp_usd_2011*mean_CDD18_db,
                                  mean_CDD18_db2_exp = ln_total_exp_usd_2011*(mean_CDD18_db^2),
                                  curr_CDD18_db2 = curr_CDD18_db^2)

# Only those with not missing values 
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ac), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ln_ely_q), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$mean_CDD18_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$mean_HDD18_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$curr_CDD18_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$curr_HDD18_db), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ln_total_exp_usd_2011), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$urban_sh), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$n_members), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$share_under18), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$edu_head_2), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ownership_d), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ely_q_impute_5_95), ] # 3648 observations
HH_OECD <- HH_OECD[complete.cases(HH_OECD$age_head), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$sex_head), ]
HH_OECD <- HH_OECD[complete.cases(HH_OECD$ln_ely_p), ]
HH_OECD <- HH_OECD %>% filter(ln_ely_q > 0)
HH_OECD <- HH_OECD %>% filter(weight > 0)

# Choosing cluster
HH_Europe <- HH_OECD %>% filter(country == "Sweden" | country == "Spain" | country == "Netherlands" | country == "Switzerland " |
                                  country == "France")
HH_NonEurope <- HH_OECD %>% filter(country == "Canada" | country == "Australia" | country == "Japan")

HH_Europe$country2 <- as.factor(HH_Europe$country)
HH_NonEurope$country2 <- as.factor(HH_NonEurope$country)

# Survey
HH_Europe_svy <- svydesign(data = HH_Europe, ids = ~adm1, weights = ~ weight)
HH_NonEurope_svy <- svydesign(data = HH_NonEurope, ids = ~adm1, weights = ~ weight)


##################################

#        Extensive margin        #

##################################

### Instead of pooling the countries: same regression but using clusters of country
# AC formula for OECD
ac_formula_oecd <- as.numeric(as.character(ac)) ~ mean_CDD18_db + mean_CDD18_db2 + curr_CDD18_db + curr_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 +
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | country2

## Logistic regression of AC on covariates
# Europe
reg_ac_eu <- feglm(ac_formula_oecd, family = binomial(link = "logit"), 
                   data = HH_Europe, weights = ~weight, cluster = c("adm1"))

# Save AME results
#margins <- margins(reg_ac_eu, design = HH_Europe_svy)
#summary(margins)
#xtable(summary(margins), display=rep('g', 8), caption = "Logit Regression for Air-conditioning Ownership - OECD Europe")
#print(xtable(summary(margins), display=rep('g', 8), caption = "Logit Regression for Air-conditioning Ownership - OECD Europe"), 
#      file= paste(output,'airconditioning/main/OECD_Europe.tex', sep=''),append=F, table.placement = "htbp",
#      caption.placement="top")

# Predicted probabilities
HH_Europe$phat0_obs <- as.numeric(predict(reg_ac_eu, type="response"))

# Non Europe
reg_ac_noneu <- feglm(ac_formula_oecd, family = binomial(link = "logit"), 
                      data = HH_NonEurope, weights = ~weight, cluster = c("adm1"))
# Save AME results
#margins <- margins(reg_ac_noneu, design = HH_NonEurope_svy)
#summary(margins)
#xtable(summary(margins), display=rep('g', 8), caption = "Logit Regression for Air-conditioning Ownership - OECD Non-Europe")
#print(xtable(summary(margins), display=rep('g', 8), caption = "Logit Regression for Air-conditioning Ownership - OECD Non-Europe"), 
#      file= paste(output,'airconditioning/main/OECD_NonEurope.tex', sep=''),append=F, table.placement = "htbp",
#      caption.placement="top")

# Predicted probabilities
HH_NonEurope$phat0_obs <- as.numeric(predict(reg_ac_noneu, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Europe$ac_obs <- ifelse(HH_Europe$phat0_obs>0.5 & !is.na(HH_Europe$phat0_obs), 1 , 0)
HH_NonEurope$ac_obs <- ifelse(HH_NonEurope$phat0_obs>0.5 & !is.na(HH_NonEurope$phat0_obs), 1 , 0)

# Selection term
HH_Europe$xb_noac = 1-HH_Europe$phat0_obs               
HH_Europe$selection = ifelse(HH_Europe$ac==1, 
                           (HH_Europe$xb_noac*log(HH_Europe$xb_noac)/HH_Europe$phat0_obs) + log(HH_Europe$phat0_obs), 
                           (HH_Europe$phat0_obs*log(HH_Europe$phat0_obs)/HH_Europe$xb_noac) + log(HH_Europe$xb_noac))

HH_NonEurope$xb_noac = 1-HH_NonEurope$phat0_obs               
HH_NonEurope$selection = ifelse(HH_NonEurope$ac==1, 
                           (HH_NonEurope$xb_noac*log(HH_NonEurope$xb_noac)/HH_NonEurope$phat0_obs) + log(HH_NonEurope$phat0_obs), 
                           (HH_NonEurope$phat0_obs*log(HH_NonEurope$phat0_obs)/HH_NonEurope$xb_noac) + log(HH_NonEurope$xb_noac))

# Survey - re-run to add new variable
HH_Europe_svy <- svydesign(data = HH_Europe, ids = ~adm1, weights = ~ weight)
HH_NonEurope_svy <- svydesign(data = HH_NonEurope, ids = ~adm1, weights = ~ weight)


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
ely_formula_oecd <- ln_ely_q ~ ac +   
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) | country2

# Without selection
model0 <- feols(ely_formula_oecd, data = HH_Europe, weights = ~weight, cluster = c("adm1")); summary(model0)

# Formula with selection
ely_formula_oecd <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head | country2

# With selection
model1 <- feols(ely_formula_oecd, data = HH_Europe, weights = ~weight, cluster = c("adm1")); summary(model1)

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country2

# With selection
model2 <- feols(ely_formula_oecd, data = HH_Europe, weights = ~weight, cluster = c("adm1")); summary(model2)

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country2

# With selection
model3 <- feols(ely_formula_oecd, data = HH_Europe, weights = ~weight, cluster = c("adm1")); summary(model3)

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) +
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection + country2

# With selection
model_eu <- svyglm(ely_formula_oecd, design = HH_Europe_svy, na.action=na.omit); summary(model_eu)

# Marginal effect of AC
ac_eff_eu <- margins(model_eu, variables = "ac", design = HH_Europe_svy)
summary(ac_eff_eu)

## Non-Europe
# Formula without selection
ely_formula_oecd <- ln_ely_q ~ ac +   
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) | country2

# Without selection
model4 <- feols(ely_formula_oecd, data = HH_NonEurope, weights = ~weight, cluster = c("adm1")); summary(model4)

# Formula with selection
ely_formula_oecd <- ln_ely_q ~ ac +   
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head  | country2

# With selection
model5 <- feols(ely_formula_oecd, data = HH_NonEurope, weights = ~weight, cluster = c("adm1")); summary(model5)

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head  + selection | country2

# With selection
model6 <- feols(ely_formula_oecd, data = HH_NonEurope, weights = ~weight, cluster = c("adm1")); summary(model6)

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) +
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head  + selection | country2

# With selection
model7 <- feols(ely_formula_oecd, data = HH_NonEurope, weights = ~weight, cluster = c("adm1")); summary(model7)

# Formula with selection and interactions
ely_formula_oecd <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) +
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection + country2

# With selection
model_noneu <- svyglm(ely_formula_oecd, design = HH_NonEurope_svy, na.action=na.omit); summary(model_noneu)

# Marginal effect of AC
ac_eff_neu <- margins(model_noneu, variables = "ac", design = HH_NonEurope_svy)
summary(ac_eff_neu)


## Compare the models
# Europe
screenreg(list(model0, model1, model2, model3), digits = 3, 
          caption = "The Effect of Air-conditioning on Residential Electricity Quantity - OECD-Europe",
          stars = c(0.1, 0.05, 0.01), 
          custom.model.names = c("No Correction Term", "Correction Term", "Interaction", "Interaction"),
          omit.coef = "(country)|(Intercept)|(selection)", 
          custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac1:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", "curr_CDD18_db" = "CDD", 
                                 "I(curr_CDD18_db^2)" = "CDD$^2$", "curr_HDD18_db" = "HDD", "I(curr_HDD18_db^2)" = "HDD$^2$",
                                 "ln_total_exp_usd_2011" = "Log(Exp)", "urban_sh" = "Urbanisation (\\%)", 
                                 "ownership_d1" = "House Ownership (Yes = 1)", "n_members" = "Household Size", 
                                 "edu_head_22" = "Secondary Edu.", "edu_head_23" = "Post Edu.", 
                                 "age_head" = "Age (Head)", "sex_head1" = "Female (Yes = 1)"),
          custom.gof.rows = list("Correction Term" = c("NO", "NO", "YES", "YES"), 
                                 "Country FE" = c("YES", "YES", "YES", "YES"))) 


screenreg(list(model4, model5, model6, model7), digits = 3, 
          caption = "The Effect of Air-conditioning on Residential Electricity Quantity - OECD-Non Europe",
          stars = c(0.1, 0.05, 0.01), 
          custom.model.names = c("No Correction Term", "Correction Term", "Interaction", "Interaction"),
          omit.coef = "(country)|(Intercept)|(selection)", 
          custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac1:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", "curr_CDD18_db" = "CDD", 
                                 "I(curr_CDD18_db^2)" = "CDD$^2$", "curr_HDD18_db" = "HDD", "I(curr_HDD18_db^2)" = "HDD$^2$",
                                 "ln_total_exp_usd_2011" = "Log(Exp)", "urban_sh" = "Urbanisation (\\%)", 
                                 "ownership_d1" = "House Ownership (Yes = 1)", "n_members" = "Household Size", 
                                 "edu_head_22" = "Secondary Edu.", "edu_head_23" = "Post Edu.", 
                                 "age_head" = "Age (Head)", "sex_head1" = "Female (Yes = 1)"),
          custom.gof.rows = list("Correction Term" = c("NO", "NO", "YES", "YES"), 
                                 "Country FE" = c("YES", "YES", "YES", "YES"))) 

# Export
texreg(list(model0, model1, model2, model3), digits = 3, 
       caption = "Air-conditioning Impact on Electricity Quantity - OECD-Europe",
       stars = c(0.1, 0.05, 0.01), 
       custom.model.names = c("No Correction Term", "Correction Term", "Interaction", "Interaction"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/heterogeneities/country/OECD_Europe.tex', sep=''), append=F,  float.pos = "htbp", label = "main: ely_afr",
       omit.coef = "(country)|(Intercept)|(selection)", 
       custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac1:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", "curr_CDD18_db" = "CDD", 
                              "I(curr_CDD18_db^2)" = "CDD$^2$", "curr_HDD18_db" = "HDD", "I(curr_HDD18_db^2)" = "HDD$^2$",
                              "ln_total_exp_usd_2011" = "Log(Exp)", "urban_sh" = "Urbanisation (\\%)", 
                              "ownership_d1" = "House Ownership (Yes = 1)", "n_members" = "Household Size", 
                              "edu_head_22" = "Secondary Edu.", "edu_head_23" = "Post Edu.", 
                              "age_head" = "Age (Head)", "sex_head1" = "Female (Yes = 1)"),
       custom.gof.rows = list("Correction Term" = c("NO", "NO", "YES", "YES"), 
                              "Country FE" = c("YES", "YES", "YES", "YES")), 
       caption.above = TRUE)
 
texreg(list(model4, model5, model6, model7), digits = 3, 
       caption = "Air-conditioning Impact on Electricity Quantity - OECD-Non Europe",
       stars = c(0.1, 0.05, 0.01), 
       custom.model.names = c("No Correction Term", "Correction Term", "Interaction", "Interaction"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/heterogeneities/country/OECD_NonEurope.tex', sep=''), append=F,  float.pos = "htbp", label = "main: ely_afr",
       omit.coef = "(country)|(Intercept)|(selection)", 
       custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac1:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", "curr_CDD18_db" = "CDD", 
                              "I(curr_CDD18_db^2)" = "CDD$^2$", "curr_HDD18_db" = "HDD", "I(curr_HDD18_db^2)" = "HDD$^2$",
                              "ln_total_exp_usd_2011" = "Log(Exp)", "urban_sh" = "Urbanisation (\\%)", 
                              "ownership_d1" = "House Ownership (Yes = 1)", "n_members" = "Household Size", 
                              "edu_head_22" = "Secondary Edu.", "edu_head_23" = "Post Edu.", 
                              "age_head" = "Age (Head)", "sex_head1" = "Female (Yes = 1)"),
       custom.gof.rows = list("Correction Term" = c("NO", "NO", "YES", "YES"), 
                              "Country FE" = c("YES", "YES", "YES", "YES")), 
       caption.above = TRUE)


# Save coefficients in data frame
dydx_ac_eu <- as.data.frame(summary(model_eu)$coefficients)
dydx_ac_eu <- tibble::rownames_to_column(dydx_ac_eu, "Variable")
dydx_ac_eu <- dydx_ac_eu %>% filter(Variable == "ac1" | Variable == "ac1:curr_CDD18_db" | Variable == "ac1:I(curr_CDD18_db^2)")
col <- colnames(dydx_ac_eu)
dydx_ac_tot_eu <- summary(ac_eff_eu)
dydx_ac_tot_eu <- dydx_ac_tot_eu %>% dplyr::select(-c(lower, upper))
colnames(dydx_ac_tot_eu) <- col
dydx_ac_tot_eu$Variable[dydx_ac_tot_eu$Variable == "ac1"] <- "ac_tot"
dydx_ac_eu <- rbind(dydx_ac_eu, dydx_ac_tot_eu)

dydx_ac_neu <- as.data.frame(summary(model_noneu)$coefficients)
dydx_ac_neu <- tibble::rownames_to_column(dydx_ac_neu, "Variable")
dydx_ac_neu <- dydx_ac_neu %>% filter(Variable == "ac1" | Variable == "ac1:curr_CDD18_db" | Variable == "ac1:I(curr_CDD18_db^2)")
col <- colnames(dydx_ac_neu)
dydx_ac_tot_neu <- summary(ac_eff_neu)
dydx_ac_tot_neu <- dydx_ac_tot_neu %>% dplyr::select(-c(lower, upper))
colnames(dydx_ac_tot_neu) <- col
dydx_ac_tot_neu$Variable[dydx_ac_tot_neu$Variable == "ac1"] <- "ac_tot"
dydx_ac_neu <- rbind(dydx_ac_neu, dydx_ac_tot_neu)

# Save the R Environment will be used for the projections
dydx_ac <- dydx_ac_eu
save(list = c("reg_ac_eu", "HH_Europe", "model3", "dydx_ac"), 
     file = paste(output,'for_projections/oecdeu_dmcf.RData', sep=''))
dydx_ac <- dydx_ac_neu
save(list = c("reg_ac_noneu", "HH_NonEurope", "model7", "dydx_ac"), 
     file = paste(output,'for_projections/oecdnoneu_dmcf.RData', sep=''))
