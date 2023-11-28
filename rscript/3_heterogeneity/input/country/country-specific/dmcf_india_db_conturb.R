
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts logit regressions for India using 2012 wave
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


house <- paste(stub,'6-Projections/data/household/India/chps/', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')
script <- paste(stub,'6-Projections/rscripts/dmcf/regressions/country/with_continuous_urbanisation/', sep='')

# Load Household data
HH_India <- readRDS(paste(house,'india_cmie.rds', sep=''))

# Add urbanisation share
source(paste0(stub, "6-Projections/rscripts/process_raw_data/add_urban/add_urban_ind.R"))
HH_India$geometry <- NULL
HH_India$adm1 <- as.character(HH_India$state)

# Interaction prices
HH_India$housing_index_lab <- as.factor(HH_India$housing_index_lab)
HH_India$mean_CDD18_db <- HH_India$meanpy_CDD18_db
HH_India$mean_hDD18_db <- HH_India$meanpy_hDD18_db
HH_India <- HH_India %>% mutate(ln_ely_p = log(ely_p_usd_2011),
                                  ln_ely_p_cdd = ln_ely_p*mean_CDD18_db,
                                  ln_ely_p_cdd2 = ln_ely_p*(mean_CDD18_db^2),
                                  ln_ely_p_own = ln_ely_p*as.numeric(as.character(ownership_d)),
                                  ln_ely_p_nme = ln_ely_p*n_members,
                                  mean_CDD18_db2 = mean_CDD18_db^2,
                                  mean_CDD18_db_exp = ln_total_exp_usd_2011*mean_CDD18_db,
                                  mean_CDD18_db2_exp = ln_total_exp_usd_2011*(mean_CDD18_db^2),
                                  curr_CDD18_db2 = curr_CDD18_db^2)

# Only those with not missing values 
HH_India <- HH_India[complete.cases(HH_India$ac), ]
HH_India <- HH_India[complete.cases(HH_India$mean_CDD18_db), ]
HH_India <- HH_India[complete.cases(HH_India$mean_HDD18_db), ]
HH_India <- HH_India[complete.cases(HH_India$curr_CDD18_db), ]
HH_India <- HH_India[complete.cases(HH_India$curr_HDD18_db), ]
HH_India <- HH_India[complete.cases(HH_India$ln_total_exp_usd_2011), ]
HH_India <- HH_India[complete.cases(HH_India$urban_sh), ]
HH_India <- HH_India[complete.cases(HH_India$n_members), ]
HH_India <- HH_India[complete.cases(HH_India$sh_under16), ]
HH_India <- HH_India[complete.cases(HH_India$ownership_d), ]
HH_India <- HH_India[complete.cases(HH_India$edu_head_2), ]
HH_India <- HH_India[complete.cases(HH_India$ln_ely_p), ]
HH_India <- HH_India[complete.cases(HH_India$age_head), ]
HH_India <- HH_India[complete.cases(HH_India$sex_head), ]
HH_India <- HH_India %>% filter(ln_ely_q > 0)
HH_India <- HH_India %>% filter(weight > 0)

# Survey - re-run to add new variable
HH_India_svy <- svydesign(data = HH_India, ids = ~ adm1, weights = ~ weight)


##################################

#        Extensive margin        #

##################################

# AC formula for India
ac_formula_ind <- as.numeric(as.character(ac)) ~ mean_CDD18_db + mean_CDD18_db2 + curr_CDD18_db + curr_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 +
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + housing_index_lab + state

# Logistic regression of AC on covariates
reg_ac <- feglm(ac_formula_ind, family = binomial(link = "logit"), 
                data = HH_India, weights = ~weight, cluster = c("adm1"))

# Save AME results
#margins <- margins(reg_ac, design = HH_India_svy)
#summary(margins)
#xtable(summary(margins), display=rep('g', 8), caption = "Logit Regression for Air-conditioning Ownership - India")
#print(xtable(summary(margins), display=rep('g', 8), caption = "Logit Regression for Air-conditioning Ownership - India"), 
#      file= paste(output,'airconditioning/main/IND.tex', sep=''),append=F, table.placement = "htbp",
#      caption.placement="top")

# Predicted probabilities
HH_India$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_India$ac_obs <- ifelse(HH_India$phat0_obs>0.5 & !is.na(HH_India$phat0_obs), 1 , 0)

table(HH_India$ac, HH_India$ac_obs)

# Selection term
HH_India$xb_noac = 1-HH_India$phat0_obs               
HH_India$selection = ifelse(HH_India$ac==1, 
                             (HH_India$xb_noac*log(HH_India$xb_noac)/HH_India$phat0_obs) + log(HH_India$phat0_obs), 
                             (HH_India$phat0_obs*log(HH_India$phat0_obs)/HH_India$xb_noac) + log(HH_India$xb_noac))

# Survey - re-run to add new variable
HH_India_svy <- svydesign(data = HH_India, ids = ~ adm1, weights = ~ weight)


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

# Formula electricity expenditure without interactions
ely_formula_ind <- ln_ely_q ~ ac +   
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) | state

# With selection
model0 <- feols(ely_formula_ind, data = HH_India, weights = ~weight, cluster = c("adm1")); summary(model0)

# Formula electricity expenditure without interactions
ely_formula_ind <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + housing_index_lab | state

# With selection
model1 <- feols(ely_formula_ind, data = HH_India, weights = ~weight, cluster = c("adm1")); summary(model1)

# Formula electricity expenditure with interactions
ely_formula_ind <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + housing_index_lab + selection | state

# With selection
model2 <- feols(ely_formula_ind, data = HH_India, weights = ~weight, cluster = c("adm1")); summary(model2)

# Formula electricity expenditure with interactions
ely_formula_ind <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) +
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + housing_index_lab + selection | state

# With selection
model3 <- feols(ely_formula_ind, data = HH_India, weights = ~weight, cluster = c("adm1")); summary(model3)

# Formula electricity expenditure with interactions
ely_formula_ind <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) +
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + housing_index_lab + selection + state

# Marginal effect of AC
model <- svyglm(ely_formula_ind, design = HH_India_svy, na.action=na.omit)
ac_eff <- margins(model, variables = "ac", design = HH_India_svy)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_India$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_India$ely_q) # Average effect 181.80 kWh


# Compare the models
screenreg(list(model0, model1, model2, model3), digits = 3, 
          caption = "The Effect of Air-conditioning on Residential Electricity Quantity - India",
          stars = c(0.1, 0.05, 0.01), 
          custom.model.names = c("No Correction Term", "Correction Term", "Interaction", "Interaction"),
          omit.coef = "(state)|(Intercept)|(selection)", 
          custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD18_db" = "AC $\\times$ CDD", 
                                 "ac1:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", "curr_CDD18_db" = "CDD", 
                                 "I(curr_CDD18_db^2)" = "CDD$^2$", "curr_HDD18_db" = "HDD", "I(curr_HDD18_db^2)" = "HDD$^2$",
                                 "ln_total_exp_usd_2011" = "Log(Exp)",
                                 "ln_ely_p" = "Log(P)",
                                 "urban_sh" = "Urbanisation (\\%)", 
                                 "ownership_d1" = "House Ownership (Yes = 1)", "n_members" = "Household Size", 
                                 "edu_head_21" = "Primary Edu.", "edu_head_22" = "Secondary Edu.", "edu_head_23" = "Post Edu.", 
                                 "age_head" = "Age (Head)", "sex_head1" = "Female (Yes = 1)", 
                                 "housing_index_lab2" = "Housing (Medium)", "housing_index_lab3" = "Housing (High)"),
          custom.gof.rows = list("Correction Term" = c("NO", "YES", "YES", "YES"), 
                                 "State FE" = c("YES", "YES", "YES", "YES")))

# Export
texreg(list(model0, model1, model2, model3), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - India",
       stars = c(0.1, 0.05, 0.01), 
       custom.model.names = c("No Correction Term", "Correction Term", "Interaction", "Interaction"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/heterogeneities/country/IND.tex', sep=''), append=F,  float.pos = "htbp", label = "main: ely_ind",
       custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac1:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", "curr_CDD18_db" = "CDD", 
                              "I(curr_CDD18_db^2)" = "CDD$^2$", "curr_HDD18_db" = "HDD", "I(curr_HDD18_db^2)" = "HDD$^2$",
                              "ln_total_exp_usd_2011" = "Log(Exp)",
                              "ln_ely_p" = "Log(P)",
                              "urban_sh" = "Urbanisation (\\%)", 
                              "ownership_d1" = "House Ownership (Yes = 1)", "n_members" = "Household Size", 
                              "edu_head_21" = "Primary Edu.", "edu_head_22" = "Secondary Edu.", "edu_head_23" = "Post Edu.", 
                              "age_head" = "Age (Head)", "sex_head1" = "Female (Yes = 1)", 
                              "housing_index_lab2" = "Housing (Medium)", "housing_index_lab3" = "Housing (High)"),
       custom.gof.rows = list("Correction Term" = c("NO", "YES", "YES", "YES"), 
                              "State FE" = c("YES", "YES", "YES", "YES")), 
       caption.above = TRUE)


# Save coefficients in data frame
dydx_ac <- as.data.frame(summary(model)$coefficients)
dydx_ac <- tibble::rownames_to_column(dydx_ac, "Variable")
dydx_ac <- dydx_ac %>% filter(Variable == "ac" | Variable == "ac:curr_CDD18_db" | Variable == "ac:I(curr_CDD18_db^2)")
col <- colnames(dydx_ac)
dydx_ac_tot <- summary(ac_eff)
dydx_ac_tot <- dydx_ac_tot %>% dplyr::select(-c(lower, upper))
colnames(dydx_ac_tot) <- col
dydx_ac_tot$Variable[dydx_ac_tot$Variable == "ac"] <- "ac_tot"
dydx_ac <- rbind(dydx_ac, dydx_ac_tot)

# Save the R Environment will be used for the projections
save(list = c("reg_ac", "HH_India", "model3", "dydx_ac"), 
     file = paste(output,'for_projections/ind_dmcf.RData', sep=''))
