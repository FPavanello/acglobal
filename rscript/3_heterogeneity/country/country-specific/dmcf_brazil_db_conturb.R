
## This R-script:
##      1) exploits CDD-dry bulb 18 deg and HDD-dry bulb 18 deg
##      2) conducts logit regressions for Brazil using 2017 wave
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
library(tibble)
library(fixest)


# Set users
user <- 'fp'
#user <- 'gf'

if (user=='fp') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

if (user=='gf') {
  stub <- "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/"
}


house <- paste(stub,'6-Projections/data/household/Fourcountries', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')
script <- paste(stub,'6-Projections/rscripts/dmcf/regressions/country/with_continuous_urbanisation/', sep='')

# Load Household data
HH_Brazil <- readRDS(paste(house,'/bra_pof.rds', sep=''))

# Add urbanisation share
source(paste0(stub, "6-Projections/rscripts/process_raw_data/add_urban/add_urban_bra.R"))
HH_Brazil$geometry <- NULL
HH_Brazil$adm1 <- as.character(HH_Brazil$state)

# Interaction prices
HH_Brazil$mean_CDD18_db <- HH_Brazil$meanpy_CDD18_db
HH_Brazil$mean_hDD18_db <- HH_Brazil$meanpy_hDD18_db
HH_Brazil <- HH_Brazil %>% mutate(ln_ely_p = log(ely_p_usd_2011),
                                  ln_ely_p_cdd = ln_ely_p*mean_CDD18_db,
                                  ln_ely_p_cdd2 = ln_ely_p*(mean_CDD18_db^2),
                                  ln_ely_p_own = ln_ely_p*as.numeric(as.character(ownership_d)),
                                  ln_ely_p_nme = ln_ely_p*n_members,
                                  mean_CDD18_db2 = mean_CDD18_db^2,
                                  mean_CDD18_db_exp = ln_total_exp_usd_2011*mean_CDD18_db,
                                  mean_CDD18_db2_exp = ln_total_exp_usd_2011*(mean_CDD18_db^2),
                                  curr_CDD18_db2 = curr_CDD18_db^2)

# Only those with not missing values 
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$ac), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$mean_CDD18_db), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$mean_HDD18_db), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$curr_CDD18_db), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$curr_HDD18_db), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$ln_total_exp_usd_2011), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$urban_sh), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$n_members), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$sh_under16), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$housing_index_lab), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$ownership_d), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$edu_head_2), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$age_head), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$sex_head), ]
HH_Brazil <- HH_Brazil[complete.cases(HH_Brazil$ln_ely_p), ]
HH_Brazil <- HH_Brazil %>% filter(ln_ely_q > 0)
HH_Brazil <- HH_Brazil %>% filter(weight > 0)

# Survey
HH_Brazil_svy <- svydesign(data = HH_Brazil, ids = ~ adm1, weights = ~ weight)


##################################

#        Extensive margin        #

##################################

# AC formula for Brazil
ac_formula_bra <- as.numeric(as.character(ac)) ~ mean_CDD18_db + mean_CDD18_db2 + curr_CDD18_db + curr_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 +
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + 
  housing_index_lab + region3

# Logistic regression of AC on covariates
reg_ac <- feglm(ac_formula_bra, family = binomial(link = "logit"), 
                data = HH_Brazil, weights = ~weight, cluster = c("adm1"))

# Save AME results
#margins <- margins(reg_ac, design = HH_Brazil_svy)
#summary(margins)
#xtable(summary(margins), display=rep('g', 8), caption = "Logit Regression for Air-conditioning Ownership - Brazil")
#print(xtable(summary(margins), display=rep('g', 8), caption = "Logit Regression for Air-conditioning Ownership - Brazil"), 
#      file= paste(output,'airconditioning/main/BRA.tex', sep=''),append=F, table.placement = "htbp",
#      caption.placement="top")

# Predicted probabilities
HH_Brazil$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Brazil$ac_obs <- ifelse(HH_Brazil$phat0_obs>0.5 & !is.na(HH_Brazil$phat0_obs), 1 , 0)

# Selection term
HH_Brazil$xb_noac = 1-HH_Brazil$phat0_obs               
HH_Brazil$selection = ifelse(HH_Brazil$ac==1, 
                             (HH_Brazil$xb_noac*log(HH_Brazil$xb_noac)/HH_Brazil$phat0_obs) + log(HH_Brazil$phat0_obs), 
                             (HH_Brazil$phat0_obs*log(HH_Brazil$phat0_obs)/HH_Brazil$xb_noac) + log(HH_Brazil$xb_noac))

# Survey - re-run to add new variable
HH_Brazil_svy <- svydesign(data = HH_Brazil, ids = ~ adm1, weights = ~ weight)


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
ely_formula_bra <- ln_ely_q ~ ac +   
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) | region3 

# Without selection
model0 <- feols(ely_formula_bra, data = HH_Brazil, weights = ~weight, cluster = c("adm1")); summary(model0)

# Formula electricity quantity without interactions
ely_formula_bra <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + housing_index_lab | region3

# With selection
model1 <- feols(ely_formula_bra, data = HH_Brazil, weights = ~weight, cluster = c("adm1")); summary(model1)

# Formula electricity quantity with interactions
ely_formula_bra <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + housing_index_lab + selection | region3

# With selection
model2 <- feols(ely_formula_bra, data = HH_Brazil, weights = ~weight, cluster = c("hhid")); summary(model2)

# Formula electricity quantity with interactions
ely_formula_bra <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) +
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + housing_index_lab + selection + region3

# With selection
model3 <- svyglm(ely_formula_bra, design = HH_Brazil_svy, na.action=na.omit); summary(model3)

# Marginal effect of AC
ac_eff <- margins(model3, variables = "ac", design = HH_Brazil_svy)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Brazil$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Brazil$ely_q) # Average effect 758.92 kWh


# Compare the models
screenreg(list(model0, model1, model2, model3), digits = 3, 
          caption = "The Effect of Air-conditioning on Residential Electricity Quantity - Brazil",
          stars = c(0.1, 0.05, 0.01), 
          custom.model.names = c("No Correction Term", "Correction Term", "Interaction", "Interaction"),
          omit.coef = "(region)|(Intercept)|(selection)", 
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
          custom.gof.rows = list("Correction Term" = c("NO", "NO", "YES", "YES"), 
                                 "Macro-Region FE" = c("YES", "YES", "YES", "YES")))

# Export
texreg(list(model0, model1, model2, model3), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - Brazil",
       stars = c(0.1, 0.05, 0.01), 
       custom.model.names = c("No Correction Term", "Correction Term", "Interaction", "Interaction"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/heterogeneities/country/BRA.tex', sep=''), append=F,  float.pos = "htbp", label = "main: ely_bra",
       omit.coef = "(region)|(Intercept)|(selection)", 
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
       custom.gof.rows = list("Correction Term" = c("NO", "NO", "YES", "YES"), 
                              "Macro-Region FE" = c("YES", "YES", "YES", "YES")), caption.above = TRUE)


# Save coefficients in data frame
dydx_ac <- as.data.frame(summary(model3)$coefficients)
dydx_ac <- tibble::rownames_to_column(dydx_ac, "Variable")
dydx_ac <- dydx_ac %>% filter(Variable == "ac1" | Variable == "ac1:curr_CDD18_db" | Variable == "ac1:I(curr_CDD18_db^2)")
col <- colnames(dydx_ac)
dydx_ac_tot <- summary(ac_eff)
dydx_ac_tot <- dydx_ac_tot %>% dplyr::select(-c(lower, upper))
colnames(dydx_ac_tot) <- col
dydx_ac_tot$Variable[dydx_ac_tot$Variable == "ac1"] <- "ac_tot"
dydx_ac <- rbind(dydx_ac, dydx_ac_tot)

# Save the R Environment will be used for the projections
save(list = c("reg_ac", "HH_Brazil", "model3", "dydx_ac"), 
     file = paste(output,'/for_projections/bra_dmcf.RData', sep=''))
