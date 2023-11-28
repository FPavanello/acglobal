
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for India using 2012 wave
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


# Set directory
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
HH_India <- readRDS(paste(house,'/ind_nss.rds', sep=''))
HH_India <- HH_India %>% filter(occupation_head != 1)

# Arrange some vars
HH_India$urban <- as.factor(HH_India$urban)
HH_India$ac <- as.factor(HH_India$ac)
HH_India$fan <- as.factor(HH_India$fan)
HH_India$refrigerator <- as.factor(HH_India$refrigerator)
HH_India$ownership <- as.factor(HH_India$ownership)
HH_India$edu_head_2 <- as.factor(HH_India$edu_head_2)
HH_India$sex_head <- as.factor(HH_India$sex_head)
HH_India$occupation_head <- as.factor(HH_India$occupation_head)

# log total exp
HH_India$ln_total_exp_usd_2011 <- log(HH_India$total_exp_usd_2011)
HH_India$ln_total_exp_usd_2011[HH_India$total_exp_usd_2011 == 0] <- NA

# Log electricity quantity
HH_India$ln_ely_q <- log(HH_India$ely_q)

# Only who has electricity access
HH_India <- HH_India %>% filter(ely_access == 1)

# Only those with not missing values 
HH_India <- HH_India[complete.cases(HH_India$ac), ]
HH_India <- HH_India[complete.cases(HH_India$mean_CDD_db), ]
HH_India <- HH_India[complete.cases(HH_India$mean_HDD_db), ]
HH_India <- HH_India[complete.cases(HH_India$curr_CDD_db), ]
HH_India <- HH_India[complete.cases(HH_India$curr_HDD_db), ]
HH_India <- HH_India[complete.cases(HH_India$ln_total_exp_usd_2011), ]
HH_India <- HH_India[complete.cases(HH_India$urban), ]
HH_India <- HH_India[complete.cases(HH_India$n_members), ]
HH_India <- HH_India[complete.cases(HH_India$sh_under16), ]
HH_India <- HH_India[complete.cases(HH_India$ownership_d), ]
HH_India <- HH_India[complete.cases(HH_India$edu_head_2), ]
HH_India <- HH_India[complete.cases(HH_India$occupation_head), ]
HH_India <- HH_India[complete.cases(HH_India$age_head), ]
HH_India <- HH_India[complete.cases(HH_India$sex_head), ]

# Check
HH_India <- HH_India[(is.finite(HH_India$ln_ely_q) & !is.na(HH_India$ln_ely_q)),]

# CDD in 100s
HH_India <- HH_India %>% mutate(mean_CDD_db = mean_CDD_db/100,
                                mean_HDD_db = mean_HDD_db/100,
                                curr_CDD_db = curr_CDD_db/100,
                                curr_HDD_db = curr_HDD_db/100)

# Filter variables and rows
HH_India <- dplyr::select(HH_India, ln_ely_q, ely_q, ac, curr_CDD_db, ln_total_exp_usd_2011, curr_HDD_db, ln_total_exp_usd_2011, 
                          urban, n_members, sh_under16, ownership_d, edu_head_2,
                          age_head , sex_head, state, mean_CDD_db, mean_HDD_db, 
                          state_district, state3, weight)

# Scale variable
HH_India <- HH_India %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD_db)),
                                std_CDD = as.numeric(scale(curr_CDD_db)),
                                std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                std_HDD = as.numeric(scale(curr_HDD_db)),
                                std_n_members = as.numeric(scale(n_members)),
                                std_age_head = as.numeric(scale(age_head)),
                                std_sh_under16 = as.numeric(scale(sh_under16)))


##################################

#        Extensive margin        #

##################################

# AC formula for India
ac_formula_ind <- ac ~ std_CDD_mean*std_texp + as.factor(urban) + std_n_members + std_sh_under16 + 
  as.factor(ownership_d) + as.factor(edu_head_2) + std_age_head + as.factor(sex_head) + 
  as.factor(state)

# Logistic regression of AC on covariates
reg_ac <- glm(ac_formula_ind, data = HH_India, family = binomial(logit), na.action=na.omit); summary(reg_ac)
vcov <- vcovCL(reg_ac, cluster = HH_India$state_district) # clusterise SEs
coeftest(reg_ac, vcov = vcovCL(reg_ac, cluster = HH_India$state_district)) # clustered standard errors

# Save AME results
margins <- margins(reg_ac, vcov = vcov)
summary(margins)
ac_margins <- summary(margins)
xtable(summary(margins), display=rep('g', 8), 
       caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - India")
print(xtable(summary(margins), display=rep('g', 8), 
             caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - India"), 
      file= paste(output,'airconditioning/standardised/IND.tex', sep=''),append=F, table.placement = "htbp",
      caption.placement="top")

# Predicted probabilities
HH_India$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 
mean(as.numeric(HH_India$ac))-1
summary(HH_India$phat0_obs)

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_India$ac_obs <- ifelse(HH_India$phat0_obs>0.5 & !is.na(HH_India$phat0_obs), 1 , 0)

table(HH_India$ac, HH_India$ac_obs)

# Selection term
HH_India$xb_noac = 1-HH_India$phat0_obs               
HH_India$selection = ifelse(HH_India$ac==1, 
                            (HH_India$xb_noac*log(HH_India$xb_noac)/HH_India$phat0_obs) + log(HH_India$phat0_obs), 
                            (HH_India$phat0_obs*log(HH_India$phat0_obs)/HH_India$xb_noac) + log(HH_India$xb_noac))


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
ely_formula_ind <- ln_ely_q ~ as.factor(ac) + std_CDD*std_texp +
  as.factor(urban) + std_n_members + std_sh_under16 + as.factor(ownership_d) + as.factor(edu_head_2) + 
  std_age_head + as.factor(sex_head) + as.factor(state)

# With selection
model0 <- lm(ely_formula_ind, data = HH_India, na.action=na.omit); summary(model0)
vcov0 <- cluster.vcov(model0, HH_India$state_district) # clusterise SEs
coeftest(model0, vcov=vcov0) # with clustered SEs
model0_cl <- coeftest(model0, vcov=vcov0) # with clustered SEs

# Save betas
betas <- coef(model0)

# Marginal effect of AC
ac_eff <- margins(model0, variables = c("ac", "std_CDD"), vcov = vcov0)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_India$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_India$ely_q) # Average effect 236.54 kWh


# Formula electricity expenditure without interactions
ely_formula_ind <- ln_ely_q ~ as.factor(ac) + std_CDD*std_texp + 
  as.factor(urban) + std_n_members + std_sh_under16 + as.factor(ownership_d) + as.factor(edu_head_2) + 
  std_age_head + as.factor(sex_head) + as.factor(state) + selection

# With selection
model1 <- lm(ely_formula_ind, data = HH_India, na.action=na.omit); summary(model1)
vcov1 <- cluster.vcov(model1, HH_India$state_district) # clusterise SEs
coeftest(model1, vcov=vcov1) # with clustered SEs
model1_cl <- coeftest(model1, vcov=vcov1) # with clustered SEs

# Save betas
betas <- coef(model1)

# Marginal effect of AC
ac_eff <- margins(model1, variables = c("ac", "std_CDD"), vcov = vcov1)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_India$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_India$ely_q) # Average effect 261.64 kWh

# Formula electricity expenditure with interactions
ely_formula_ind <- ln_ely_q ~ as.factor(ac) + as.factor(ac)*std_CDD +
  std_CDD*std_texp + as.factor(urban) + std_n_members + std_sh_under16 + 
  as.factor(ownership_d) + as.factor(edu_head_2) + 
  std_age_head + as.factor(sex_head) + as.factor(state) + selection

# With selection
model2 <- lm(ely_formula_ind, data = HH_India, na.action=na.omit); summary(model2)
vcov2 <- cluster.vcov(model2, HH_India$state_district) # clusterise SEs
coeftest(model2, vcov=vcov2) # with clustered SEs
model2_cl <- coeftest(model2, vcov=vcov2) # with clustered SEs

# Save betas
betas <- coef(model2)

# Marginal effect of AC
ac_eff <- margins(model2, variables = c("ac", "std_CDD"), vcov = vcov2)
summary(ac_eff)
ely_margins <- summary(margins(model2, vcov = vcov2))
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_India$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_India$ely_q) # Average effect 181.80 kWh


# Compare the models
screenreg(list(model0_cl, model1_cl, model2_cl), digits = 3, 
          caption = "Air-conditioning Impact on Electricity Quantity using Standardised Covariates - India",
          stars = c(0.1, 0.05, 0.01), custom.model.names = c("No Correction Term", "Correction Term", "Interaction"),
          omit.coef = "(state)|(Intercept)", reorder.coef = c(1, 15, 2, 3, 13, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14), 
          custom.coef.names = c("AC", "CDD", "Log(Exp)", "Urban (Yes = 1)", "Household Size", "Share of Minors", 
                                "House Ownership (Yes = 1)", "Primary Edu.", "Secondary Edu.", "Post Edu.",
                                "Age (Head)", "Gender Head (Yes = 1)", 
                                "CDD $\\times$ Log(Exp)", "Correction Term", "AC $\\times$ CDD"),
          custom.gof.rows = list("State FE" = c("YES", "YES", "YES")))
texreg(list(model0_cl, model1_cl, model2_cl), digits = 3, 
       caption = "Air-conditioning Impact on Electricity Quantity using Standardised Covariates - India",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("No Correction Term", "Correction Term", "Interaction"),
       omit.coef = "(state)|(Intercept)", reorder.coef = c(1, 15, 2, 3, 13, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14), 
       custom.coef.names = c("AC", "CDD", "Log(Exp)", "Urban (Yes = 1)", "Household Size", "Share of Minors", 
                             "House Ownership (Yes = 1)", "Primary Edu.", "Secondary Edu.", "Post Edu.",
                             "Age (Head)", "Gender Head (Yes = 1)", 
                             "CDD $\\times$ Log(Exp)", "Correction Term", "AC $\\times$ CDD"),
       custom.gof.rows = list("State FE" = c("YES", "YES", "YES")))

# Export
texreg(list(model0_cl, model1_cl, model2_cl), digits = 3, 
       caption = "Air-conditioning Impact on Electricity Quantity using Standardised Covariates - India",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("No Correction Term", "Correction Term", "Interaction"),
       custom.note = "Clustered std. errors at the district level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/main/IND.tex', sep=''), append=F,  float.pos = "htbp", label = "main: ely_ind",
       omit.coef = "(state)|(Intercept)", reorder.coef = c(1, 15, 2, 3, 13, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14), 
       custom.coef.names = c("AC", "CDD", "Log(Exp)", "Urban (Yes = 1)", "Household Size", "Share of Minors", 
                             "House Ownership (Yes = 1)", "Primary Edu.", "Secondary Edu.", "Post Edu.",
                             "Age (Head)", "Gender Head (Yes = 1)", 
                             "CDD $\\times$ Log(Exp)", "Selection", "AC $\\times$ CDD"),
       custom.gof.rows = list("State FE" = c("YES", "YES", "YES")), caption.above = TRUE)

# Export
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/ind_dmcf.RData', sep=''))
