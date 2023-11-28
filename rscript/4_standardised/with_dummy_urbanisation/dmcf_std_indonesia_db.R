
## This R-script:
##      1) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      2) conducts STANDARDISED logit regressions for Indonesia using 2016 wave
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
HH_Indonesia <- readRDS(paste(house,'/idn_susenas.rds', sep=''))

# log total exp
HH_Indonesia$ln_total_exp_usd_2011 <- log(HH_Indonesia$total_exp_usd_2011)
HH_Indonesia$ln_total_exp_usd_2011[HH_Indonesia$total_exp_usd_2011 == 0] <- NA

# Log electricity quantity
HH_Indonesia$ln_ely_q <- log(HH_Indonesia$ely_q)

# Only those with not missing values 
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$ac), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$mean_CDD_db), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$mean_HDD_db), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$curr_CDD_db), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$curr_HDD_db), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$ln_total_exp_usd_2011), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$urban), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$n_members), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$sh_under16), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$housing_index_lab), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$ownership_d), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$edu_head_2), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$occupation_head), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$age_head), ]
HH_Indonesia <- HH_Indonesia[complete.cases(HH_Indonesia$sex_head), ]

# Arrange some vars
HH_Indonesia$state <- as.factor(HH_Indonesia$state)
HH_Indonesia$urban <- as.factor(HH_Indonesia$urban)
HH_Indonesia$ac <- as.factor(HH_Indonesia$ac)
HH_Indonesia$refrigerator <- as.factor(HH_Indonesia$refrigerator)
HH_Indonesia$ownership <- as.factor(HH_Indonesia$ownership)
HH_Indonesia$edu_head_2 <- as.factor(HH_Indonesia$edu_head_2)
HH_Indonesia$sex_head <- as.factor(HH_Indonesia$sex_head)
HH_Indonesia$housing_index_lab <- as.factor(HH_Indonesia$housing_index_lab)
HH_Indonesia$occupation_head <- as.factor(HH_Indonesia$occupation_head)

# Only who has electricity access
HH_Indonesia <- HH_Indonesia %>% filter(ely_access == 1)

# CDD in 100s
HH_Indonesia <- HH_Indonesia %>% mutate(mean_CDD_db = mean_CDD_db/100)
HH_Indonesia <- HH_Indonesia %>% mutate(curr_CDD_db = curr_CDD_db/100)
HH_Indonesia <- HH_Indonesia %>% mutate(mean_HDD_db = mean_HDD_db/100)
HH_Indonesia <- HH_Indonesia %>% mutate(curr_HDD_db = curr_HDD_db/100)

# Check
HH_Indonesia <- HH_Indonesia[(is.finite(HH_Indonesia$ln_ely_q) & !is.na(HH_Indonesia$ln_ely_q)),]


# Scale variable
HH_Indonesia <- HH_Indonesia %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD_db)),
                                        std_CDD = as.numeric(scale(curr_CDD_db)),
                                        std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                        std_HDD = as.numeric(scale(curr_HDD_db)),
                                        std_n_members = as.numeric(scale(n_members)),
                                        std_age_head = as.numeric(scale(age_head)),
                                        std_sh_under16 = as.numeric(scale(sh_under16)))


##################################

#        Extensive margin        #

##################################

# AC formula for Indonesia
ac_formula_idn <- ac ~ std_CDD_mean*std_texp + 
  as.factor(urban) + std_n_members + std_sh_under16 + 
  as.factor(ownership_d) + as.factor(edu_head_2) + std_age_head + as.factor(sex_head) + 
  as.factor(housing_index_lab) + as.factor(state)

# Logistic regression of AC on covariates
reg_ac <- glm(ac_formula_idn, data = HH_Indonesia, family = binomial(logit), na.action=na.omit); summary(reg_ac)
vcov <- vcovCL(reg_ac, cluster = HH_Indonesia$state_district) # clusterise SEs
coeftest(reg_ac, vcov = vcovCL(reg_ac, cluster = HH_Indonesia$state_district)) # clustered standard errors

# Save AME results
margins <- margins(reg_ac, vcov = vcov)
summary(margins)
ac_margins <- summary(margins)
xtable(summary(margins), display=rep('g', 8), 
       caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Indonesia")
print(xtable(summary(margins), display=rep('g', 8), 
             caption = "Logit Regression for Air-conditioning Ownership using Standardised Covariates - Indonesia"), 
      file= paste(output,'airconditioning/standardised/IDN.tex', sep=''),append=F, table.placement = "htbp",
      caption.placement="top")

# Predicted probabilities
HH_Indonesia$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Indonesia$ac_obs <- ifelse(HH_Indonesia$phat0_obs>0.5 & !is.na(HH_Indonesia$phat0_obs), 1 , 0)

table(HH_Indonesia$ac, HH_Indonesia$ac_obs)

# Selection term
HH_Indonesia$xb_noac = 1-HH_Indonesia$phat0_obs               
HH_Indonesia$selection = ifelse(HH_Indonesia$ac==1, 
                                (HH_Indonesia$xb_noac*log(HH_Indonesia$xb_noac)/HH_Indonesia$phat0_obs) + log(HH_Indonesia$phat0_obs), 
                                (HH_Indonesia$phat0_obs*log(HH_Indonesia$phat0_obs)/HH_Indonesia$xb_noac) + log(HH_Indonesia$xb_noac))


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
ely_formula_idn <- ln_ely_q ~ ac + std_CDD*std_texp + 
  urban + std_n_members + std_sh_under16 + ownership_d + edu_head_2 + 
  housing_index_lab + std_age_head + sex_head + state

# With selection
model0 <- lm(ely_formula_idn, data = HH_Indonesia, na.action=na.omit); summary(model0)
vcov0 <- cluster.vcov(model0, HH_Indonesia$state_district) # clusterise SEs
coeftest(model0, vcov=vcov0) # with clustered SEs
model0_cl <- coeftest(model0, vcov=vcov0) # with clustered SEs

# Save betas
betas <- coef(model0)

# Marginal effect of AC
ac_eff <- margins(model0, variables = "ac", vcov = vcov0)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Indonesia$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Indonesia$ely_q) # Average effect 907.70 kWh


# Formula electricity quantity without interactions
ely_formula_idn <- ln_ely_q ~ ac + std_CDD*std_texp + 
  urban + std_n_members + std_sh_under16 + ownership_d + edu_head_2 + 
  housing_index_lab + std_age_head + sex_head + state + selection

# With selection
model1 <- lm(ely_formula_idn, data = HH_Indonesia, na.action=na.omit); summary(model)
vcov1 <- cluster.vcov(model1, HH_Indonesia$state_district) # clusterise SEs
coeftest(model1, vcov=vcov1) # with clustered SEs
model1_cl <- coeftest(model1, vcov=vcov1) # with clustered SEs

# Save betas
betas <- coef(model1)

# Marginal effect of AC
ac_eff <- margins(model1, variables = "ac", vcov = vcov1)
summary(ac_eff)
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Indonesia$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Indonesia$ely_q) # Average effect 954.21 kWh


# Formula electricity quantity with interactions
ely_formula_idn <- ln_ely_q ~ ac + ac*std_CDD + std_CDD*std_texp + 
  urban + std_n_members + std_sh_under16 + ownership_d + edu_head_2 + 
  housing_index_lab + std_age_head + sex_head + state + selection

# With selection
model2 <- lm(ely_formula_idn, data = HH_Indonesia, na.action=na.omit); summary(model2)
vcov2 <- cluster.vcov(model2, HH_Indonesia$state_district) # clusterise SEs
coeftest(model2, vcov=vcov2) # with clustered SEs
model2_cl <- coeftest(model2, vcov=vcov2) # with clustered SEs

# Save betas
betas <- coef(model2)

# Marginal effect of AC
ac_eff <- margins(model2, variables = "ac", vcov = vcov2)
summary(ac_eff)
ely_margins <- summary(margins(model2, vcov = vcov2))
mean(ac_eff$dydx_ac1) # same as #213
mean(HH_Indonesia$ely_q)*exp(mean(ac_eff$dydx_ac1)) - mean(HH_Indonesia$ely_q) # Average effect 777.58 kWh


# Compare the models
screenreg(list(model0_cl, model1_cl, model2_cl), digits = 3,
          caption = "Air-conditioning Impact on Electricity Quantity using Standardised Covariates - Indonesia",
          stars = c(0.1, 0.05, 0.01), custom.model.names = c("No Correction Term", "Correction Term", "Interaction"),
          omit.coef = "(state)|(Intercept)", reorder.coef = c(1, 17, 2, 3, 15, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16), 
          custom.coef.names = c("AC", "CDD", "Log(Exp)", "Urban (Yes = 1)", "Household Size", "Share of Minors", 
                                "House Ownership (Yes = 1)", "Primary Edu.", "Secondary Edu.", "Post Edu.", "Housing Quality (Med.)", 
                                "Housing Quality (High)", "Age (Head)", "Gender Head (Yes = 1)", 
                                "CDD $\\times$ Log(Exp)", "Correction Term", "AC $\\times$ CDD"),
          custom.gof.rows = list("State FE" = c("YES", "YES", "YES")))
texreg(list(model0_cl, model1_cl, model2_cl), digits = 3, 
       caption = "Air-conditioning Impact on Electricity Quantity using Standardised Covariates - Indonesia",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("No Correction Term", "Correction Term", "Interaction"),
       omit.coef = "(state)|(Intercept)", reorder.coef = c(1, 17, 2, 3, 15, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16), 
       custom.coef.names = c("AC", "CDD", "Log(Exp)", "Urban (Yes = 1)", "Household Size", "Share of Minors", 
                             "House Ownership (Yes = 1)", "Primary Edu.", "Secondary Edu.", "Post Edu.", "Housing Quality (Med.)", 
                             "Housing Quality (High)", "Age (Head)", "Gender Head (Yes = 1)", 
                             "CDD $\\times$ Log(Exp)", "Correction Term", "AC $\\times$ CDD"),
       custom.gof.rows = list("State FE" = c("YES", "YES", "YES")))

# Export
texreg(list(model0_cl, model1_cl, model2_cl), digits = 3, 
       caption = "Air-conditioning Impact on Electricity Quantity using Standardised Covariates - Indonesia",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("No Correction Term", "Correction Term", "Interaction"),
       custom.note = "Clustered std. errors at the district level in parentheses. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/standardised/IDN.tex', sep=''), append=F,  float.pos = "htbp", label = "main: ely_mex",
       omit.coef = "(state)|(Intercept)", reorder.coef = c(1, 17, 2, 3, 15, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16), 
       custom.coef.names = c("AC", "CDD", "Log(Exp)", "Urban (Yes = 1)", "Household Size", "Share of Minors", 
                             "House Ownership (Yes = 1)", "Primary Edu.", "Secondary Edu.", "Post Edu.", "Housing Quality (Med.)", 
                             "Housing Quality (High)", "Age (Head)", "Gender Head (Yes = 1)", 
                             "CDD $\\times$ Log(Exp)", "Correction Term", "AC $\\times$ CDD"),
       custom.gof.rows = list("State FE" = c("YES", "YES", "YES")), caption.above = TRUE)

# Export
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/idn_dmcf.RData', sep=''))

