
## This R-script:
##      1) exploits CDD-dry bulb 18 deg
##      2) conducts logit regressions for the global data set
##      3) run intensive margin regressions: electricity expenditure on climate + covariates
##         using Dubin and McFadden (1984) approach
##      4) focus on heterogeneities by income country-specific decile. We do so by creating sub-samples
##      5) t-tests on the coefficients

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
library(purrr)
library(StatMeasures)
library(survey)
library(reshape)
library(reshape2)
library(fixest)
library(marginaleffects)


# Set directory
user <- 'fp'

if (user=='fp') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

house <- paste(stub,'6-Projections/data/household/', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')

# Load global data
global <- readRDS(paste(house,'global.rds', sep=''))

# Interaction prices
global$mean_CDD18_db <- global$meanpy_CDD18_db
global$mean_hDD18_db <- global$meanpy_hDD18_db
global <- global %>% mutate(ln_ely_p_cdd = ln_ely_p*mean_CDD18_db,
                            ln_ely_p_cdd2 = ln_ely_p*(mean_CDD18_db^2),
                            ln_ely_p_own = ln_ely_p*as.numeric(as.character(ownership_d)),
                            ln_ely_p_nme = ln_ely_p*n_members,
                            mean_CDD18_db2 = mean_CDD18_db^2,
                            mean_CDD18_db_exp = ln_total_exp_usd_2011*mean_CDD18_db,
                            mean_CDD18_db2_exp = ln_total_exp_usd_2011*(mean_CDD18_db^2),
                            curr_CDD18_db2 = curr_CDD18_db^2)

# Check
global <- global[complete.cases(global$ln_ely_q), ]
global <- global[complete.cases(global$ac), ]
global <- global[complete.cases(global$ln_total_exp_usd_2011), ]
global <- global[complete.cases(global$mean_CDD18_db), ]
global <- global[complete.cases(global$urban), ]
global <- global[complete.cases(global$ownership_d), ]
global <- global[complete.cases(global$n_members), ]
global <- global[complete.cases(global$age_head), ]
global <- global[complete.cases(global$country), ]
global <- global[complete.cases(global$weight), ]
global <- global[complete.cases(global$sex_head), ]
global <- global[complete.cases(global$urban_sh), ]
global <- global[complete.cases(global$ln_ely_p), ]
global <- global[complete.cases(global$curr_CDD18_db), ]
global <- global[complete.cases(global$curr_HDD18_db), ]
global <- global[complete.cases(global$adm1), ]
global <- global %>% filter(ln_ely_q > 0)
global <- global %>% filter(weight > 0)

# Decile of total expenditure
#global <- global %>% mutate(exp_cap_usd_2011 = total_exp_usd_2011/n_members)
setDT(global)[,qnt_inc := cut(total_exp_usd_2011, breaks = quantile(total_exp_usd_2011, probs = seq(0, 1, 0.2)),
                              labels = c(1, 2, 3, 4, 5), include.lowest = TRUE), by = country]


################################

#        INCOME QUINTILES      #

################################

# AC formula for global
ac_formula <- as.numeric(as.character(ac)) ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | country

# Split income sample
global_inc <- split(global, global$qnt_inc)

# Logit - extensive margin on subsamples
results_ext <- lapply(global_inc, function(data){
  
  # Logit
  feglm(ac_formula, family = binomial(link = "logit"), data = data, weights = ~weight, cluster = c("adm1"))
  
} 
)

lapply(results_ext, summary)

# Find the selection term
selection_inc <- lapply(global_inc, function(data){

  # Logit
  results <- feglm(ac_formula, family = binomial(link = "logit"), data = data, weights = ~weight, cluster = c("adm1"))
  
  # Save data set for which there are obs both in first and second stage
  data <- data[obs(results),]
  
  # Predicted probabilities
  data$phat0_obs <- as.numeric(predict(results, type="response")) 
  
  # Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
  data$ac_obs <- ifelse(data$phat0_obs>0.5 & !is.na(data$phat0_obs), 1 , 0)
  
  # Selection term for intensive margin part
  data$xb_noac = 1-data$phat0_obs               
  data$selection = ifelse(data$ac==1, 
                          (data$xb_noac*log(data$xb_noac)/data$phat0_obs) + log(data$phat0_obs), 
                          (data$phat0_obs*log(data$phat0_obs)/data$xb_noac) + log(data$xb_noac))
  
  return(data)
  
} 
)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | country

## Results intensive margin
results <- lapply(selection_inc, function(data){
  
  # Regression
  model <- feols(ely_formula, data = data, weights = ~weight, cluster = c("adm1"))
  
} 
)

lapply(results, summary) # Drivers across quintiles

# Save the models
model1 <- results[[1]]
model2 <- results[[2]]
model3 <- results[[3]]
model4 <- results[[4]]
model5 <- results[[5]]

## Effect of AC - At the averages
ac <- lapply(selection_inc, function(data){
  
  # Formula
  ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
    curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
    urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection + country
  
  # Survey
  survey <- svydesign(data = data, ids = ~ adm1, weights = ~ weight)
  
  # Regression
  model <- svyglm(ely_formula, design = survey, na.action=na.omit)
  ac_eff <- summary(margins(model, variables = "ac", design = survey))
} 
)

lapply(ac, summary) # Effect of AC across quintiles

# Save the effects
ac_eff1 <- ac[[1]]
ac_eff1$qnt_inc <- "1st"
ac_eff2 <- ac[[2]]
ac_eff2$qnt_inc <- "2nd"
ac_eff3 <- ac[[3]]
ac_eff3$qnt_inc <- "3rd"
ac_eff4 <- ac[[4]]
ac_eff4$qnt_inc <- "4th"
ac_eff5 <- ac[[5]]
ac_eff5$qnt_inc <- "5th"

ac_eff <- rbind(ac_eff1, ac_eff2, ac_eff3, ac_eff4, ac_eff5)

# T-tests: Z = (b1-b2)/sqrt((SEb1)^2+(SEb2)^2) > z0.99 = 2.576
( ac_eff[1,2] - ac_eff[2,2] )/( sqrt( (ac_eff[1,3])^2 + (ac_eff[2,3])^2 ) ) # 1v2
( ac_eff[1,2] - ac_eff[3,2] )/( sqrt( (ac_eff[1,3])^2 + (ac_eff[3,3])^2 ) ) # 1v3
( ac_eff[1,2] - ac_eff[4,2] )/( sqrt( (ac_eff[1,3])^2 + (ac_eff[4,3])^2 ) ) # 1v4
( ac_eff[1,2] - ac_eff[5,2] )/( sqrt( (ac_eff[1,3])^2 + (ac_eff[5,3])^2 ) ) # 1v5
( ac_eff[2,2] - ac_eff[3,2] )/( sqrt( (ac_eff[2,3])^2 + (ac_eff[3,3])^2 ) ) # 2v3 
( ac_eff[2,2] - ac_eff[4,2] )/( sqrt( (ac_eff[2,3])^2 + (ac_eff[4,3])^2 ) ) # 2v4 
( ac_eff[2,2] - ac_eff[5,2] )/( sqrt( (ac_eff[2,3])^2 + (ac_eff[5,3])^2 ) ) # 2v5
( ac_eff[3,2] - ac_eff[4,2] )/( sqrt( (ac_eff[3,3])^2 + (ac_eff[4,3])^2 ) ) # 3v4
( ac_eff[3,2] - ac_eff[5,2] )/( sqrt( (ac_eff[3,3])^2 + (ac_eff[5,3])^2 ) ) # 3v5
( ac_eff[4,2] - ac_eff[5,2] )/( sqrt( (ac_eff[3,3])^2 + (ac_eff[5,3])^2 ) ) # 4v5

# No differences!

## Effect of AC - At different level of CDD
ac_cdd <- lapply(selection_inc, function(data){
  
          # Regression
          model <- feols(ely_formula, data = data, weights = ~weight, cluster = c("adm1"))
#          ac_eff <- margins(model, at = list(curr_CDD_db = c(0, 5, 10, 15)), variables = "ac", design = data) # it doesn't work
          
          # Dataframes for saving growth rates and SEs
          AMEs_cdd <- data.frame(Var = "AME")
          SEs_cdd <- data.frame(Var = "SE")
          
          # Coeff and var-cov
          vcov <- vcov(model)
          betas <- coef(model)
          
          # Compute manually the AMEs
          cdd <- seq(0, 30, 1)
          for (val in cdd){
            
            form_ses <- sprintf(paste0("~ x1 + x17*", val, " + x18*(", val, ")^2"))
            AMEs_cdd <- AMEs_cdd %>% dplyr::mutate(!!paste0("cdd", val) := betas[1] + betas[17]*val + betas[18]*(val)^2)
            SEs_cdd <- SEs_cdd %>% dplyr::mutate(!!paste0("cdd", val) := msm::deltamethod(as.formula(form_ses), 
                                                                              coef(model), cov = vcov, ses = TRUE))
            
          }
          
          ac_eff_cdd <- rbind(AMEs_cdd, SEs_cdd)
          
          return(ac_eff_cdd)
  
} 
)

# Save the effects
ac_cdd1 <- ac_cdd[[1]]
ac_cdd1$qnt_inc <- "1st"
ac_cdd2 <- ac_cdd[[2]]
ac_cdd2$qnt_inc <- "2nd"
ac_cdd3 <- ac_cdd[[3]]
ac_cdd3$qnt_inc <- "3rd"
ac_cdd4 <- ac_cdd[[4]]
ac_cdd4$qnt_inc <- "4th"
ac_cdd5 <- ac_cdd[[5]]
ac_cdd5$qnt_inc <- "5th"

ac_cdd <- rbind(ac_cdd1, ac_cdd2, ac_cdd3, ac_cdd4, ac_cdd5)
ac_cdd <- melt(setDT(ac_cdd), id = c("qnt_inc", "Var"), 
               value.name = c("value"), variable.name = "cdd")
ac_cdd <- cast(ac_cdd, qnt_inc + cdd ~ Var)
ac_cdd <- ac_cdd %>% mutate(CI_lb = AME - 1.96*SE,
                            CI_ub = AME + 1.96*SE)

# Compare the models
screenreg(list(model1, model2, model3, model4, model5), digits = 3, 
          caption = "Impact of Air-conditioning on Electricity Quantity by Income Quintile - Global",
          stars = c(0.1, 0.05, 0.01), custom.model.names = c("1st Quintile", "2nd Quintile", "3rd Quintile",
                                                             "4th Quintile", "5th Quintile"),
          custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD18_db" = "AC $\\times$ CDD",
                                 "ac1:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$"),
          custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES"), 
                                 "Country FE" = c("YES", "YES", "YES", "YES", "YES"), 
                                 "Countries" = c("25", "25", "25", "25", "25")))

# Export
texreg(list(model1, model2, model3, model4, model5), digits = 3, 
       caption = "Impact of Air-conditioning on Electricity Quantity by Income Quintile - Global",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("1st Quintile", "2nd Quintile", "3rd Quintile",
                                                          "4th Quintile", "5th Quintile"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'electricity/heterogeneities/subsamples/Global_qinc.tex', sep=''), append=F,  float.pos = "htbp", label = "main: ely_afr",
       custom.coef.map = list("ac1"= "AC", "ac1:curr_CDD18_db" = "AC $\\times$ CDD",
                              "ac1:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$"),
       custom.gof.rows = list("Controls" = c("YES", "YES", "YES", "YES", "YES"),
                              "Correction Term" = c("YES", "YES", "YES", "YES", "YES"), 
                              "Country FE" = c("YES", "YES", "YES", "YES", "YES"), 
                              "Countries" = c("25", "25", "25", "25", "25")), caption.above = TRUE)

# Save the R Environment will be used for the projections
save(list = c("ac_eff", "ac_cdd"), 
     file = paste(output,'for_graphs/subsamples/global_csquintiles.RData', sep=''))
