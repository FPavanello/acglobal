
# .rs.restartR()
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Load packages
library(data.table)
library(plyr)
library(dplyr)
library(FSA)
library(haven)
library(stringr)
library(tidyverse)
library(sandwich)
library(lmtest)
library(ResourceSelection)
library(multiwayvcov)
library(texreg)
library(xtable)
library(stargazer)
library(effects)
library(fixest)
library(marginaleffects)

# Set users
user <- 'user'

if (user=='user') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

house <- paste(stub,'data/household/', sep='')
interm <- paste(stub,'results/regressions/for_graphs/subsamples/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/'
output <- paste(stub,'output/figures/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/figures/'
script <- 'C:/Users/Standard/Documents/Github/acglobal/rscript/3_heterogeneity/input/country/country-specific/'


# Load global data
global <- readRDS(paste(house,'global.rds', sep=''))

# Interaction prices
global <- global %>% mutate(ln_ely_p_cdd = ln_ely_p*mean_CDD18_db,
                            ln_ely_p_cdd2 = ln_ely_p*(mean_CDD18_db^2),
                            ln_ely_p_own = ln_ely_p*ownership_d,
                            ln_ely_p_nme = ln_ely_p*n_members,
                            mean_CDD18_db2 = mean_CDD18_db^2,
                            mean_CDD18_db_exp = ln_total_exp_usd_2011*mean_CDD18_db,
                            mean_CDD18_db2_exp = ln_total_exp_usd_2011*(mean_CDD18_db^2),
                            curr_CDD18_db2 = curr_CDD18_db^2,
                            edu_head_2 = as.factor(edu_head_2),
                            housing_index_lab = as.factor(housing_index_lab))

# Check
global <- global[complete.cases(global$ln_ely_q), ]
global <- global[complete.cases(global$ac), ]
global <- global[complete.cases(global$ln_total_exp_usd_2011), ]
global <- global[complete.cases(global$mean_CDD18_db), ]
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

# Select countries
HH_China <- dplyr::filter(global, country == "China")

# Macroarea
east <-c("Beijing", "Fujian", "Guangdong", "Hainan", "Hebei", "Jiangsu", "Liaoning", "Guangxi", 
         "Shandong", "Shanghai", "Tianjin", "Zhejiang")
west <- c("Sichuan", "Chongqing", "Gansu", "Guizhou", "Shanxi", "Sichuan", "Yunnan", "Xinjiang Uygur",
          "Ningxia Hui")
central <- c("Anhui", "Heilongjiang", "Henan", "Hubei", "Hunan", "Jiangxi", "Jilin", "Shaanxi")

HH_China$macroarea <- ifelse(HH_China$adm1 %in% east, "East", 
                           ifelse(HH_China$adm1 %in% west, "West", 
                                  ifelse(HH_China$adm1 %in% central, "Central", NA)))


# AC formula for China
ac_formula_chn <- ac ~ mean_CDD18_db + mean_CDD18_db2 + curr_CDD18_db + curr_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 +
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head + 
  housing_index_lab | macroarea

# Logistic regression of AC on covariates
reg_ac <- feglm(ac_formula_chn, family = binomial(link = "logit"), 
                data = HH_China, weights = ~weight, cluster = c("adm1"))

# Predicted probabilities
HH_China$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_China$ac_obs <- ifelse(HH_China$phat0_obs>0.5 & !is.na(HH_China$phat0_obs), 1 , 0)

# Selection term for intensive margin part
HH_China$xb_noac = 1-HH_China$phat0_obs               
HH_China$selection = ifelse(HH_China$ac==1, 
                             (HH_China$xb_noac*log(HH_China$xb_noac)/HH_China$phat0_obs) + log(HH_China$phat0_obs), 
                             (HH_China$phat0_obs*log(HH_China$phat0_obs)/HH_China$xb_noac) + log(HH_China$xb_noac))

# Formula electricity expenditure without selection
ely_formula_chn <- ln_ely_q ~ ac +   
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) | macroarea

# With selection
model0 <- feols(ely_formula_chn, data = HH_China, weights = ~weight, cluster = c("adm1")); summary(model0)

# Formula electricity expenditure without interactions
ely_formula_chn <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + housing_index_lab | macroarea

# With selection
model1 <- feols(ely_formula_chn, data = HH_China, weights = ~weight, cluster = c("adm1")); summary(model1)

# Formula electricity expenditure with interactions
ely_formula_chn <- ln_ely_q ~ ac + 
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + housing_index_lab + selection | macroarea

# With selection
model2 <- feols(ely_formula_chn, data = HH_China, weights = ~weight, cluster = c("adm1")); summary(model2)

# Formula electricity expenditure with interactions
ely_formula_chn <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) +
  curr_CDD18_db + I(curr_CDD18_db^2) + 
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + housing_index_lab + selection | macroarea

# With selection
model3 <- feols(ely_formula_chn, data = HH_China, weights = ~weight, cluster = c("adm1")); summary(model3)

# Marginal effect of AC
ac_eff <- avg_slopes(model3, variables = "ac", slope = "dydx", wts = HH_China$weight)
summary(ac_eff)

# Save coefficients in data frame
dydx_ac <- summary(ac_eff)

# Save the R Environment will be used for the projections
save(list = c("reg_ac", "HH_China", "model3", "dydx_ac"), 
     file = paste(interm,'chn_dmcf.RData', sep=''))
