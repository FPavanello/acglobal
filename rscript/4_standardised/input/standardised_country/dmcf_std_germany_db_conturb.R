
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
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/'
}

house <- paste(stub,'6-Projections/repo/household/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/'
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/figures/'
script <- 'C:/Users/Standard/Documents/Github/acglobal/rscript/4_standardised/input/standardised_country/'



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
                            edu_head_2 = as.factor(edu_head_2))

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
HH_Germany <- dplyr::filter(global, country == "Germany")

# Macro-region
HH_Germany$macroarea <- ifelse((HH_Germany$adm1_code==15 | HH_Germany$adm1_code==8 | HH_Germany$adm1_code==6| HH_Germany$adm1_code==5| HH_Germany$adm1_code==9| HH_Germany$adm1_code==4 | HH_Germany$adm1_code==3| HH_Germany$adm1_code==14), 1, 
                               ifelse((HH_Germany$adm1_code==10 | HH_Germany$adm1_code==7 | HH_Germany$adm1_code==16 | HH_Germany$adm1_code==13), 2, 3))

# Ref category for education
HH_Germany$edu_head_2 <- relevel(HH_Germany$edu_head_2,"1") # we do so since there are no edu = 0

# Scale variable
HH_Germany <- HH_Germany %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD18_db)),
                                    std_CDD = as.numeric(scale(curr_CDD18_db)),
                                    std_elyp = as.numeric(scale(ln_ely_p)),
                                    std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                    std_HDD = as.numeric(scale(curr_HDD18_db)),
                                    std_urban_sh = as.numeric(scale(urban_sh)),
                                    std_n_members = as.numeric(scale(n_members)),
                                    std_age_head = as.numeric(scale(age_head)))

# AC formula for Germany - no ownership_id at the moment
ac_formula_deu <- ac ~ std_CDD_mean + I(std_CDD_mean^2) + std_CDD_mean*std_texp + I(std_CDD_mean^2)*std_texp + std_texp + std_CDD + I(std_CDD^2) + 
  std_elyp + std_elyp*std_CDD_mean + std_elyp*I(std_CDD_mean^2) + std_HDD + I(std_HDD^2) + std_elyp*std_n_members + 
   std_urban_sh + std_n_members + edu_head_2 + std_age_head + sex_head | macroarea

# Logistic regression of AC on covariates
reg_ac <- feglm(ac_formula_deu, family = binomial(link = "logit"), 
                data = HH_Germany, weights = ~weight, cluster = c("adm1"))

# Average marginal effects (AMEs)
ac_margins <- avg_slopes(reg_ac, wts = HH_Germany$weight)
gc()

# Predicted probabilities
HH_Germany$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Germany$ac_obs <- ifelse(HH_Germany$phat0_obs>0.5 & !is.na(HH_Germany$phat0_obs), 1 , 0)

# Selection term for intensive margin part
HH_Germany$xb_noac = 1-HH_Germany$phat0_obs               
HH_Germany$selection = ifelse(HH_Germany$ac==1, 
                              (HH_Germany$xb_noac*log(HH_Germany$xb_noac)/HH_Germany$phat0_obs) + log(HH_Germany$phat0_obs), 
                              (HH_Germany$phat0_obs*log(HH_Germany$phat0_obs)/HH_Germany$xb_noac) + log(HH_Germany$xb_noac))

# Formula electricity expenditure without selection
ely_formula_deu <- ln_ely_q ~ ac + ac*std_CDD + ac*I(std_CDD^2) + std_CDD + I(std_CDD^2) + 
  std_texp + std_HDD + I(std_HDD^2) + std_elyp + 
  std_urban_sh + std_n_members + edu_head_2 +
  std_age_head + sex_head + selection | macroarea

# With selection
model <- feols(ely_formula_deu, data = HH_Germany, weights = ~weight, cluster = c("adm1")); summary(model)

# Marginal effect of AC
ely_margins <- avg_slopes(model, slope = "dydx", wts = HH_Germany$weight)

# Export
save(list = c("ely_margins", "ac_margins"), file = paste(interm,'standardised/deu_dmcf.RData', sep=''))
