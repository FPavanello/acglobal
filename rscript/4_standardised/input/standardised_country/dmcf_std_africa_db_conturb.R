
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

house <- paste(stub,'data/household/', sep='')
interm <- paste(stub,'results/regressions/for_graphs/subsamples/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/'
output <- paste(stub,'output/figures/', sep='')
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
HH_Africa <- dplyr::filter(global, country == "Burkina Faso" | country == "Ghana" | country == "Kenya" | 
                             country == "Malawi" | country == "Niger" | country == "Nigeria")

# Scale variable
HH_Africa <- HH_Africa %>% mutate(std_CDD_mean = as.numeric(scale(mean_CDD18_db)),
                                  std_elyp = as.numeric(scale(ln_ely_p)),
                                  std_CDD = as.numeric(scale(curr_CDD18_db)),
                                  std_texp = as.numeric(scale(ln_total_exp_usd_2011)),
                                  std_HDD = as.numeric(scale(curr_HDD18_db)),
                                  std_urban_sh = as.numeric(scale(urban_sh)),
                                  std_n_members = as.numeric(scale(n_members)),
                                  std_age_head = as.numeric(scale(age_head)))

# AC formula for Africa
ac_formula_afr <- ac ~ std_CDD_mean + I(std_CDD_mean^2) + std_CDD_mean*std_texp + I(std_CDD_mean^2)*std_texp + std_texp +
  std_elyp + std_elyp*std_CDD_mean + std_elyp*I(std_CDD_mean^2) + std_CDD + I(std_CDD^2) + std_elyp*ownership_d + std_elyp*std_n_members + 
  std_urban_sh + std_n_members + ownership_d + edu_head_2 + std_age_head + sex_head + 
  housing_index_lab | country

# Logistic regression of AC on covariates
reg_ac <- feglm(ac_formula_afr, family = binomial(link = "logit"), data = HH_Africa, weights = ~weight, cluster = c("adm1"))

# Average marginal effects (AMEs)
ac_margins <- summary(avg_slopes(reg_ac, wts = HH_Africa$weight))
gc()

# Predicted probabilities
HH_Africa$phat0_obs <- as.numeric(predict(reg_ac, type="response")) 

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
HH_Africa$ac_obs <- ifelse(HH_Africa$phat0_obs>0.5 & !is.na(HH_Africa$phat0_obs), 1 , 0)

# Selection term
HH_Africa$xb_noac = 1-HH_Africa$phat0_obs               
HH_Africa$selection = ifelse(HH_Africa$ac==1, 
                             (HH_Africa$xb_noac*log(HH_Africa$xb_noac)/HH_Africa$phat0_obs) + log(HH_Africa$phat0_obs), 
                             (HH_Africa$phat0_obs*log(HH_Africa$phat0_obs)/HH_Africa$xb_noac) + log(HH_Africa$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula_afr  <- ln_ely_q ~ ac + ac*std_CDD + ac*I(std_CDD^2) + std_CDD + I(std_CDD^2) + 
  std_texp + std_HDD + I(std_HDD^2) + std_elyp +
  std_urban_sh + ownership_d + std_n_members + edu_head_2 + std_age_head + sex_head + 
  housing_index_lab + selection | adm1

# With selection
model <- feols(ely_formula_afr, data = HH_Africa, weights = ~weight, cluster = c("adm1")); summary(model)

# Marginal effect of AC
ely_margins <- summary(avg_slopes(model, slope = "dydx", wts = HH_Africa$weight))

# Save the R Environment will be used for the projections
save(list = c("ely_margins", "ac_margins"), file = paste(output,'/for_graphs/standardised/afr_dmcf.RData', sep=''))
