
##########################################

#           Table S17, Table S18

##########################################

# Free memory
.rs.restartR()
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
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'data/household/', sep='')
interm <- paste(stub,'results/regressions/for_graphs/subsamples/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/'
output <- paste(stub,'output/supplementary/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/supplementary/'


# Load global data
global <- readRDS(paste(house,'global.rds', sep=''))

# Interaction prices
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

# Shares
ivshare <- global %>% dplyr::group_by(country, subnat) %>% dplyr::summarise(cap_solar = mean(cap_solar, na.rm = TRUE), cap_wind = mean(cap_wind, na.rm = TRUE), 
                                                                            cap_hydro = mean(cap_hydro, na.rm = TRUE), cap_nuclear  = mean(cap_nuclear, na.rm = TRUE), 
                                                                            cap_otherres = mean(cap_otherres, na.rm = TRUE), cap_coal = mean(cap_coal, na.rm = TRUE),
                                                                            cap_gas = mean(cap_gas, na.rm = TRUE), cap_oil = mean(cap_oil, na.rm = TRUE), cap_other = mean(cap_other, na.rm = TRUE))
ivshare <- ivshare %>% mutate(cap_tot = cap_solar + cap_wind + cap_hydro + cap_nuclear + cap_otherres + cap_coal + cap_gas + cap_oil + cap_other)
ivshare <- ivshare %>% mutate(sh_solar = ifelse(cap_solar == 0, 0, cap_solar/cap_tot),
                              sh_wind = ifelse(cap_wind == 0, 0, cap_wind/cap_tot),
                              sh_hydro = ifelse(cap_hydro == 0, 0, cap_hydro/cap_tot),
                              sh_nuclear = ifelse(cap_nuclear == 0, 0, cap_nuclear/cap_tot),
                              sh_otherres = ifelse(cap_otherres == 0, 0, cap_otherres/cap_tot),
                              sh_coal = ifelse(cap_coal == 0, 0, cap_coal/cap_tot),
                              sh_gas = ifelse(cap_gas == 0, 0, cap_gas/cap_tot),
                              sh_oil = ifelse(cap_oil == 0, 0, cap_oil/cap_tot),
                              sh_other = ifelse(cap_other == 0, 0, cap_other/cap_tot),
                              sh_ren = sh_solar + sh_wind + sh_otherres)
ivshare <- ivshare %>% dplyr::select(-c(cap_solar, cap_wind, cap_nuclear, cap_otherres, cap_hydro, cap_coal, cap_gas, cap_oil, cap_other))
gc()

# Merge
global <- merge(global, ivshare, by = c("country", "subnat"))


## 1) Electricity capacity by fuel
# Formula first-stage
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),]

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity q without selection
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + ln_ely_p + selection

# Without selection
model1 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model1)


# Formula electricity q without selection 
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | ln_ely_p ~ sh_coal + sh_gas + sh_hydro + sh_otherres + sh_nuclear + sh_other

# Without selection
model_iv1 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_iv1)

# First-stage
elyp_formula  <- ln_ely_p ~ ac + sh_coal + sh_gas + sh_hydro + sh_otherres + sh_nuclear + sh_other + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection

model_ivfs <- feols(elyp_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_ivfs)
fitstat(model_iv1, ~ ivwald + ivwaldall + wh, cluster = "adm1")
stat1 <- fitstat(model_iv1, ~ ivwald + ivwaldall + wh, cluster = "adm1")
kp1 <- stat1$`ivwald1::ln_ely_p`$stat
gc()

# Mean electricity quantity
meaniv1 <- weighted.mean(exp(sec[obs(model_iv1),]$ln_ely_q), sec[obs(model_iv1),]$weight)
meaniv1


# Formula first-stage
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),]

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity q without selection
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + ln_ely_p + selection | adm1

# Without selection
model2 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model1)


# Formula electricity q without selection 
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1 | ln_ely_p ~ sh_coal + sh_gas + sh_hydro + sh_otherres + sh_nuclear + sh_other

# Without selection
model_iv2 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_iv1)

# First-stage
elyp_formula  <- ln_ely_p ~ ac + sh_coal + sh_gas + sh_hydro + sh_otherres + sh_nuclear + sh_other + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

model_ivfs <- feols(elyp_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_ivfs)
fitstat(model_iv2, ~ ivwald + ivwaldall + wh, cluster = "adm1")
stat2 <- fitstat(model_iv2, ~ ivwald + ivwaldall + wh, cluster = "adm1")
kp2 <- stat2$`ivwald1::ln_ely_p`$stat
gc()

# Mean electricity quantity
meaniv2 <- weighted.mean(exp(sec[obs(model_iv2),]$ln_ely_q), sec[obs(model_iv2),]$weight)
meaniv2

# Country
cntry <- length(unique.default(sec$country))
cntry

# Clean
gc()


# Export
texreg(list(model1, model_iv1, model2, model_iv2), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - Instrumenting Electricity Prices",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF","DMF","DMF","DMF"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'Table_S17.tex', sep=''), append=F,  
       float.pos = "htbp", label = "si: tableS17",
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "ln_ely_p" = "Log(P)",
                              "fit_ln_ely_p" = "Log(P)"),
       custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES"),
                              "ADM-1 FE" = c("NO", "NO", "YES", "YES"),
                              "Mean Outcome" = c(meaniv1, meaniv1, meaniv2, meaniv2),
                              "Kleibergen-Paap Wald test" = c("", round(kp1, digits = 3), "", round(kp2, digits = 3)),
                              "Countries" = c(cntry, cntry, cntry, cntry)), 
       caption.above = TRUE)

# Clean
rm(model_iv1, model_iv2, model_ivfs, model1, model2, stat1, stat2, kp1, kp2, sec, fs, meaniv1, meaniv2)
gc()


## 2) ADM-1/Country as instruments for air-conditioning
# Formula
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  ln_ely_p + urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head

# OLS
modelac0 <- feols(ac_formula, data = global, weights = ~weight, cluster = c("adm1")); summary(modelac0)

# Formula
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | ln_ely_p ~ as.factor(country)

# IV
modelac1 <- feols(ac_formula, data = global, weights = ~weight, cluster = c("adm1")); summary(modelac1)

# First-stage
fitstat(modelac1, ~ ivwald + ivwaldall + wh + sargan, cluster = "adm1")
statac1 <- fitstat(modelac1, ~ ivwald + ivwaldall + wh + sargan, cluster = "adm1")
kpac1 <- statac1$`ivwald1::ln_ely_p`$stat
gc()

# Formula
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 + 
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 + 
  urban_sh + ownership_d + 
  n_members + edu_head_2 + age_head + sex_head | ln_ely_p ~ as.factor(adm1)

# IV
modelac2 <- feols(ac_formula, data = global, weights = ~weight, cluster = c("adm1")); summary(modelac2)

# First-stage
fitstat(modelac2, ~ ivwald + ivwaldall + wh + sargan, cluster = "adm1")
statac2 <- fitstat(modelac2, ~ ivwald + ivwaldall + wh + sargan, cluster = "adm1")
kpac2 <- statac2$`ivwald1::ln_ely_p`$stat
gc()

# Mean electricity quantity
mean <- weighted.mean(exp(global$ln_ely_q), global$weight)
mean


# Export
# Only AC
texreg(list(modelac0, modelac1, modelac2), digits = 3, 
       caption = "Air-conditioning ownership - Instrumenting Electricity Prices",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("LPM","2SLS", "2SLS"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'Table_S18.tex', sep=''), append=F,  
       float.pos = "htbp", label = "si: tableS18",
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.coef.map = list(
         "ln_ely_p" = "Log(P)",
         "fit_ln_ely_p" = "Log(P)"),
       custom.gof.rows = list("Controls" = c("YES", "YES", "YES"),
                              "Mean Outcome" = c(mean, mean, mean),
                              "Kleibergen-Paap Wald test" = c("", round(kpac1, digits = 3), round(kpac2, digits = 3)),
                              "Countries" = c("25", "25", "25")), 
       caption.above = TRUE)

# Clean
rm(modelac0, modelac1, modelac2, kpac1, kpac2, mean, statac1, statac2)
gc()
