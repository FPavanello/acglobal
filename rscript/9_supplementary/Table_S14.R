
#############################################################

#                         Table S14

#############################################################

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
library(DescTools)

# Set users
user <- 'user'

if (user=='user') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

house <- paste(stub,'6-Projections/repo/household/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/supplementary/'


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


# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  curr_HDD18_db + I(curr_HDD18_db^2) +
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fs <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fs)

# Save data set for which there are obs both in first and second stage
sec <- global[obs(fs),] # 52k observations lost

# Predicted probabilities
sec$phat0_obs <- as.numeric(predict(fs, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1 , 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev0 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_lev0)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev1 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_lev1)


# Winsorize based on electricity consumption, income and CDD
global$ely_qw <- DescTools::Winsorize(global$ely_q, quantile(global$ely_q, probs = c(0.05, 0.95), na.rm = TRUE))
global$ln_total_exp_usd_2011w <- Winsorize(global$total_exp_usd_2011, quantile(global$total_exp_usd_2011, probs = c(0.05, 0.95), na.rm = TRUE))
global$curr_CDD18_dbw <- Winsorize(global$curr_CDD18_db, quantile(global$curr_CDD18_db, probs = c(0.05, 0.95), na.rm = TRUE))
global$mean_CDD18_dbw <- Winsorize(global$mean_CDD18_db, quantile(global$mean_CDD18_db, probs = c(0.05, 0.95), na.rm = TRUE))
gc()

# AC formula for global
ac_formula <- ac ~ mean_CDD18_dbw + I(mean_CDD18_dbw^2) +
  ln_total_exp_usd_2011w*mean_CDD18_dbw + ln_total_exp_usd_2011w*I(mean_CDD18_dbw^2) + ln_total_exp_usd_2011w + curr_CDD18_dbw + I(curr_CDD18_dbw^2) +  
  curr_HDD18_db + I(curr_HDD18_db^2) + ln_ely_p + ln_ely_p*mean_CDD18_dbw + ln_ely_p*I(mean_CDD18_dbw^2) + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head  | adm1

# Logistic regression
fsw <- feglm(ac_formula, family = binomial(link = "logit"), data = global, weights = ~weight, cluster = c("adm1")); summary(fsw)

# Save data set for which there are obs both in first and second stage
secw <- global[obs(fsw),] # 52k observations lost

# Predicted probabilities
secw$phat0_obs <- as.numeric(predict(fsw, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
secw$ac_obs <- ifelse(secw$phat0_obs>0.5 & !is.na(secw$phat0_obs), 1 , 0)

# Selection term
secw$xb_noac = 1-secw$phat0_obs               
secw$selection = ifelse(secw$ac==1, 
                        (secw$xb_noac*log(secw$xb_noac)/secw$phat0_obs) + log(secw$phat0_obs), 
                        (secw$phat0_obs*log(secw$phat0_obs)/secw$xb_noac) + log(secw$xb_noac))
gc()

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_qw ~ ac +
  curr_CDD18_dbw + I(curr_CDD18_dbw^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011w + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev2 <- feols(ely_formula, data = secw, weights = ~weight, cluster = c("adm1")); summary(model_lev2)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_qw ~ ac + ac*curr_CDD18_dbw + ac*I(curr_CDD18_dbw^2) + 
  curr_CDD18_dbw + I(curr_CDD18_dbw^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011w + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev3 <- feols(ely_formula, data = secw, weights = ~weight, cluster = c("adm1")); summary(model_lev3)


# Trim based on electricity consumption
cut_point_top <- quantile(global$ely_q, 0.95)
cut_point_bottom <- quantile(global$ely_q, 0.05)

# Filter
global_trim <- filter(global, ely_q >= cut_point_bottom & ely_q <= cut_point_top)

# AC formula for global
ac_formula <- ac ~ mean_CDD18_db + mean_CDD18_db2 +
  mean_CDD18_db_exp + mean_CDD18_db2_exp + ln_total_exp_usd_2011 + curr_CDD18_db + curr_CDD18_db2 +  
  curr_HDD18_db + I(curr_HDD18_db^2) +
  ln_ely_p + ln_ely_p_cdd + ln_ely_p_cdd2 + ln_ely_p_nme + ln_ely_p_own + 
  urban_sh + ownership_d +
  n_members + edu_head_2 + age_head + sex_head | adm1

# Logistic regression
fst <- feglm(ac_formula, family = binomial(link = "logit"), data = global_trim, weights = ~weight, cluster = c("adm1")); summary(fst)

# Save data set for which there are obs both in first and second stage
sect <- global_trim[obs(fst),] # 52k observations lost

# Predicted probabilities
sect$phat0_obs <- as.numeric(predict(fst, type="response"))

# Find old HHs classified as owning an AC (NB: using ROC curve we have seen that we are GOOD at predicting those who have AC)
sect$ac_obs <- ifelse(sect$phat0_obs>0.5 & !is.na(sect$phat0_obs), 1 , 0)

# Selection term
sect$xb_noac = 1-sect$phat0_obs               
sect$selection = ifelse(sect$ac==1, 
                        (sect$xb_noac*log(sect$xb_noac)/sect$phat0_obs) + log(sect$phat0_obs), 
                        (sect$phat0_obs*log(sect$phat0_obs)/sect$xb_noac) + log(sect$xb_noac))
gc()

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev4 <- feols(ely_formula, data = sect, weights = ~weight, cluster = c("adm1")); summary(model_lev4)

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model_lev5 <- feols(ely_formula, data = sect, weights = ~weight, cluster = c("adm1")); summary(model_lev5)


# Mean electricity quantity
mean_lev1 <- weighted.mean(sec$ely_q, sec$weight)
mean_lev1
mean_lev2 <- weighted.mean(secw$ely_qw, secw$weight)
mean_lev2
mean_lev3 <- weighted.mean(sect$ely_qw, sect$weight)
mean_lev3

# Number of countries
cwin1 <- length(unique.default(sec$country))
cwin2 <- length(unique.default(secw$country))
cwin3 <- length(unique.default(sect$country))

# Clean
rm(fs, sec, fsw, secw, global_trim, sect, fst)
gc()


# Export
texreg(list(model_lev0, model_lev1, model_lev2, model_lev3, model_lev4, model_lev5), digits = 3, 
       caption = "Robustness Checks - Electricity in Levels",
       stars = c(0.1, 0.05, 0.01), 
       custom.model.names = c("Full", "Full", 
                              "Winsorized", "Winsorized", "Trimmed", "Trimmed"),
       custom.note = "'Subnational' means at the most disaggregated geographical information for each country. Regressions are conducted using survey weights.
          Standard errors are clustered at the ADM1 level. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$.", 
       file = paste(output,'TableS14.tex', sep=''), append=F,  float.pos = "H", label = "si: tableS14",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:curr_CDD18_dbw" = "AC $\\times$ CDD",
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$", 
                              "ac:I(curr_CDD18_dbw^2)" = "AC $\\times$ CDD$^2$"),
       custom.gof.rows = list("Correction Term" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "ADM-1 FE" = c("YES", "YES", "YES", "YES", "YES", "YES"), 
                              "Mean Outcome (kWh)" = c(mean_lev1, mean_lev1, mean_lev2, mean_lev2, mean_lev3, mean_lev3), 
                              "Countries" = c(cwin1, cwin1, cwin2, cwin2, cwin3, cwin3)), 
       caption.above = TRUE)
