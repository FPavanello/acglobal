
##########################################

#                 Table A4

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
output <- paste(stub,'output/tables/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/tables/'


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
model0 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model0)


# Formula electricity q without selection
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | ln_ely_p ~ as.factor(country)

# Without selection
model_iv00 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_iv00)

# First-stage
fitstat(model_iv00, ~ ivwald + ivwaldall + wh + sargan, cluster = "adm1")
stat00 <- fitstat(model_iv00, ~ ivwald + ivwaldall + wh + sargan, cluster = "adm1")
kp00 <- stat00$`ivwald1::ln_ely_p`$stat
gc()


# Formula electricity q without selection
ely_formula  <- ln_ely_q ~ ac +
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | ln_ely_p ~ as.factor(adm1)

# Without selection
model_iv0 <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model_iv0)

# First-stage
fitstat(model_iv0, ~ ivwald + ivwaldall + wh + sargan, cluster = "adm1")
stat0 <- fitstat(model_iv0, ~ ivwald + ivwaldall + wh + sargan, cluster = "adm1")
kp0 <- stat0$`ivwald1::ln_ely_p`$stat
gc()

# Mean electricity quantity
mean <- weighted.mean(exp(sec$ln_ely_q), sec$weight)
mean

# Country
cntry <- length(unique.default(sec$country))
cntry


# Export
texreg(list(model0, model_iv00, model_iv0), digits = 3, 
       caption = "The Effect of Air-conditioning on Residential Electricity Quantity - Instrumenting Electricity Prices",
       stars = c(0.1, 0.05, 0.01), custom.model.names = c("DMF","2SLS", "2SLS"),
       custom.note = "Clustered standard errors at the ADM1 level in parentheses. Regressions are conducted using survey weights. $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$", 
       file = paste(output,'TableA4.tex', sep=''), append=F,  
       float.pos = "htbp", label = "main: tableA4",
       omit.coef = "(country)|(Intercept)|(selection)",
       custom.coef.map = list("ac"= "AC", "ac:curr_CDD18_db" = "AC $\\times$ CDD", 
                              "ac:I(curr_CDD18_db^2)" = "AC $\\times$ CDD$^2$",
                              "ln_ely_p" = "Log(P)",
                              "fit_ln_ely_p" = "Log(P)"),
       custom.gof.rows = list("Controls" = c("YES", "YES", "YES"),
                              "Correction Term" = c("YES", "YES", "YES"),
                              "Mean Outcome" = c(mean, mean, mean),
                              "Kleibergen-Paap Wald test" = c("", round(kp0, digits = 3), round(kp00, digits = 3)),
                              "Countries" = c(cntry, cntry, cntry)), 
       caption.above = TRUE)

# Clean
rm(model_iv0, model_iv00, model0, stat0, stat00, sec, fs)
gc()
