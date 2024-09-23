
##########################################

#               Figure 2

##########################################

# Free memory
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
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
library(msm)
library(margins)
library(texreg)
library(xtable)
library(stargazer)
library(effects)
library(survey)
library(purrr)
library(patchwork)
library(ggsci)
library(prediction)
library(reshape)
library(reshape2)
library(mob)
library(fixest)

# Set users
user <- 'user'

if (user=='user') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'data/household/', sep='')
output <- paste(stub,'output/figures/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/figures/'


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



## Run the model
# AC formula for global
ac_formula <- as.numeric(as.character(ac)) ~ mean_CDD18_db + mean_CDD18_db2 + 
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
sec$ac_obs <- ifelse(sec$phat0_obs>0.5 & !is.na(sec$phat0_obs), 1, 0)

# Selection term
sec$xb_noac = 1-sec$phat0_obs               
sec$selection = ifelse(sec$ac==1, 
                       (sec$xb_noac*log(sec$xb_noac)/sec$phat0_obs) + log(sec$phat0_obs), 
                       (sec$phat0_obs*log(sec$phat0_obs)/sec$xb_noac) + log(sec$xb_noac))

# Formula electricity expenditure for electricity expenditure
ely_formula  <- ln_ely_q ~ ac + ac*curr_CDD18_db + ac*I(curr_CDD18_db^2) + 
  curr_CDD18_db + I(curr_CDD18_db^2) + curr_HDD18_db + I(curr_HDD18_db^2) + ln_total_exp_usd_2011 + ln_ely_p +
  urban_sh + ownership_d + n_members + edu_head_2 + age_head + sex_head + selection | adm1

# With selection
model <- feols(ely_formula, data = sec, weights = ~weight, cluster = c("adm1")); summary(model)

## Marginal Effect
# Dataframes for saving growth rates and SEs
AMEs_cdd <- data.frame(Var = "AME")
SEs_cdd <- data.frame(Var = "SE")

# Coeff and var-cov
vcov <- vcov(model)
betas <- coef(model)

# Compute manually the AMEs
cdd <- seq(0, 40, 1)
for (val in cdd){
  
  form_ses <- sprintf(paste0("~ x1 + x17*", val, " + x18*(", val, ")^2"))
  AMEs_cdd <- AMEs_cdd %>% dplyr::mutate(!!paste0("cdd", val) := betas[1] + betas[17]*val + betas[18]*(val)^2)
  SEs_cdd <- SEs_cdd %>% dplyr::mutate(!!paste0("cdd", val) := msm::deltamethod(as.formula(form_ses), 
                                                                                coef(model), cov = vcov, ses = TRUE))
  
}

ac_eff_cdd <- rbind(AMEs_cdd, SEs_cdd)
ac_cdd <- melt(setDT(ac_eff_cdd), id = c("Var"), 
               value.name = c("value"), variable.name = "cdd")
ac_cdd <- cast(ac_cdd, cdd ~ Var)
ac_cdd <- ac_cdd %>% mutate(CI_lb = AME - 1.96*SE,
                            CI_ub = AME + 1.96*SE)

# Adjust
ac_cdd$cdd <- as.character(ac_cdd$cdd)
ac_cdd$cdd <- substr(ac_cdd$cdd, 4, 5)
ac_cdd$cdd <- as.numeric(ac_cdd$cdd)

# Create histogram bars
for (i in seq(0,60,1)) {
  
  global$bin[global$curr_CDD18_db >= i & global$curr_CDD18_db < i+1] <- i
  
}

global$value <- 1 
hist <- global %>% dplyr::group_by(bin) %>% dplyr::select(bin, weight, value) %>% dplyr::summarise(Freq = sum(value*weight))
hist <- hist %>% filter(bin <= 40)
#hist <- global %>% dplyr::select(curr_CDD18_db) %>%
#  do(data.frame(table(cut(.$curr_CDD18_db, breaks=seq(0, 31, 1), include.lowest=T))))
ac_cdd <- cbind(ac_cdd, hist)
ac_cdd <- ac_cdd %>% mutate(cdd = cdd*100,
                            AME = AME*100,
                            CI_lb = CI_lb*100,
                            CI_ub = CI_ub*100)

# Average Effect
ac_c <- ac_cdd %>% 
  ggplot(aes(x= cdd, y=AME)) + geom_point() + 
  geom_errorbar(aes(ymin=CI_lb, ymax=CI_ub), alpha=0.6, size=0.6, width=0.4) +
  geom_bar(aes(x=cdd, y=Freq/10000000),stat="identity", fill = "steelblue", alpha = 0.3) + # scale for plotting reasons
  geom_hline(aes(yintercept=35.1), linetype="dashed", color = "red") +
  xlab("Cooling Degree Days") +
  ylab("Additional Electricity Demand due to Air-conditioning \n(Marginal Effect, %)") +
  theme_classic() +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color="lightgray"),
        panel.grid.minor=element_blank(),
        plot.title = element_text(size = 20, family = "Tahoma", hjust = 0.5), 
        text = element_text(size = 20, family = "Tahoma"), 
        #                  axis.title.x = element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_blank())

ac_c

# Combine plots
ggsave(paste(output, 'Figure2.png', sep = ''), last_plot(), scale=2.5, width = 4.5, height = 4)
