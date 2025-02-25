rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

## 1) Load libraries and data ##
library(sandwich)
library(lmtest)
library(foreign)
library(ResourceSelection)
library(optmatch)
library(tidyverse)
library(haven)
library(psych)
library(rnaturalearthdata)
library(sf)
library(gdata)
library(nngeo)
library(caret)
library(MatchIt)
library(ggsci)
library(gdata)
library(jtools)
library(glm2)
library(reshape2)
library(cobalt)
library(relaimpo)
library(domir)
library(pscl)
library(margins)

# Set users
user <- 'gf'

if (user=='gf') {
  stub <- "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

#######################

rdss = list.files(path=paste0(stub, "data/household"), recursive = T, full.names = T, pattern = "rds", ignore.case = T)

rdss = rdss[!grepl("global", rdss)]
rdss = rdss[!grepl("CDD", rdss)]
rdss = rdss[!grepl("us.Rds", rdss)]
rdss = rdss[!grepl("nss", rdss)]
rdss = rdss[!grepl("nss", rdss)]
rdss = rdss[!grepl("submission", rdss)]
rdss <- rdss[c(1,2,3,4,5,6,7,8,10,11,12,13,15,16,18,19,20, 21)]

rdss_l = lapply(rdss, read_rds)

for (i in 1:length(rdss_l)){
  
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ln_ely_q), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ac), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ln_total_exp_usd_2011), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$mean_CDD_db), ]
  #rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$urban), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ownership_d), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$n_members), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$age_head), ]
  #rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$country), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$weight), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$sex_head), ]
  #rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$urban_sh), ]
  
  rdss_l[[i]]$hhid = as.character(rdss_l[[i]]$hhid)
  
  rdss_l[[i]] = dplyr::select(rdss_l[[i]], weight, ln_total_exp_usd_2011, ln_ely_q, starts_with("ely_p"), hhid)
  rdss_l[[i]]$country <- basename(rdss[[i]])
}

rdss_l <- bind_rows(rdss_l)

#tapply(rdss_l$ely_p_usd_2011, rdss_l$country, summary)

########

load(paste0(stub, "results/regressions/for_projections/global_wgt_dmcf.Rdata"))

future_ac_adoption = read_rds(paste0(stub, "results/household_level/", "global_ac.Rds"))
output_ac = read_rds(paste0(stub, "results/household_level/", "global_ely_tot.Rds"))
output_impact_ac = read_rds(paste0(stub, "results/household_level/", "global_ely_due_to_ac.Rds"))

colnames(future_ac_adoption) = paste0("future_ac_adoption_", colnames(future_ac_adoption))
colnames(output_ac) = paste0("output_ac_", colnames(output_ac))
colnames(output_impact_ac) = paste0("output_impact_ac_", colnames(output_impact_ac))

###

weights = read_rds(paste0(stub, "results/household_level/", "weights_global_ac.Rds"))

data_orig = read_rds(paste0(stub, "results/household_level/", "data_global_ac.Rds"))

data_orig = bind_cols(data_orig, future_ac_adoption, output_ac, output_impact_ac, weights)

rdss_l$weight = as.numeric(rdss_l$weight)

rdss_l$ln_ely_q = round(rdss_l$ln_ely_q, 3)
data_orig$ln_ely_q = round(data_orig$ln_ely_q, 3)
data_orig$ln_total_exp_usd_2011 = NULL

###

rdss_l$ely_p_usd_2011 <- ifelse(is.na(rdss_l$ely_p_usd_2011), rdss_l$ely_p, rdss_l$ely_p_usd_2011)
rdss_l$rdss_l$ely_p <-NULL
rdss_l$ely_p_impute_5_95 <-NULL

###

rdss_l = merge(data_orig, rdss_l, by=c("ln_ely_q", "weight"))

rdss_l$total_exp_usd_2011 = exp(rdss_l$ln_total_exp_usd_2011)
rdss_l$ely_p_usd_2011 = rdss_l$ely_p_usd_2011.x
rdss_l$ely_q = exp(rdss_l$ln_ely_q)

rdss_l =  rdss_l %>% group_by(country.x) %>% dplyr::mutate(q_hist=ntile(exp(ln_total_exp_usd_2011), 5), q_ssp2_2050 =ntile(exp(exp_cap_usd_SSP2_2050), 5), q_ssp5_2050 =ntile(exp(exp_cap_usd_SSP5_2050), 5))

rdss_l = dplyr::select(rdss_l, country.x, starts_with("q_"), ely_p_usd_2011, starts_with("exp_cap"), starts_with("future_ac_adoption"),starts_with("output_ac"),starts_with("output_impact_ac"),starts_with("weight_"), hhid.x, total_exp_usd_2011)

pr1 = as.numeric(predict(reg_ely))

load(paste0(stub, "results/regressions/for_projections/global_wgt_dmcf.RData"))

sec$ac[sec$ac == '1'] <- '0'
sec$ac <- as.numeric(sec$ac)
pr2 = as.numeric(predict(reg_ely, newdata=sec))
diff_exp = exp(pr1)

diff_exp = data.frame(diff_exp=diff_exp, hhid=as.character(sec$hhid), country.x=sec$country)
diff_exp$diff_exp = ifelse(diff_exp$diff_exp<=0, NA, diff_exp$diff_exp)

#rdss_l_bk = rdss_l

rdss_l$hhid = rdss_l$hhid.x

rdss_l = merge(rdss_l, diff_exp, c("country.x", "hhid"))

rdss_l$output_impact_ac_SSP2.2020 = rdss_l$diff_exp
rdss_l$output_impact_ac_SSP5.2020 = rdss_l$diff_exp
rdss_l$exp_cap_usd_SSP2_2020 = log(rdss_l$total_exp_usd_2011+1)
rdss_l$exp_cap_usd_SSP5_2020 = log(rdss_l$total_exp_usd_2011+1)

rdss_l_o = data.frame(country=c(rdss_l$country.x))

rdss_l_o$wq_hist=rdss_l$q_hist
rdss_l_o$wq_ssp2_2050=rdss_l$q_ssp2_2050
rdss_l_o$wq_ssp5_2050=rdss_l$q_ssp5_2050

for (ssp in c("SSP2", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  for (year in seq(2020, 2050, 10)){
    
    rdss_l_o[,paste0("ely_share_", ssp, "_", (year))] = (rdss_l[,paste0("output_impact_ac_", ssp, ".", (year))] * rdss_l$ely_p_usd_2011) / exp(rdss_l[,paste0("exp_cap_usd_", ssp, "_", (year))])
    
  }}

#

rdss_l_o = split(rdss_l_o, rdss_l_o$country)
rdss_l = split(rdss_l, rdss_l$country.x)

for (kappa in 1:length(rdss_l)){
  print(kappa)
  rdss_l_o[[kappa]][,paste0("ely_share_", ssp, "_", (year))] = rdss_l_o[[kappa]][,paste0("ely_share_", ssp, "_", (year))] * (rdss_l[[kappa]][,paste0("weight_", year, "_", rcp, "_", tolower(ssp))]/mean(rdss_l[[kappa]][,paste0("weight_", year, "_", rcp, "_", tolower(ssp))], na.rm=T))
  
}

rdss_l_o = bind_rows(rdss_l_o)


#

#rdss_l_o = na.omit(rdss_l_o)

rdss_l_o_bk = rdss_l_o

rdss_l_o = reshape2::melt(rdss_l_o, c(1:4))

boundone_m <- filter(rdss_l_o, variable=="ely_share_SSP2_2020")

boundone_m$variable <- as.character(boundone_m$variable)
boundone_m$variable[boundone_m$variable =="ely_share_SSP2_2020" ] <- " Survey year"
boundone_m$variable[boundone_m$variable =="ely_share_SSP2_2050" ] <- "2050, SSP2"
boundone_m$variable[boundone_m$variable =="ely_share_SSP5_2050" ] <- "2050, SSP5"

unique(boundone_m$country)
unique(boundone_m$variable)
boundone_m$country = as.character(boundone_m$country)

boundone_m$country[boundone_m$country=="Ghana" | boundone_m$country=="Niger" | boundone_m$country=="Nigeria" | boundone_m$country=="Tanzania"| boundone_m$country=="Malawi"| boundone_m$country=="Burkina Faso"| boundone_m$country=="Kenya"] = "Africa"
boundone_m$country[boundone_m$country=="Australia" | boundone_m$country=="Japan" | boundone_m$country=="Canada"] = "OECD-NONEU"
boundone_m$country[boundone_m$country=="France" | boundone_m$country=="Germany"| boundone_m$country=="Netherlands"| boundone_m$country=="Spain"| boundone_m$country=="Sweden"] = "OECD-EU"

##############

calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

#boundone_m$wq = ifelse(boundone_m$variable==" Survey year", boundone_m$wq_hist, ifelse(boundone_m$variable=="2050, SSP2", boundone_m$wq_ssp2_2050, boundone_m$wq_ssp5_2050))

boundone_m$wq = boundone_m$wq_hist

boundone_m = filter(boundone_m, value<1)

ggplot(data=boundone_m, aes(x=as.factor(wq), y=value))+
  theme_classic()+
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", alpha=0.8, position = "dodge", fill="lightblue") + 
  facet_wrap(vars(country), scales = "free", ncol=3)+
  ylab("Proportion of household expenditure for total electricity consumption")+
  xlab("Total household expenditure quintile in the survey year")+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_npg(name="")+
  theme(legend.position = "bottom", legend.direction = "horizontal", strip.background = element_blank())

ggsave(paste0(stub, "results/energy_poverty/ely_hist_share.png"), width = 6, height = 6, scale=1.2)

boundone_m_bk2 = boundone_m

#######################

rdss = list.files(path=paste0(stub, "data/household"), recursive = T, full.names = T, pattern = "rds", ignore.case = T)

rdss = rdss[!grepl("global", rdss)]
rdss = rdss[!grepl("CDD", rdss)]
rdss = rdss[!grepl("us.Rds", rdss)]
rdss = rdss[!grepl("nss", rdss)]
rdss = rdss[!grepl("nss", rdss)]
rdss = rdss[!grepl("submission", rdss)]
rdss <- rdss[c(1,2,3,4,5,6,7,8,10,11,12,13,15,16,18,19,20, 21)]

rdss_l = lapply(rdss, read_rds)

for (i in 1:length(rdss_l)){
  
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ln_ely_q), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ac), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ln_total_exp_usd_2011), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$mean_CDD_db), ]
  #rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$urban), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ownership_d), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$n_members), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$age_head), ]
  #rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$country), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$weight), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$sex_head), ]
  #rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$urban_sh), ]
  
  rdss_l[[i]]$hhid = as.character(rdss_l[[i]]$hhid)
  
  rdss_l[[i]] = dplyr::select(rdss_l[[i]], weight, ln_total_exp_usd_2011, ln_ely_q, starts_with("ely_p"), hhid)
  rdss_l[[i]]$country <- basename(rdss[[i]])
}

rdss_l <- bind_rows(rdss_l)

#tapply(rdss_l$ely_p_usd_2011, rdss_l$country, summary)

########

load(paste0(stub, "results/regressions/for_projections/global_wgt_dmcf.Rdata"))

future_ac_adoption = read_rds(paste0(stub, "results/household_level/", "global_ac.Rds"))
output_ac = read_rds(paste0(stub, "results/household_level/", "global_ely_tot.Rds"))
output_impact_ac = read_rds(paste0(stub, "results/household_level/", "global_ely_due_to_ac.Rds"))

colnames(future_ac_adoption) = paste0("future_ac_adoption_", colnames(future_ac_adoption))
colnames(output_ac) = paste0("output_ac_", colnames(output_ac))
colnames(output_impact_ac) = paste0("output_impact_ac_", colnames(output_impact_ac))

###

weights = read_rds(paste0(stub, "results/household_level/", "weights_global_ac.Rds"))

data_orig = read_rds(paste0(stub, "results/household_level/", "data_global_ac.Rds"))

data_orig = bind_cols(data_orig, future_ac_adoption, output_ac, output_impact_ac, weights)

rdss_l$weight = as.numeric(rdss_l$weight)

rdss_l$ln_ely_q = round(rdss_l$ln_ely_q, 3)
data_orig$ln_ely_q = round(data_orig$ln_ely_q, 3)
data_orig$ln_total_exp_usd_2011 = NULL

###

rdss_l$ely_p_usd_2011 <- ifelse(is.na(rdss_l$ely_p_usd_2011), rdss_l$ely_p, rdss_l$ely_p_usd_2011)
rdss_l$rdss_l$ely_p <-NULL
rdss_l$ely_p_impute_5_95 <-NULL

###

rdss_l = merge(data_orig, rdss_l, by=c("ln_ely_q", "weight"))

rdss_l$total_exp_usd_2011 = exp(rdss_l$ln_total_exp_usd_2011)
rdss_l$ely_p_usd_2011 = rdss_l$ely_p_usd_2011.x
rdss_l$ely_q = exp(rdss_l$ln_ely_q)

rdss_l =  rdss_l %>% group_by(country.x) %>% dplyr::mutate(q_hist=ntile(exp(ln_total_exp_usd_2011), 5), q_ssp2_2050 =ntile(exp(exp_cap_usd_SSP2_2050), 5), q_ssp5_2050 =ntile(exp(exp_cap_usd_SSP5_2050), 5))

rdss_l = dplyr::select(rdss_l, country.x, starts_with("q_"), ely_p_usd_2011, starts_with("exp_cap"), starts_with("future_ac_adoption"),starts_with("output_ac"),starts_with("output_impact_ac"),starts_with("weight_"), hhid.x, total_exp_usd_2011)

pr1 = as.numeric(predict(reg_ely))
load(paste0(stub, "results/regressions/for_projections/global_wgt_dmcf.RData"))

sec$ac <- sec$ac[1]
sec$ac <- as.numeric(sec$ac)
pr2 = as.numeric(predict(reg_ely, newdata=sec))
diff_exp = exp(pr1) - exp(pr2)
load(paste0(stub, "results/regressions/for_projections/global_wgt_dmcf.RData"))

diff_exp = ifelse(diff_exp<=0 | sec$ac =="0" , NA, diff_exp)

diff_exp = data.frame(diff_exp=diff_exp, hhid=as.character(sec$hhid), country.x=sec$country)

#rdss_l_bk = rdss_l

rdss_l$hhid = rdss_l$hhid.x

rdss_l = merge(rdss_l, diff_exp, c("country.x", "hhid"))

rdss_l$output_impact_ac_SSP2.2020 = rdss_l$diff_exp
rdss_l$output_impact_ac_SSP5.2020 = rdss_l$diff_exp
rdss_l$exp_cap_usd_SSP2_2020 = log(rdss_l$total_exp_usd_2011+1)
rdss_l$exp_cap_usd_SSP5_2020 = log(rdss_l$total_exp_usd_2011+1)

rdss_l_o = data.frame(country=c(rdss_l$country.x))

rdss_l_o$wq_hist=rdss_l$q_hist
rdss_l_o$wq_ssp2_2050=rdss_l$q_ssp2_2050
rdss_l_o$wq_ssp5_2050=rdss_l$q_ssp5_2050

for (ssp in c("SSP2", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  for (year in seq(2020, 2050, 10)){
    
    rdss_l_o[,paste0("ely_share_", ssp, "_", (year))] = (rdss_l[,paste0("output_impact_ac_", ssp, ".", (year))] * rdss_l$ely_p_usd_2011) / exp(rdss_l[,paste0("exp_cap_usd_", ssp, "_", (year))])
    
  }}

#

rdss_l_o = split(rdss_l_o, rdss_l_o$country)
rdss_l = split(rdss_l, rdss_l$country.x)

for (kappa in 1:length(rdss_l)){
  print(kappa)
  rdss_l_o[[kappa]][,paste0("ely_share_", ssp, "_", (year))] = rdss_l_o[[kappa]][,paste0("ely_share_", ssp, "_", (year))] * (rdss_l[[kappa]][,paste0("weight_", year, "_", rcp, "_", tolower(ssp))]/mean(rdss_l[[kappa]][,paste0("weight_", year, "_", rcp, "_", tolower(ssp))], na.rm=T))
  
}

rdss_l_o = bind_rows(rdss_l_o)


#

#rdss_l_o = na.omit(rdss_l_o)

rdss_l_o_bk = rdss_l_o

rdss_l_o = reshape2::melt(rdss_l_o, c(1:4))

boundone_m <- filter(rdss_l_o, variable=="ely_share_SSP2_2020")

boundone_m$variable <- as.character(boundone_m$variable)
boundone_m$variable[boundone_m$variable =="ely_share_SSP2_2020" ] <- " Survey year"
boundone_m$variable[boundone_m$variable =="ely_share_SSP2_2050" ] <- "2050, SSP2"
boundone_m$variable[boundone_m$variable =="ely_share_SSP5_2050" ] <- "2050, SSP5"

unique(boundone_m$country)
unique(boundone_m$variable)
boundone_m$country = as.character(boundone_m$country)

boundone_m$country[boundone_m$country=="Ghana" | boundone_m$country=="Niger" | boundone_m$country=="Nigeria" | boundone_m$country=="Tanzania"| boundone_m$country=="Malawi"| boundone_m$country=="Burkina Faso"| boundone_m$country=="Kenya"] = "Africa"
boundone_m$country[boundone_m$country=="Australia" | boundone_m$country=="Japan" | boundone_m$country=="Canada"] = "OECD-NONEU"
boundone_m$country[boundone_m$country=="France" | boundone_m$country=="Germany"| boundone_m$country=="Netherlands"| boundone_m$country=="Spain"| boundone_m$country=="Sweden"] = "OECD-EU"

##############

calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

#boundone_m$wq = ifelse(boundone_m$variable==" Survey year", boundone_m$wq_hist, ifelse(boundone_m$variable=="2050, SSP2", boundone_m$wq_ssp2_2050, boundone_m$wq_ssp5_2050))

boundone_m$wq = boundone_m$wq_hist

boundone_m = filter(boundone_m, value<1)

ggplot(data=boundone_m, aes(x=as.factor(wq), y=value))+
  theme_classic()+
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", alpha=0.8, position = "dodge", fill="lightblue") + 
  facet_wrap(vars(country), scales = "free", ncol=3)+
  ylab("Proportion of household expenditure for AC electricity consumption")+
  xlab("Total household expenditure quintile in the survey year")+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_npg(name="")+
  theme(legend.position = "bottom", legend.direction = "horizontal", strip.background = element_blank())

ggsave(paste0(stub, "results/energy_poverty/ely_hist_share_wac_ACZOOM.png"), width = 6, height = 6, scale=1.2)

boundone_m_bk = boundone_m

####

rdss = list.files(path=paste0(stub, "data/household"), recursive = T, full.names = T, pattern = "rds", ignore.case = T)

rdss = rdss[!grepl("global", rdss)]
rdss = rdss[!grepl("CDD", rdss)]
rdss = rdss[!grepl("us.Rds", rdss)]
rdss = rdss[!grepl("nss", rdss)]
rdss = rdss[!grepl("nss", rdss)]
rdss = rdss[!grepl("submission", rdss)]
rdss <- rdss[c(1,2,3,4,5,6,7,8,10,11,12,13,15,16,18,19,20, 21)]

rdss_l = lapply(rdss, read_rds)

for (i in 1:length(rdss_l)){
  
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ln_ely_q), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ac), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ln_total_exp_usd_2011), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$mean_CDD_db), ]
  #rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$urban), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$ownership_d), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$n_members), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$age_head), ]
  #rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$country), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$weight), ]
  rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$sex_head), ]
  #rdss_l[[i]] <- rdss_l[[i]][complete.cases(rdss_l[[i]]$urban_sh), ]
  
  rdss_l[[i]]$hhid = as.character(rdss_l[[i]]$hhid)
  
  rdss_l[[i]] = dplyr::select(rdss_l[[i]], weight, ln_total_exp_usd_2011, ln_ely_q, starts_with("ely_p"), hhid)
  rdss_l[[i]]$country <- basename(rdss[[i]])
}

rdss_l <- bind_rows(rdss_l)

#tapply(rdss_l$ely_p_usd_2011, rdss_l$country, summary)

########

load(paste0(stub, "results/regressions/for_projections/global_wgt_dmcf.Rdata"))

future_ac_adoption = read_rds(paste0(stub, "results/household_level/", "global_ac.Rds"))
output_ac = read_rds(paste0(stub, "results/household_level/", "global_ely_tot.Rds"))
output_impact_ac = read_rds(paste0(stub, "results/household_level/", "global_ely_due_to_ac.Rds"))

colnames(future_ac_adoption) = paste0("future_ac_adoption_", colnames(future_ac_adoption))
colnames(output_ac) = paste0("output_ac_", colnames(output_ac))
colnames(output_impact_ac) = paste0("output_impact_ac_", colnames(output_impact_ac))

###

weights = read_rds(paste0(stub, "results/household_level/", "weights_global_ac.Rds"))

data_orig = read_rds(paste0(stub, "results/household_level/", "data_global_ac.Rds"))

data_orig = bind_cols(data_orig, future_ac_adoption, output_ac, output_impact_ac, weights)

rdss_l$weight = as.numeric(rdss_l$weight)

rdss_l$ln_ely_q = round(rdss_l$ln_ely_q, 3)
data_orig$ln_ely_q = round(data_orig$ln_ely_q, 3)
data_orig$ln_total_exp_usd_2011 = NULL

###

rdss_l$ely_p_usd_2011 <- ifelse(is.na(rdss_l$ely_p_usd_2011), rdss_l$ely_p, rdss_l$ely_p_usd_2011)
rdss_l$rdss_l$ely_p <-NULL
rdss_l$ely_p_impute_5_95 <-NULL

###

rdss_l = merge(data_orig, rdss_l, by=c("ln_ely_q", "weight"))

rdss_l$total_exp_usd_2011 = exp(rdss_l$ln_total_exp_usd_2011)
rdss_l$ely_p_usd_2011 = rdss_l$ely_p_usd_2011.x
rdss_l$ely_q = exp(rdss_l$ln_ely_q)

rdss_l =  rdss_l %>% group_by(country.x) %>% dplyr::mutate(q_hist=ntile(exp(ln_total_exp_usd_2011), 5), q_ssp2_2050 =ntile(exp(exp_cap_usd_SSP2_2050), 5), q_ssp5_2050 =ntile(exp(exp_cap_usd_SSP5_2050), 5))

rdss_l = dplyr::select(rdss_l, country.x, starts_with("q_"), ely_p_usd_2011, starts_with("exp_cap"), starts_with("future_ac_adoption"),starts_with("output_ac"),starts_with("output_impact_ac"),starts_with("weight_"), hhid.x, total_exp_usd_2011)

pr1 = as.numeric(predict(reg_ely))
load(paste0(stub, "results/regressions/for_projections/global_wgt_dmcf.RData"))

sec$ac <- sec$ac[1]
sec$ac <- as.numeric(sec$ac)
pr2 = as.numeric(predict(reg_ely, newdata=sec))
diff_exp = exp(pr1)
load(paste0(stub, "results/regressions/for_projections/global_wgt_dmcf.RData"))
diff_exp = ifelse(diff_exp<=0 | sec$ac==0 , NA, diff_exp)

diff_exp = data.frame(diff_exp=diff_exp, hhid=as.character(sec$hhid), country.x=sec$country)

#rdss_l_bk = rdss_l

rdss_l$hhid = rdss_l$hhid.x

rdss_l = merge(rdss_l, diff_exp, c("country.x", "hhid"))

rdss_l$output_impact_ac_SSP2.2020 = rdss_l$diff_exp
rdss_l$output_impact_ac_SSP5.2020 = rdss_l$diff_exp
rdss_l$exp_cap_usd_SSP2_2020 = log(rdss_l$total_exp_usd_2011+1)
rdss_l$exp_cap_usd_SSP5_2020 = log(rdss_l$total_exp_usd_2011+1)

rdss_l_o = data.frame(country=c(rdss_l$country.x))

rdss_l_o$wq_hist=rdss_l$q_hist
rdss_l_o$wq_ssp2_2050=rdss_l$q_ssp2_2050
rdss_l_o$wq_ssp5_2050=rdss_l$q_ssp5_2050

for (ssp in c("SSP2", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  for (year in seq(2020, 2050, 10)){
    
    rdss_l_o[,paste0("ely_share_", ssp, "_", (year))] = (rdss_l[,paste0("output_impact_ac_", ssp, ".", (year))] * rdss_l$ely_p_usd_2011) / exp(rdss_l[,paste0("exp_cap_usd_", ssp, "_", (year))])
    
  }}

#

rdss_l_o = split(rdss_l_o, rdss_l_o$country)
rdss_l = split(rdss_l, rdss_l$country.x)

for (kappa in 1:length(rdss_l)){
  print(kappa)
  rdss_l_o[[kappa]][,paste0("ely_share_", ssp, "_", (year))] = rdss_l_o[[kappa]][,paste0("ely_share_", ssp, "_", (year))] * (rdss_l[[kappa]][,paste0("weight_", year, "_", rcp, "_", tolower(ssp))]/mean(rdss_l[[kappa]][,paste0("weight_", year, "_", rcp, "_", tolower(ssp))], na.rm=T))
  
}

rdss_l_o = bind_rows(rdss_l_o)


#

#rdss_l_o = na.omit(rdss_l_o)

rdss_l_o_bk = rdss_l_o

rdss_l_o = reshape2::melt(rdss_l_o, c(1:4))

boundone_m <- filter(rdss_l_o, variable=="ely_share_SSP2_2020")

boundone_m$variable <- as.character(boundone_m$variable)
boundone_m$variable[boundone_m$variable =="ely_share_SSP2_2020" ] <- " Survey year"
boundone_m$variable[boundone_m$variable =="ely_share_SSP2_2050" ] <- "2050, SSP2"
boundone_m$variable[boundone_m$variable =="ely_share_SSP5_2050" ] <- "2050, SSP5"

unique(boundone_m$country)
unique(boundone_m$variable)
boundone_m$country = as.character(boundone_m$country)

boundone_m$country[boundone_m$country=="Ghana" | boundone_m$country=="Niger" | boundone_m$country=="Nigeria" | boundone_m$country=="Tanzania"| boundone_m$country=="Malawi"| boundone_m$country=="Burkina Faso"| boundone_m$country=="Kenya"] = "Africa"
boundone_m$country[boundone_m$country=="Australia" | boundone_m$country=="Japan" | boundone_m$country=="Canada"] = "OECD-NONEU"
boundone_m$country[boundone_m$country=="France" | boundone_m$country=="Germany"| boundone_m$country=="Netherlands"| boundone_m$country=="Spain"| boundone_m$country=="Sweden"] = "OECD-EU"

##############

calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

#boundone_m$wq = ifelse(boundone_m$variable==" Survey year", boundone_m$wq_hist, ifelse(boundone_m$variable=="2050, SSP2", boundone_m$wq_ssp2_2050, boundone_m$wq_ssp5_2050))

boundone_m$wq = boundone_m$wq_hist

boundone_m = filter(boundone_m, value<1)

boundone_m_bk$type = "AC electricity"
boundone_m_bk2$type = " Total electricity (all HHs)"
boundone_m$type = "  Total electricity (HHs with AC)"

boundone_m_b = bind_rows(boundone_m,boundone_m_bk, boundone_m_bk2)

boundone_m_b = filter(boundone_m_b, value<1)

ggplot(data=boundone_m_b, aes(x=as.factor(wq), y=value, fill=type))+
  theme_classic()+
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", alpha=0.8, position = "dodge") + 
  facet_wrap(vars(country), scales = "free", ncol=3)+
  ylab("Proportion of household total expenditure for (AC) electricity consumption")+
  xlab("Total household expenditure quintile in the survey year")+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_npg(name="")+
  theme(legend.position = "bottom", legend.direction = "horizontal", strip.background = element_blank())

ggsave(paste0(stub, "results/energy_poverty/ely_hist_share_wac.png"), width = 6, height = 6, scale=1.2)
