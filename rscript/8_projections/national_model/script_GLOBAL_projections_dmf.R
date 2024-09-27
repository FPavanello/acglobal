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
library(raster)
library(rnaturalearthdata)
library(sf)
library(gdata)
library(exactextractr)
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


# Use CDDs 18 threshold?
CDD18=T
CDD_type = ifelse(CDD18==T, "cdd_18", "cdd_24"); HDD_type = ifelse(CDD18==T, "hdd_18", "hdd_15")




# country specific data / parameters

countryname <- "GLOBAL"
countryiso3 <- "GLOBAL"

# 1 import DMF environment with trained models and data

ll <- list.files(paste0(stub, "results/drivers_evolution"), pattern = "Rdata", full.names = T)[!grepl("2", list.files(paste0(stub, "results/drivers_evolution"), pattern = "Rdata"))]

listone <- list()

for (i in 1:length(ll)){
  print(i)
  load(ll[i])
  data_c_sp_export$ac <- as.factor(as.character(data_c_sp_export$ac))
  data_c_sp_export$hhid <- as.character(data_c_sp_export$hhid)
  data_c_sp_export$ownership_d <- as.factor(as.character(data_c_sp_export$ownership_d))
  data_c_sp_export$sex_head <- as.numeric(as.character(data_c_sp_export$sex_head))
  data_c_sp_export$sex_head <- as.factor(as.character(ifelse(data_c_sp_export$sex_head==2, 1, data_c_sp_export$sex_head)))
  listone[[i]] <- data_c_sp_export
  rm(data_c_sp_export)
  rm(orig_data)
}

listone[[2]]$country <- listone[[2]]$country.x
listone[[3]]$country <- listone[[3]]$country.x
listone[[4]]$country <- listone[[4]]$country.x
listone[[5]]$country <- listone[[5]]$country.x
listone[[6]]$country <- "Indonesia"
listone[[7]]$country <- "India"
listone[[8]]$country <- listone[[8]]$country.x
listone[[9]]$country <- listone[[9]]$country.x
listone[[10]]$country <- listone[[10]]$country.x
listone[[11]]$country <- listone[[11]]$country.x
listone[[12]]$country <- listone[[12]]$country.x
listone[[13]]$country <- listone[[13]]$country.x

ll_cols <- Reduce(intersect, lapply(listone, colnames))

for (i in 1:length(ll)){
  
  listone[[i]] <- dplyr::select(listone[[i]] , all_of(ll_cols))

}

listone <- bind_rows(listone)

CDD18=T
if(CDD18==T){colnames(listone) = gsub("mean_CDD18_db", "mean_CDD18_db", colnames(listone) )}

# filter some countries
#listone = filter(listone, country!="Kenya")

load(paste0(stub, "results/regressions/for_projections/global_wgt_dmcf.Rdata"))

## 3) Make projections based on trained models and extracted data ##
# 3.1) AC adoption projections

listone$age_SSP1_2010 <- NULL
listone$age_SSP2_2010 <- NULL
listone$age_SSP3_2010 <- NULL
listone$age_SSP4_2010 <- NULL
listone$age_SSP5_2010 <- NULL

listone$weighted_mean.URB_gr_SSP1_2010 <- NULL
listone$weighted_mean.URB_gr_SSP2_2010 <- NULL
listone$weighted_mean.URB_gr_SSP3_2010 <- NULL
listone$weighted_mean.URB_gr_SSP4_2010 <- NULL
listone$weighted_mean.URB_gr_SSP5_2010 <- NULL

listone <- listone %>% drop_na(starts_with("exp_cap_usd_"), starts_with("mean_CDD_"), starts_with("mean_HDD_"), starts_with("edu_"), starts_with("age_"), starts_with("weighted_mean.URB_"))

orig_data <- listone
orig_data$ac <- NULL

orig_data_bk <- orig_data
data_c_sp <- listone

output <- list()

# loop for all ssps and time-steps

for (ssp in c("SSP2", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  orig_data_bk <- orig_data
  
  output2 <- list()
  
  for (year in seq(2020, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$country <- as.factor(data_c_sp$country)
    
    orig_data$ln_total_exp_usd_2011 = data_c_sp[,paste0("exp_cap_usd_", ssp, "_", (year))]
    
    orig_data$mean_CDD18_db  = data_c_sp[,paste0("mean_CDD_", year, "_", rcp, "_", tolower(ssp))] / 100 # 10-year average bins
    
    orig_data$mean_CDD18_db2 = orig_data$mean_CDD18_db^2
    
    orig_data$mean_CDD18_db_exp = orig_data$mean_CDD18_db*orig_data$ln_total_exp_usd_2011
    
    orig_data$mean_CDD18_db2_exp = (orig_data$mean_CDD18_db^2)*orig_data$ln_total_exp_usd_2011
    
    orig_data$curr_CDD18_db = orig_data$mean_CDD18_db 
    
    orig_data$curr_CDD18_db2 = orig_data$curr_CDD18_db^2
    
    orig_data$ln_ely_p = log(orig_data$ely_p_usd_2011)
    
    orig_data$ln_ely_p_cdd = orig_data$mean_CDD18_db*orig_data$ln_ely_p 
    
    orig_data$ln_ely_p_cdd2 = (orig_data$mean_CDD18_db^2)*orig_data$ln_ely_p 
    
    orig_data$ln_ely_p_nme = orig_data$ln_ely_p * orig_data$n_members
      
    orig_data$ln_ely_p_own = orig_data$ln_ely_p * as.numeric(as.character(orig_data$n_members))
  
    orig_data$mean_HDD18_db  = data_c_sp[,paste0("mean_HDD_", year, "_", rcp, "_", tolower(ssp))] / 100  # 10-year average bins
    
    orig_data$edu_head_2 <- as.factor(data_c_sp[,paste0("edu_", year, "_", ssp)])
    
    orig_data$age_head <- round(orig_data$age_head * (1 + data_c_sp[,paste0("age_", ssp, "_", year)]), 0)
  
    #orig_data$ownership_d = as.factor(data_c_sp[,paste0("ownership_d_", ssp, "_", (year))])
    
    #orig_data$fan = as.factor(data_c_sp[,paste0("fan_", ssp, "_", (year))])
    
    #orig_data$housing_index_lab = as.factor(data_c_sp[,paste0("housing_index_lab_", year, "_s1")])
    
    orig_data$urban_sh = data_c_sp[,paste0("weighted_mean.URB_", ssp, "_", (year))]    
    #
    projected <- predict(reg_ac, orig_data, type="response")
    # projected <- ifelse(as.numeric(projected)>0.5, 1, 0)
    # 
    # if (year>2020){
    #   projected <- ifelse(output2[[as.character(year-10)]]==1, 1, projected)
    # }
    
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output[[as.character(ssp)]] <- output2
  
}

###

# bind results together
future_ac_adoption <- unlist(output, recursive = FALSE)
future_ac_adoption <- as.data.frame(do.call("cbind", future_ac_adoption))

data_c_sp <- bind_cols(data_c_sp, future_ac_adoption)

data_c_sp$geometry.x <- NULL
data_c_sp$geometry.y <- NULL

#############

# project weights (dynamic loops for names)

for (ssp in c("SSP2", "SSP5")){
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  for (year in seq(2020, 2050, 10)){
    
    if (year == 2020){
      
      data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))] <- data_c_sp$weight 
      
    } else {
      
      
      data_c_sp[,paste0("pop_gr", year, "_", rcp, "_", tolower(ssp))] = (data_c_sp[,paste0("sum.POP_", ssp, "_", year)] / data_c_sp[,paste0("sum.POP_", ssp, "_", (year-10))]) - 1
      
      
      data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))] = (1 + data_c_sp[,paste0("pop_gr", year, "_", rcp, "_", tolower(ssp))]) *  data_c_sp[,paste0("weight_", (year - 10), "_", rcp, "_", tolower(ssp))]
      
    }}}

    

future_ac_adoption_g <- split(future_ac_adoption, data_c_sp$country)

for (ctry in 1:length(future_ac_adoption_g)){
  
  for (ssp in c("SSP2", "SSP5")){
    
    rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
    
    for (year in seq(2020, 2050, 10)){
      
      
      future_ac_adoption_g[[ctry]][,paste0(ssp, ".", year, "_wgt")] =  future_ac_adoption_g[[ctry]][,paste0(ssp, ".", year)] * (data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))][data_c_sp$country==names(future_ac_adoption_g)[ctry]] / sum(data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))][data_c_sp$country==names(future_ac_adoption_g)[ctry]], na.rm=T))
      

    }}
  
  future_ac_adoption_g[[ctry]] <- dplyr::select(future_ac_adoption_g[[ctry]], 9:16)
  
  }


national_summary_ac_g <- list()

for (ctry in 1:length(future_ac_adoption_g)){

  national_summary_ac_g[[ctry]] <- future_ac_adoption_g[[ctry]] %>%
  dplyr::summarise_at(vars(contains("wgt")), sum, na.rm=T) %>%
  pivot_longer(cols = 1:8, names_to = c('ssp', 'year'), names_sep = ".") %>%
  mutate(ssp=rep(c("SSP245","SSP585"), each=4), year=rep(seq(2020, 2050, 10), 2))

}


###########################################

# 3.2) Electricity consumption projections
# 3.2.1) Predict consumption without AC

orig_data <- listone
orig_data$ely_q <- NULL

orig_data_bk <- orig_data

output <- list()

# loop for all ssps and time-steps

for (ssp in c("SSP2", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  output2 <- list()
  
  for (year in seq(2020, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$country <- as.factor(data_c_sp$country)
    
    orig_data$ac = as.factor(0)
    
    orig_data$ln_ely_p <- log(data_c_sp$ely_p_usd_2011)
    
    orig_data$ln_total_exp_usd_2011 = data_c_sp[,paste0("exp_cap_usd_", ssp, "_", (year))]
    
    orig_data$curr_CDD18_db  = data_c_sp[,paste0("mean_CDD_", year, "_", rcp, "_", tolower(ssp))] / 100 # 10-year average bins
    
    orig_data$curr_HDD18_db  = data_c_sp[,paste0("mean_HDD_", year, "_", rcp, "_", tolower(ssp))] / 100  # 10-year average bins
    
  orig_data$edu_head_2 <- as.factor(data_c_sp[,paste0("edu_", year, "_", ssp)])
    
    orig_data$age_head <- round(orig_data$age_head * (1 + data_c_sp[,paste0("age_", ssp, "_", year)]), 0)
    
    orig_data$ownership_d = as.factor(orig_data$ownership_d)
    
    #orig_data$fan = as.factor(data_c_sp[,paste0("fan_", ssp, "_", (year))])
    
    #orig_data$housing_index_lab = as.factor(data_c_sp[,paste0("housing_index_lab_", year, "_s1")])
    
    orig_data$urban_sh = data_c_sp[,paste0("weighted_mean.URB_", ssp, "_", (year))]    
    #
    projected <- as.numeric(predict(reg_ely, orig_data))
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output[[as.character(ssp)]] <- output2
  
}

###

# bind results together
output <- unlist(output, recursive = FALSE)
output_noac <- as.data.frame(do.call("cbind", output))


# 3.2.2) Predict consumption with AC

orig_data <- listone
orig_data$ely_q <- NULL

orig_data_bk <- orig_data

output <- list()

# loop for all ssps and time-steps

for (ssp in c("SSP2", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  output2 <- list()
  
  for (year in seq(2020, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$country <- as.factor(data_c_sp$country)
    
    orig_data$ac = as.factor(ifelse(data_c_sp[,paste0(ssp, ".", (year))]>0.5, 1, 0))
    
    orig_data$ln_ely_p <- log(data_c_sp$ely_p_usd_2011)
    
    orig_data$ln_total_exp_usd_2011 = data_c_sp[,paste0("exp_cap_usd_", ssp, "_", (year))]
    
    orig_data$curr_CDD18_db  = data_c_sp[,paste0("mean_CDD_", year, "_", rcp, "_", tolower(ssp))] / 100 # 10-year average bins
    
    orig_data$curr_HDD18_db  = data_c_sp[,paste0("mean_HDD_", year, "_", rcp, "_", tolower(ssp))] / 100  # 10-year average bins
    
    orig_data$edu_head_2 <- as.factor(data_c_sp[,paste0("edu_", year, "_", ssp)])
    
    orig_data$age_head <- round(orig_data$age_head * (1 + data_c_sp[,paste0("age_", ssp, "_", year)]), 0)
    
    orig_data$ownership_d = as.factor(orig_data$ownership_d)
    
    #orig_data$fan = as.factor(data_c_sp[,paste0("fan_", ssp, "_", (year))])
    
    #orig_data$housing_index_lab = as.factor(data_c_sp[,paste0("housing_index_lab_", year, "_s1")])
    
    orig_data$urban_sh = data_c_sp[,paste0("weighted_mean.URB_", ssp, "_", (year))]    
    #
    projected <- as.numeric(predict(reg_ely, orig_data))
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output[[as.character(ssp)]] <- output2
  
}


###

# bind results together
output <- unlist(output, recursive = FALSE)
output_ac <- as.data.frame(do.call("cbind", output))

#####
# 3.2.3) Estimate impact of AC ownership on electricity consumption

#### weighted statistics!!!

output_impact_ac <- exp(output_ac) - exp(output_noac)
output_impact_ac[output_impact_ac<=0] <- NA

output_impact_ac_g <- split(output_impact_ac, data_c_sp$country)

for (ctry in 1:length(output_impact_ac_g)){

for (ssp in c("SSP2", "SSP5")){
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  for (year in seq(2020, 2050, 10)){
    
    output_impact_ac_g[[ctry]][,paste0(ssp, ".", year)] =  output_impact_ac_g[[ctry]][,paste0(ssp, ".", year)] * (data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))][data_c_sp$country==names(output_impact_ac_g)[ctry]] / mean(data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))][data_c_sp$country==names(output_impact_ac_g)[ctry]][!is.na(output_impact_ac_g[[ctry]][,paste0(ssp, ".", year)])], na.rm=T))
    
  }}}

national_summary_cons_g <- list()

for (ctry in 1:length(output_impact_ac_g)){

  national_summary_cons_g[[ctry]] <- output_impact_ac_g[[ctry]] %>%
    dplyr::summarise_all(mean, na.rm=T) %>%
    pivot_longer(cols = 1:8, names_to = c('ssp', 'year'), names_sep = ".") %>%
    mutate(ssp=rep(c("SSP245","SSP585"), each=4), year=rep(seq(2020, 2050, 10), 2))
  
  national_summary_cons_g[[ctry]][is.na(national_summary_cons_g[[ctry]])] <- 0
  
 }

#

national_summary_total <- list()

for (ctry in 1:length(output_impact_ac_g)){
  
national_summary_total[[ctry]] <- output_impact_ac_g[[ctry]] %>%
  dplyr::summarise_all(mean, na.rm=T) %>%
  pivot_longer(cols = 1:8, names_to = c('ssp', 'year'), names_sep = ".") %>%
  dplyr::mutate(ssp=rep(c("SSP245","SSP585"), each=4), year=rep(seq(2020, 2050, 10), 2))

national_summary_total[[ctry]][is.na(national_summary_total[[ctry]])] <- 0

}

#

# elaborate to get total consumption

custom_shape <- data_c_sp[!duplicated(data_c_sp[,c('sum.POP_SSP1_2010')]),]

pop_long <- custom_shape %>% dplyr::select(country, starts_with("sum.POP")) %>% group_by(country) %>%  dplyr::summarise_all(., "sum") %>% pivot_longer(2:51)
pop_long$name <- gsub("sum.POP_", "", pop_long$name )
pop_long$ssp <- substr(pop_long$name, 1, 4)
pop_long$year <- substr(pop_long$name, 6, 9)
pop_long <- filter(pop_long, year>2010 & year<=2050)

for (country in unique(custom_shape$country)){

pop_long$value[pop_long$country==country] <- pop_long$value[pop_long$country==country] / weighted.mean(data_c_sp$n_members[data_c_sp$country==country], data_c_sp$weight[data_c_sp$country==country], na.rm=T) # https://globaldatalab.org/areadata/hhsize

}

pop_long <- filter(pop_long, grepl("SSP2", name) | grepl("SSP5", name))

names(national_summary_ac_g) <-  unique(custom_shape$country)

national_summary_ac <- bind_rows(national_summary_ac_g, .id = "country")

pop_long$value_orig <- pop_long$value

# multiply by AC ownership
pop_long$value =  pop_long$value * national_summary_ac$value

names(national_summary_ac_g) <-  unique(custom_shape$country)

national_summary_cons <- bind_rows(national_summary_cons_g, .id = "country")

# multiply by average consumption due to AC
national_summary_cons$value_tot <- pop_long$value * national_summary_cons$value

national_summary_cons <- group_by(national_summary_cons, ssp, year) %>% dplyr::summarise(value_tot = sum(value_tot, na.rm=T))

# aggregate to global scale

divider <- bind_rows(national_summary_ac_g, .id = "country")
divider$value <- divider$value * pop_long$value
divider <- group_by(divider, ssp, year) %>% dplyr::summarise(value = sum(value, na.rm=T))

national_summary_cons$value <- national_summary_total$value

national_summary_ac$ll1 <- pop_long$value
national_summary_ac$ll2 <- pop_long$value_orig

national_summary_ac <- group_by(national_summary_ac, ssp, year) %>% dplyr::summarise(value = sum(ll1)/sum(ll2))

national_summary_total <- bind_rows(national_summary_total, .id = "country")

national_summary_total$ll1 <- pop_long$value

national_summary_total <- group_by(national_summary_total, ssp, year) %>% dplyr::summarise(value = weighted.mean(value, ll1))

#######

weights = dplyr::select(data_c_sp, starts_with("weight_"))

write_rds(orig_data, paste0(stub, "results/household_level/", "data_global_ac.Rds"))
write_rds(weights, paste0(stub, "results/household_level/", "weights_global_ac.Rds"))

write_rds(future_ac_adoption,  paste0(stub, "results/household_level/", "global_ac.Rds"))
write_rds(output_ac,  paste0(stub, "results/household_level/", "global_ely_tot.Rds"))
write_rds(output_impact_ac,  paste0(stub, "results/household_level/", "global_ely_due_to_ac.Rds"))

###############
# export projections data

write.csv(national_summary_ac, paste0(stub, "results/household_level/", countryiso3, "_national_ac_penetration.csv"))
write.csv(national_summary_cons,  paste0(stub, "results/household_level/", countryiso3, "_national_ac_consumption.csv"))
write.csv(national_summary_total,  paste0(stub, "results/household_level/", countryiso3, "_national_ac_consumption_total.csv"))

###########

output_ac <- exp(output_ac)
output_ac[output_ac<=0] <- NA

output_ac_g <- split(output_ac, data_c_sp$country)

for (ctry in 1:length(output_ac_g)){
  
  for (ssp in c("SSP2", "SSP5")){
    
    rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
    
    for (year in seq(2020, 2050, 10)){
      
      output_ac_g[[ctry]][,paste0(ssp, ".", year)] =  output_ac_g[[ctry]][,paste0(ssp, ".", year)] * (data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))][data_c_sp$country==names(output_ac_g)[ctry]] / mean(data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))][data_c_sp$country==names(output_ac_g)[ctry]][!is.na(output_ac_g[[ctry]][,paste0(ssp, ".", year)])], na.rm=T))
      
    }}}

national_summary_cons_g <- list()

for (ctry in 1:length(output_ac_g)){
  
  national_summary_cons_g[[ctry]] <- output_ac_g[[ctry]] %>%
    dplyr::summarise_all(mean, na.rm=T) %>%
    pivot_longer(cols = 1:8, names_to = c('ssp', 'year'), names_sep = ".") %>%
    mutate(ssp=rep(c("SSP245","SSP585"), each=4), year=rep(seq(2020, 2050, 10), 2))
  
  national_summary_cons_g[[ctry]][is.na(national_summary_cons_g[[ctry]])] <- 0
  
}

#

national_summary_total <- list()

for (ctry in 1:length(output_ac_g)){
  
  national_summary_total[[ctry]] <- output_ac_g[[ctry]] %>%
    dplyr::summarise_all(mean, na.rm=T) %>%
    pivot_longer(cols = 1:8, names_to = c('ssp', 'year'), names_sep = ".") %>%
    dplyr::mutate(ssp=rep(c("SSP245","SSP585"), each=4), year=rep(seq(2020, 2050, 10), 2))
  
  national_summary_total[[ctry]][is.na(national_summary_total[[ctry]])] <- 0
  
}

#

# elaborate to get total consumption

custom_shape <- data_c_sp[!duplicated(data_c_sp[,c('sum.POP_SSP1_2010')]),]

pop_long <- custom_shape %>% dplyr::select(country, starts_with("sum.POP")) %>% group_by(country) %>%  dplyr::summarise_all(., "sum") %>% pivot_longer(2:51)
pop_long$name <- gsub("sum.POP_", "", pop_long$name )
pop_long$ssp <- substr(pop_long$name, 1, 4)
pop_long$year <- substr(pop_long$name, 6, 9)
pop_long <- filter(pop_long, year>2010 & year<=2050)

for (country in unique(custom_shape$country)){
  
  pop_long$value[pop_long$country==country] <- pop_long$value[pop_long$country==country] / weighted.mean(data_c_sp$n_members[data_c_sp$country==country], data_c_sp$weight[data_c_sp$country==country], na.rm=T) # https://globaldatalab.org/areadata/hhsize
  
}

pop_long <- filter(pop_long, grepl("SSP2", name) | grepl("SSP5", name))

national_summary_ac <- bind_rows(national_summary_ac_g, .id = "country")

pop_long$value_orig <- pop_long$value

# multiply by AC ownership
pop_long$value =  pop_long$value * national_summary_ac$value

national_summary_cons <- bind_rows(national_summary_cons_g, .id = "country")

# multiply by average consumption due to AC
national_summary_cons$value_tot <- pop_long$value * national_summary_cons$value

national_summary_cons <- group_by(national_summary_cons, ssp, year) %>% dplyr::summarise(value_tot = sum(value_tot, na.rm=T))

# aggregate to global scale

divider <- bind_rows(national_summary_ac_g, .id = "country")
divider$value <- divider$value * pop_long$value
divider <- group_by(divider, ssp, year) %>% dplyr::summarise(value = sum(value, na.rm=T))

national_summary_cons$value <- national_summary_total$value

national_summary_ac$ll1 <- pop_long$value
national_summary_ac$ll2 <- pop_long$value_orig

national_summary_ac <- group_by(national_summary_ac, ssp, year) %>% dplyr::summarise(value = sum(ll1)/sum(ll2))

national_summary_total <- bind_rows(national_summary_total, .id = "country")

national_summary_total$ll1 <- pop_long$value

national_summary_total <- group_by(national_summary_total, ssp, year) %>% dplyr::summarise(value = weighted.mean(value, ll1))

national_summary_cons$value <- national_summary_total$value

write.csv(national_summary_cons,  paste0(stub, "results/household_level/", countryiso3, "_national_ely_consumption.csv"))
