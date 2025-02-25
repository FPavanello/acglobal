
##########################################

#               Figure A4

##########################################

rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
library(tidyverse)
library(data.table)
library(patchwork)
library(ggsci)

# Set users
user <- 'user'

if (user=='user') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/"
}

output <- paste(stub,'output/figures/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/figures/'

setwd(stub)

# infer observed value from data

global <- readRDS(paste(stub,'6-Projections/data/household/global.rds', sep=''))

global_sum_stat <- global

global_sum_stat$country_macro= NA
global_sum_stat$country_macro[global_sum_stat$country=="Burkina Faso" | global_sum_stat$country=="Ghana"| global_sum_stat$country=="Kenya"| global_sum_stat$country=="Malawi"| global_sum_stat$country=="Niger"| global_sum_stat$country=="Nigeria"| global_sum_stat$country=="Tanzania"] <- "Africa"
global_sum_stat$country_macro[global_sum_stat$country=="Australia" | global_sum_stat$country=="Canada"| global_sum_stat$country=="Japan"] <- "OECD_NONEU"
global_sum_stat$country_macro[global_sum_stat$country=="France" | global_sum_stat$country=="Netherlands"| global_sum_stat$country=="Spain"| global_sum_stat$country=="Sweden" | global_sum_stat$country=="Switzerland" | global_sum_stat$country=="Germany"] <- "OECD_EU"

global_sum_stat <- group_by(global_sum_stat, country_macro) %>% dplyr::mutate(weight = ifelse(is.na(weight), mean(weight, na.rm=T), weight)) %>%  dplyr::summarise(ac=weighted.mean(as.numeric(as.character(ac)), weight, na.rm=T))

global_sum_stat_gp <- group_by(global) %>% dplyr::mutate(weight = ifelse(is.na(weight), mean(weight, na.rm=T), weight)) %>%  dplyr::summarise(ac=weighted.mean(as.numeric(as.character(ac)), weight, na.rm=T))

global_sum_stat_gp$country_macro <- "Global pool"

global_sum_stat <- bind_rows(global_sum_stat, global_sum_stat_gp)
rm(global_sum_stat_gp)

setwd(paste0(stub, "6-Projections/results/household_level"))

#########

l1 <- list.files(pattern="national_ac_penetration_glomod.csv", full.names = T, recursive = F)
l1_wsd <- list.files(pattern="national_ac_penetration_glomod_wsd.csv", full.names = T, recursive = F)

l2 <- list.files(pattern="_national_ac_consumption_glomod.csv", full.names = T, recursive = F)
l2_wsd <- list.files(pattern="_national_ac_consumption_glomod_wsd.csv", full.names = T, recursive = F)

l3 <- list.files(pattern="_national_ac_consumption_total_glomod.csv", full.names = T, recursive = F)
l3_wsd <- list.files(pattern="_national_ac_consumption_total_glomod_wsd.csv", full.names = T, recursive = F)

l1 <- l1[!grepl("DEU", l1)]
l1_wsd <- l1_wsd[!grepl("DEU", l1_wsd)]
l2 <- l2[!grepl("DEU", l2)]
l2_wsd <- l2_wsd[!grepl("DEU", l2_wsd)]
l3 <- l3[!grepl("DEU", l3)]
l3_wsd <- l3_wsd[!grepl("DEU", l3_wsd)]

l1_d <- lapply(l1, read.csv)
l1_wsd_d <- lapply(l1_wsd, read.csv)

l2_d <- lapply(l2, read.csv)
l2_wsd_d <- lapply(l2_wsd, read.csv)

l3_d <- lapply(l3, read.csv)
l3_wsd_d <- lapply(l3_wsd, read.csv)

###############

for (i in 1:length(l1)){
  
  l1_d[[i]]$country <- gsub("_national_ac_penetration_glomod.csv", "", basename(l1)[i])
  l1_d[[i]]$type <- "All drivers"
  l1_d[[i]]$country <- ifelse(l1_d[[i]]$country=="GLOBAL", "Global pool", l1_d[[i]]$country)
  
  l1_wsd_d[[i]]$country <- gsub("_national_ac_penetration_glomod.csv", "", basename(l1)[i])
  l1_wsd_d[[i]]$type <- "Only macro drivers"
  l1_wsd_d[[i]]$country <- ifelse(l1_wsd_d[[i]]$country=="GLOBAL", "Global pool", l1_wsd_d[[i]]$country)
  
  
}

l1_d <- bind_rows(l1_d) %>%
  group_by(ssp, country, type) %>% 
  group_modify(~ add_row(.x, X=0, year=2010 ,.before=0))


for (ctry in unique(l1_d$country)){
l1_d$value[l1_d$year==2010 & l1_d$country==ctry] <- rep(global_sum_stat$ac[global_sum_stat$country_macro==ctry], 2)
}

l1_d$value[l1_d$year==2010 & l1_d$country=="Global pool"] <- predict(lm("value ~ year", data = as.data.frame(l1_d[l1_d$country=="Global pool" & l1_d$year!=2010,])), data.frame(year=2010))


l1_wsd_d <- bind_rows(l1_wsd_d) %>%
  group_by(ssp, country, type) %>% 
  group_modify(~ add_row(.x, X=0, year=2010 ,.before=0))


for (ctry in unique(l1_wsd_d$country)){
  l1_wsd_d$value[l1_wsd_d$year==2010 & l1_wsd_d$country==ctry] <- rep(global_sum_stat$ac[global_sum_stat$country_macro==ctry], 2)
}

l1_wsd_d$value[l1_wsd_d$year==2010 & l1_wsd_d$country=="Global pool"] <- predict(lm("value ~ year", data = as.data.frame(l1_wsd_d[l1_wsd_d$country=="Global pool" & l1_wsd_d$year!=2010,])), data.frame(year=2010))

  
l1 <- bind_rows(l1_d, l1_wsd_d)

l1 <- l1 %>%
  group_by(ssp, country) %>% 
  mutate(value = ifelse((value < lead(value)) | year == 2050, value, lead(value)))

g1 <- ggplot(l1 %>% filter(country!="Global pool"))+
  geom_line(aes(x=year, y=value*100, colour=ssp, linetype=type))+
  facet_wrap(vars(country), scales = "free_y")+
  xlab("Year")+
  ylab("AC penetration rate (%)")

g1_g <- ggplot(l1 %>% filter(country=="Global pool"))+
  geom_line(aes(x=year, y=value*100, colour=ssp, linetype=type))+
  facet_wrap(vars(country), scales = "free_y")+
  xlab("Year")+
  ylab("AC penetration rate (%)")

#ggsave("6-Projections/results/graphs/ac_penetration_sensitivity_drivers.png", scale=1)

#

for (i in 1:length(l2)){
  
  l2_d[[i]]$country <- gsub("_national_ac_consumption_glomod.csv", "", basename(l2)[i])
  l2_d[[i]]$type <- "All drivers"
  l2_d[[i]]$country <- ifelse(l2_d[[i]]$country=="GLOBAL", "Global pool", l2_d[[i]]$country)
  
  l2_wsd_d[[i]]$country <- sub("_national_ac_consumption_glomod.csv", "", basename(l2)[i])
  l2_wsd_d[[i]]$type <- "Only macro drivers"
  l2_wsd_d[[i]]$country <- ifelse(l2_wsd_d[[i]]$country=="GLOBAL", "Global pool", l2_wsd_d[[i]]$country)
  
  
}

l2 <- bind_rows(bind_rows(l2_d), bind_rows(l2_wsd_d))

g2 <- ggplot(l2 %>% filter(country!="Global pool"))+
  geom_line(aes(x=year, y=value_tot/1e9, colour=ssp, linetype=type))+
  facet_wrap(vars(country), scales="free_y")+
  xlab("Year")+
  ylab("National AC-induced ely. cons. (TWh/yr)")

g2_g <- ggplot(l2 %>% filter(country=="Global pool"))+
  geom_line(aes(x=year, y=value_tot/1e9, colour=ssp, linetype=type))+
  facet_wrap(vars(country), scales="free_y")+
  xlab("Year")+
  ylab("National AC-induced ely. cons. (TWh/yr)")

#

for (i in 1:length(l3)){
  
  l3_d[[i]]$country <- gsub("_national_ac_consumption_total_glomod.csv", "", basename(l3)[i])
  l3_d[[i]]$type <- "All drivers"
  l3_d[[i]]$country <- ifelse(l3_d[[i]]$country=="GLOBAL", "Global pool", l3_d[[i]]$country)
  
  l3_wsd_d[[i]]$country <- gsub("_national_ac_consumption_total_glomod.csv", "", basename(l3)[i])
  l3_wsd_d[[i]]$type <- "Only macro drivers"
  l3_wsd_d[[i]]$country <- ifelse(l3_wsd_d[[i]]$country=="GLOBAL", "Global pool", l3_wsd_d[[i]]$country)
  
  
}

l3 <- bind_rows(bind_rows(l3_d), bind_rows(l3_wsd_d))

g3 <- ggplot(l3 %>% filter(country!="Global pool"))+
  geom_line(aes(x=year, y=value, colour=ssp, linetype=type))+
  facet_wrap(vars(country), scales="free_y")+
  xlab("Year")+
  ylab("Avg. per-capita AC-induced ely. cons. (kWh/hh/yr)")

g3_g <- ggplot(l3 %>% filter(country=="Global pool"))+
  geom_line(aes(x=year, y=value, colour=ssp, linetype=type))+
  facet_wrap(vars(country), scales="free_y")+
  xlab("Year")+
  ylab("Avg. per-capita AC-induced ely. cons. (kWh/hh/yr)")

(g1_g + g2_g + g1 + g2)  + plot_layout(guides = "collect")+ plot_annotation(tag_levels = "A") & theme_classic() & theme(legend.position = "bottom", legend.direction = "horizontal", axis.text.x = element_text(angle = 90))

# Save
ggsave(paste0(output, "FigureA4.png"), scale=1.8, width = 6, height = 5)

