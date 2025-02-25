
##########################################

#               Figure S2

##########################################

rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()

# Load packages
library(data.table)
library(plyr)
library(dplyr)
library(FSA)
library(stringr)
library(tidyverse)
library(wbstats)
library(ggrepel)
#library(sjPlot)

# Set users
user <- 'user'

if (user=='user') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'data/household/', sep='')
proj <- paste(stub,'results/household_level/', sep='')
output <- paste(stub,'output/figures/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/supplementary/'


# Load AC projections at 2020
setwd("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/household_level")

proj_f = list.files(pattern="national_ac_penetration.csv")
proj = lapply(proj_f, read.csv)

for (i in 1:length(proj)){
  
  proj[[i]]$country <- gsub("_national_ac_penetration.csv", "", basename(proj_f)[i])
  proj[[i]]$type <- "All drivers"
  

}

proj <- bind_rows(proj)

#

proj_gm_f = list.files(pattern="_national_ac_penetration_glomod.csv")
proj_gm = lapply(proj_gm_f, read.csv)

for (i in 1:length(proj_gm)){
  
  proj_gm[[i]]$country <- gsub("_national_ac_penetration_glomod.csv", "", basename(proj_gm_f)[i])
  proj_gm[[i]]$type <- "All drivers"
  
  
}

proj_gm <- bind_rows(proj_gm)

# Load household data
global <- readRDS(paste(house,'global.rds', sep=''))

# Check
global <- global[complete.cases(global$ln_ely_q), ]
global <- global[complete.cases(global$ac), ]
global <- global[complete.cases(global$ln_total_exp_usd_2011), ]
global <- global[complete.cases(global$mean_CDD_db), ]
global <- global[complete.cases(global$ownership_d), ]
global <- global[complete.cases(global$n_members), ]
global <- global[complete.cases(global$age_head), ]
global <- global[complete.cases(global$country), ]
global <- global[complete.cases(global$weight), ]

# Obtain current statistics
stat = global %>% group_by(country) %>%  dplyr::summarise(stat=weighted.mean(as.numeric(as.character(ac)), weight,  na.rm = T))

stat$country_c <- c("ARG", "OECD_NONEU","BRA", "Africa",  "OECD_NONEU", "CHN", "OECD_EU", "OECD_EU", "Africa","IND", "IDN", "ITA", "OECD_NONEU", "Africa", "Africa", "MEX", "OECD_EU", "Africa", "Africa", "PAK", "OECD_EU", "OECD_EU", "OECD_EU", "Africa", "USA")

# pop and weighting 
pop = wbstats::wb_data(indicator = "SP.POP.TOTL", mrv=1)

stat$iso3c <- countrycode::countrycode(stat$country, 'country.name', 'iso3c')

stat = merge(stat, pop, by.x="iso3c", by.y="iso3c")

stat = stat %>% group_by(country_c) %>%  dplyr::summarise(stat=weighted.mean(stat, SP.POP.TOTL))

#

proj = filter(proj_gm, year==2020)

bound = merge(stat, proj, by.x="country_c", by.y="country")

bound <- group_by(bound, country_c) %>% dplyr::summarise(stat=mean(stat), value=mean(value))

bound$bias = round(bound$value - bound$stat, 3)

# Figure S2
g1 <- ggplot(bound)+
    geom_abline()+
    geom_point(aes(x=value*100, y=stat*100))+
  geom_label_repel(aes(x=value*100, y=stat*100, label=paste0(country_c, "\n", bias*100, " %")))+
  xlab("Modeled AC prevalence, 2020 (SSP245)")+
  ylab("Surveyed AC penetration rate (survey year)")

g1

# Save
ggsave(paste(output, 'FigureS2.png', sep = ''), g1, scale=.75, height = 10, width = 10)
