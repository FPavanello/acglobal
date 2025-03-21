
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
library(spatstat)
library(patchwork)
#library(sjPlot)

# Set users
user <- 'user'

if (user=='user') {
  stub <- "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste0(stub, "repo/household/")
output <- 'C:/Users/Utente/Downloads/acglobal/output/supplementary/'

# Load ely projections at 2020
setwd(paste0(stub, "repo/interm/projections"))

proj_f = list.files(pattern="_national_ely_consumption.csv")
proj = lapply(proj_f, read.csv)

for (i in 1:length(proj)){
  
  proj[[i]]$country <- gsub("_national_ely_consumption.csv", "", basename(proj_f)[i])
  proj[[i]]$type <- "All drivers"
  

}

proj <- bind_rows(proj)

proj_gm_f = list.files(pattern="_national_ely_consumption_glomod.csv")
proj_gm = lapply(proj_gm_f, read.csv)

for (i in 1:length(proj_gm)){
  
  proj_gm[[i]]$country <- gsub("_national_ely_consumption_glomod.csv", "", basename(proj_gm_f)[i])
  proj_gm[[i]]$type <- "All drivers"
  
  
}

proj_gm <- bind_rows(proj_gm)

# Load global data
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
gl <- data.frame(country=unique(global$country), country_c=c("Africa", "Africa", "Africa", "Africa", "Africa", "Africa", "Africa", "ARG", "BRA", "CHN", "OECD_EU", "IND", "IDN", "ITA", "MEX", "PAK", "OECD_EU", "OECD_EU", "OECD_EU", "OECD_NONEU", "OECD_NONEU", "OECD_NONEU", "OECD_EU", "OECD_EU", "USA"))
  
global <- merge(global, gl, 'country')

stat = global %>% group_by(country_c) %>%  dplyr::summarise(stat=weighted.mean(ely_q, weight, na.rm = T))

proj = filter(proj, year==2020)

bound = merge(stat, proj, by.x="country_c", by.y="country")

bound <- group_by(bound, country_c) %>% dplyr::summarise(stat=mean(stat), value=mean(value))

bound$bias = round(((bound$value - bound$stat)/bound$stat)*100, 0)

# Figure S2 - country specific models
g1 <- ggplot(bound)+
    geom_abline()+
    geom_point(aes(x=value, y=stat))+
  geom_label_repel(aes(x=value, y=stat, label=paste0(country_c, "\n", bias)))+
  xlab("Modeled electricity consumption, kWh/HH/yr, 2020")+
  ylab("Suveyed electricity consumption, kWh/HH/yr (survey year)")+
  scale_x_log10()+
  scale_y_log10()

#

proj_gm = filter(proj_gm, year==2020)

bound = merge(stat, proj_gm, by.x="country_c", by.y="country")

bound <- group_by(bound, country_c) %>% dplyr::summarise(stat=mean(stat), value=mean(value))

bound$bias = round(((bound$value - bound$stat)/bound$stat)*100, 0)

# Figure S2 - global model
g2 <- ggplot(bound)+
  geom_abline()+
  geom_point(aes(x=value, y=stat))+
  geom_label_repel(aes(x=value, y=stat, label=paste0(country_c, "\n", bias)))+
  xlab("Modeled electricity consumption, kWh/HH/yr, 2020")+
  ylab("Surveyed electricity consumption, kWh/HH/yr (survey year)")+
  scale_x_log10()+
  scale_y_log10()

g1+g2

# Save
ggsave(paste(output, 'FigureS3.png', sep = ''), last_plot(), scale=1.4, height = 5.5, width = 10)
