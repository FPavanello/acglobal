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

output <- 'C:/Users/Utente/Downloads/acglobal/output/tables/'


ll_list <- list.files(path=paste0(stub, "repo/interm/projections/"), pattern="national_ac_penetration_glomod.csv", full.names = T)

ll_list <- ll_list[!grepl("DEU", ll_list)]

ll <- lapply(ll_list, read.csv)

for (i in 1:length(ll)){
  ll[[i]]$country <- gsub("_national_ac_penetration_glomod.csv", "", basename(ll_list[i]))
}

ll <- as.data.frame(do.call(rbind, ll))
ll$X <- NULL
ll$stat <- "AC penetration rate"

#

ll2_list <- list.files(path=paste0(stub, "repo/interm/projections/"), pattern="_national_ac_consumption_glomod.csv", full.names = T)

ll2_list <- ll2_list[!grepl("DEU", ll2_list)]

ll2 <- lapply(ll2_list, read.csv)

for (i in 1:length(ll2)){
  ll2[[i]]$country <- gsub("_national_ac_consumption_glomod.csv", "", basename(ll2_list[i]))
  ll2[[i]]$value <- ll2[[i]]$value_tot
  ll2[[i]]$value_tot <- NULL
}

ll2 <- as.data.frame(do.call(rbind, ll2))
ll2$X <- NULL
ll2$stat <- "Total AC-induced consumption"

#

ll3_list <- list.files(path=paste0(stub, "repo/interm/projections/"), pattern="_national_ac_consumption_total_glomod.csv", full.names = T)

ll3_list <- ll3_list[!grepl("DEU", ll3_list)]

ll3 <- lapply(ll3_list, read.csv)

for (i in 1:length(ll3)){
  ll3[[i]]$country <- gsub("_national_ac_consumption_total_glomod.csv", "", basename(ll3_list[i]))
}

ll3 <- as.data.frame(do.call(rbind, ll3))
ll3$X <- NULL
ll3$stat <- "Per-capita AC consumption"

##

ll <- bind_rows(ll, ll2, ll3)
ll <- dplyr::filter(ll, year==2020 | year==2050)
ll$ssp <- as.character(ll$ssp)
ll$ssp <- ifelse(ll$year==2020, "2020", ll$ssp)
ll$ssp <- ifelse(ll$ssp=="SSP245", "SSP245 (2050)", ll$ssp)
ll$ssp <- ifelse(ll$ssp=="SSP585", "SSP585 (2050)", ll$ssp)

ll <- group_by(ll, ssp, year, country, stat) %>% dplyr::summarise(value=median(value, na.rm=T))

ll$country <- stringr::str_remove(ll$country, "projections")

ll$country <- ifelse(ll$country=="GLOBAL", " GLOBAL pool", ll$country)

#

############

# Calculate carbon emissions

require(data.table)
ci <- fread(paste0(stub, "rscripts/projections/carbon_intensity/AR6_Scenarios_Database_ISO3_v1.0.csv"), header = T)

emis_lab <- unique(ci$Variable)[grep("Electricity", unique(ci$Variable))][88]
ci_emis <- filter(ci, Variable %in% emis_lab)
ci_emis <- filter(ci_emis, Scenario %in% c("SSP2_BASE", "SSP5-baseline"))
ci_emis <- dplyr::select(ci_emis, Scenario, Region, Unit, `2020`, `2050`)
ci_emis$`2050` <- ci_emis$`2050` * 1000000000000  # convert to grams
ci_emis$`2020` <- ci_emis$`2020` * 1000000000000  # convert to grams

ely_lab <- unique(ci$Variable)[grep("Electricity", unique(ci$Variable))][136:152]
ci_ely <- filter(ci, Variable %in% ely_lab[c(1:5, 9:10, 12, 15:17)])
ci_ely <- filter(ci_ely, Scenario %in% c("SSP2_BASE", "SSP5-baseline"))
ci_ely <- group_by(ci_ely, Scenario, Region, Unit) %>% dplyr::summarise(`2020` = sum(`2020`, na.rm=T), `2050`=sum(`2050`, na.rm=T))
ci_ely$`2050` <- ci_ely$`2050` * 277777777777.78 # convert to kWh
ci_ely$`2020` <- ci_ely$`2020` * 277777777777.78 # convert to kWh

ci_emis <- ci_emis[with(ci_emis, order(Region, Scenario)), ]
ci_ely <- ci_ely[with(ci_ely, order(Region, Scenario)), ]

ci_emis$ci_intensity_2020 <- ci_emis$`2020` / ci_ely$`2020`
ci_emis$ci_intensity_2050 <- ci_emis$`2050` / ci_ely$`2050`

ll_emis <- ll
ll_emis <- filter(ll_emis, stat=="Total AC-induced consumption")

# parse regions

ll_emis$country_bk <- ll_emis$country
ll_emis$country[ll_emis$country=="ARG"] <- "BRA"
ll_emis$country[ll_emis$country=="Africa"] <- "ZAF"
ll_emis$country[ll_emis$country=="DEU"] <- "EU"
ll_emis$country[ll_emis$country=="ITA"] <- "EU"
ll_emis$country[ll_emis$country=="OECD_EU"] <- "EU"
ll_emis$country[ll_emis$country=="OECD_NONEU"] <- "JPN"
ll_emis$country[ll_emis$country=="PAK"] <- "IND"

# merge

ci_emis <- pivot_longer(ci_emis, 6:7)
ci_emis$name <- gsub("ci_intensity_", "", ci_emis$name)
ci_emis$Scenario <- ifelse(ci_emis$Scenario=="SSP5-baseline", "SSP585 (2050)", "SSP245 (2050)")
ci_emis$Scenario <- ifelse(ci_emis$name=="2020", "2020", ci_emis$Scenario)

colnames(ci_emis)[7] <- "emission_factor"
colnames(ll_emis)[5] <- "ely_ac"

merger <- merge(ll_emis, ci_emis, by.x=c("ssp", "country"), by.y=c("Scenario", "Region"), all.x=T)

merger <- group_by(merger, year, country, country_bk) %>% mutate(`2020`=ifelse(is.na(`2020`), mean(`2020`, na.rm=T), `2020`), `2050`=ifelse(is.na(`2050`), mean(`2050`, na.rm=T), `2050`), emission_factor=ifelse(is.na(emission_factor), mean(emission_factor, na.rm=T), emission_factor)) %>% ungroup()

merger$emission_factor[merger$country==" GLOBAL pool"] <- weighted.mean(merger$emission_factor, merger$`2050`, na.rm=T)

merger$mtco2 <- merger$ely_ac * merger$emission_factor / 1000000000000 # convert to Mt

merger <- group_by(merger, ssp, country, country_bk) %>% dplyr::summarise(mtco2=mean(mtco2, na.rm=T))

merger$country <- NULL
colnames(merger)[2] <- "country"
merger$stat <- "CO2 emissions"
merger$year <- ifelse(merger$ssp==2020, 2020, 2050)
colnames(merger)[3] <- "value"

ll_emis <- bind_rows(ll, merger)
ll_emis <- filter(ll_emis, stat=="CO2 emissions")

ll_emis$stat <- "Total CO2 emissions (Mt/yr)"

ll2 <- ll_emis %>% 
  mutate_if(is.numeric, round, 1)

ll2$year <- NULL 

library(modelsummary)

emptycol <- function(x) " "

datasummary(value * country ~ stat * ssp * (Mean),
            data = ll2, output = paste0(output, "TableA15.tex"))
#############

ll$year <- NULL 

library(modelsummary)

ll$value <- ifelse(ll$stat=="Total AC-induced consumption", ll$value/1e9, ll$value)
ll$value <- ifelse(ll$stat=="AC penetration rate", ll$value*100, ll$value)

ll$stat <- ifelse(ll$stat=="Total AC-induced consumption", "Total AC electr. (TWh)", ifelse(ll$stat=="Per-capita AC consumption", "Per cap. AC electr. (avg. kWh/hh/yr)", ifelse(ll$stat=="AC penetration rate", "AC penetr.rate (%)", "Total CO2 emissions (Mt/yr)")))

ll2 <- ll %>% 
  mutate_if(is.numeric, round, 1)

emptycol <- function(x) " "

datasummary(value * country ~ stat * ssp * (Mean),
            data = ll2, output = paste0(output, "Table4.tex"))
