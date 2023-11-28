
## This R-script:
##      1) descriptive statistics whole data set (Table 2)


rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

## 1) Load libraries and data ##
library(data.table)
library(xtable)
#library(foreign)
#library(tidyverse)
#library(haven)
#library(reshape2)
#library(survey)
#library(stargazer)
library(reldist)
library(Hmisc)

# Set users
user <- 'fp'
#user <- 'gf'

if (user=='fp') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

if (user=='gf') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'data/household/', sep='') 
output <- paste(stub,'results/graphs/tables/', sep='') 


###################################

# TABLE 2: Descriptive Statistics #

###################################

# Load global data
#load(paste0(stub, "results/regressions/for_projections/global_dmcf.RData"))
global <- readRDS(paste(house,'global.rds', sep=''))

# Check
global <- global[complete.cases(global$ln_ely_q), ]
global <- global[complete.cases(global$ac), ]
global <- global[complete.cases(global$ln_total_exp_usd_2011), ]
global <- global[complete.cases(global$mean_CDD18_db), ]
global <- global[complete.cases(global$urban), ]
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

# Number of households per country
summary(global$country)

# Transform log data
global <- dplyr::mutate(global, ely = exp(ln_ely_q), total_exp = exp(ln_total_exp_usd_2011))

# Correlation for omitted category CDD18 and CDD24
global <- dplyr::mutate(global, omitted = curr_CDD18_db - curr_CDD_db)
cor(as.numeric(global$ac), global$omitted, use = "complete.obs") # -0.3

# Survey
global_svy <- svydesign(data = global, ids = ~1, weights = ~ weight)


## Summary global data set
# Vector of variables
variables <- c("ely", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "edu_head_2", "age_head", "sex_head")

# Load function
devtools::source_gist("c4d1089a501d3567be9fb784b1c5a6ab")

# Define 
datasets <- list("Global" = global)
variable_names <- list(variables)
labels <- list(variables)
colnames <- c("Mean", "SD", "10th", "25th", "Median", "75th", "90th")

# Function to get descriptives
myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, global$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, global$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, global$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, global$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, global$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, global$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, global$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

## Full sample
# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics",
                       label = "tab:desc_all",
                       file = paste(output,'desc_global.tex', sep=''))


## Africa
# Filter
africa <- dplyr::filter(global, country == "Burkina Faso" | country == "Ghana" | country == "Kenya" | 
                                country == "Malawi" | country == "Niger" | country == "Nigeria")

# Function to get descriptives
datasets <- list("Africa" = africa)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, africa$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, africa$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, africa$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, africa$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, africa$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, africa$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, africa$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- Africa",
                       label = "tab:desc_afr",
                       file = paste(output,'desc_africa.tex', sep=''))

## Argentina
# Filter
argentina <- dplyr::filter(global, country == "Argentina")

# Function to get descriptives
datasets <- list("Argentina" = argentina)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, argentina$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, argentina$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, argentina$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, argentina$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, argentina$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, argentina$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, argentina$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- Argentina",
                       label = "tab:desc_arg",
                       file = paste(output,'desc_arg.tex', sep=''))

## Brazil
# Filter
brazil <- dplyr::filter(global, country == "Brazil")

# Function to get descriptives
datasets <- list("Brazil" = brazil)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, brazil$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, brazil$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, brazil$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, brazil$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, brazil$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, brazil$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, brazil$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- Brazil",
                       label = "tab:desc_bra",
                       file = paste(output,'desc_bra.tex', sep=''))

## China
# Filter
china <- dplyr::filter(global, country == "China")

# Function to get descriptives
datasets <- list("China" = china)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, china$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, china$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, china$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, china$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, china$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, china$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, china$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- China",
                       label = "tab:desc_chn",
                       file = paste(output,'desc_chn.tex', sep=''))

## Germany
# Filter
germany <- dplyr::filter(global, country == "Germany")

# Function to get descriptives
datasets <- list("Germany" = germany)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, germany$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, germany$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, germany$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, germany$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, germany$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, germany$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, germany$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- Germany",
                       label = "tab:desc_deu",
                       file = paste(output,'desc_deu.tex', sep=''))

## Indonesia
# Filter
indonesia <- dplyr::filter(global, country == "Indonesia")

# Function to get descriptives
datasets <- list("Indonesia" = indonesia)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, indonesia$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, indonesia$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, indonesia$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, indonesia$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, indonesia$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, indonesia$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, indonesia$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- Indonesia",
                       label = "tab:desc_idn",
                       file = paste(output,'desc_idn.tex', sep=''))

## India
# Filter
india <- dplyr::filter(global, country == "India")

# Function to get descriptives
datasets <- list("India" = india)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, india$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, india$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, india$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, india$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, india$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, india$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, india$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- India",
                       label = "tab:desc_ind",
                       file = paste(output,'desc_ind.tex', sep=''))

## Italy
# Filter
italy <- dplyr::filter(global, country == "Italy")

# Function to get descriptives
datasets <- list("Italy" = italy)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, italy$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, italy$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, italy$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, italy$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, italy$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, italy$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, italy$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- Italy",
                       label = "tab:desc_ita",
                       file = paste(output,'desc_ita.tex', sep=''))

## Mexico
# Filter
mexico <- dplyr::filter(global, country == "Mexico")

# Function to get descriptives
datasets <- list("Mexico" = mexico)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, mexico$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, mexico$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, mexico$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, mexico$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, mexico$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, mexico$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, mexico$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- Mexico",
                       label = "tab:desc_mex",
                       file = paste(output,'desc_mex.tex', sep=''))

## EPIC - EU
# Filter
eu <- dplyr::filter(global, country == "France" | country == "Netherlands" | country == "Switzerland" |
                      country == "Sweden" | country == "Spain")

# Function to get descriptives
datasets <- list("EU" = eu)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, eu$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, eu$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, eu$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, eu$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, eu$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, eu$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, eu$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- EPIC-EU",
                       label = "tab:desc_eu",
                       file = paste(output,'desc_eu.tex', sep=''))

## EPIC - NONEU
# Filter
noneu <- dplyr::filter(global, country == "Australia" | country == "Canada" | country == "Japan")

# Function to get descriptives
datasets <- list("NONEU" = noneu)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, noneu$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, noneu$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, noneu$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, noneu$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, noneu$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, noneu$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, noneu$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- EPIC-NONEU",
                       label = "tab:desc_noneu",
                       file = paste(output,'desc_noneu.tex', sep=''))

## Pakistan
# Filter
pak <- dplyr::filter(global, country == "Pakistan")

# Function to get descriptives
datasets <- list("Pakistan" = pak)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, pak$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, pak$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, pak$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, pak$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, pak$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, pak$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, pak$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- Pakistan",
                       label = "tab:desc_pak",
                       file = paste(output,'desc_pak.tex', sep=''))


## United States
# Filter
usa <- dplyr::filter(global, country == "United States")

# Function to get descriptives
datasets <- list("United States" = usa)

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, usa$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, usa$weight, na.rm = TRUE))
  df <- reldist::wtd.quantile (x, q=0.10, usa$weight, na.rm = TRUE)
  tf <- reldist::wtd.quantile (x, q=0.25, usa$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.50, usa$weight, na.rm = TRUE)
  sf <- reldist::wtd.quantile (x, q=0.75, usa$weight, na.rm = TRUE)
  nf <- reldist::wtd.quantile (x, q=0.90, usa$weight, na.rm = TRUE)
  return(c(m, sd, df, tf, md, sf, nf))
}

# Table
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics --- United States",
                       label = "tab:desc_usa",
                       file = paste(output,'desc_usa.tex', sep=''))








#############

as.num <- function(X){as.numeric(as.character(X))}

list_country_files <- list.files(path=paste0(stub, "results/regressions/for_projections"), full.names = T, pattern = "RData")[-c(6,7)]

load(list_country_files[1])
HH_Africa <- HH_Africa %>% mutate_if(is.factor,as.num)
HH_Africa <- dplyr::select(HH_Africa, -country, -state, -phat0_obs, -ac_obs, -xb_noac, -selection, -hhid)
HH_Africa <- filter(HH_Africa, ln_ely_q>0 & ln_ely_q<13)
stargazer(HH_Africa, type="latex", out="desc_afr.tex", digits = 2)

load(list_country_files[2])
HH_Argentina <- HH_Argentina %>% mutate_if(is.factor,as.num)
HH_Argentina <- dplyr::select(HH_Argentina, -country, -state, -subregion, -phat0_obs, -ac_obs, -xb_noac, -selection, -hhid)
HH_Argentina <- dplyr::select(HH_Argentina, colnames(HH_Africa))
HH_Argentina <- filter(HH_Argentina, ln_ely_q>0 & ln_ely_q<13)
stargazer(HH_Argentina, type="latex", out="desc_arg.tex", digits = 2)

load(list_country_files[3])
HH_Brazil <- HH_Brazil %>% mutate_if(is.factor,as.num)
HH_Brazil <- dplyr::select(HH_Brazil, -state, -region3, -phat0_obs, -ac_obs, -xb_noac, -selection, -hhid)
HH_Brazil <- filter(HH_Brazil, ln_ely_q>0 & ln_ely_q<13)
HH_Brazil <- dplyr::select(HH_Brazil, colnames(HH_Africa))
stargazer(HH_Brazil, type="latex", out="desc_bra.tex", digits = 2)

load(list_country_files[4])
HH_China <- HH_China %>% mutate_if(is.factor,as.num)
HH_China <- dplyr::select(HH_China, -state, -state3, -macroarea, -phat0_obs, -ac_obs, -xb_noac, -selection, -hhid)
HH_China <- filter(HH_China, ln_ely_q>0 & ln_ely_q<13)
HH_China <- dplyr::select(HH_China, colnames(HH_Africa))
stargazer(HH_China, type="latex", out="desc_chn.tex", digits = 2)

load(list_country_files[5])
HH_Germany <- HH_Germany %>% mutate_if(is.factor,as.num)
HH_Germany <- dplyr::select(HH_Germany, -state, -country, -phat0_obs, -ac_obs, -xb_noac, -selection, -hhid)
HH_Germany <- filter(HH_Germany, ln_ely_q>0 & ln_ely_q<13)
HH_Germany <- dplyr::select(HH_Germany, colnames(HH_Africa)[-13])
stargazer(HH_Germany, type="latex", out="desc_deu.tex", digits = 2)

load(list_country_files[6])
HH_Indonesia <- HH_Indonesia %>% mutate_if(is.factor,as.num)
HH_Indonesia <- dplyr::select(HH_Indonesia, -state, -phat0_obs, -ac_obs, -xb_noac, -selection)
HH_Indonesia <- filter(HH_Indonesia, ln_ely_q>0 & ln_ely_q<13)
HH_Indonesia <- dplyr::select(HH_Indonesia, colnames(HH_Africa))
stargazer(HH_Indonesia, type="latex", out="desc_idn.tex", digits = 2)

load(list_country_files[7])
HH_India <- HH_India %>% mutate_if(is.factor,as.num)
HH_India <- dplyr::select(HH_India, -state, -phat0_obs, -ac_obs, -xb_noac, -selection)
HH_India <- filter(HH_India, ln_ely_q>0 & ln_ely_q<13)
HH_India <- dplyr::select(HH_India, colnames(HH_Africa)[-13])
stargazer(HH_India, type="latex", out="desc_ind.tex", digits = 2)

load(list_country_files[8])
HH_Italy <- HH_Italy %>% mutate_if(is.factor,as.num)
HH_Italy <- dplyr::select(HH_Italy, -state, -home_tenure, -n_rooms, -house_type, -macroarea, -selection)
HH_Italy <- filter(HH_Italy, ln_ely_q>0 & ln_ely_q<13)
HH_Italy <- dplyr::select(HH_Italy, colnames(HH_Africa)[-c(6,13,8,9)])
stargazer(HH_Italy, type="latex", out="desc_ita.tex", digits = 2)

load(list_country_files[9])
HH_Mexico <- HH_Mexico %>% mutate_if(is.factor,as.num)
HH_Mexico <- dplyr::select(HH_Mexico, -state, -phat0_obs, -ac_obs, -xb_noac, -selection)
HH_Mexico <- filter(HH_Mexico, ln_ely_q>0 & ln_ely_q<13)
HH_Mexico <- dplyr::select(HH_Mexico, colnames(HH_Africa))
stargazer(HH_Mexico, type="latex", out="desc_mex.tex", digits = 2)

load(list_country_files[10])
HH_Europe <- HH_Europe %>% mutate_if(is.factor,as.num)
HH_Europe <- dplyr::select(HH_Europe, -country2, -phat0_obs, -ac_obs, -xb_noac, -selection)
HH_Europe <- filter(HH_Europe, ln_ely_q>0 & ln_ely_q<13)
HH_Europe <- dplyr::select(HH_Europe, colnames(HH_Africa)[-c(8,13)])
stargazer(HH_Europe, type="latex", out="desc_eu.tex", digits = 2)

load(list_country_files[11])
HH_NonEurope <- HH_NonEurope %>% mutate_if(is.factor,as.num)
HH_NonEurope <- dplyr::select(HH_NonEurope, -country2, -phat0_obs, -ac_obs, -xb_noac, -selection)
HH_NonEurope <- filter(HH_NonEurope, ln_ely_q>0 & ln_ely_q<13)
HH_NonEurope <- dplyr::select(HH_NonEurope, colnames(HH_Africa)[-c(8,13)])
stargazer(HH_NonEurope, type="latex", out="desc_noneu.tex", digits = 2)

load(list_country_files[12])
HH_Pakistan <- HH_Pakistan %>% mutate_if(is.factor,as.num)
HH_Pakistan <- dplyr::select(HH_Pakistan, -state, -phat0_obs, -ac_obs, -xb_noac, -selection, -hhid)
HH_Pakistan <- filter(HH_Pakistan, ln_ely_q>0 & ln_ely_q<13)
HH_Pakistan <- dplyr::select(HH_Pakistan, colnames(HH_Africa))
stargazer(HH_Pakistan, type="latex", out="desc_pak.tex", digits = 2)

load(list_country_files[13])
HH_USA <- HH_USA %>% mutate_if(is.factor,as.num)
HH_USA <- dplyr::select(HH_USA, -state, -macroarea, -phat0_obs, -ac_obs, -xb_noac, -selection, -hhid)
HH_USA <- filter(HH_USA, ln_ely_q>0 & ln_ely_q<13)
HH_USA <- dplyr::select(HH_USA, colnames(HH_Africa)[-13])
stargazer(HH_USA, type="latex", out="desc_usa.tex", digits = 2)



##########################################################################################################
# DO NOT USE!

descriptive <- function(variable, design){
  
  formula <- make.formula(variable)
  m <- svymean(formula, design, na.rm=TRUE)
  c(mean = coef(m), se = SE(m))
  
}

sapply(variables, descriptive, design = global_svy)

# Load function
devtools::source_gist("c4d1089a501d3567be9fb784b1c5a6ab")

# Define 
datasets <- list("Global" = global)
variable_names <- list(variables)
labels <- list(variables)
colnames <- c("Mean", "Median", "SD")

# We can define a descriptive function:
myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- coef(svymean(~ x, global_svy, na.rm=TRUE))
  sd <- SE(svymean(~ x, global_svy, na.rm=TRUE))
  return(c(m, sd))
}

myDescriptives = function(x) {
  x <- as.numeric(x)
  m <- weighted.mean(x, global$weight, na.rm = TRUE)
  md <- reldist::wtd.quantile (x, q=0.5, global$weight, na.rm = TRUE)
  sd <- sqrt(wtd.var(x, global$weight, na.rm = TRUE))
  return(c(m, md, sd))
}

createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Descriptive statistics")

stargazer(global)
#stargazer(data_c, type="latex", out="desc_global.tex", digits = 2)








