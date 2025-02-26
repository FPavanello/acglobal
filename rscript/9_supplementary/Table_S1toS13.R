
##########################################

#             Table S1 to S13

##########################################

# Free memory
.rs.restartR()
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
library(data.table)
library(xtable)
library(reldist)
library(Hmisc)

# Set users
user <- 'user'

if (user=='user') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

# Set directory
house <- paste(stub,'6-Projections/repo/household/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/supplementary/'


# Load
global <- readRDS(paste(house,'global.rds', sep=''))

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
global <- global[complete.cases(global$urban_sh), ]
global <- global[complete.cases(global$n_members), ]
global <- global[complete.cases(global$age_head), ]
global <- global[complete.cases(global$edu_head_2), ]

global <- dplyr::filter(global, ln_ely_q > 0)
global <- dplyr::filter(global, weight > 0)

# Education
global <- dplyr::mutate(global, 
                        noedu = ifelse(edu_head_2 == 0, 1, 0),
                        prim = ifelse(edu_head_2 == 1, 1, 0),
                        sec = ifelse(edu_head_2 == 2, 1, 0),
                        post = ifelse(edu_head_2 == 3, 1, 0))

# Load function
devtools::source_gist("c4d1089a501d3567be9fb784b1c5a6ab")

# Define 
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


## Africa
# Filter
africa <- dplyr::filter(global, country == "Burkina Faso" | country == "Ghana" | country == "Kenya" | 
                                country == "Malawi" | country == "Niger" | country == "Nigeria")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "tv", "pc", "wshm")

# Function to get descriptives
datasets <- list("Africa" = africa)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS1",
                       file = paste(output,'TableS1.tex', sep=''))

## Argentina
# Filter
argentina <- dplyr::filter(global, country == "Argentina")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "tv", "wshm")

# Function to get descriptives
datasets <- list("Argentina" = argentina)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS2",
                       file = paste(output,'TableS2.tex', sep=''))

## Brazil
# Filter
brazil <- dplyr::filter(global, country == "Brazil")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "tv", "pc", "wshm")

# Function to get descriptives
datasets <- list("Brazil" = brazil)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS3",
                       file = paste(output,'TableS3.tex', sep=''))

## China
# Filter
china <- dplyr::filter(global, country == "China")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "tv", "pc", "wshm")

# Function to get descriptives
datasets <- list("China" = china)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS4",
                       file = paste(output,'TableS4.tex', sep=''))

## Germany
# Filter
germany <- dplyr::filter(global, country == "Germany")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head")

# Function to get descriptives
datasets <- list("Germany" = germany)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS5",
                       file = paste(output,'TableS5.tex', sep=''))

## Indonesia
# Filter
indonesia <- dplyr::filter(global, country == "Indonesia")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "tv", "pc")

# Function to get descriptives
datasets <- list("Indonesia" = indonesia)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS6",
                       file = paste(output,'TableS6.tex', sep=''))

## India
# Filter
india <- dplyr::filter(global, country == "India")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "tv", "wshm")

# Function to get descriptives
datasets <- list("India" = india)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS7",
                       file = paste(output,'TableS7.tex', sep=''))

## Italy
# Filter
italy <- dplyr::filter(global, country == "Italy")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "tv", "pc", "wshm")

# Function to get descriptives
datasets <- list("Italy" = italy)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS8",
                       file = paste(output,'TableS8.tex', sep=''))

## Mexico
# Filter
mexico <- dplyr::filter(global, country == "Mexico")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "tv", "pc", "wshm")

# Function to get descriptives
datasets <- list("Mexico" = mexico)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS9",
                       file = paste(output,'TableS9.tex', sep=''))

## EPIC - EU
# Filter
eu <- dplyr::filter(global, country == "France" | country == "Netherlands" | country == "Switzerland" |
                      country == "Sweden" | country == "Spain")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "tv", "pc", "wshm")

# Function to get descriptives
datasets <- list("EU" = eu)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS10",
                       file = paste(output,'TableS10.tex', sep=''))

## EPIC - NONEU
# Filter
noneu <- dplyr::filter(global, country == "Australia" | country == "Canada" | country == "Japan")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "tv", "pc", "wshm")

# Function to get descriptives
datasets <- list("NONEU" = noneu)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS11",
                       file = paste(output,'TableS11.tex', sep=''))

## Pakistan
# Filter
pak <- dplyr::filter(global, country == "Pakistan")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "tv", "wshm")

# Function to get descriptives
datasets <- list("Pakistan" = pak)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS12",
                       file = paste(output,'TableS12.tex', sep=''))


## United States
# Filter
usa <- dplyr::filter(global, country == "United States")

# Vector of variables
variables <- c("ely_q", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp_usd_2011", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "noedu", "prim", "sec", "post", "age_head", "sex_head", "ref", "wshm")

# Function to get descriptives
datasets <- list("United States" = usa)
variable_names <- list(variables)
labels <- list(variables)

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
                       label = "si:tableS13",
                       file = paste(output,'TableS13.tex', sep=''))

# Clean
rm(africa, argentina, brazil, china, datasets, eu, germany, global, india, indonesia, italy, labels, mexico, noneu, pak, usa, variable_names)
gc()
