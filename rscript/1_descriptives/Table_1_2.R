
## This R-script:
##      1) descriptive statistics whole data set (Table 1-2)


rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

## 1) Load libraries and data ##
library(data.table)
library(xtable)
library(reldist)
library(Hmisc)

# Set users
user <- 'user'

if (user=='user') {
  stub <- "add your repository"
}

house <- paste(stub,'data/household/', sep='') 
output <- paste(stub,'output/tables/', sep='') 


#############################################################################

# TABLE 1-2: Number of households for each country + Descriptive Statistics #

#############################################################################

# Load global data
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


# Table 1: Number of households per country
summary(global$country) 

# Transform log data
global <- dplyr::mutate(global, ely = exp(ln_ely_q), total_exp = exp(ln_total_exp_usd_2011))

# Survey
global_svy <- svydesign(data = global, ids = ~1, weights = ~ weight)


## Summary global data set
# Vector of variables
variables <- c("ely", "ac", "mean_CDD18_db", "curr_CDD18_db", "curr_HDD18_db", "total_exp", "ely_p_usd_2011", 
               "urban_sh", "ownership_d", "n_members", "edu_head_2", "age_head", "sex_head")

# Load function
devtools::source_gist("c4d1089a501d3567be9fb784b1c5a6ab") # Github repo containing the function to compute weighted descriptives in Latex format

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

# Table 2: Descriptive Statistics
createDescriptiveTable(datasets,
                       summary_function = myDescriptives,
                       column_names = colnames,
                       variable_names = variables,
                       variable_labels = labels,
                       arraystretch = 1.3,
                       title = "Weighted Descriptive statistics",
                       label = "tab:desc_all",
                       file = paste(output,'desc_global.tex', sep=''))
