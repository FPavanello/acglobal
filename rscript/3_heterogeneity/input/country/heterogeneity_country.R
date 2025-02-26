
# Run country-specific models

.rs.restartR()
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Set users
user <- 'user'

if (user=='user') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'6-Projections/repo/household/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/'
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/figures/'
script <- 'C:/Users/Standard/Documents/Github/acglobal/rscript/3_heterogeneity/input/country/country-specific/'


# Run
source(paste(script,'dmcf_africa_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_argentina_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_brazil_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_china_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_germany_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_india_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_indonesia_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_italy_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_mexico_db_conturb.R', sep=''))# Done
source(paste(script,'dmcf_oecd_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_pakistan_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_united_states_ahs_db_conturb.R', sep='')) # Done

