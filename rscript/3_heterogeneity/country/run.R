
## This R-script:
##      1) run all country-specific scripts
##      2) exploits CDD-dry bulb 24 deg and HDD-dry bulb 15 deg
##      3) conducts logit regressions for Africa
##      4) run intensive margin regressions: electricity expenditure on climate + covariates
##         using Dubin and McFadden (1984) approach

.rs.restartR()
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Set directory
user <- 'fp'
#user <- 'gf'

if (user=='fp') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

house <- paste(stub,'6-Projections/data/household/', sep='')
output <- paste(stub,'6-Projections/results/regressions/', sep='')
script <- paste(stub,'6-Projections/rscripts/dmcf/regressions/country/with_continuous_urbanisation/', sep='')

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

