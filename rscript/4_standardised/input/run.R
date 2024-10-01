
# Run country-specific models

.rs.restartR()
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Set users
user <- 'user'

if (user=='user') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'data/household/', sep='')
interm <- paste(stub,'results/regressions/for_graphs/subsamples/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/'
output <- paste(stub,'output/figures/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/figures/'
script <- 'C:/Users/Standard/Documents/Github/acglobal/rscript/4_standardised/input/standardised_country/'


# Run
source(paste(script,'dmcf_std_africa_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_std_argentina_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_std_brazil_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_std_china_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_std_germany_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_std_global_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_std_india_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_std_indonesia_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_std_italy_db_cont_urb.R', sep='')) # Done
source(paste(script,'dmcf_std_mexico_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_std_oecd_db_conturb.R', sep='')) # Done
source(paste(script,'dmcf_std_pakistan_db_conturb.R', sep='')) # Done 
source(paste(script,'dmcf_std_usa_ahs_db_conturb.R', sep='')) # Done

