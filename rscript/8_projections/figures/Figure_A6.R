
## This R-script:
##      1) Figure 3: expenditure distribution analysis

# Free memory
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
library(tidyverse)
library(patchwork)
#remotes::install_github('rpkgs/gg.layers')
#library(gg.layers)
library(ggsci)

# Set users
user <- 'fp'
user <- 'gf'

if (user=='fp') {
  stub <- 'G:/Il mio Drive/'
}

if (user=='gf') {
  stub <- 'F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

setwd(stub)

############################
############################

# A' la Burke analysis

lll <- list.files(path="F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/energy_poverty/", full.names = T, pattern = ".Rdata")[-6]

###

ff <- function(X){
 
  print(X)
  
load(X)

colnames(output_impact_ac) <- paste0("cons_", colnames(output_impact_ac))

rr <- bind_cols(data_c_sp, output_impact_ac)

rr <- dplyr::select(rr, ely_p_usd_2011, weight, total_exp_usd_2011, starts_with("cons_"))

rr$country <- basename(X)

return(rr)
 
}

library(pbapply)

llll <- pblapply(lll, ff)

llll <- bind_rows(llll)

rr <- reshape2::melt(llll, c(1:3, 12))

rr$share <- (rr$ely_p_usd_2011 * rr$value) / rr$total_exp_usd_2011

library(stringr)

rr$year <- substr(rr$variable, 11, 14)
rr$year <- as.numeric(rr$year)

rr$scenario <- substr(rr$variable, 6, 9)

rr$country <- gsub("\\.Rdata", "", rr$country)

rr <- na.omit(rr)

r <- ggplot(rr %>% filter(year %in% c(2020, 2050) & country!="DEU"))+
  theme_classic()+
  gg.layers::geom_boxplot2(aes(y=share*100, fill=as.factor(year), x=scenario, weight=weight), width.errorbar = 0.1)+
  ylab("% of household electricity consumption owing to AC")+
  xlab("Scenario")+
  scale_fill_discrete(name="Year")+
  facet_wrap(vars(country), scales = "free_y")

ggsave("6-Projections/results/graphs/ceteris_paribus_shares.png", r, height = 3.75, width = 3.75, scale=1.5, bg="white")

