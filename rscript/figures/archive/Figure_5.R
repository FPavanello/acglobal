
## This R-script:
##      1) plots AC elasticities in each country sample by income tercile - coefficient plot (Figure 5)

# Free memory
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
library(tidyverse)
library(data.table)
library(patchwork)
library(ggsci)

# Set users
user <- 'fp'
#user <- 'gf'

if (user=='fp') {
  stub <- 'G:/Il mio Drive/'
}

if (user=='gf') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

results <- paste(stub,'6-Projections/results/regressions/', sep='')
graphs <- paste(stub,'6-Projections/results/graphs/', sep='')


########################################################

#    PLOT OF MFXs BY COUNTRY AND BY INCOME TERCILES    #

########################################################

# List of all files
list_p <- list.files(path= paste(results, '/for_graphs/subsamples/', sep = ''), full.names = T, pattern = "terciles")

# Function to merge all files
merger <- function(X){
  
  load(list_p[X])
  
  ac_effect <-  ac_eff
  ac_effect$country <- gsub("_terciles.RData", "", basename(list_p[X]))
  
  return(ac_effect)
  
}

# Load AC effects' terciles
ac_effect_trc <- lapply(1:length(list_p), merger)
ac_effect_trc <- as.data.frame(do.call(rbind, ac_effect_trc))

# Adjustments
#ac_effect_trc$country[ac_effect_trc$country=="deu"] <- "oecd_eu"
#ac_effect_trc <- group_by(ac_effect_trc, country, ter_inc) %>% summarise_all(., "mean")

ac_effect_trc$ter_inc[ac_effect_trc$ter_inc=="Low"] <- "Poor" 
ac_effect_trc$ter_inc[ac_effect_trc$ter_inc=="Medium"] <-  "Med"
ac_effect_trc$ter_inc[ac_effect_trc$ter_inc=="High"] <- "Rich"

ac_effect_trc$country[ac_effect_trc$country=="afr"] <- "Africa" 
ac_effect_trc$country[ac_effect_trc$country=="arg"] <-  "ARG"
ac_effect_trc$country[ac_effect_trc$country=="bra"] <- "BRA"
ac_effect_trc$country[ac_effect_trc$country=="chn"] <- "CHN"
ac_effect_trc$country[ac_effect_trc$country=="deu"] <- "DEU"
ac_effect_trc$country[ac_effect_trc$country=="oecd_eu"] <- "EPIC - EU"
ac_effect_trc$country[ac_effect_trc$country=="oecd_non_eu"] <- "EPIC - NonEU"
ac_effect_trc$country[ac_effect_trc$country=="idn"] <- "IDN"
ac_effect_trc$country[ac_effect_trc$country=="ind"] <- "IND"
ac_effect_trc$country[ac_effect_trc$country=="ita"] <- "ITA"
ac_effect_trc$country[ac_effect_trc$country=="mex"] <- "MEX"
ac_effect_trc$country[ac_effect_trc$country=="pak"] <- "PAK"
ac_effect_trc$country[ac_effect_trc$country=="usa"] <- "USA"
ac_effect_trc <- filter(ac_effect_trc, country!="global")

# Significance
ac_effect_trc$sign <- ifelse(ac_effect_trc$p <0.05, 1, 0)

# Coefficients Plot
ac_trc <- ac_effect_trc %>% ggplot(aes(x= factor(ter_inc, levels = c("Poor", "Med", "Rich")), y=AME*100)) + 
          geom_point() + 
          geom_errorbar(aes(ymin=lower*100, ymax=upper*100, colour=as.factor(sign)), alpha=0.6, size=0.6, width=0.4) +
          geom_hline(yintercept = 0, colour="black", linetype="dashed") +
          scale_colour_manual(name="95% significance", values = c("#bf0000", "#4DBBD5B2"), labels=c("No", "Yes")) + 
          ylab("Additional Electricity Demand due to AC (Marginal Effect, %)")+
          xlab("Expenditure Tercile")+
          theme_classic() +
          theme(panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_line(color="lightgray"),
                panel.grid.minor=element_blank(),
                plot.title = element_text(size = 12, family = "Tahoma", hjust = 0.5),
                text = element_text(size = 12, family = "Tahoma"), 
#                axis.title.x = element_blank(),
                axis.ticks=element_blank()) +
          facet_wrap(vars(country), scales="free")

ac_trc

# Save
ggsave(paste(graphs, 'figure_5.png', sep = ''), last_plot(), scale=2.5, width = 3.5, height = 2)

#png(paste0(graphs, 'figure_5.png',sep=''), height=2, width=3.5, scale=2.5, units='in', res=1200)
#print(ac_trc)
#dev.off()
