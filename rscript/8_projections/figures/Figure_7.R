
##########################################

#               Figure 7

##########################################

rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
library(tidyverse)
library(patchwork)
#remotes::install_github('rpkgs/gg.layers')
#library(gg.layers)
library(ggsci)

# Set users
user <- 'user'

if (user=='user') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'data/household/', sep='')
output <- paste(stub,'output/figures/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/figures/'
interm <- "C:/Users/Standard/Documents/Github/acglobal/interm/"
interm <- paste(stub,'6-Projections/results/graphs/', sep='')


# Load data
list_p <- list.files(path=output, full.names = T, pattern = "distribution")
list_p <- list_p[grepl("Rdata", list_p)][]

merger <- function(X){
  
  load(list_p[[X]])
  
  output_impact_ac2$country <- gsub("_distribution.Rdata", "", gsub(output, "", list_p[X]))
  
  return(output_impact_ac2)
  
}

pp <- lapply(1:length(list_p), merger)

pp <- as.data.frame(do.call(rbind, pp))

# Labels
pp$country[pp$country=="OECD-EU"] <- "OECD - EU"
pp$country[pp$country=="OECD-NONEU"] <- "OECD - non EU"

pp$country[pp$country=="ARG"] <- "Argentina"
pp$country[pp$country=="BRA"] <- "Brazil"
pp$country[pp$country=="CHN"] <- "China"
pp$country[pp$country=="IDN"] <- "Indonesia"
pp$country[pp$country=="IND"] <- "India"
pp$country[pp$country=="ITA"] <- "Italy"
pp$country[pp$country=="MEX"] <- "Mexico"
pp$country[pp$country=="PAK"] <- "Pakistan"
pp$country[pp$country=="USA"] <- "United States"

# Arrange
pp <- filter(pp, year==2020 | year==2050)
pp <- filter(pp, country!="DEU")
pp <- filter(pp, value>0 & value<20000)
pp = filter(pp, ely_p < 1)

pp$exp = pp$value*pp$ely_p

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.1, .9), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

pp = pp %>%
  group_by(ssp, year, country) %>%
  mutate(exp = remove_outliers(exp))

pp <- pp %>% group_by(ssp, year, country) %>% dplyr::mutate(mean=mean(exp[exp>15], na.rm=T)) %>% ungroup()

pp <- pp %>% filter(!(ssp=="SSP585" & year==2020))

pp$ssp <- ifelse(pp$year==2020 & pp$ssp=="SSP245", "2020", pp$ssp)

# Distribution
distribution_plots <- ggplot(pp %>% filter(value>15), aes(x = exp, y = after_stat(count), group=ssp, colour=as.factor(ssp)))+ theme_classic() +   geom_density(adjust=.7)+ facet_wrap(vars(country), ncol=3, scales="free")+
  geom_vline(aes(xintercept = mean, group=ssp, colour=as.factor(ssp)), linetype="dashed", size=.3)+
  ylab("Count of households")+
  xlab("Expenditure in AC electricity consumption, distribution among households with AC, 2011 PPP USD/HH/yr.")+
  scale_color_manual(name="Year", values=c("grey", "orange", "darkred"))+
  labs(caption = "Solid line: density curve; dashed line: mean value")+
  theme(legend.position = "bottom", legend.direction = "horizontal", strip.background = element_blank())+guides(colour=guide_legend(ncol=4,byrow=TRUE))

distribution_plots

# Save
ggsave(paste0(output, "Figure7.png"), distribution_plots, scale=1.2, width = 7, height = 10*(2/3))

