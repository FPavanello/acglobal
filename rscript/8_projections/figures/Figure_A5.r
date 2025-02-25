
## This R-script:
##      1) Figure 2: Decomposition analysis

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

data <- paste(stub,'6-Projections/results/regressions/for_graphs/standardised', sep='')
output <- paste(stub,'6-Projections/results/graphs/', sep='')

list_p <- list.files(path=output, full.names = T, pattern = "decompose")
list_p <- list_p[grepl("Rdata", list_p)]
list_p <- list_p[grepl("SSP", list_p)]
list_p = list_p[!grepl("DEU", list_p)]

merger <- function(X){
  
  load(list_p[[X]])
  
  all$country <- sub("\\_decompose.*", "", gsub(output, "", list_p[X])) 
  all$ssp <- gsub(".Rdata", "", sub(".*\\_decompose_", "", gsub(output, "", list_p[X])))
  
  return(all)
  
}

pp <- lapply(1:length(list_p), merger)

pp <- as.data.frame(do.call(rbind, pp))

pp$country[pp$country=="DEU"] <- "OECD-EU"
pp$country[pp$country=="OECD_EU"] <- "OECD-EU"
pp$country[pp$country=="OECD_NONEU"] <- "OECD-nonEU"


pp$country[pp$country=="ARG"] <- "Argentina"
pp$country[pp$country=="BRA"] <- "Brazil"
pp$country[pp$country=="CHN"] <- "China"
pp$country[pp$country=="IDN"] <- "Indonesia"
pp$country[pp$country=="IND"] <- "India"
pp$country[pp$country=="ITA"] <- "Italy"
pp$country[pp$country=="MEX"] <- "Mexico"
pp$country[pp$country=="PAK"] <- "Pakistan"
pp$country[pp$country=="USA"] <- "United States"

pp <- group_by(pp, type, time, x, type_d, country, ssp) %>% dplyr::summarise(value=mean(value, na.rm=T))

###########

pp$type[pp$type== "Future AC change"] <- "Future AC ownership change"
pp$type[pp$type== "AC"] <- "AC ownership"
pp$type[pp$type== "Social drivers"] <- "Socio-demographic drivers"
pp$type[pp$type== "Future social drivers change"] <- "Future socio-demographic change"

pp$type_d <- factor(pp$type, ordered = TRUE,
                     levels = c("Socio-demographic drivers", "Urbanisation", "AC ownership", "CDDs/HDDs","Expenditure", "Future socio-demographic change", "Future urbanisation", "Future AC ownership change", "Future CDDs/HDDs change", "Future expenditure change"))


col_palette <- c("#fcba03", "#0f1f61", "#700c0c", "#3ed2f0", "#02ab4b", "#fcdf8d", "#7591ff", "#a86c6c", "#dbf9ff", "#90fcbf")


decomposition_plots <- ggplot(pp)+
  theme_classic()+
  geom_bar(aes(x=interaction(x, ssp), y=value, fill= (type_d), group= forcats::fct_rev(type_d)), position="stack", stat="identity", colour="black") + scale_fill_manual(name="Driver", values=col_palette) + 
  scale_x_discrete(labels=rep(c("SSP245", "SSP585"), length(unique(pp$x))))+
  geom_hline(yintercept = 1, size=2, linetype="dashed") + scale_y_continuous(labels = scales::label_percent()) + xlab("") + ylab("Relative contribution; 100% = current consumption") + ggtitle(paste0("Drivers of HH (per-capita) electricity consumption, current and growth until 2050"))+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) + guides(fill=guide_legend(ncol=1))+
  facet_wrap(vars(country), scales = "free")+
  labs(caption = "Contributions above the dashed line represent the projected growth until 2050")+
  theme(legend.position = "bottom", legend.direction = "horizontal")+guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
decomposition_plots

ggsave(paste0(output, "fig2.png"), last_plot(), scale=2.1, width = 4.5, height = 4)
