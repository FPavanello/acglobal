
## This R-script:
##      1) create maps with the countries covered in the data set

rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Load packages
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(readxl)
library(countrycode)
#devtools::install_github("yutannihilation/ggsflabel")
library(ggsflabel)

# Set users
user <- 'fp'
#user <- 'gf'

if (user=='fp') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

if (user=='gf') {
  stub <- 'H:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

data <- paste(stub,'6-Projections/', sep='')
output <- paste(stub,'6-Projections/results/graphs/Paper1/', sep='')

# Load global map
sf <- st_as_sf(rnaturalearthdata::countries110)

# Load list of countries
countries <- read_xlsx(paste(data,'list_of_countries_final.xlsx', sep=''))
countries$iso3c <- countrycode(countries$Country, 'country.name', 'iso3c')

# Merge
sf <- merge(sf, countries, by.x="iso_a3", by.y="iso3c", all.x=T)
sf <- filter(sf, sovereignt != "Antarctica")
sf$included <- ifelse(!is.na(sf$Country), 1, 0)
sf$sovereignt <- ifelse(sf$included==1, sf$sovereignt, NA)
sf_only <- filter(sf, !is.na(Country))

# Map
map <-  ggplot()+
        geom_sf(data=sf, aes(fill=as.factor(included)))+
        scale_fill_manual(values = c("white", "orange"))+
        theme_classic()+
        theme(legend.position = "none")+ 
        geom_sf_label_repel(data=sf_only, aes(label = paste0(sovereignt, "\n", `Time Span`)), colour = "black") +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())

map

# Save
ggsave(paste(output, 'FigureS1_2.png', sep = ''), last_plot(), scale=2.5, width = 4.5, height = 4)

