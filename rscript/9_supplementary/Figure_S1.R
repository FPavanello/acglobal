
##########################################

#               Figure S1

##########################################

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
user <- 'user'

if (user=='user') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'6-Projections/repo/household/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/supplementary/'


# Load global map
sf <- st_as_sf(rnaturalearthdata::countries110)

# Load pooled data
global <- readRDS(paste(house,'global.rds', sep='')) %>% dplyr::select(country)

# Summary statistics
all_merge_s <- dplyr::group_by(global, country) %>%  dplyr::summarise(n = n())

# Merge
sf$name[sf$name == "United States of America"] <- "United States"
sf <- merge(sf, all_merge_s, by.x="name", by.y="country", all.x=T) 
sf <- filter(sf, sovereignt != "Antarctica")

sf$iso_a3 <- ifelse(is.na(sf$n), NA, sf$iso_a3)

# Coordinates
sf <- st_transform(sf, "ESRI:54009")
gc() # clean

# Included dummy
sf$included <- ifelse(!is.na(sf$n), 1, 0)
sf$sovereignt <- ifelse(sf$included==1, sf$sovereignt, NA)
sf_only <- filter(sf, !is.na(n))

# Map
map <-  ggplot()+
  geom_sf(data=sf, aes(fill=as.factor(included)))+
  scale_fill_manual(values = c("white", "orange"))+
  theme_classic()+
  theme(legend.position = "none")+ 
  geom_sf_label_repel(data=sf_only, aes(label = paste0(sovereignt)), colour = "black") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

map

# Save
ggsave(paste(output, 'FigureS1.png', sep = ''), last_plot(), scale=2.5, width = 4.5, height = 4)

