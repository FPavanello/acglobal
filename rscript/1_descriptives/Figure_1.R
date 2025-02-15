
##########################################

#               Figure 1

##########################################

# Free memory
.rs.restartR()
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
library(tidyverse)
library(rnaturalearthdata)
library(sf)
library(reshape2)
library(spatstat)
#devtools::install_github("yutannihilation/ggsflabel")
library(ggsflabel)
library(patchwork)

# Set users
user <- 'user'

if (user=='user') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'data/household/', sep='')
interm <- paste(stub,'results/regressions/for_graphs/subsamples/', sep='')
interm <- "C:/Users/Standard/Documents/Github/acglobal/interm/"
output <- paste(stub,'output/figures/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/figures/'


# Load pooled data
global <- readRDS(paste(house,'global.rds', sep='')) %>% dplyr::select(country, ac, ely_q, contains("CDD"), total_exp_usd_2011, weight)
global <- dplyr::filter(global, weight > 0) # filter

# Summary statistics
all_merge_s <- dplyr::group_by(global, country) %>%  dplyr::summarise(ac=weighted.mean(ac, weight, na.rm = T), 
                                                                      ely_q=weighted.median(ely_q, weight, na.rm = T), 
                                                                      mean_CDD_db=weighted.mean(mean_CDD18_db, weight, na.rm = T), 
                                                                      total_exp_usd_2011=weighted.median(total_exp_usd_2011, weight, 
                                                                                                         na.rm = T))

# Prepare data for the figures
all_merge_s_m <- reshape2::melt(all_merge_s, 1)

# Check
e <- ggplot(all_merge_s_m)+
  geom_col(aes(x=country, y=value))+
  facet_wrap(vars(variable), scales = "free", ncol=4)
e # check

# Load world map
sf <- st_as_sf(rnaturalearthdata::countries110)

# Merge
sf$name[sf$name == "United States of America"] <- "United States"
sf <- merge(sf, all_merge_s, by.x="name", by.y="country", all.x=T) 
sf <- filter(sf, sovereignt != "Antarctica")

sf$iso_a3 <- ifelse(is.na(sf$ac), NA, sf$iso_a3)

# Coordinates
sf <- st_transform(sf, "ESRI:54009")
gc() # clean

# Panel A
a <- ggplot()+
  geom_sf(data=sf, aes(fill=ac*100))+
  scale_fill_stepsn(name="AC ownership (%)", colours=rev(c('#00429d', '#2e59a8', '#4771b2', '#5d8abd', '#73a2c6', '#8abccf', '#a5d5d8', '#c5eddf', '#ffffe0')), n.breaks=6)+
  theme_void()+
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.title=element_text(size=13))+ 
  geom_sf_label_repel(data=sf, aes(label = iso_a3), colour = "black")+
  guides(fill = guide_colourbar(barwidth = 12, barheight = 1))
a

# Panel B
b <- ggplot()+
  geom_sf(data=sf, aes(fill=ely_q))+
  scale_fill_stepsn(name="Median electricity consumption (kWh/hh/yr)", colours=rev(c('#005a74', '#3c6c7c', '#5e8084', '#7c938d', '#98a897', '#b4bda2', '#cfd2af', '#eae8c0', '#ffffe0')), n.breaks=6)+
  theme_void()+
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.title=element_text(size=13))+ 
  geom_sf_label_repel(data=sf, aes(label = iso_a3), colour = "black")+
  guides(fill = guide_colourbar(barwidth = 12, barheight = 1))
b

# Panel C
c <- ggplot()+
  geom_sf(data=sf, aes(fill=mean_CDD_db*100))+
  scale_fill_stepsn(name="Median historical CDDs / yr", colours=rev(c('#00004d', '#2c1d5f', '#4b3971', '#685783', '#847795', '#a197a8', '#bfb9ba', '#dedccd', '#ffffe0')), n.breaks=6)+
  theme_void()+
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.title=element_text(size=13))+ 
  geom_sf_label_repel(data=sf, aes(label = iso_a3), colour = "black")+
  guides(fill = guide_colourbar(barwidth = 12, barheight = 1))
c

# Panel D
d <- ggplot()+
  geom_sf(data=sf, aes(fill=total_exp_usd_2011))+
  scale_fill_stepsn(name="Median expenditure/hh/yr (2011 PPP USD)", colours=rev(c('#784441', '#8c594a', '#a06e54', '#b48460', '#c79b6f', '#dab382', '#eacb99', '#f8e4b7', '#ffffe0')), n.breaks=6)+
  theme_void()+
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.title=element_text(size=13))+ 
  geom_sf_label_repel(data=sf, aes(label = iso_a3), colour = "black")+
  guides(fill = guide_colourbar(barwidth = 12, barheight = 1))
d

saver1 <- a + b + c + d + plot_annotation(tag_levels = "A")
saver1

# Figure 1
ggsave(paste(output,'Figure1.png', sep=''), saver1, width = 12, height = 8, scale=1.2)

