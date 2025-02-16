
##########################################

#               Figure A1

##########################################

# Free memory
.rs.restartR()
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
library(tidyverse)
library(reshape2)
library(spatstat)
#devtools::install_github("yutannihilation/ggsflabel")
library(ggsflabel)
library(patchwork)
library(gtools)

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

# Filter
all_merge <- filter(global, ely_q<15000)

# Quintiles
all_merge <- all_merge %>% group_by(country) %>% mutate(total_exp_usd_2011_d=gtools::quantcut(total_exp_usd_2011, q = 5, na.rm = TRUE)) %>% ungroup()
all_merge <- all_merge %>% group_by(country) %>% mutate(mean_CDD_db_d=gtools::quantcut(mean_CDD_db*100, q = 5, na.rm = TRUE)) %>% ungroup()

# Tile
all_merge_tile <- all_merge %>% group_by(country, total_exp_usd_2011_d, mean_CDD_db_d) %>% summarise(ac=mean(ac, na.rm=T), ely_q=mean(ely_q, na.rm=T))

# Plotting heat maps
tile_ac <- ggplot(na.exclude(all_merge_tile)) + geom_tile(aes(x = total_exp_usd_2011_d, y = mean_CDD_db_d, fill = ac*100)) + 
  scale_fill_stepsn(name="AC ownership (%)", colours=rev(c('#00429d', '#2e59a8', '#4771b2', '#5d8abd', '#73a2c6', '#8abccf', '#a5d5d8', '#c5eddf', '#ffffe0')), n.breaks=6) + 
  facet_wrap(vars(country), scales = "free", nrow=5) + xlab("Expenditure/hh/yr") + ylab("Historical CDDs / yr") + scale_x_discrete(labels=c("Q1","Q2","Q3","Q4","Q5")) + scale_y_discrete(labels=c("Q1","Q2","Q3","Q4","Q5")) +   theme(legend.position = "bottom", legend.direction = "horizontal", axis.title=element_text(size=14),legend.title=element_text(size=14)) +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1))
tile_ac

tile_ely <- ggplot(na.exclude(all_merge_tile)) + geom_tile(aes(x = total_exp_usd_2011_d, y = mean_CDD_db_d, fill = ely_q)) + 
  scale_fill_stepsn(name="Electricity consumption (kWh/hh/yr)", colours=rev(c('#005a74', '#3c6c7c', '#5e8084', '#7c938d', '#98a897', '#b4bda2', '#cfd2af', '#eae8c0', '#ffffe0')), n.breaks=6) + 
  facet_wrap(vars(country), scales = "free", nrow=5) + xlab("Expenditure/hh/yr") + ylab("Historical CDDs / yr") + scale_x_discrete(labels=c("Q1","Q2","Q3","Q4","Q5")) + scale_y_discrete(labels=c("Q1","Q2","Q3","Q4","Q5")) +   theme(legend.position = "bottom", legend.direction = "horizontal", axis.title=element_text(size=14),legend.title=element_text(size=14)) +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1))
tile_ely

saver2 <- tile_ac + tile_ely + plot_annotation(tag_levels = "A") + plot_layout(ncol=2)
saver2

# Figure A1
ggsave(paste(output, 'FigureA1.png', sep=''), saver2, width = 10, height = 7, scale=1.4)

