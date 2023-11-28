## 1) Load libraries and data ##
library(sandwich)
library(lmtest)
library(foreign)
library(ResourceSelection)
library(optmatch)
library(tidyverse)
library(haven)
library(psych)
library(raster)
library(rnaturalearthdata)
library(sf)
library(gdata)
library(exactextractr)
library(nngeo)
library(caret)
library(MatchIt)
library(ggsci)
library(gdata)
library(jtools)
library(glm2)
library(reshape2)
library(cobalt)

ARG <- readRDS("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Argentina/ENGH/argentina_engh.rds") %>% dplyr::select(ac, ely_q, contains("CDD"), total_exp_usd_2011, weight)

CHN <- readRDS("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/China/china.rds") %>% dplyr::select(ac, ely_q, contains("CDD"), contains("_total_exp_usd_2011"), weight)
CHN$total_exp_usd_2011 <- exp(CHN$ln_total_exp_usd_2011)
CHN$ln_total_exp_usd_2011 <- NULL

DEU <- readRDS("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Germany/germany.Rds")  %>% dplyr::select(ac, ely_q, contains("CDD"), contains("_total_exp_usd_2011"), weight)
DEU$total_exp_usd_2011 <- exp(DEU$ln_total_exp_usd_2011)
DEU$ln_total_exp_usd_2011 <- NULL

GHA <- readRDS("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Ghana/ghana.rds")  %>% dplyr::select(ac, ely_q, contains("CDD"), contains("_total_exp_usd_2011"), weight)
GHA$total_exp_usd_2011 <- exp(GHA$ln_total_exp_usd_2011)
GHA$ln_total_exp_usd_2011 <- NULL

ITA <- readRDS("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Italy/italy_hbs.rds") %>% dplyr::select(ac, ely_q, contains("CDD"), total_exp_usd_2011)

#GRE <- readRDS("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Italy/")

#BGA <- readRDS("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Italy/")

KEN <- readRDS("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Kenya/IHBS/kenya_ihbs.rds") %>% dplyr::select(ac, ely_q, contains("CDD"), contains("_total_exp_usd_2011"), weight)
KEN$total_exp_usd_2011 <- exp(KEN$ln_total_exp_usd_2011)
KEN$ln_total_exp_usd_2011 <- NULL

NGA <- read_rds("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Nigeria/NGA/nigeria_ghs.rds") %>% dplyr::select(ac, ely_q, contains("CDD"), contains("_total_exp_usd_2011"), weight)
NGA$total_exp_usd_2011 <- exp(NGA$ln_total_exp_usd_2011)
NGA$ln_total_exp_usd_2011 <- NULL


PAK <- read_rds("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Pakistan/LSM-IHS/pakistan_lsmihs.rds") %>% dplyr::select(ac, ely_q, contains("CDD"), contains("_total_exp_usd_2011"), weight)
PAK$total_exp_usd_2011 <- exp(PAK$ln_total_exp_usd_2011)
PAK$ln_total_exp_usd_2011 <- NULL

TZA <- readRDS("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Tanzania/HBS/tanzania_hbs.rds") %>% dplyr::select(ac, ely_q, contains("CDD"), contains("_total_exp_usd_2011"), weight)
TZA$total_exp_usd_2011 <- exp(TZA$ln_total_exp_usd_2011)
TZA$ln_total_exp_usd_2011 <- NULL

USA <- readRDS("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/United States/us_ahs.Rds") %>% dplyr::select(ac, ely_q, contains("CDD"), contains("_total_exp_usd_2011"), weight)
USA$total_exp_usd_2011 <- exp(USA$ln_total_exp_usd_2011)
USA$ln_total_exp_usd_2011 <- NULL

#

BRA <- read_rds("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Fourcountries/bra_pof.rds")
IDN <- read_rds("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Fourcountries/idn_susenas.rds")
IND <- read_rds("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Fourcountries/ind_nss.rds")
MEX <- read_rds("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/Fourcountries/mex_enigh.rds")

IND <- IND %>% filter(country=="India")  %>% dplyr::select(ac, ely_q, mean_CDD_1970_2011_db , total_exp_usd_2011, weight)
IDN <-  IDN %>% filter(country=="Indonesia")  %>% dplyr::select(ac, ely_q, mean_CDD_1970_2016_db , total_exp_usd_2011, weight) 
MEX <-  MEX %>% filter(country=="Mexico")  %>% dplyr::select(ac, ely_q, mean_CDD_1970_2016_db , total_exp_usd_2011, weight) 
BRA <-  BRA %>% filter(country=="Brazil")  %>% dplyr::select(ac, ely_q, mean_CDD_1970_2016 , total_exp_usd_2011, weight) 

DTA <- read_dta("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/OECD/EPIC/Household_OECD_mod.dta") %>% dplyr::select(ac, ely_q, CDD_mean_1986_2011,income, country) %>% rename(total_exp_usd_2011=income)

SWE <- DTA %>% filter(country=="Sweden")
ESP <- DTA %>% filter(country=="Spain")
NLD <- DTA %>% filter(country=="Netherlands")
FRA <- DTA %>% filter(country=="France")
rm(DTA)

DTA <- read_dta("G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/data/household/OECD/EPIC/Household_OECD_mod.dta") %>% dplyr::select(ac, ely_q, CDD_mean_1986_2011,income, country) %>% rename(total_exp_usd_2011=income)

CAN <- DTA %>% filter(country=="Canada")
AUS <- DTA %>% filter(country=="Australia")
JPN <- DTA %>% filter(country=="Japan")
rm(DTA)

#

all <- mget(ls(), envir = globalenv())
rm(list=setdiff(ls(), "all"))

nomi <- names(all)

all_bk <- all

#

for (i in 1:length(all)){

  all[[i]]$ac<- as.numeric(as.character( all[[i]]$ac))
  all[[i]]$weight <- ifelse(is.null(all[[i]]$weight), 1, all[[i]]$weight)
  all[[i]]$weight<- as.numeric(as.character(all[[i]]$weight))
  
  if (is.null(all[[i]]$mean_CDD_db)){
    
    all[[i]]$mean_CDD_db <- pull(all[[i]] %>% dplyr::select(contains("CDD")))
    
  }
  

  all[[i]] <- all[[i]] %>% dplyr::select(ac, ely_q, mean_CDD_db, total_exp_usd_2011, weight)
  all[[i]] <- all[[i]] %>% dplyr::select(sort(names(.)))
}

library(data.table)

all_merge <- rbindlist(all, idcol=nomi, fill=T)

library(spatstat)

all_merge_s <- dplyr::group_by(all_merge, ARG) %>%  dplyr::summarise(ac=weighted.mean(ac, weight, na.rm = T), ely_q=weighted.median(ely_q, weight, na.rm = T), mean_CDD_db=weighted.mean(mean_CDD_db, weight, na.rm = T), total_exp_usd_2011=weighted.median(total_exp_usd_2011, weight, na.rm = T))

#


all_merge_s_m <- reshape2::melt(all_merge_s, 1)

e <- ggplot(all_merge_s_m)+
  geom_col(aes(x=ARG, y=value))+
  facet_wrap(vars(variable), scales = "free", ncol=4)

library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ggsflabel)

sf <- st_as_sf(rnaturalearthdata::countries110)

sf <- merge(sf, all_merge_s, by.x="iso_a3", by.y="ARG", all.x=T)
sf <- filter(sf, sovereignt != "Antarctica")

sf$iso_a3 <- ifelse(is.na(sf$ac), NA, sf$iso_a3)

#

sf <- st_transform(sf, "ESRI:54009")


a <- ggplot()+
  geom_sf(data=sf, aes(fill=ac))+
  scale_fill_stepsn(name="AC ownership (%)", colours=rev(c('#00429d', '#2e59a8', '#4771b2', '#5d8abd', '#73a2c6', '#8abccf', '#a5d5d8', '#c5eddf', '#ffffe0')), n.breaks=6)+
  theme_void()+
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.title=element_text(size=13))+ 
  geom_sf_label_repel(data=sf, aes(label = iso_a3), colour = "black")+
  guides(fill = guide_colourbar(barwidth = 12, barheight = 1))

b <- ggplot()+
  geom_sf(data=sf, aes(fill=ely_q))+
  scale_fill_stepsn(name="Median electricity consumption (kWh/hh/yr)", colours=rev(c('#005a74', '#3c6c7c', '#5e8084', '#7c938d', '#98a897', '#b4bda2', '#cfd2af', '#eae8c0', '#ffffe0')), n.breaks=6)+
  theme_void()+
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.title=element_text(size=13))+ 
  geom_sf_label_repel(data=sf, aes(label = iso_a3), colour = "black")+
  guides(fill = guide_colourbar(barwidth = 12, barheight = 1))

c <- ggplot()+
  geom_sf(data=sf, aes(fill=mean_CDD_db))+
  scale_fill_stepsn(name="Median historical CDDs / yr", colours=rev(c('#00004d', '#2c1d5f', '#4b3971', '#685783', '#847795', '#a197a8', '#bfb9ba', '#dedccd', '#ffffe0')), n.breaks=6)+
  theme_void()+
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.title=element_text(size=13))+ 
  geom_sf_label_repel(data=sf, aes(label = iso_a3), colour = "black")+
  guides(fill = guide_colourbar(barwidth = 12, barheight = 1))

d <- ggplot()+
  geom_sf(data=sf, aes(fill=total_exp_usd_2011))+
  scale_fill_stepsn(name="Median expenditure/hh/yr (2011 PPP USD)", colours=rev(c('#784441', '#8c594a', '#a06e54', '#b48460', '#c79b6f', '#dab382', '#eacb99', '#f8e4b7', '#ffffe0')), n.breaks=6)+
  theme_void()+
  theme(legend.position = "bottom", legend.direction = "horizontal",legend.title=element_text(size=13))+ 
  geom_sf_label_repel(data=sf, aes(label = iso_a3), colour = "black")+
  guides(fill = guide_colourbar(barwidth = 12, barheight = 1))

library(patchwork)

saver1 <- a + b + c + d + plot_annotation(tag_levels = "A")

ggsave("G:/.shortcut-targets-by-id/13znqeVDfPULc4J_lQLbyW_Kmfa03o63F/3-Research/GIACOMO/inequality paper/fig1.png", saver1, width = 12, height = 8, scale=1.2)

#

all_merge <- filter(all_merge, ely_q<15000)

all_merge <- all_merge %>% group_by(ARG) %>% mutate(total_exp_usd_2011_d=gtools::quantcut(total_exp_usd_2011, q = 5, na.rm = TRUE)) %>% ungroup()

all_merge <- all_merge %>% group_by(ARG) %>% mutate(mean_CDD_db_d=gtools::quantcut(mean_CDD_db, q = 5, na.rm = TRUE)) %>% ungroup()

all_merge_tile <- all_merge %>% group_by(ARG, total_exp_usd_2011_d, mean_CDD_db_d) %>% summarise(ac=mean(as.numeric(as.character(ac)), na.rm=T), ely_q=mean(as.numeric(as.character(ely_q)), na.rm=T))

tile_ac <- ggplot(na.exclude(all_merge_tile)) + geom_tile(aes(x = total_exp_usd_2011_d, y = mean_CDD_db_d, fill = ac)) + scale_fill_stepsn(name="AC ownership (%)", colours=rev(c('#00429d', '#2e59a8', '#4771b2', '#5d8abd', '#73a2c6', '#8abccf', '#a5d5d8', '#c5eddf', '#ffffe0')), n.breaks=6) + facet_wrap(vars(ARG), scales = "free", nrow=7) + xlab("Expenditure/hh/yr") + ylab("Historical CDDs / yr") + scale_x_discrete(labels=c("Q1","Q2","Q3","Q4","Q5")) + scale_y_discrete(labels=c("Q1","Q2","Q3","Q4","Q5")) +   theme(legend.position = "bottom", legend.direction = "horizontal", axis.title=element_text(size=14),legend.title=element_text(size=14)) +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1))

tile_ely <- ggplot(na.exclude(all_merge_tile)) + geom_tile(aes(x = total_exp_usd_2011_d, y = mean_CDD_db_d, fill = ely_q)) + scale_fill_stepsn(name="Electricity consumption (kWh/hh/yr)", colours=rev(c('#005a74', '#3c6c7c', '#5e8084', '#7c938d', '#98a897', '#b4bda2', '#cfd2af', '#eae8c0', '#ffffe0')), n.breaks=6) + facet_wrap(vars(ARG), scales = "free", nrow=7) + xlab("Expenditure/hh/yr") + ylab("Historical CDDs / yr") + scale_x_discrete(labels=c("Q1","Q2","Q3","Q4","Q5")) + scale_y_discrete(labels=c("Q1","Q2","Q3","Q4","Q5")) +   theme(legend.position = "bottom", legend.direction = "horizontal", axis.title=element_text(size=14),legend.title=element_text(size=14)) +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1))

saver2 <- tile_ac + tile_ely + plot_annotation(tag_levels = "A") + plot_layout(ncol=2)
ggsave("G:/.shortcut-targets-by-id/13znqeVDfPULc4J_lQLbyW_Kmfa03o63F/3-Research/GIACOMO/inequality paper/fig2.png", saver2, width = 10, height = 8, scale=1.4)

