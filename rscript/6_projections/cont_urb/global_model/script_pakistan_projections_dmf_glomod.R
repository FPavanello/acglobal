rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

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
library(relaimpo)
library(domir)
library(pscl)
library(margins)

# Set users
user <- 'fp'
user <- 'gf'

if (user=='fp') {
  stub <- 'F:/Il mio Drive/'
}

if (user=='gf') {
  stub <- "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}


# Use CDDs 18 threshold?
CDD18=T
CDD_type = ifelse(CDD18==T, "cdd_18", "cdd_24"); HDD_type = ifelse(CDD18==T, "hdd_18", "hdd_15")



# country specific data / parameters

countryname <- "Pakistan"
countryiso3 <- "PAK"

# 1 import DMF environment with trained models and data

load(paste0(stub, "results/regressions/for_projections/pak_dmcf.RData"))

data_c <- HH_Pakistan
#rm(HHPakistan)



# shapefile

gadm <- read_sf(paste0(stub, "data/shapefiles/Pakistan/gadm40_PAK_2.shp"))

# assumed parameters

rate_improvement_housing_1 <- 0.00375
rate_improvement_housing_2 <- 0.01

#######################


## 2) future drivers data for projections ##
# 2.1) import data#

gadm_0 <- st_as_sf(readRDS(paste0(stub, "rscripts/projections/ac_ely/gadm36_PAK_0_sp.rds")))
gadm_1 <- st_as_sf(readRDS(paste0(stub, "rscripts/projections/ac_ely/gadm36_PAK_1_sp.rds")))

custom_shape <- gadm %>% dplyr::select(geometry)

###

setwd(paste0(stub, "data/projections/new_data_jan_2022"))

# gdp downscaled (SSPS)
gdp_ssps <- list.files(path="gdp_downscaled_ssps", recursive = T, pattern="tif", full.names = T)

gdp_ssps_data <- lapply(gdp_ssps, raster)
gdp_ssps_data <- split(gdp_ssps_data, rep(1:5, each=26))
gdp_ssps_data <- lapply(gdp_ssps_data, stack)
names(gdp_ssps_data[[1]]) <- paste0("GDP_SSP1_", seq(1850, 2100, by=10))
names(gdp_ssps_data[[2]]) <- paste0("GDP_SSP2_", seq(1850, 2100, by=10))
names(gdp_ssps_data[[3]]) <- paste0("GDP_SSP3_", seq(1850, 2100, by=10))
names(gdp_ssps_data[[4]]) <- paste0("GDP_SSP4_", seq(1850, 2100, by=10))
names(gdp_ssps_data[[5]]) <- paste0("GDP_SSP5_", seq(1850, 2100, by=10))

for (i in 1:5){
  gdp_ssps_data[[i]] <- raster::subset(gdp_ssps_data[[i]], 17:26)
}

gdp_ssps_data <- stack(gdp_ssps_data[[1]], gdp_ssps_data[[2]], gdp_ssps_data[[3]],gdp_ssps_data[[4]],gdp_ssps_data[[5]])

# pop downscaled (SSPS)
pop_ssps <- list.files(path="pop_downscaled_spps", recursive = T, pattern="nc", full.names = T)
pop_ssps_data <- lapply(pop_ssps, stack)

for (i in 1:5){
  pop_ssps_data[[i]] <- raster::subset(pop_ssps_data[[i]], c(5+c(10*c(0:9))))
  
  pop_ssps_data[[i]] <- crop(pop_ssps_data[[i]], extent(gadm))
  
}

names(pop_ssps_data[[1]]) <- paste0("POP_SSP1_", seq(2010, 2100, by=10))
names(pop_ssps_data[[2]]) <- paste0("POP_SSP2_", seq(2010, 2100, by=10))
names(pop_ssps_data[[3]]) <- paste0("POP_SSP3_", seq(2010, 2100, by=10))
names(pop_ssps_data[[4]]) <- paste0("POP_SSP4_", seq(2010, 2100, by=10))
names(pop_ssps_data[[5]]) <- paste0("POP_SSP5_", seq(2010, 2100, by=10))

pop_ssps_data <- stack(pop_ssps_data[[1]], pop_ssps_data[[2]], pop_ssps_data[[3]],pop_ssps_data[[4]],pop_ssps_data[[5]])

# urban share downscaled (SSPS)
urban_ssps <- list.files(path="urban_share_downscaled_ssps/UrbanFraction_1_8_dgr_NETCDF_Projections_SSPs1-5_2010-2100_v1/", recursive = T, pattern="nc", full.names = T)

urban_ssps_data <- lapply(urban_ssps, raster)
urban_ssps_data <- lapply(urban_ssps_data, crop, extent(gadm))
urban_ssps_data <- split(urban_ssps_data, rep(1:5, each=10))
urban_ssps_data <- lapply(urban_ssps_data, stack)

names(urban_ssps_data[[1]]) <- paste0("URB_SSP1_", seq(2010, 2100, by=10))
names(urban_ssps_data[[2]]) <- paste0("URB_SSP2_", seq(2010, 2100, by=10))
names(urban_ssps_data[[3]]) <- paste0("URB_SSP3_", seq(2010, 2100, by=10))
names(urban_ssps_data[[4]]) <- paste0("URB_SSP4_", seq(2010, 2100, by=10))
names(urban_ssps_data[[5]]) <- paste0("URB_SSP5_", seq(2010, 2100, by=10))

#

urban_ssps_data <- stack(urban_ssps_data[[1]], urban_ssps_data[[2]], urban_ssps_data[[3]],urban_ssps_data[[4]],urban_ssps_data[[5]])

# pop features (SSPS)
pop_features <- read.csv("pop_features_ssps/SspDb_country_data_2013-06-12.csv")

pop_features <- filter(pop_features, REGION==countryiso3)

pop_features$SCENARIO <- substr(pop_features$SCENARIO, 1, 4)

pop_features <- pop_features %>% group_by(SCENARIO, REGION, VARIABLE) %>%  mutate_at(vars(contains('X')), funs(median))

pop_features <- dplyr::select(pop_features, 1,2,3,4,18, 20, 22, 24, 26, 28, 30, 32, 34, 36)

# 2.2) extract data into districts shapefile #

# extract the data
gdp_extracted <- exactextractr::exact_extract(gdp_ssps_data, custom_shape, fun="sum", max_cells_in_memory = 93355200 )
pop_extracted <- exactextractr::exact_extract(pop_ssps_data, custom_shape, fun="sum", max_cells_in_memory = 93355200 )

# weight_raster <- clamp(disaggregate(pop_ssps_data[[15]], fact=4), lower=1, useValues=T)
# 
# weight_raster <- projectRaster(weight_raster, urban_ssps_data[[1]])

urban_ssps_data <- raster::aggregate(urban_ssps_data, fact=4, fun="mean")

urban_extracted <- list()

for (i in 1:nlayers(urban_ssps_data)){
  weight_raster <- pop_ssps_data[[i]]
  weight_raster <- projectRaster(weight_raster, urban_ssps_data[[i]])
  urban_extracted[[i]] <- exactextractr::exact_extract(urban_ssps_data[[i]], custom_shape, fun="weighted_mean", weights=weight_raster, max_cells_in_memory = 93355200)
}

urban_extracted <- bind_cols(urban_extracted)

urban_extracted <- exactextractr::exact_extract(urban_ssps_data, custom_shape, fun="weighted_mean", weights=weight_raster, max_cells_in_memory = 93355200 )

#gini_extracted <- filter(gini, gini$variable==countryiso3) %>% dplyr::select(-1)

custom_shape <- bind_cols(custom_shape, gdp_extracted, pop_extracted, urban_extracted)

geo_bk <- custom_shape$geometry

custom_shape$geometry<-NULL

#############

# calculate gdp per capita

gdp_capita <- gdp_extracted / pop_extracted
gdp_capita[sapply(gdp_capita, is.infinite)] <- NA

colnames(gdp_capita) <- gsub("GDP_", "GDP_capita_", colnames(custom_shape)[1:50])

# adjust GDP per-capita growth to expenditure growth
gdp_capita = gdp_capita  - (mean(gdp_capita$sum.GDP_capita_SSP1_2010, na.rm=T) - mean(exp(data_c$ln_total_exp_usd_2011), na.rm=T))

# merge
custom_shape <- bind_cols(custom_shape,gdp_capita)

# calculate gdp per capita growth
growth_rate <- function(x)(x/lag(x)-1)*100

gdp_capita_growth <- gdp_capita

for (i in 1:nrow(gdp_capita)){
  gdp_capita_growth[i,c(1:10)] <- growth_rate(as.numeric(gdp_capita[i,c(1:10)])) 
  gdp_capita_growth[i,c(11:20)] <- growth_rate(as.numeric(gdp_capita[i,c(11:20)])) 
  gdp_capita_growth[i,c(21:30)] <- growth_rate(as.numeric(gdp_capita[i,c(21:30)])) 
  gdp_capita_growth[i,c(31:40)] <- growth_rate(as.numeric(gdp_capita[i,c(31:40)])) 
  gdp_capita_growth[i,c(41:50)] <- growth_rate(as.numeric(gdp_capita[i,c(41:50)])) 
}

colnames(gdp_capita_growth) <- gsub("GDP_", "GDP_capita_yearly_avg_growth_", colnames(custom_shape)[1:50])

custom_shape <- bind_cols(custom_shape,gdp_capita_growth)

# calculate urbanisation rate absolute change

urban_rate <- as.data.frame(custom_shape[,c(101:150)])

growth_rate_abs <- function(x)(x - lag(x))

urban_rate_growth <- urban_rate

for (i in 1:nrow(gdp_capita)){
  urban_rate_growth[i,c(1:10)] <- growth_rate_abs(as.numeric(urban_rate[i,c(1:10)])) 
  urban_rate_growth[i,c(11:20)] <- growth_rate_abs(as.numeric(urban_rate[i,c(11:20)])) 
  urban_rate_growth[i,c(21:30)] <- growth_rate_abs(as.numeric(urban_rate[i,c(21:30)])) 
  urban_rate_growth[i,c(31:40)] <- growth_rate_abs(as.numeric(urban_rate[i,c(31:40)])) 
  urban_rate_growth[i,c(41:50)] <- growth_rate_abs(as.numeric(urban_rate[i,c(41:50)]))
}

colnames(urban_rate_growth) <- gsub("sum.GDP_capita_yearly_avg_growth", "weighted_mean.URB_gr", colnames(gdp_capita_growth))

colnames(custom_shape)[c(101:150)] <- gsub("weighted_mean.URB_gr", "weighted_mean.URB", colnames(urban_rate_growth))

custom_shape <- bind_cols(custom_shape,urban_rate_growth)

custom_shape$geometry <- geo_bk

####
# other vars from formula (through scenario data, or, in some cases, assumptions)
####

# check what we had included in the ac drivers formula
#formula

# Gender head (national SSPs)

pop_features_all <- filter(pop_features, VARIABLE=="Population") %>% ungroup()
pop_features_all <- dplyr::select(pop_features_all, c(2, 5:14))

pop_features_all <- reshape2::melt(pop_features_all, c(1))
pop_features_all$variable <- as.numeric(as.character(gsub("X", "", pop_features_all$variable )))

pop_features_all <- group_by(pop_features_all, SCENARIO, variable) %>% dplyr::summarise(value=first(value))

pop_features_all_before_melt <- pop_features_all

pop_features_all <- pivot_wider(
  pop_features_all,
  names_from = c(SCENARIO, variable),
  values_from = value,
  names_prefix = "pop_",
  values_fn = first
)

pop_features_gender <- filter(pop_features, VARIABLE=="Population|Female")  %>% ungroup()
pop_features_gender <- dplyr::select(pop_features_gender, c(2, 5:14))

pop_features_gender <- reshape2::melt(pop_features_gender, c(1))
pop_features_gender$variable <- as.numeric(as.character(gsub("X", "", pop_features_gender$variable )))

pop_features_gender <- pop_features_gender[order( pop_features_gender$SCENARIO, pop_features_gender$variable ),]

pop_features_gender$value <- pop_features_gender$value / pop_features_all_before_melt$value

pop_features_gender <- pop_features_gender %>% group_by(SCENARIO) %>% mutate(value=growth_rate(value) / 100)

pop_features_gender <- pivot_wider(
  pop_features_gender,
  names_from = c(SCENARIO, variable),
  values_from = value,
  names_prefix = "women_"
)

custom_shape <- bind_cols(custom_shape, pop_features_gender)

# Age  (national SSPs)
pop_features_age <-filter(pop_features, grepl("Aged",VARIABLE) & !grepl("Education",VARIABLE))
pop_features_age <- dplyr::select(pop_features_age, c(3:4, 5:14))

pop_features_age$REGION<-NULL

pop_features_age$VARIABLE <- gsub("Population|Female|Aged", "", pop_features_age$VARIABLE)
pop_features_age$VARIABLE <- gsub("Population|Male|Aged", "", pop_features_age$VARIABLE)
pop_features_age$VARIABLE  <- gsub("\\|", "", pop_features_age$VARIABLE)
pop_features_age$VARIABLE  <- gsub("\\+", "", pop_features_age$VARIABLE)
pop_features_age$VARIABLE <- sapply(strsplit(pop_features_age$VARIABLE, split = "-", fixed = TRUE), function(k) mean(as.numeric(k)))

pop_features_age <- dplyr::group_by(pop_features_age, SCENARIO) %>% dplyr::mutate((X2010=X2010*VARIABLE)/sum(X2010), (X2020=X2020*VARIABLE)/sum(X2020), (X2030=X2030*VARIABLE)/sum(X2030), (X2040=X2040*VARIABLE)/sum(X2040), (X2050=X2050*VARIABLE)/sum(X2050))

pop_features_age <- dplyr::select(pop_features_age, c(1, 2, 13:17))

colnames(pop_features_age)[3:7] <- paste0("X", seq(2010, 2050, 10))

pop_features_age <- pop_features_age %>% group_by(SCENARIO) %>% dplyr::summarise(X2010=sum(X2010*VARIABLE), X2020=sum(X2020*VARIABLE), X2030=sum(X2030*VARIABLE), X2040=sum(X2040*VARIABLE), X2050=sum(X2050*VARIABLE))

pop_features_age <- reshape2::melt(pop_features_age, c(1))
pop_features_age$variable <- as.numeric(as.character(gsub("X", "", pop_features_age$variable )))

pop_features_age <- pop_features_age %>% group_by(SCENARIO) %>% mutate(value=growth_rate(value) / 100)

pop_features_age <- pivot_wider(
  pop_features_age,
  names_from = c(SCENARIO, variable),
  values_from = value,
  names_prefix = "age_"
)


custom_shape <- bind_cols(custom_shape, pop_features_age)

# Edu (national SSPs)

pop_features_edu <-filter(pop_features, grepl("Education",VARIABLE))
pop_features_edu <- dplyr::select(pop_features_edu, c(3:4, 5:14))
pop_features_edu$REGION<-NULL
pop_features_edu$VARIABLE <- gsub('.*\\|', "", pop_features_edu$VARIABLE)

pop_features_edu <- group_by(pop_features_edu, SCENARIO) %>% dplyr::mutate((X2010=X2010)/sum(X2010), (X2020=X2020)/sum(X2020), (X2030=X2030)/sum(X2030), (X2040=X2040)/sum(X2040), (X2050=X2050)/sum(X2050))

pop_features_edu <- dplyr::select(pop_features_edu, c(1, 2, 13:17))

colnames(pop_features_edu)[3:7] <- paste0("X", seq(2010, 2050, 10))

pop_features_edu <- group_by(pop_features_edu, SCENARIO, VARIABLE) %>% dplyr::summarise_all(., "sum")

pop_features_edu$VARIABLE[pop_features_edu$VARIABLE=="No Education"] <- "0"
pop_features_edu$VARIABLE[pop_features_edu$VARIABLE=="Primary Education"] <- "1"
pop_features_edu$VARIABLE[pop_features_edu$VARIABLE=="Secondary Education"] <-"2"
pop_features_edu$VARIABLE[pop_features_edu$VARIABLE=="Tertiary Education"] <- "3"
pop_features_edu$VARIABLE <- as.numeric(pop_features_edu$VARIABLE)

pop_features_edu <- split(pop_features_edu, pop_features_edu$SCENARIO)

for (i in 1:length(pop_features_edu)){
  
  pop_features_edu[[i]]$SCENARIO<- NULL
  pop_features_edu[[i]] <-  reshape(melt(pop_features_edu[[i]], c(1)), idvar = "variable", timevar = "VARIABLE", direction = "wide")
  colnames(pop_features_edu[[i]]) <- gsub("\\.", "\\_", colnames(pop_features_edu[[i]]))
  colnames(pop_features_edu[[i]]) <- gsub("variable", "year", colnames(pop_features_edu[[i]]))
  colnames(pop_features_edu[[i]]) <- gsub("value", "edu", colnames(pop_features_edu[[i]]))
  pop_features_edu[[i]]$year <- as.numeric(as.character(gsub("X", "",  pop_features_edu[[i]]$year )))
}


#

for (scen in 1:5){
  
  edu_ssp1 <- pop_features_edu[[scen]]
  row.names(edu_ssp1) <- 1:5
  edu_ssp1_hhs <- data.frame(edu_2010=as.numeric(as.character(data_c$edu_head_2)))
  edu_ssp1_hhs <- na.omit(edu_ssp1_hhs)
  
  edu_ssp1_hhs$edu_2020 <- NA
  edu_ssp1_hhs$edu_2030 <- NA
  edu_ssp1_hhs$edu_2040 <- NA
  edu_ssp1_hhs$edu_2050 <- NA
  
  
  for (year in seq(2020, 2050, 10)){
    
    for (i in c(2:0)){
      
      edu_ssp1_hhs[,paste0("edu_", year)] <- ifelse(edu_ssp1_hhs[,paste0("edu_", year-10)]==3, 3, edu_ssp1_hhs[,paste0("edu_", year)])
      
      if(as.numeric(unlist(edu_ssp1[paste0("edu_", i)])[2]) == 0){
        
        edu_ssp1_hhs[,paste0("edu_", year)][edu_ssp1_hhs[,paste0("edu_", year-10)] ==i] <- rep(i+1,sum(edu_ssp1_hhs[,paste0("edu_", year-10)]==i, na.rm=T))
        
      } else{
        
        edu_ssp1_hhs[,paste0("edu_", year)][edu_ssp1_hhs[,paste0("edu_", year-10)]==i] <- as.vector(sample(c(i, i+1), prob = c(edu_ssp1[,paste0("edu_", i)][edu_ssp1$year==year], edu_ssp1[,paste0("edu_", i+1)][edu_ssp1$year==year])/sum(c(edu_ssp1[,paste0("edu_", i)][edu_ssp1$year==year], edu_ssp1[,paste0("edu_", i+1)][edu_ssp1$year==year])), replace = T, size=sum(edu_ssp1_hhs[,paste0("edu_", year-10)]==i, na.rm=T)))
        
      }}}
  
  pop_features_edu[[scen]] <- edu_ssp1_hhs
  
}

for (i in 1:5){
  
  ssps <- c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")
  
  colnames(pop_features_edu[[i]]) <- paste0(colnames(pop_features_edu[[i]]), "_", ssps[i])
  
}

names(pop_features_edu) <- NULL
pop_features_edu <- as.data.frame(do.call(cbind, pop_features_edu))

data_c <- data_c[!is.na(data_c$edu_head_2),]

pop_features_edu <- dplyr::select(pop_features_edu, -contains("2010"))

data_c <- bind_cols(data_c, pop_features_edu)

########################################
########################################

# Housing index (shift)

# Scenario 1: design it

housing_index_lab_s1 <- data.frame(year=seq(2010, 2050, by=10), housing_index_lab_0=NA , housing_index_lab_1=NA , housing_index_lab_2=NA , housing_index_lab_3=NA , housing_index_lab_4=NA)

housing_index_lab_s1[1,c(2:6)] <-  table(data_c$housing_index) / sum(table(data_c$housing_index))

for (i in 2:5){
  for (j in 3:6){

    if(housing_index_lab_s1[i-1,j-1]<=rate_improvement_housing_1) {
      housing_index_lab_s1[i,j] <- housing_index_lab_s1[i-1,j]+housing_index_lab_s1[i-1,j-1]
      housing_index_lab_s1[i,j-1] <- 0

    } else{

      housing_index_lab_s1[i,j] <- housing_index_lab_s1[i-1,j] + rate_improvement_housing_1
      housing_index_lab_s1[i,j-1] <- housing_index_lab_s1[i-1,j-1] - rate_improvement_housing_1
      
    }

  }}

housing_index_lab_s1[2:5,2:6] <- housing_index_lab_s1[2:5,2:6] + (1- rowSums(housing_index_lab_s1[2:5,2:6])) / 5

# Scenario 2: design it

housing_index_lab_s2 <- data.frame(year=seq(2010, 2050, by=10), housing_index_lab_0=NA , housing_index_lab_1=NA , housing_index_lab_2=NA , housing_index_lab_3=NA , housing_index_lab_4=NA)

housing_index_lab_s2[1,c(2:6)] <-  table(data_c$housing_index) / sum(table(data_c$housing_index))

for (i in 2:5){
  for (j in 3:6){

    if(housing_index_lab_s2[i-1,j-1]<=rate_improvement_housing_2) {
      housing_index_lab_s2[i,j] <- housing_index_lab_s2[i-1,j]+housing_index_lab_s2[i-1,j-1]
      housing_index_lab_s2[i,j-1] <- 0

    } else{

      housing_index_lab_s2[i,j] <- housing_index_lab_s2[i-1,j] + rate_improvement_housing_1
      housing_index_lab_s2[i,j-1] <- housing_index_lab_s2[i-1,j-1] - rate_improvement_housing_1
      
    }

  }}

housing_index_lab_s2[2:5,2:6] <- housing_index_lab_s2[2:5,2:6] + (1- rowSums(housing_index_lab_s2[2:5,2:6])) / 5


# Scenario 1: apply it to households

housing_index_lab_s1_hhs <- data.frame(housing_index_lab_2010=as.numeric(as.character(data_c$housing_index)))
housing_index_lab_s1_hhs <- na.omit(housing_index_lab_s1_hhs)

housing_index_lab_s1_hhs$housing_index_lab_2020 <- NA
housing_index_lab_s1_hhs$housing_index_lab_2030 <- NA
housing_index_lab_s1_hhs$housing_index_lab_2040 <- NA
housing_index_lab_s1_hhs$housing_index_lab_2050 <- NA


for (year in seq(2020, 2050, 10)){

  for (i in c(3:0)){

    housing_index_lab_s1_hhs[,paste0("housing_index_lab_", year)] <- ifelse(housing_index_lab_s1_hhs[,paste0("housing_index_lab_", year-10)]==4, 4, housing_index_lab_s1_hhs[,paste0("housing_index_lab_", year)])

    if(as.numeric(unlist(housing_index_lab_s1[paste0("housing_index_lab_", i)])[2]) == 0){

      housing_index_lab_s1_hhs[,paste0("housing_index_lab_", year)][housing_index_lab_s1_hhs[,paste0("housing_index_lab_", year-10)] ==i] <- rep(i+1,sum(housing_index_lab_s1_hhs[,paste0("housing_index_lab_", year-10)]==i, na.rm=T))

    } else{

      housing_index_lab_s1_hhs[,paste0("housing_index_lab_", year)][housing_index_lab_s1_hhs[,paste0("housing_index_lab_", year-10)]==i] <- as.vector(sample(c(i, i+1), prob = c(housing_index_lab_s1[,paste0("housing_index_lab_", i)][housing_index_lab_s1$year==year], housing_index_lab_s1[,paste0("housing_index_lab_", i+1)][housing_index_lab_s1$year==year])/sum(c(housing_index_lab_s1[,paste0("housing_index_lab_", i)][housing_index_lab_s1$year==year], housing_index_lab_s1[,paste0("housing_index_lab_", i+1)][housing_index_lab_s1$year==year])), replace = T, size=sum(housing_index_lab_s1_hhs[,paste0("housing_index_lab_", year-10)]==i, na.rm=T)))

    }}}


# Scenario 2: apply it to households

housing_index_lab_s2_hhs <- data.frame(housing_index_lab_2010=as.numeric(as.character(data_c$housing_index)))
housing_index_lab_s2_hhs <- na.omit(housing_index_lab_s2_hhs)

housing_index_lab_s2_hhs$housing_index_lab_2020 <- NA
housing_index_lab_s2_hhs$housing_index_lab_2030 <- NA
housing_index_lab_s2_hhs$housing_index_lab_2040 <- NA
housing_index_lab_s2_hhs$housing_index_lab_2050 <- NA


for (year in seq(2020, 2050, 10)){

  for (i in c(3:0)){

    housing_index_lab_s2_hhs[,paste0("housing_index_lab_", year)] <- ifelse(housing_index_lab_s2_hhs[,paste0("housing_index_lab_", year-10)]==4, 4, housing_index_lab_s2_hhs[,paste0("housing_index_lab_", year)])

    if(as.numeric(unlist(housing_index_lab_s2[paste0("housing_index_lab_", i)])[2]) == 0){

      housing_index_lab_s2_hhs[,paste0("housing_index_lab_", year)][housing_index_lab_s2_hhs[,paste0("housing_index_lab_", year-10)] ==i] <- rep(i+1,sum(housing_index_lab_s2_hhs[,paste0("housing_index_lab_", year-10)]==i, na.rm=T))

    } else{

      housing_index_lab_s2_hhs[,paste0("housing_index_lab_", year)][housing_index_lab_s2_hhs[,paste0("housing_index_lab_", year-10)]==i] <- as.vector(sample(c(i, i+1), prob = c(housing_index_lab_s2[,paste0("housing_index_lab_", i)][housing_index_lab_s2$year==year], housing_index_lab_s2[,paste0("housing_index_lab_", i+1)][housing_index_lab_s2$year==year])/sum(c(housing_index_lab_s2[,paste0("housing_index_lab_", i)][housing_index_lab_s2$year==year], housing_index_lab_s2[,paste0("housing_index_lab_", i+1)][housing_index_lab_s2$year==year])), replace = T, size=sum(housing_index_lab_s2_hhs[,paste0("housing_index_lab_", year-10)]==i, na.rm=T)))

    }}}


# Scenarios: merge with data

data_c <- data_c[!is.na(data_c$housing_index),]

# housing_index_lab_s1_hhs <- as.data.frame(ifelse(housing_index_lab_s1_hhs==0 | housing_index_lab_s1_hhs==1, 1, ifelse(housing_index_lab_s1_hhs==2, 2, ifelse(housing_index_lab_s1_hhs==3 | housing_index_lab_s1_hhs==4, 3, NA))))
# 
# housing_index_lab_s2_hhs <- as.data.frame(ifelse(housing_index_lab_s2_hhs==0 | housing_index_lab_s2_hhs==1, 1, ifelse(housing_index_lab_s2_hhs==2, 2, ifelse(housing_index_lab_s2_hhs==3 | housing_index_lab_s2_hhs==4, 3, NA))))

housing_index_lab_s1_hhs$housing_index_lab_2010<-NULL
colnames(housing_index_lab_s1_hhs) <- paste0(colnames(housing_index_lab_s1_hhs), "_s1")

housing_index_lab_s2_hhs$housing_index_lab_2010<-NULL
colnames(housing_index_lab_s2_hhs) <- paste0(colnames(housing_index_lab_s2_hhs), "_s2")

data_c <- bind_cols(data_c, housing_index_lab_s1_hhs, housing_index_lab_s2_hhs)

###

# Merge spatial projection data with household survey data

library(fuzzyjoin); library(dplyr);data_c$NAME_2.y<-NULL;

data_c_map <- data_c
data_c_map$NAME_2 <- data_c$district
data_c_map <- data_c_map[!duplicated(data_c_map[ , c("NAME_2")]), ] 

data_c_sp <- stringdist_join(data_c_map, gadm, 
                             by = "NAME_2",
                             mode = "left",
                             ignore_case = TRUE, 
                             method = "jw", 
                             max_dist = 99, 
                             distance_col = "dist") %>%
  group_by(NAME_2.x) %>%
  slice_min(order_by = dist, n = 1)

data_c_sp <- dplyr::select(data_c_sp, NAME_2.x, NAME_2.y, geometry)

data_c_sp <- merge(data_c, data_c_sp, by.x="district", by.y="NAME_2.x")

custom_shape$NAME_2 <- gadm$NAME_2

data_c_sp <- merge(data_c_sp, custom_shape, by.x="NAME_2.y", by.y="NAME_2")


#
# Project future expenditure

data_c_sp$exp_cap_usd_SSP1_2010 <- exp(data_c_sp$ln_total_exp_usd_2011)
data_c_sp$exp_cap_usd_SSP2_2010 <- exp(data_c_sp$ln_total_exp_usd_2011)
data_c_sp$exp_cap_usd_SSP3_2010 <- exp(data_c_sp$ln_total_exp_usd_2011)
data_c_sp$exp_cap_usd_SSP4_2010 <- exp(data_c_sp$ln_total_exp_usd_2011)
data_c_sp$exp_cap_usd_SSP5_2010 <- exp(data_c_sp$ln_total_exp_usd_2011)

data_c_sp$exp_cap_usd_SSP1_2020 <- NA
data_c_sp$exp_cap_usd_SSP2_2020 <- NA
data_c_sp$exp_cap_usd_SSP3_2020 <- NA
data_c_sp$exp_cap_usd_SSP4_2020 <- NA
data_c_sp$exp_cap_usd_SSP5_2020 <- NA

data_c_sp$exp_cap_usd_SSP1_2030 <- NA
data_c_sp$exp_cap_usd_SSP2_2030 <- NA
data_c_sp$exp_cap_usd_SSP3_2030 <- NA
data_c_sp$exp_cap_usd_SSP4_2030 <- NA
data_c_sp$exp_cap_usd_SSP5_2030 <- NA

data_c_sp$exp_cap_usd_SSP1_2040 <- NA
data_c_sp$exp_cap_usd_SSP2_2040 <- NA
data_c_sp$exp_cap_usd_SSP3_2040 <- NA
data_c_sp$exp_cap_usd_SSP4_2040 <- NA
data_c_sp$exp_cap_usd_SSP5_2040 <- NA

data_c_sp$exp_cap_usd_SSP1_2050 <- NA
data_c_sp$exp_cap_usd_SSP2_2050 <- NA
data_c_sp$exp_cap_usd_SSP3_2050 <- NA
data_c_sp$exp_cap_usd_SSP4_2050 <- NA
data_c_sp$exp_cap_usd_SSP5_2050 <- NA


load(paste0(stub, "rscripts/global_spline/supporting_data/adj_factors.Rds"))

ss <- filter(ss, ISO3==countryiso3)

#

for (ssp in c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")){
  
  for (year in seq(2020, 2050, 10)){
    
    
    data_c_sp[,paste0("exp_cap_usd_", ssp, "_", year)] <- data_c_sp[,paste0("exp_cap_usd_", ssp, "_", year-10)] * (1+((data_c_sp[,paste0("sum.GDP_capita_yearly_avg_growth_", ssp, "_", year)]/100)))
    
  }}


for (ssp in c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")){
  
  for (year in seq(2010, 2050, 10)){
    
    
    data_c_sp[,paste0("exp_cap_usd_", ssp, "_", year)] <- log( data_c_sp[,paste0("exp_cap_usd_", ssp, "_", year)])
    
  }}

#


# # House ownership (shift based on growth of expenditure)
# 
# table(data_c_sp$ownership_d)
# 
# model_ownership_d <- glm2("ownership_d ~ exp_cap_usd_SSP2_2010",data = data_c_sp, family = binomial(logit), na.action=na.omit, weights=weight)
# 
# output <- list()
# 
# for (ssp in c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")){
# 
#   for (year in seq(2020, 2050, 10)){
# 
# 
#     output[[paste0("ownership_d_", ssp, "_", year)]] <- ifelse(exp(model_ownership_d$coefficients[1] + model_ownership_d$coefficients[2]*data_c_sp[,paste0("exp_cap_usd_", ssp, "_", year)]) / (1 + exp(model_ownership_d$coefficients[1] + model_ownership_d$coefficients[2]*data_c_sp[,paste0("exp_cap_usd_", ssp, "_", year)]))>0.75, 1, 0)
# 
#   }}
# 
# output <- as.data.frame(do.call(cbind, output))
# 
# data_c_sp <- bind_cols(data_c_sp, output)
# 
# # Fan ownership (shift based on growth of expenditure)
# 
# table(data_c_sp$fan)
# 
# model_fan <- glm2("fan ~ exp_cap_usd_SSP2_2010",data = data_c_sp, family = binomial(logit), na.action=na.omit, weights=weight)
# 
# #
# 
# data_c_sp$fan_SSP1_2010 <- data_c_sp$fan
# data_c_sp$fan_SSP2_2010 <- data_c_sp$fan
# data_c_sp$fan_SSP3_2010 <- data_c_sp$fan
# data_c_sp$fan_SSP4_2010 <- data_c_sp$fan
# data_c_sp$fan_SSP5_2010 <- data_c_sp$fan
# 
# output <- list()
# 
# for (ssp in c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5")){
# 
#   for (year in seq(2020, 2050, 10)){
# 
# 
#     output[[paste0("fan_", ssp, "_", year)]] <- ifelse(exp(model_fan$coefficients[1] + model_fan$coefficients[2]*data_c_sp[,paste0("exp_cap_usd_", ssp, "_", year)]) / (1 + exp(model_fan$coefficients[1] + model_fan$coefficients[2]*data_c_sp[,paste0("exp_cap_usd_", ssp, "_", year)]))>0.75, 1, 0)
# 
#   }}
# 
# output <- as.data.frame(do.call(cbind, output))
# 
# data_c_sp <- bind_cols(data_c_sp, output)

gc()

####### import CDDs and HDDs projections and merge them to the survey data

# import CMIP6 data

cdd_hist_cmip6 <- readRDS(paste0(stub, "data/projections/climate/processed/", CDD_type, "/", countryiso3, "_CDD_model_ensemble_median_hist.rds"))

cdd_245_cmip6 <-  readRDS(paste0(stub, "data/projections/climate/processed/", CDD_type, "/", countryiso3, "_CDD_model_ensemble_median_rcp45_ssp2.rds"))

cdd_585_cmip6 <-  readRDS(paste0(stub, "data/projections/climate/processed/", CDD_type, "/", countryiso3, "_CDD_model_ensemble_median_rcp85_ssp5.rds"))

#

hdd_hist_cmip6 <-  readRDS(paste0(stub, "data/projections/climate/processed/", HDD_type, "/", countryiso3, "_HDD_model_ensemble_median_hist.rds"))

hdd_245_cmip6 <- readRDS(paste0(stub, "data/projections/climate/processed/", HDD_type, "/", countryiso3, "_HDD_model_ensemble_median_rcp45_ssp2.rds"))

hdd_585_cmip6 <-  readRDS(paste0(stub, "data/projections/climate/processed/", HDD_type, "/", countryiso3, "_HDD_model_ensemble_median_rcp85_ssp5.rds"))

#

cmip6_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c("country", "state2", "district2"), all.x = TRUE),
                       list(cdd_hist_cmip6, cdd_245_cmip6, cdd_585_cmip6, hdd_hist_cmip6, hdd_245_cmip6, hdd_585_cmip6))


data_c_map <- dplyr::select(data_c_map, district)
data_c_map$district2 <- data_c_map$district

cmip6_merged <- stringdist_join(cmip6_merged, data_c_map, 
                                by = "district2",
                                mode = "left",
                                ignore_case = TRUE, 
                                method = "jw", 
                                max_dist = 99, 
                                distance_col = "dist") %>%
  group_by(district2.y) %>%
  slice_min(order_by = dist, n = 1)


data_c_sp <- merge(data_c_sp, cmip6_merged, by.x="district", by.y="district2.y")


for (year in seq(2020, 2050, 10)){
  
  for (scen in c("rcp45_ssp2", "rcp85_ssp5")){
    
    if (year <2050){
      data_c_sp[,paste0("delta_", scen, "_", year)] <- (data_c_sp[,paste0("mean_CDD_", year-5, "_", scen)] + data_c_sp[,paste0("mean_CDD_", year+5, "_", scen)])/2 - data_c_sp$mean_CDD_modelens_median_sy
    }
    
    if (year == 2050){
      data_c_sp[,paste0("delta_", scen, "_", year)] <- (data_c_sp[,paste0("mean_CDD_", year, "_", scen)]) - data_c_sp$mean_CDD_modelens_median_sy
    }
    
  }}

for (year in seq(2020, 2050, 10)){
  
  for (scen in c("rcp45_ssp2", "rcp85_ssp5")){
    
    data_c_sp[,paste0("mean_CDD_", year, "_", scen)] <- data_c_sp[,paste0("delta_", scen, "_", year)] + (data_c_sp$mean_CDD18_db * 100)
    
    
  }} 


for (year in seq(2020, 2050, 10)){
  
  for (scen in c("rcp45_ssp2", "rcp85_ssp5")){
    #
    
    if (year <2050){
      data_c_sp[,paste0("delta_", scen, "_", year)] <- (data_c_sp[,paste0("mean_HDD_", year-5, "_", scen)] + data_c_sp[,paste0("mean_HDD_", year+5, "_", scen)])/2 - data_c_sp$mean_HDD_modelens_median_sy
    }
    
    if (year == 2050){
      data_c_sp[,paste0("delta_", scen, "_", year)] <- (data_c_sp[,paste0("mean_HDD_", year, "_", scen)]) - data_c_sp$mean_HDD_modelens_median_sy
    }
    
  }}    

for (year in seq(2020, 2050, 10)){
  
  for (scen in c("rcp45_ssp2", "rcp85_ssp5")){
    
    data_c_sp[,paste0("mean_HDD_", year, "_", scen)] <- data_c_sp[,paste0("delta_", scen, "_", year)] + (data_c_sp$mean_HDD18_db * 100)
    
  }}


# set negative to 0

for (year in seq(2020, 2050, 10)){
  
  for (scen in c("rcp45_ssp2", "rcp85_ssp5")){
    
    data_c_sp[,paste0("mean_CDD_", year, "_", scen)] <- ifelse(data_c_sp[,paste0("mean_CDD_", year, "_", scen)]<0, 0, data_c_sp[,paste0("mean_CDD_", year, "_", scen)])
    
    #
    
    data_c_sp[,paste0("mean_HDD_", year, "_", scen)] <- ifelse(data_c_sp[,paste0("mean_HDD_", year, "_", scen)]<0, 0, data_c_sp[,paste0("mean_HDD_", year, "_", scen)])
    
  }}

########
# Plot future data#

# custom_shape_plotting <- custom_shape
# custom_shape_plotting$state <- gadm$CVE_ENT
# 
# gdp_trend <- dplyr::select(custom_shape_plotting, starts_with("sum.GDP_S"), state)
# pop_trend <- dplyr::select(custom_shape_plotting, starts_with("sum.POP_"), state)
# gdp_capita_trend <- dplyr::select(custom_shape_plotting, starts_with("sum.GDP_capita_S"), state)
# gdp_capita_growth_trend <- dplyr::select(custom_shape_plotting, starts_with("sum.GDP_capita_yearly_avg_growth_S"), state)
# urb_trend <- dplyr::select(custom_shape_plotting, starts_with("weighted_mean.URB_"), state)
# 
# #
# 
# gdp_trend <-reshape2::melt(gdp_trend, "state")
# gdp_trend$variable <- gsub("sum.GDP_", "", gdp_trend$variable)
# gdp_trend$SSP <- substr(gdp_trend$variable , 1, 4)
# gdp_trend$year <- substr(gdp_trend$variable , 6, 9)
# gdp_trend$variable <- NULL
# 
# gdp_trend_national <- gdp_trend %>% group_by(SSP, year) %>% dplyr::summarise(value=sum(value, na.rm=T))
# 
# plot_gdp_1 <- ggplot(gdp_trend_national)+
#   geom_line(aes(x=year, y=value, colour=SSP, group=SSP))+
#   scale_colour_npg()
# 
# ggsave(paste0(stub, "results/graphs/plot_gdp_1_", countryiso3, ".png"), plot_gdp_1)
# 
# gdp_trend_regional <- gdp_trend %>% group_by(SSP, year, state) %>% dplyr::summarise(value=sum(value, na.rm=T))
# 
# plot_gdp_2 <- ggplot(gdp_trend_regional)+
#   geom_line(aes(x=year, y=value, colour=SSP, group=SSP))+
#   facet_wrap(~ state, ncol=10)
# 
# ggsave(paste0(stub, "results/graphs/plot_gdp_2_", countryiso3, ".png"), plot_gdp_2)
# 
# #
# 
# pop_trend <-reshape2::melt(pop_trend, "state")
# pop_trend$variable <- gsub("sum.POP_", "", pop_trend$variable)
# pop_trend$SSP <- substr(pop_trend$variable , 1, 4)
# pop_trend$year <- substr(pop_trend$variable , 6, 9)
# pop_trend$variable <- NULL
# 
# pop_trend_national <- pop_trend %>% group_by(SSP, year) %>% dplyr::summarise(value=sum(value, na.rm=T))
# 
# plot_pop_1 <- ggplot(pop_trend_national)+
#   geom_line(aes(x=year, y=value, colour=SSP, group=SSP))+
#   scale_colour_npg()
# 
# ggsave(paste0(stub, "results/graphs/plot_pop_1_", countryiso3, ".png"), plot_pop_1)
# 
# pop_trend_regional <- pop_trend %>% group_by(SSP, year, state) %>% dplyr::summarise(value=sum(value, na.rm=T))
# 
# plot_pop_2 <- ggplot(pop_trend_regional)+
#   geom_line(aes(x=year, y=value, colour=SSP, group=SSP))+
#   facet_wrap(~ state, ncol=10)
# 
# ggsave(paste0(stub, "results/graphs/plot_pop_2_", countryiso3, ".png"), plot_pop_2)
# 
# #
# gdp_capita_trend <-reshape2::melt(gdp_capita_trend, "state")
# gdp_capita_trend$variable <- gsub("sum.GDP_capita_", "", gdp_capita_trend$variable)
# gdp_capita_trend$SSP <- substr(gdp_capita_trend$variable , 1, 4)
# gdp_capita_trend$year <- substr(gdp_capita_trend$variable , 6, 9)
# gdp_capita_trend$variable <- NULL
# 
# gdp_capita_trend_national <- gdp_capita_trend %>% group_by(SSP, year) %>% dplyr::summarise(value=mean(value, na.rm=T))
# 
# plot_gdp_capita_1 <- ggplot(gdp_capita_trend_national)+
#   geom_line(aes(x=year, y=value, colour=SSP, group=SSP))+
#   scale_colour_npg()
# 
# ggsave(paste0(stub, "results/graphs/plot_gdp_capita_1_", countryiso3, ".png"), plot_gdp_capita_1)
# 
# 
# gdp_capita_trend_regional <- gdp_capita_trend %>% group_by(SSP, year, state) %>% dplyr::summarise(value=mean(value, na.rm=T))
# 
# plot_gdp_capita_2 <- ggplot(gdp_capita_trend_regional)+
#   geom_line(aes(x=year, y=value, colour=SSP, group=SSP))+
#   facet_wrap(~ state, ncol=10)
# 
# ggsave(paste0(stub, "results/graphs/plot_gdp_capita_2_", countryiso3, ".png"), plot_gdp_capita_2)
# 
# #
# 
# urb_trend <-reshape2::melt(urb_trend, "state")
# urb_trend$variable <- gsub("weighted_mean.URB_", "", urb_trend$variable)
# urb_trend$SSP <- substr(urb_trend$variable , 1, 4)
# urb_trend$year <- substr(urb_trend$variable , 6, 9)
# urb_trend$variable <- NULL
# 
# urb_trend_national <- urb_trend %>% group_by(SSP, year) %>% dplyr::summarise(value=mean(value, na.rm=T))
# 
# plot_urb_1 <- ggplot(urb_trend_national)+
#   geom_line(aes(x=year, y=value, colour=SSP, group=SSP))+
#   scale_colour_npg()
# 
# ggsave(paste0(stub, "results/graphs/plot_urb_1_", countryiso3, ".png"), plot_urb_1)
# 
# 
# urb_trend_regional <- urb_trend %>% group_by(SSP, year, state) %>% dplyr::summarise(value=mean(value, na.rm=T))
# 
# plot_urb_2 <- ggplot(urb_trend_regional)+
#   geom_line(aes(x=year, y=value, colour=SSP, group=SSP))+
#   facet_wrap(~ state, ncol=10)
# 
# ggsave(paste0(stub, "results/graphs/plot_urb_2_", countryiso3, ".png"), plot_urb_2)
# 
# 
# ##
# 
# housing_index_lab_s1_plot <- pivot_longer(
#   housing_index_lab_s1,
#   cols = 2:6,
#   names_prefix = "housing_index_lab_")
# 
# housing_index_lab_plot <- ggplot(housing_index_lab_s1_plot)+
#   geom_col(aes(x=year, y=value, fill=name))+
#   scale_fill_npg()
# 
# ggsave(paste0(stub, "results/graphs/housing_index_lab_improvement_s1_", countryiso3, ".png"), housing_index_lab_plot)
# 
# #
# 
# housing_index_lab_s2_plot <- pivot_longer(
#   housing_index_lab_s2,
#   cols = 2:6,
#   names_prefix = "housing_index_lab_")
# 
# housing_index_lab_plot <- ggplot(housing_index_lab_s2_plot)+
#   geom_col(aes(x=year, y=value, fill=name))+
#   scale_fill_npg()
# 
# ggsave(paste0(stub, "results/graphs/housing_index_lab_improvement_s2_", countryiso3, ".png"), housing_index_lab_plot)

#

data_c_sp_export <- data_c_sp

data_c_sp_export$geometry.x <- NULL
data_c_sp_export$geometry.y <- NULL

orig_data <- data_c_sp

load(paste0(stub, "results/regressions/for_projections/global_wgt_dmcf.Rdata")); data_c_sp$country = data_c_sp$country.x

#

data_c_sp_export <- data_c_sp

data_c_sp_export$geometry.x <- NULL
data_c_sp_export$geometry.y <- NULL

orig_data <- data_c_sp

load(paste0(stub, "results/regressions/for_projections/global_wgt_dmcf.Rdata")); data_c_sp$country = data_c_sp$country.x

## 3) Make projections based on trained models and extracted data ##
# 3.1) AC adoption projections

orig_data <- data_c_sp
orig_data$ac <- NULL

orig_data_bk <- orig_data

library(imputeTS)
data_c_sp <- na_mean(data_c_sp)

output <- list()

# loop for all ssps and time-steps

for (ssp in c("SSP2", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  orig_data_bk <- data_c_sp
  
  output2 <- list()
  
  for (year in seq(2020, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$ln_total_exp_usd_2011 = data_c_sp[,paste0("exp_cap_usd_", ssp, "_", (year))]
    
    orig_data$mean_CDD18_db  = data_c_sp[,paste0("mean_CDD_", year, "_", rcp, "_", tolower(ssp))] / 100 # 10-year average bins
    
    orig_data$mean_CDD18_db2 = orig_data$mean_CDD18_db^2
    
    orig_data$mean_CDD18_db_exp = orig_data$mean_CDD18_db*orig_data$ln_total_exp_usd_2011
    
    orig_data$mean_CDD18_db2_exp = (orig_data$mean_CDD18_db^2)*orig_data$ln_total_exp_usd_2011
	
	orig_data$curr_CDD18_db = orig_data$mean_CDD18_db 
    
    orig_data$curr_CDD18_db2 = orig_data$curr_CDD18_db^2
	
    orig_data$mean_HDD18_db  = data_c_sp[,paste0("mean_HDD_", year, "_", rcp, "_", tolower(ssp))] / 100  # 10-year average bins
    
    orig_data$edu_head_2 <- as.factor(data_c_sp[,paste0("edu_", year, "_", ssp)])
    
    orig_data$age_head <- round(orig_data$age_head * (1 + data_c_sp[,paste0("age_", ssp, "_", year)]), 0)
    
    orig_data$country = data_c$country; orig_data$ownership_d = as.factor(orig_data$ownership_d)
    
    #orig_data$fan = as.factor(data_c_sp[,paste0("fan_", ssp, "_", (year))])
    
    orig_data$housing_index = as.factor(data_c_sp[,paste0("housing_index_lab_", year, "_s1")])
    
	orig_data$urban_sh = data_c_sp[,paste0("weighted_mean.URB_", ssp, "_", (year))]
    
    #
    projected <- predict(reg_ac, orig_data, type="response")
    # projected <- ifelse(as.numeric(projected)>0.5, 1, 0)
    
    # if (year>2020){
      # projected <- ifelse(output2[[as.character(year-10)]]==1, 1, projected)
    # }
    
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output[[as.character(ssp)]] <- output2
  
}

###

# bind results together
future_ac_adoption <- unlist(output, recursive = FALSE)
future_ac_adoption <- as.data.frame(do.call("cbind", future_ac_adoption))

data_c_sp <- bind_cols(data_c_sp, future_ac_adoption)

data_c_sp$geometry.x <- NULL
data_c_sp$geometry.y <- NULL

#############

# project weights (dynamic loops for names)

for (ssp in c("SSP2", "SSP5")){
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  for (year in seq(2020, 2050, 10)){
    
    if (year == 2020){
      
      data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))] <- data_c_sp$weight 
      
    } else {
      
      
      data_c_sp[,paste0("pop_gr", year, "_", rcp, "_", tolower(ssp))] = (data_c_sp[,paste0("sum.POP_", ssp, "_", year)] / data_c_sp[,paste0("sum.POP_", ssp, "_", (year-10))]) - 1
      
      
      data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))] = (1 + data_c_sp[,paste0("pop_gr", year, "_", rcp, "_", tolower(ssp))]) *  data_c_sp[,paste0("weight_", (year - 10), "_", rcp, "_", tolower(ssp))]
      
    }
    
    future_ac_adoption[,paste0(ssp, ".", year, "_wgt")] =  future_ac_adoption[,paste0(ssp, ".", year)] * (data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))] / sum(data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))], na.rm=T))
    
  }}

national_summary_ac <- future_ac_adoption %>%
  dplyr::summarise_at(vars(contains("wgt")), sum, na.rm=T) %>%
  pivot_longer(cols = 1:8, names_to = c('ssp', 'year'), names_sep = ".") %>%
  mutate(ssp=rep(c("SSP245","SSP585"), each=4), year=rep(seq(2020, 2050, 10), 2))


#### add bias correction?

future_ac_adoption$state3 <- data_c_sp$district

regional_summary_ac <- future_ac_adoption %>%
  group_by(state3) %>%
  dplyr::summarise_all(mean, na.rm=T) %>%
  pivot_longer(cols = 2:9, names_to = c('ssp', 'year'), names_sep = ".") %>%
  mutate(ssp=rep(rep(c("SSP245", "SSP585"), each=4), length(unique(data_c_sp$district))), year=rep(rep(seq(2020, 2050, 10), 2), length(unique(data_c_sp$district))))

# plot projections

line_country_ac <- ggplot(national_summary_ac)+
  geom_line(aes(x=year, y=value*100, colour=ssp, group=ssp))+
  scale_colour_npg(name="Scenario")+
  xlab("Year")+
  ylab("AC penetration rate (%)")

ggsave(paste0(stub, "results/graphs/line_country_ac_", countryiso3, ".png"), line_country_ac)

line_region_ac <- ggplot(regional_summary_ac)+
  geom_line(aes(x=year, y=value*100, colour=ssp, group=ssp))+
  facet_wrap(vars(state3))+
  scale_colour_npg(name="Scenario")+
  xlab("Year")+
  ylab("AC penetration rate (%)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=.1))

ggsave(paste0(stub, "results/graphs/line_region_ac_", countryiso3, ".png"), line_region_ac)


regional_summary_ac$NAME_2 <- regional_summary_ac$state3

regional_summary_ac <- stringdist_join(regional_summary_ac, gadm, 
                             by = "NAME_2",
                             mode = "left",
                             ignore_case = TRUE, 
                             method = "jw", 
                             max_dist = 99, 
                             distance_col = "dist") %>%
  group_by(NAME_2.x) %>%
  slice_min(order_by = dist, n = 1)

regional_summary_ac <- st_as_sf(regional_summary_ac)

map_ac <- ggplot(regional_summary_ac)+
  geom_sf(aes(fill=value*100), size = 0.2)+
  facet_wrap(vars(ssp, year), ncol = 4)+
  scale_fill_viridis_c(name="AC penetration (%)")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

ggsave(paste0(stub, "results/graphs/map_ac_", countryiso3, ".png"), map_ac, scale=1.75)

####

# Run decomposition analysis

# historical importance

data_c_bk <- data_c

ely_formula <- ln_ely_q ~ ac + curr_CDD18_db + curr_HDD18_db + ln_total_exp_usd_2011 + n_members + 
  sh_under16 + as.factor(ownership_d) + edu_head_2 + housing_index_lab + 
  age_head + sex_head + urban_sh

lm1 <- lm(ely_formula, data = data_c, na.action=na.omit)

metrics2 <- calc.relimp(lm1, type = c("lmg"), rela=T)

metrics2 <- metrics2@lmg

metrics2_n <- melt(metrics2)
metrics2_n$var <- rownames(metrics2_n)
rownames(metrics2_n) <- NULL

metrics2_n$type <- ifelse(metrics2_n$var=="ac", "AC", "Social drivers")
metrics2_n$type <- ifelse(metrics2_n$var=="urban_sh", "Urbanisation",  metrics2_n$type)
metrics2_n$type <- ifelse(metrics2_n$var=="ln_total_exp_usd_2011", "Expenditure", metrics2_n$type)
metrics2_n$type <- ifelse(metrics2_n$var=="curr_CDD18_db", "CDDs/HDDs", metrics2_n$type)
metrics2_n$type <- ifelse(metrics2_n$var=="curr_HDD18_db", "CDDs/HDDs", metrics2_n$type)
metrics2_n$x = 1

hist_dec <- ggplot(metrics2_n)+
  geom_bar(aes(x=x, y=value, fill=type, group=type), position="stack", stat="identity") + ggsci::scale_fill_npg(name="Driver") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Relative importance")

#ggsave(paste0(paste0(stub, "results/graphs/historical_decomp.png")))

#

year <- 2050

for (ssp in c("SSP2", "SSP5")){
  
rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")

orig_data_bk <- data_c_sp

baseline_hist <- as.numeric(predict(lm1, type="response"))

#

orig_data <- orig_data_bk

orig_data$ac = as.factor(ifelse(data_c_sp[,paste0(ssp, ".", (year))]>0.5, 1, 0))

orig_data$ln_total_exp_usd_2011 = data_c_sp[,paste0("exp_cap_usd_", ssp, "_", (year))]

orig_data$curr_CDD18_db  = data_c_sp[,paste0("mean_CDD_", year, "_", rcp, "_", tolower(ssp))] / 100 # survey-year 

orig_data$curr_HDD18_db  = data_c_sp[,paste0("mean_HDD_", year, "_", rcp, "_", tolower(ssp))] / 100 # survey-year 

orig_data$edu_head_2 <- as.factor(data_c_sp[,paste0("edu_", year, "_", ssp)])

orig_data$age_head <- round(orig_data$age_head * (1 + data_c_sp[,paste0("age_", ssp, "_", year)]), 0)

orig_data$country = data_c$country; orig_data$ownership_d = as.factor(orig_data$ownership_d)

#orig_data$fan = as.factor(data_c_sp[,paste0("fan_", ssp, "_", (year))])

orig_data$housing_index = as.factor(as.character(data_c_sp[,paste0("housing_index_lab_", year, "_s1")]))

orig_data$urban_sh = data_c_sp[,paste0("weighted_mean.URB_", ssp, "_", (year))]

#

total <- as.numeric(predict(lm1, orig_data, type="response"))

# 

orig_data <- orig_data_bk

orig_data$ac = as.factor(ifelse(data_c_sp[,paste0(ssp, ".", (year))]>0.5, 1, 0))

decomp_ac <- as.numeric(predict(lm1, orig_data, type="response"))

#


orig_data <- orig_data_bk

orig_data$ln_total_exp_usd_2011 =  data_c_sp[,paste0("exp_cap_usd_", ssp, "_", (year))]

decomp_exp <- as.numeric(predict(lm1, orig_data, type="response"))

#


orig_data <- orig_data_bk

orig_data$curr_CDD18_db = data_c_sp[,paste0("mean_CDD_", year, "_", rcp, "_", tolower(ssp))] / 100 # year specific

decomp_cdds <- as.numeric(predict(lm1, orig_data, type="response"))

#

orig_data <- orig_data_bk

orig_data$curr_HDD18_db = data_c_sp[,paste0("mean_HDD_", year, "_", rcp, "_", tolower(ssp))] / 100 # year specific

decomp_hdds <- as.numeric(predict(lm1, orig_data, type="response"))

#

orig_data <- orig_data_bk

orig_data$urban_sh = data_c_sp[,paste0("weighted_mean.URB_", ssp, "_", (year))]

decomp_urb <- as.numeric(predict(lm1, orig_data, type="response"))


baseline_hist <- weighted.mean(exp(baseline_hist), data_c$weight, na.rm=T)
decomp_ac <- weighted.mean(exp(decomp_ac), data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))], na.rm=T)
decomp_exp <- weighted.mean(exp(decomp_exp), data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))], na.rm=T)
decomp_cdds <- weighted.mean(exp(decomp_cdds), data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))], na.rm=T)
decomp_hdds <- weighted.mean(exp(decomp_hdds), data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))], na.rm=T)
decomp_urb <- weighted.mean(exp(decomp_urb), data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))], na.rm=T)
decomp_total <- weighted.mean(exp(total), data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))], na.rm=T)


decomp_fut <- as.data.frame(rbind(baseline_hist, decomp_ac - baseline_hist, decomp_exp - baseline_hist, decomp_cdds - baseline_hist, decomp_hdds - baseline_hist,  decomp_urb - baseline_hist, decomp_total - baseline_hist))

decomp_fut$V1  <- ifelse(decomp_fut$V1<0, 0, decomp_fut$V1)

rownames(decomp_fut) <- 1:nrow(decomp_fut)

decomp_fut$V1[7] <- decomp_fut$V1[7] - decomp_fut$V1[6] - decomp_fut$V1[5] - decomp_fut$V1[4] - decomp_fut$V1[3] - decomp_fut$V1[2]

decomp_fut$V1  <- ifelse(decomp_fut$V1<0, 0, decomp_fut$V1)

decomp_fut$type <- NA
decomp_fut$type <- ifelse(rownames(decomp_fut)==1, "Historical", decomp_fut$type)
decomp_fut$type <- ifelse(rownames(decomp_fut)==2, "Future AC change", decomp_fut$type)
decomp_fut$type <- ifelse(rownames(decomp_fut)==3, "Future expenditure change", decomp_fut$type)
decomp_fut$type <- ifelse(rownames(decomp_fut)==4, "Future CDDs/HDDs change", decomp_fut$type)
decomp_fut$type <- ifelse(rownames(decomp_fut)==5, "Future CDDs/HDDs change", decomp_fut$type)
decomp_fut$type <- ifelse(rownames(decomp_fut)==6, "Future urbanisation", decomp_fut$type)
decomp_fut$type <- ifelse(rownames(decomp_fut)==7, "Future social drivers change", decomp_fut$type)

decomp_fut$x = 1

for (i in 2:7){
  decomp_fut$V1[i] <- (decomp_fut$V1[i] + decomp_fut$V1[1]) / decomp_fut$V1[1] - 1
}

decomp_fut <- filter(decomp_fut, type!="Historical")

# ggplot(decomp_fut)+
#   geom_bar(aes(x=x, y=V1, fill=type, group=type), position="stack", stat="identity") + ggsci::scale_fill_npg()

decomp_fut$time <- "Future"

colnames(decomp_fut)[1] <- "value"


#

metrics2_n_s <- group_by(metrics2_n, type) %>% dplyr::summarise(value=sum(value))

metrics2_n_s$time <- "Historical"

metrics2_n_s$x=1

#

all <- bind_rows(metrics2_n_s, decomp_fut)


#



all$type_d <- factor(all$type, ordered = TRUE,
                     levels = c("Fixed effects", "Social drivers", "AC", "CDDs/HDDs","Expenditure", "Urbanisation", "Future social drivers change", "Future AC change", "Future CDDs/HDDs change", "Future urbanisation", "Future expenditure change"))


# 
# library("scales")
# show_col(pal_npg("nrc")(10))
# pal_npg("nrc")(10)

decomposition_plot <- ggplot(all)+
  geom_bar(aes(x=x, y=value, fill= forcats::fct_rev(type_d), group= forcats::fct_rev(type_d)), position="stack", stat="identity", colour="black") + scale_fill_discrete(name="Driver") + 
  geom_hline(yintercept = 1, size=2, linetype="dashed") + scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 2.5, 0.25)) + xlab("") + geom_text(aes(x=0.5, y=0.5, label="Historical"),angle=90) + geom_text(aes(x=0.5, y=1.15, label="Future"),angle=90) + ylab("Relative contribution") + ggtitle(paste0("Drivers of electricity consumption, current and 2050 (SSP2-RCP45)"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + guides(fill=guide_legend(ncol=1))

ggsave(paste0(stub, "results/graphs/", countryiso3, "_decompose_", ssp, ".png"), last_plot())

save(all, decomposition_plot, file=paste0(stub, "results/graphs/", countryiso3, "_decompose_", ssp, ".Rdata"))

}

###########################################

# 3.2) Electricity consumption projections
# 3.2.1) Predict consumption without AC

orig_data <- data_c_sp
orig_data$ely_q <- NULL

orig_data_bk <- orig_data

output <- list()

# loop for all ssps and time-steps

for (ssp in c("SSP2", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  output2 <- list()
  
  for (year in seq(2020, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$ac = as.factor(0)
    
    orig_data$ln_total_exp_usd_2011 = data_c_sp[,paste0("exp_cap_usd_", ssp, "_", (year))]
    
    orig_data$curr_CDD18_db  = data_c_sp[,paste0("mean_CDD_", year, "_", rcp, "_", tolower(ssp))] / 100 # 10-year average bins
    
    orig_data$curr_HDD18_db  = data_c_sp[,paste0("mean_HDD_", year, "_", rcp, "_", tolower(ssp))] / 100  # 10-year average bins
    
    orig_data$edu_head_2 <- as.factor(data_c_sp[,paste0("edu_", year, "_", ssp)])
    
    orig_data$age_head <- round(orig_data$age_head * (1 + data_c_sp[,paste0("age_", ssp, "_", year)]), 0)
    
    orig_data$country = data_c$country; orig_data$ownership_d = as.factor(orig_data$ownership_d)
    
    #orig_data$fan = as.factor(data_c_sp[,paste0("fan_", ssp, "_", (year))])
    
    orig_data$housing_index = as.factor(data_c_sp[,paste0("housing_index_lab_", year, "_s1")])
    
	orig_data$urban_sh = data_c_sp[,paste0("weighted_mean.URB_", ssp, "_", (year))]
    
    #
    projected <- predict(model3, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output[[as.character(ssp)]] <- output2
  
}

###

# bind results together
output <- unlist(output, recursive = FALSE)
output_noac <- as.data.frame(do.call("cbind", output))

# 3.2.2) Predict consumption with AC

orig_data <- data_c_sp
orig_data$ely_q <- NULL

orig_data_bk <- orig_data

output <- list()

# loop for all ssps and time-steps

for (ssp in c("SSP2", "SSP5")){
  
  print(ssp)
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  output2 <- list()
  
  for (year in seq(2020, 2050, 10)){
    
    orig_data <- orig_data_bk
    
    print(year)
    
    orig_data$ac = as.factor(ifelse(data_c_sp[,paste0(ssp, ".", (year))]>0.5, 1, 0))
    
    orig_data$ln_total_exp_usd_2011 = data_c_sp[,paste0("exp_cap_usd_", ssp, "_", (year))]
    
    orig_data$curr_CDD18_db  = data_c_sp[,paste0("mean_CDD_", year, "_", rcp, "_", tolower(ssp))] / 100 # 10-year average bins
    
    orig_data$curr_HDD18_db  = data_c_sp[,paste0("mean_HDD_", year, "_", rcp, "_", tolower(ssp))] / 100  # 10-year average bins
    
    orig_data$edu_head_2 <- as.factor(data_c_sp[,paste0("edu_", year, "_", ssp)])
    
    orig_data$age_head <- round(orig_data$age_head * (1 + data_c_sp[,paste0("age_", ssp, "_", year)]), 0)
    
    orig_data$country = data_c$country; orig_data$ownership_d = as.factor(orig_data$ownership_d)
    
    #orig_data$fan = as.factor(data_c_sp[,paste0("fan_", ssp, "_", (year))])
    
    orig_data$housing_index = as.factor(data_c_sp[,paste0("housing_index_lab_", year, "_s1")])
    
	orig_data$urban_sh = data_c_sp[,paste0("weighted_mean.URB_", ssp, "_", (year))]
    
    #
    projected <- predict(model3, orig_data)
    
    output2[[as.character(year)]] <- projected
    
  }
  
  output[[as.character(ssp)]] <- output2
  
}


###

# bind results together
output <- unlist(output, recursive = FALSE)
output_ac <- as.data.frame(do.call("cbind", output))

#####
# 3.2.3) Estimate impact of AC ownership on electricity consumption

output_impact_ac <- exp(output_ac) - exp(output_noac)

output_impact_ac[output_impact_ac<=0] <- NA

save(output_ac, output_impact_ac, file = paste0(stub, "results/energy_poverty/", countryiso3, ".Rdata"))


##

for (ssp in c("SSP2", "SSP5")){
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  for (year in seq(2020, 2050, 10)){
    
    output_impact_ac[,paste0(ssp, ".", year)] =  output_impact_ac[,paste0(ssp, ".", year)] * (data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))] / mean(data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))][!is.na(output_impact_ac[,paste0(ssp, ".", year)])], na.rm=T))
    
  }}

national_summary_cons <- output_impact_ac %>%
  dplyr::summarise_all(mean, na.rm=T) %>%
  pivot_longer(cols = 1:8, names_to = c('ssp', 'year'), names_sep = ".") %>%
  mutate(ssp=rep(c("SSP245","SSP585"), each=4), year=rep(seq(2020, 2050, 10), 2))
  
line_country_ely_q <- ggplot(national_summary_cons)+
  geom_line(aes(x=year, y=value, colour=ssp, group=ssp))+
  xlab("Year")+
  ylab("Mean additional electricity consumption due to AC \n(kWh/hh/year)")+
  scale_colour_npg(name="Scenario")

ggsave(paste0(stub, "results/graphs/line_country_ely_q_", countryiso3, ".png"), line_country_ely_q)

#

# output_share_ac <- output_impact_ac / output_ac
#
# national_summary_cons <- output_share_ac %>%
#   dplyr::summarise_all("mean", na.rm=T) %>%
#   pivot_longer(cols = 1:20, names_to = c('ssp', 'year'), names_sep = ".") %>%
#   mutate(ssp=rep(c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5"), each=4), year=rep(seq(2020, 2050, 10), 5))
#
# line_country_ely_q_share <- ggplot(national_summary_cons)+
#   geom_line(aes(x=year, y=value, colour=ssp, group=ssp))+
#   xlab("Year")+
#   ylab("Average % electricity consumption due to AC")+
#   scale_colour_npg()
#
# ggsave(paste0(stub, "results/graphs/line_country_ely_q_share_", countryiso3, ".png"), line_country_ely_q_share, scale=1.5)

#

national_summary_total <- output_impact_ac %>%
  dplyr::summarise_all("mean", na.rm=T) %>%
  pivot_longer(cols = 1:8, names_to = c('ssp', 'year'), names_sep = ".") %>%
  mutate(ssp=rep(c("SSP245","SSP585"), each=4), year=rep(seq(2020, 2050, 10), 2))

#

# elaborate to get total consumption

pop_long <- custom_shape %>% dplyr::select(starts_with("sum.POP")) %>% dplyr::summarise_all(., "sum") %>% pivot_longer(1:50)
pop_long$name <- gsub("sum.POP_", "", pop_long$name )
pop_long$ssp <- substr(pop_long$name, 1, 4)
pop_long$year <- substr(pop_long$name, 6, 9)
pop_long <- filter(pop_long, year>2010 & year<=2050)

hhsize <- 3


pop_long$value <- pop_long$value / weighted.mean(data_c$n_members, data_c$weight) # https://globaldatalab.org/areadata/hhsize

pop_long <- filter(pop_long, grepl("SSP2", name) | grepl("SSP5", name))

# multiply by AC ownership
pop_long$value =  pop_long$value * national_summary_ac$value

# multiply by average consumption due to AC
national_summary_cons$value_tot <- pop_long$value * national_summary_cons$value

line_country_ely_q_total <- ggplot(national_summary_cons)+
  geom_line(aes(x=year, y=value_tot/1e9, colour=ssp, group=ssp))+
  xlab("Year")+
  ylab("Total electricity consumption due to AC (TWh), \ninclusive of population growth")+
  scale_colour_npg(name="Scenario")

ggsave(paste0(stub, "results/graphs/line_country_ely_q_tot_", countryiso3, ".png"), line_country_ely_q_total)

#

output_impact_ac$ely_p = data_c$ely_p_usd_2011

output_impact_ac2 <- output_impact_ac %>%
  pivot_longer(cols = 1:8, names_to = c('ssp', 'year'), names_sep = ".") %>%
  mutate(ssp=rep(rep(c("SSP245","SSP585"), each=4), nrow(output_impact_ac)), year=rep(rep(seq(2020, 2050, 10), 2), nrow(output_impact_ac)))

output_impact_ac2 <- output_impact_ac2 %>% group_by(ssp, year) %>% mutate(mean=median(value, na.rm=T))

distr_country_ely_q <- ggplot(output_impact_ac2, aes(x = value, y=..density.., group=year, colour=as.factor(year))) +   geom_density(adjust=.5)+ facet_wrap(.~ssp, nrow=2)+
  xlim(c(0, 2500))+
  geom_vline(aes(xintercept = mean, group=year, colour=as.factor(year)), linetype="dashed", size=.5)+
  ylab("Density")+
  xlab("Electricity consumption due to AC, distribution among households")+
  scale_colour_npg(name="Year")+
  labs(caption = "Solid line: density curve; dashed line: median value")

ggsave(paste0(stub, "results/graphs/distr_country_ely_q_", countryiso3, ".png"), distr_country_ely_q)

save(output_impact_ac2, distr_country_ely_q, file=paste0(stub, "results/graphs/", countryiso3, "_distribution.Rdata"))

###############
# export projections data

write.csv(national_summary_ac, paste0(stub, "results/household_level/", countryiso3, "_national_ac_penetration_glomod.csv"))
write.csv(national_summary_cons,  paste0(stub, "results/household_level/", countryiso3, "_national_ac_consumption_glomod.csv"))
write.csv(national_summary_total,  paste0(stub, "results/household_level/", countryiso3, "_national_ac_consumption_total_glomod.csv"))

regional_summary_ac$geometry<-NULL
write.csv(as.data.frame(regional_summary_ac),  paste0(stub, "results/household_level/", countryiso3, "_regional_ac_penetration_glomod.csv"))
# write.csv(regional_summary_cons,  paste0(stub, "results/household_level/", countryiso3, "_regional_consumption_glomod.csv"))
# write.csv(regional_summary_total,  paste0(stub, "results/household_level/", countryiso3, "_regional_consumption_total_glomod.csv"))

###

for (ssp in c("SSP2", "SSP5")){
  
  rcp <- ifelse(ssp=="SSP2", "rcp45", "rcp85")
  
  for (year in seq(2020, 2050, 10)){
    
    output_ac[,paste0(ssp, ".", year)] =  exp(output_ac[,paste0(ssp, ".", year)]) * (data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))] / mean(data_c_sp[,paste0("weight_", year, "_", rcp, "_", tolower(ssp))][!is.na(output_ac[,paste0(ssp, ".", year)])], na.rm=T))
    
  }}

national_summary_cons_totelectr <- output_ac %>%
  dplyr::summarise_all(mean, na.rm=T) %>%
  pivot_longer(cols = 1:8, names_to = c('ssp', 'year'), names_sep = ".") %>%
  mutate(ssp=rep(c("SSP245","SSP585"), each=4), year=rep(seq(2020, 2050, 10), 2))

# elaborate to get total consumption

pop_long <- custom_shape %>% dplyr::select(starts_with("sum.POP")) %>% dplyr::summarise_all(., "sum") %>% pivot_longer(1:50)
pop_long$name <- gsub("sum.POP_", "", pop_long$name )
pop_long$ssp <- substr(pop_long$name, 1, 4)
pop_long$year <- substr(pop_long$name, 6, 9)
pop_long <- filter(pop_long, year>2010 & year<=2050)

pop_long$value <- pop_long$value / weighted.mean(data_c$n_members, data_c$weight) # https://globaldatalab.org/areadata/hhsize

pop_long <- filter(pop_long, grepl("SSP2", name) | grepl("SSP5", name))

# multiply by average consumption due to AC
national_summary_cons_totelectr$value_tot <- pop_long$value * national_summary_cons_totelectr$value

write.csv(national_summary_cons_totelectr,  paste0(stub, "results/household_level/", countryiso3, "_national_ely_consumption_glomod.csv"))