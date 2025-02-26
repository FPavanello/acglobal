
##########################################

#                 Figure 5

##########################################

# Free memory
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
#remotes::install_github('rpkgs/gg.layers')
library(gg.layers)
library(tidyverse)
library(patchwork)
library(ggsci)

# Set users
user <- 'user'

if (user=='user') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'6-Projections/repo/household/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/standardised/'
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/figures/'


## Figure 5
# Set
list_p <- list.files(path=paste(interm), full.names = T, pattern = "RData")

# Function
merger <- function(X){
  
  load(list_p[[X]])
  
  ely_margins$term <- ifelse(ely_margins$term=="edu_head_2" & ely_margins$contrast=="mean(3) - mean(0)", "edu_head_thigher", ely_margins$term)
  
  ely_margins$term <- ifelse(ely_margins$term=="edu_head_2" & ely_margins$contrast=="mean(2) - mean(0)", "edu_head_secondary", ely_margins$term)
  
  ely_margins$term <- ifelse(ely_margins$term=="edu_head_2" & ely_margins$contrast=="mean(1) - mean(0)", "edu_head_primary", ely_margins$term)
  
  ely_margins <- filter(ely_margins, term!="edu_head_2")
  
  ely_margins$term <- ifelse(ely_margins$term=="housing_index_lab" & ely_margins$contrast=="mean(3) - mean(1)", "housing_index_3", ely_margins$term)
  
  ely_margins$term <- ifelse(ely_margins$term=="housing_index_lab" & ely_margins$contrast=="mean(2) - mean(1)", "housing_index_2", ely_margins$term)
  
  ely_margins$term <- ifelse(ely_margins$term=="sex_head", "sex_head", ely_margins$term)
  ely_margins$term <- ifelse(ely_margins$term=="ac", "ac", ely_margins$term)
  ely_margins$term <- ifelse(ely_margins$term=="ownership_d", "ownership_d", ely_margins$term)
  
  
  ely_margins$country <- toupper(unlist(qdapRegex::ex_between(list_p[X], "C:/Users/Standard/Documents/Github/acglobal/interm/standardised/", "_dmcf.RData")))
  
  ely_margins <- dplyr::filter(ely_margins, !grepl("state|country|region|macro|selection", term))
  
  ely_margins$driver_type <- ifelse(grepl("CDD|exp|HDD|elyp", ely_margins$term, ignore.case = T), "Income and climate drivers", ifelse(grepl("^ac", ely_margins$term), "AC", "Socio-economic drivers"))
  
  ely_margins$glob <- as.factor(ifelse(ely_margins$country=="GLOBAL", 1, 0))
  
  return(ely_margins)
  
}

# Call
pp <- lapply(1:length(list_p), merger)

# Append
pp <- as.data.frame(do.call(rbind, pp))

# Significance
pp$sign <- ifelse(pp$p.value<0.05, 1, 0)

# Plot
p5<- ggplot(pp %>% filter(sign==1)) +
  theme_classic()+
  geom_hline(yintercept = 0, colour="black")+
  geom_boxplot2(aes(x= term, y=estimate, fill=driver_type))+
  scale_fill_manual(name="Driver type", values = c("#bf0000", "#4DBBD5B2", "#00A087B2")) + 
  geom_point(aes(x= term, y=estimate, shape=glob), size=2.5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Average marginal effect on \nelectricity consumption (%)")+
  scale_y_continuous(limits = c(-0.4, 0.5), labels = scales::label_percent())+
  xlab("Drivers of electricity consumption")+
  scale_shape_manual(name="Models", values = c(1, 17), labels=c("Countries", "Global"))+
  scale_x_discrete(labels = c("AC ownership", "Primary education", "Secondary education", "Higher education", 
                              "Housing index - 2", "Housing index - 3", "Home ownership", "Female HH head", 
                              "Age HH head", "Yearly CDDs", "Electricity price", "Yearly HDDs", "HH size", 
                              "Total HH expenditure", "Urbanisation"))

# Save
ggsave(paste(output, 'Figure5.png', sep = ''), last_plot(), scale=2.5, width = 3.5, height = 4/2)


## Figure A2
# Load
list_p <- list.files(path=paste(interm), full.names = T, pattern = "RData")

# Function
merger <- function(X){
  
  load(list_p[[X]])
  
  ely_margins$term <- ifelse(ely_margins$term=="edu_head_2" & ely_margins$contrast=="mean(3) - mean(0)", "edu_head_thigher", ely_margins$term)
  
  ely_margins$term <- ifelse(ely_margins$term=="edu_head_2" & ely_margins$contrast=="mean(2) - mean(0)", "edu_head_secondary", ely_margins$term)
  
  ely_margins$term <- ifelse(ely_margins$term=="edu_head_2" & ely_margins$contrast=="mean(1) - mean(0)", "edu_head_primary", ely_margins$term)
  
  ely_margins <- filter(ely_margins, term!="edu_head_2")
  
  ely_margins$term <- ifelse(ely_margins$term=="housing_index_lab" & ely_margins$contrast=="mean(3) - mean(1)", "housing_index_3", ely_margins$term)
  
  ely_margins$term <- ifelse(ely_margins$term=="housing_index_lab" & ely_margins$contrast=="mean(2) - mean(1)", "housing_index_2", ely_margins$term)
  
  ely_margins$term <- ifelse(ely_margins$term=="sex_head", "sex_head", ely_margins$term)
  ely_margins$term <- ifelse(ely_margins$term=="ac", "ac", ely_margins$term)
  ely_margins$term <- ifelse(ely_margins$term=="ownership_d", "ownership_d", ely_margins$term)
  
  ely_margins$country <- toupper(unlist(qdapRegex::ex_between(list_p[X], paste(interm), "_dmcf.RData")))
  
  ely_margins <- dplyr::filter(ely_margins, !grepl("state|country|region|macro|selection", term))
  
  ely_margins$driver_type <- ifelse(grepl("CDD|exp|HDD|elyp", ely_margins$term, ignore.case = T), "Income and climate drivers", ifelse(grepl("^ac", ely_margins$term), "AC", "Socio-economic drivers"))
  
  ely_margins$oecd = ifelse(ely_margins$country %in% toupper(c("deu", "ita", "oecdeu", "oecdnoneu", "usa", "mex")), "OECD", "Non-OECD")
  
  return(ely_margins)
  
}

# Call
pp <- lapply(1:length(list_p), merger)

# Append
pp <- as.data.frame(do.call(rbind, pp))

# Significance
pp$sign <- ifelse(pp$p.value<0.05, 1, 0)

# Plot
p5<- ggplot(pp %>% filter(sign==1)) +
  theme_classic()+
  facet_wrap(vars(as.factor(oecd)), ncol=1)+
  geom_hline(yintercept = 0, colour="black")+
  geom_boxplot2(aes(x= term, y=estimate, fill=driver_type))+
  scale_fill_manual(name="Driver type", values = c("#bf0000", "#4DBBD5B2", "#00A087B2")) + 
  geom_point(aes(x= term, y=estimate), size=1.25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Average marginal effect on \nelectricity consumption (%)")+
  scale_y_continuous(limits = c(-0.4, 0.5), labels = scales::label_percent())+
  xlab("Drivers of electricity consumption")+
  scale_shape_manual(name="Models", values = c(1, 17), labels=c("Countries", "Global"))+
  scale_x_discrete(labels = c("AC ownership", "Primary education", "Secondary education", "Higher education", "Housing index - 2", "Housing index - 3", "Home ownership", "Female HH head", "Age HH head", "Yearly CDDs", "Electricity price", "Yearly HDDs", "HH size", "Total HH expenditure", "Urbanisation"))

p5

# Save
ggsave(paste(output, 'FigureA2.png', sep = ''), last_plot(), scale=2.5, width = 3.5, height = 4/2)
