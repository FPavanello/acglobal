
## This R-script:
##      1) Figure 1: Boxplot of social and macro-drivers of air-conditoning adoption and electricity quantity

# Free memory
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
library(gg.layers)
library(tidyverse)
library(patchwork)
#remotes::install_github('rpkgs/gg.layers')
library(ggsci)

# Set users
user <- 'fp'
user <- 'gf'

if (user=='fp') {
  stub <- 'F:/Il mio Drive/'
}

if (user=='gf') {
  stub <- 'F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

data <- paste(stub,'6-Projections/results/regressions/for_graphs/standardised', sep='')
output <- paste(stub,'6-Projections/results/graphs/', sep='')

list_p <- list.files(path="F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/regressions/for_graphs/standardised", full.names = T, pattern = "RData")

merger <- function(X){
  
  load(list_p[[X]])
  
    ac_margins$factor <- ifelse(ac_margins$factor=="ownership_d1", "ownership_d", ac_margins$factor)
    ac_margins$factor <- ifelse(ac_margins$factor=="sex_head1", "sex_head", ac_margins$factor)
        ac_margins$factor <- ifelse(ac_margins$factor=="sex_head2", "sex_head", ac_margins$factor)
    ac_margins$country <- toupper(unlist(qdapRegex::ex_between(list_p[X], "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/regressions/for_graphs/standardised/", "_dmcf.RData")))
    ac_margins$country <- toupper(unlist(qdapRegex::ex_between(list_p[X], "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/regressions/for_graphs/standardised/", "_dmcf.RData")))
    
  ac_margins <- dplyr::filter(ac_margins, !grepl("state|country|region|macro", factor))
  
  ac_margins$driver_type <- ifelse(grepl("CDD|HDD|exp|elyp", ac_margins$factor, ignore.case = T), "Income and climate drivers", "Socio-economic drivers")
  
  ac_margins$glob <- as.factor(ifelse(ac_margins$country=="GLOBAL", 1, 0))
  
  return(ac_margins)
  
}
  
pp <- lapply(1:length(list_p), merger)

pp <- as.data.frame(do.call(rbind, pp))

pp <- filter(pp, factor!="std_year_postschool")

pp$sign <- ifelse(pp$p<0.05, 1, 0)
 
pp <- filter(pp, factor!="std_elyp")

p3<- ggplot(pp %>% filter(sign==1)) +
  theme_classic()+
  geom_hline(yintercept = 0, colour="black", linetype="dashed")+
  geom_boxplot2(aes(x= factor, y=AME, fill=driver_type))+
  scale_fill_manual(name="Driver type", values = c("#4DBBD5B2", "#00A087B2")) + 
  geom_point(aes(x= factor, y=AME, shape=glob), size=2.5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")+
  ylab("Average marginal effect on AC \nadoption probability (%)")+
  scale_y_continuous(labels = scales::label_percent())+
  xlab("Drivers of AC adoption")+
  scale_shape_manual(name="Type", values = c(1, 17))+
  scale_x_discrete(labels = c("Primary education", "Secondary education", "Higer education", "Housing index - 2", "Housing index - 3", "Home ownership", "Female HH head", "Age HH head", "Contempor. CDDs", "Historical CDDs", "Contempor. HDDs", "HH size",  "Total HH expenditure", "Urbanisation"))
  

#ggsave("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/graphs/p2.png", p3, scale=2.5, width = 7)

###

# redo for ely

list_p <- list.files(path="F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/regressions/for_graphs/standardised", full.names = T, pattern = "RData")

merger <- function(X){
  
  load(list_p[[X]])

  ely_margins$factor <- ifelse(ely_margins$factor=="sex_head1", "sex_head", ely_margins$factor)
    ely_margins$factor <- ifelse(ely_margins$factor=="ac1", "ac", ely_margins$factor)
  ely_margins$factor <- ifelse(ely_margins$factor=="ownership_d1", "ownership_d", ely_margins$factor)
  ely_margins$factor <- ifelse(ely_margins$factor=="sex_head2", "sex_head", ely_margins$factor)
  ely_margins$country <- toupper(unlist(qdapRegex::ex_between(list_p[X], "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/regressions/for_graphs/standardised/", "_dmcf.RData")))
  
  ely_margins <- dplyr::filter(ely_margins, !grepl("state|country|region|macro|selection", factor))
  
  ely_margins$driver_type <- ifelse(grepl("CDD|exp|HDD|elyp", ely_margins$factor, ignore.case = T), "Income and climate drivers", ifelse(grepl("^ac", ely_margins$factor), "AC", "Socio-economic drivers"))
  
  ely_margins$glob <- as.factor(ifelse(ely_margins$country=="GLOBAL", 1, 0))
  
  return(ely_margins)
  
}

pp <- lapply(1:length(list_p), merger)

pp <- as.data.frame(do.call(rbind, pp))

#remotes::install_github('rpkgs/gg.layers')

# p4<- ggplot(pp, aes(x= factor, y=AME, fill=driver_type)) + 
#   geom_hline(yintercept = 0, colour="black")+
#   geom_violin()+
#   ggsci::scale_fill_npg(name="Driver type") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggsave("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/graphs/p3.png", p4, scale=2.5, width = 7)

pp$sign <- ifelse(pp$p<0.05, 1, 0)

#View(pp %>% filter(sign==1))

p5<- ggplot(pp %>% filter(sign==1)) +
  theme_classic()+
  geom_hline(yintercept = 0, colour="black")+
  geom_boxplot2(aes(x= factor, y=AME, fill=driver_type))+
  scale_fill_manual(name="Driver type", values = c("#bf0000", "#4DBBD5B2", "#00A087B2")) + 
  geom_point(aes(x= factor, y=AME, shape=glob), size=2.5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Average marginal effect on \nelectricity consumption (%)")+
  scale_y_continuous(limits = c(-0.4, 0.5), labels = scales::label_percent())+
  xlab("Drivers of electricity consumption")+
  scale_shape_manual(name="Models", values = c(1, 17), labels=c("Countries", "Global"))+
scale_x_discrete(labels = c("AC ownership", "Primary education", "Secondary education", "Higer education", "Housing index - 2", "Housing index - 3", "Home ownership", "Female HH head", "Age HH head", "Yearly CDDs", "Electricity price", "Yearly HDDs", "HH size", "Total HH expenditure", "Urbanisation"))

#ggsave("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/graphs/p4.png", p5, scale=2.5, width = 7)

###

p_final <- p3 + p5 + plot_layout(guides = "collect", ncol=1) + plot_annotation(tag_levels = 'A')

ggsave("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/graphs/fig1.png", last_plot(), scale=2.5, width = 3.5, height = 4)

p5

ggsave("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/graphs/fig1_b_only.png", last_plot(), scale=2.5, width = 3.5, height = 4/2)


#########


## This R-script:
##      1) Figure 1: Boxplot of social and macro-drivers of air-conditoning adoption and electricity quantity

# Free memory
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Set users
user <- 'fp'
user <- 'gf'

if (user=='fp') {
  stub <- 'F:/Il mio Drive/'
}

if (user=='gf') {
  stub <- 'F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

data <- paste(stub,'6-Projections/results/regressions/for_graphs/standardised', sep='')
output <- paste(stub,'6-Projections/results/graphs/', sep='')

list_p <- list.files(path="F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/regressions/for_graphs/standardised", full.names = T, pattern = "RData")[-6]

merger <- function(X){
  
  load(list_p[[X]])
  
  ac_margins$factor <- ifelse(ac_margins$factor=="ownership_d1", "ownership_d", ac_margins$factor)
  ac_margins$factor <- ifelse(ac_margins$factor=="sex_head1", "sex_head", ac_margins$factor)
  ac_margins$factor <- ifelse(ac_margins$factor=="sex_head2", "sex_head", ac_margins$factor)
  ac_margins$country <- toupper(unlist(qdapRegex::ex_between(list_p[X], "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/regressions/for_graphs/standardised/", "_dmcf.RData")))
  ac_margins$country <- toupper(unlist(qdapRegex::ex_between(list_p[X], "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/regressions/for_graphs/standardised/", "_dmcf.RData")))
  
  ac_margins <- dplyr::filter(ac_margins, !grepl("state|country|region|macro", factor))
  
  ac_margins$driver_type <- ifelse(grepl("CDD|HDD|exp|elyp", ac_margins$factor, ignore.case = T), "Income and climate drivers", "Socio-economic drivers")
  
  ac_margins$glob <- as.factor(ifelse(ac_margins$country=="GLOBAL", 1, 0))
  
  ac_margins$oecd = ifelse(ac_margins$country %in% toupper(c("deu", "ita", "oecdeu", "oecdnoneu", "usa", "mex")), "OECD", "Non-OECD")
  
  return(ac_margins)
  
}

pp <- lapply(1:length(list_p), merger)

pp <- as.data.frame(do.call(rbind, pp))

pp <- filter(pp, factor!="std_year_postschool")

pp$sign <- ifelse(pp$p<0.05, 1, 0)

pp <- filter(pp, factor!="std_elyp")

p3<- ggplot(pp %>% filter(sign==1)) +
  theme_classic()+
  facet_wrap(vars(as.factor(oecd)), ncol=1)+
  geom_hline(yintercept = 0, colour="black", linetype="dashed")+
  geom_boxplot2(aes(x= factor, y=AME, fill=driver_type))+
  scale_fill_manual(name="Driver type", values = c("#4DBBD5B2", "#00A087B2")) + 
  geom_point(aes(x= factor, y=AME), size=1.25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")+
  ylab("Average marginal effect on AC \nadoption probability (%)")+
  scale_y_continuous(labels = scales::label_percent())+
  xlab("Drivers of AC adoption")+
  scale_x_discrete(labels = c("Primary education", "Secondary education", "Higer education", "Housing index - 2", "Housing index - 3", "Home ownership", "Female HH head", "Age HH head", "Contempor. CDDs", "Historical CDDs", "Contempor. HDDs", "HH size",  "Total HH expenditure", "Urbanisation"))


#ggsave("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/graphs/p2.png", p3, scale=2.5, width = 7)

###

# redo for ely

list_p <- list.files(path="F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/regressions/for_graphs/standardised", full.names = T, pattern = "RData")[-6]

merger <- function(X){
  
  load(list_p[[X]])
  
  ely_margins$factor <- ifelse(ely_margins$factor=="sex_head1", "sex_head", ely_margins$factor)
  ely_margins$factor <- ifelse(ely_margins$factor=="ac1", "ac", ely_margins$factor)
  ely_margins$factor <- ifelse(ely_margins$factor=="ownership_d1", "ownership_d", ely_margins$factor)
  ely_margins$factor <- ifelse(ely_margins$factor=="sex_head2", "sex_head", ely_margins$factor)
  ely_margins$country <- toupper(unlist(qdapRegex::ex_between(list_p[X], "F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/regressions/for_graphs/standardised/", "_dmcf.RData")))
  
  ely_margins <- dplyr::filter(ely_margins, !grepl("state|country|region|macro|selection", factor))
  
  ely_margins$driver_type <- ifelse(grepl("CDD|exp|HDD|elyp", ely_margins$factor, ignore.case = T), "Income and climate drivers", ifelse(grepl("^ac", ely_margins$factor), "AC", "Socio-economic drivers"))
  
  ely_margins$oecd = ifelse(ely_margins$country %in% toupper(c("deu", "ita", "oecdeu", "oecdnoneu", "usa", "mex")), "OECD", "Non-OECD")
  
  #ely_margins$glob <- as.factor(ifelse(ely_margins$country=="GLOBAL", 1, 0))
  
  return(ely_margins)
  
}

pp <- lapply(1:length(list_p), merger)

pp <- as.data.frame(do.call(rbind, pp))

#remotes::install_github('rpkgs/gg.layers')

# p4<- ggplot(pp, aes(x= factor, y=AME, fill=driver_type)) + 
#   geom_hline(yintercept = 0, colour="black")+
#   geom_violin()+
#   ggsci::scale_fill_npg(name="Driver type") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggsave("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/graphs/p3.png", p4, scale=2.5, width = 7)

pp$sign <- ifelse(pp$p<0.05, 1, 0)

#View(pp %>% filter(sign==1))

p5<- ggplot(pp %>% filter(sign==1)) +
  theme_classic()+
  facet_wrap(vars(as.factor(oecd)), ncol=1)+
  geom_hline(yintercept = 0, colour="black")+
  geom_boxplot2(aes(x= factor, y=AME, fill=driver_type))+
  scale_fill_manual(name="Driver type", values = c("#bf0000", "#4DBBD5B2", "#00A087B2")) + 
  geom_point(aes(x= factor, y=AME), size=1.25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Average marginal effect on \nelectricity consumption (%)")+
  scale_y_continuous(limits = c(-0.4, 0.5), labels = scales::label_percent())+
  xlab("Drivers of electricity consumption")+
  scale_shape_manual(name="Models", values = c(1, 17), labels=c("Countries", "Global"))+
  scale_x_discrete(labels = c("AC ownership", "Primary education", "Secondary education", "Higer education", "Housing index - 2", "Housing index - 3", "Home ownership", "Female HH head", "Age HH head", "Yearly CDDs", "Electricity price", "Yearly HDDs", "HH size", "Total HH expenditure", "Urbanisation"))

p_final <- p3 + p5 + plot_layout(guides = "collect", ncol=1) + plot_annotation(tag_levels = 'A')

p5

ggsave("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/results/graphs/fig1_oecd.png", last_plot(), scale=2.5, width = 3.5, height = 4/2)
