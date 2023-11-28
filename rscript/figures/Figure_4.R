
## This R-script:
##      1) plots AC elasticities by country - coefficient plot (Figure 3)

# Free memory
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
library(tidyverse)
library(data.table)
library(patchwork)
library(ggsci)
library(ggpmisc)
library(ggtext)

# Set users
user <- 'fp'
#user <- 'gf'

if (user=='fp') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

if (user=='gf') {
  stub <- 'F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

house <- paste(stub,'6-Projections/data/household/', sep='')
results <- paste(stub,'6-Projections/results/regressions/', sep='')
graphs <- paste(stub,'6-Projections/results/graphs/Paper1/', sep='')


#################################

#    PLOT OF MFXs BY COUNTRY    #

#################################

## Load regression environments where marginal effects are saved
# Global
#load(paste(results, 'for_projections/global_dmcf.RData', sep=''))
#rm(global, model1, model2, reg_ac)
#glo <- dydx_ac %>% mutate(country = "Global")
#rm(dydx_ac)

# Africa
load(paste(results, 'for_projections/afr_dmcf.RData', sep=''))
rm(HH_Africa, model1, model2, model3, reg_ac)
afr <- dydx_ac %>% mutate(country = "Africa")
rm(dydx_ac)

# Argentina
load(paste(results, 'for_projections/arg_dmcf.RData', sep=''))
rm(HH_Argentina, model1, model2, model3, reg_ac)
arg <- dydx_ac %>% mutate(country = "Argentina")
rm(dydx_ac)

# Brazil
load(paste(results, 'for_projections/bra_dmcf.RData', sep=''))
rm(HH_Brazil, model1, model2, model3, reg_ac)
bra <- dydx_ac %>% mutate(country = "Brazil")
rm(dydx_ac)

# China
load(paste(results, 'for_projections/chn_dmcf.RData', sep=''))
rm(HH_China, model1, model2, model3, reg_ac)
chn <- dydx_ac %>% mutate(country = "China")
rm(dydx_ac)

# Germany
load(paste(results, 'for_projections/deu_dmcf.RData', sep=''))
rm(HH_Germany, model1, model2, model3, reg_ac)
deu <- dydx_ac %>% mutate(country = "Germany")
rm(dydx_ac)

# Indonesia
load(paste(results, 'for_projections/idn_dmcf.RData', sep=''))
rm(HH_Indonesia, model1, model2, model3, reg_ac)
idn <- dydx_ac %>% mutate(country = "Indonesia")
rm(dydx_ac)

# India
load(paste(results, 'for_projections/ind_dmcf.RData', sep=''))
rm(HH_India, model1, model2, model3, reg_ac)
ind <- dydx_ac %>% mutate(country = "India")
rm(dydx_ac)

# Italy
load(paste(results, 'for_projections/ita_dmcf.RData', sep=''))
rm(HH_Italy, model1, model2, reg_ac)
ita <- dydx_ac %>% mutate(country = "Italy")
rm(dydx_ac)

# Mexico
load(paste(results, 'for_projections/mex_dmcf.RData', sep=''))
rm(HH_Mexico, model1, model2, model3, reg_ac)
mex <- dydx_ac %>% mutate(country = "Mexico")
rm(dydx_ac)

# OECD-EU
load(paste(results, 'for_projections/oecdeu_dmcf.RData', sep=''))
rm(HH_Europe, model1, model2, model3, reg_ac_eu)
oeu <- dydx_ac %>% mutate(country = "OECD-EU")
rm(dydx_ac)

# OECD-NonEU
load(paste(results, 'for_projections/oecdnoneu_dmcf.RData', sep=''))
rm(HH_NonEurope, model5, model6, model7, reg_ac_noneu)
neu <- dydx_ac %>% mutate(country = "OECD-NonEU")
rm(dydx_ac)

# Pakistan
load(paste(results, 'for_projections/pak_dmcf.RData', sep=''))
rm(HH_Pakistan, model1, model2, model3, reg_ac)
pak <- dydx_ac %>% mutate(country = "Pakistan")
rm(dydx_ac)

# USA
load(paste(results, 'for_projections/usa_dmcf.RData', sep=''))
rm(HH_USA, model1, model2, model3, reg_ac)
usa <- dydx_ac %>% mutate(country = "USA")
rm(dydx_ac)

# Rbind
all <- rbind(afr, arg, bra, chn, deu, idn, ind, ita, mex, neu, oeu, pak, usa)

# Significance
all$sign <- ifelse(all$`Pr(>|t|)` <0.05, 1, 0)

# Order country based on expenditure
global <- readRDS(paste(house,'global.rds', sep=''))

global$country2 <- as.character(global$country)
global$country2[global$country == "United States"] <- "USA"
global$country2[global$country == "Ghana"] <- "Africa"
global$country2[global$country == "Nigeria"] <- "Africa"
global$country2[global$country == "Tanzania"] <- "Africa"
global$country2[global$country == "Niger"] <- "Africa"
global$country2[global$country == "Malawi"] <- "Africa"
global$country2[global$country == "Burkina Faso"] <- "Africa"
global$country2[global$country == "Kenya"] <- "Africa"
global$country2[global$country == "Spain"] <- "OECD-EU"
global$country2[global$country == "France"] <- "OECD-EU"
global$country2[global$country == "Netherlands"] <- "OECD-EU"
global$country2[global$country == "Switzerland"] <- "OECD-EU"
global$country2[global$country == "Germany"] <- "OECD-EU"
global$country2[global$country == "Sweden"] <- "OECD-EU"
global$country2[global$country == "Japan"] <- "OECD-NonEU"
global$country2[global$country == "Australia"] <- "OECD-NonEU"
global$country2[global$country == "Canada"] <- "OECD-NonEU"
global <- global[complete.cases(global$weight), ]
global <- global[complete.cases(global$ln_total_exp_usd_2011), ]
global <- global %>% mutate(total_exp_usd_2011 = exp(ln_total_exp_usd_2011))

global <- global %>% 
  dplyr::group_by(country2) %>% 
  dplyr::summarise(mean = weighted.mean(total_exp_usd_2011/n_members, weight, na.rm = TRUE)) %>% # expenditure per capita
  dplyr::arrange(mean)

colnames(global)[1] <- "country"
order <- global$country

# Weight Germany inside OECD-EU
germany_pop <- 83000000
germany_ac <- 0.01

epic_eu_pop <- (67e6 + 17e6 + 47e6 + 10e6)
epic_eu_ac <- (0.1267 * 67e6 + 0.1302 * 17e6 + 0.5444 * 47e6  + 0.1688 * 10e6) / ( 67e6 + 17e6 + 47e6  +  10e6)

weight1 <- (germany_pop*germany_ac) 
weight2 <- (epic_eu_pop*epic_eu_ac) 
weight1_sh <- (weight1) / (weight1+weight2)
weight2_sh <- (weight2) / (weight1+weight2)

all[all$country=="OECD-EU",c(2:5)] <- all[all$country=="OECD-EU",c(2:5)]*weight2_sh + all[all$country=="Germany",c(2:5)]*weight1_sh
all <- filter(all, country!="Germany")
all[all$country=="OECD-EU",5] <- pnorm((all[all$country=="OECD-EU",2] / all[all$country=="OECD-EU",3]), lower.tail=FALSE)
all[all$country=="OECD-EU",7] <- ifelse(all[all$country=="OECD-EU",5] <0.05, 1, 0)

# Adjust confidence intervals
all <- all %>% mutate(CI_lb = Estimate - 1.96*`Std. Error`,
                      CI_ub = Estimate + 1.96*`Std. Error`)

# Factor order
all <- all %>% filter(country != "Global")
all$country <- factor(all$country, levels = c(order))

# Test LM
all <- all %>% dplyr::filter(Variable == "ac_tot")
reg <- lm(Estimate*100 ~ as.numeric(country), data = all); summary(reg) # p-value 0.01318

## Coefficient plot
ac_country <- all %>% filter(Variable == "ac_tot") %>%
              ggplot(aes(x= country, y=Estimate*100)) + 
              geom_point(position = position_dodge(0.4)) + 
              geom_smooth(aes(as.numeric(country), Estimate*100), method='lm', se=T) + 
              geom_errorbar(aes(ymin=CI_lb*100, ymax=CI_ub*100), # colour=as.factor(sign)
                            alpha=0.6, size=0.6, width=0.4, position = "dodge") +
              geom_hline(aes(yintercept=0), linetype="dashed", color = "black") + 
#              scale_colour_manual(name="95% significance", 
#                                  values = c("#bf0000", "#4DBBD5B2"), labels=c("No", "Yes")) +
              ylab("Additional Electricity Demand due to AC (Marginal Effect, %)") +
              xlab("Country") +
              theme_classic() +
              theme(panel.background=element_blank(),
                    panel.border=element_blank(),
                    panel.grid.major=element_line(color="lightgray"),
                    panel.grid.minor=element_blank(),
                    plot.title = element_text(size = 12, family = "Tahoma", hjust = 0.5),
                    axis.text.x = element_text(size = 8, family = "Tahoma"), 
                    axis.text.y = element_text(size = 12, family = "Tahoma"),
#                    axis.title.x = element_blank(),
                    axis.ticks=element_blank()) +
              annotate("text", x = "OECD-NonEU", family = "mono", y=90, label = "P-value = 0.013", fontface = "bold")

ac_country

# Save
ggsave(paste(graphs, 'Figure4.png', sep = ''), last_plot(), scale=2.5, width = 3.5, height = 2)
