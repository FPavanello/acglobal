
##########################################

#                 Figure 4

##########################################

# Free memory
.rs.restartR()
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
user <- 'user'

if (user=='user') {
  stub <- "G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/"
}

house <- paste(stub,'data/household/', sep='')
interm <- paste(stub,'results/regressions/for_graphs/subsamples/', sep='')
interm <- 'C:/Users/Standard/Documents/Github/acglobal/interm/'
output <- paste(stub,'output/figures/', sep='')
output <- 'C:/Users/Standard/Documents/Github/acglobal/output/figures/'


## Load
# Africa
load(paste(interm, 'afr_dmcf.RData', sep=''))
rm(HH_Africa, model3, reg_ac)
afr <- dydx_ac %>% mutate(country = "Africa")
rm(dydx_ac)

# Argentina
load(paste(interm, 'arg_dmcf.RData', sep=''))
rm(HH_Argentina, model3, reg_ac)
arg <- dydx_ac %>% mutate(country = "Argentina")
rm(dydx_ac)

# Brazil
load(paste(interm, 'bra_dmcf.RData', sep=''))
rm(HH_Brazil, model3, reg_ac)
bra <- dydx_ac %>% mutate(country = "Brazil")
rm(dydx_ac)

# China
load(paste(interm, 'chn_dmcf.RData', sep=''))
rm(HH_China, model3, reg_ac)
chn <- dydx_ac %>% mutate(country = "China")
rm(dydx_ac)

# Germany
load(paste(interm, 'deu_dmcf.RData', sep=''))
rm(HH_Germany, model3, reg_ac)
deu <- dydx_ac %>% mutate(country = "Germany")
rm(dydx_ac)

# Indonesia
load(paste(interm, 'idn_dmcf.RData', sep=''))
rm(HH_Indonesia, model3, reg_ac)
idn <- dydx_ac %>% mutate(country = "Indonesia")
rm(dydx_ac)

# India
load(paste(interm, 'ind_dmcf.RData', sep=''))
rm(HH_India, model3, reg_ac)
ind <- dydx_ac %>% mutate(country = "India")
rm(dydx_ac)

# Italy
load(paste(interm, 'ita_dmcf.RData', sep=''))
rm(HH_Italy, model3, reg_ac)
ita <- dydx_ac %>% mutate(country = "Italy")
rm(dydx_ac)

# Mexico
load(paste(interm, 'mex_dmcf.RData', sep=''))
rm(HH_Mexico, model3, reg_ac)
mex <- dydx_ac %>% mutate(country = "Mexico")
rm(dydx_ac)

# OECD-EU
load(paste(interm, 'oecdeu_dmcf.RData', sep=''))
rm(HH_Europe, model3, reg_ac_eu)
oeu <- dydx_ac %>% mutate(country = "OECD-EU")
rm(dydx_ac)

# OECD-NonEU
load(paste(interm, 'oecdnoneu_dmcf.RData', sep=''))
rm(HH_NonEurope, model7, reg_ac_noneu)
neu <- dydx_ac %>% mutate(country = "OECD-NonEU")
rm(dydx_ac)

# Pakistan
load(paste(interm, 'pak_dmcf.RData', sep=''))
rm(HH_Pakistan, model3, reg_ac)
pak <- dydx_ac %>% mutate(country = "Pakistan")
rm(dydx_ac)

# USA
load(paste(interm, 'usa_dmcf.RData', sep=''))
rm(HH_USA, model3, reg_ac)
usa <- dydx_ac %>% mutate(country = "USA")
rm(dydx_ac)

# Rbind
all <- rbind(afr, arg, bra, chn, deu, idn, ind, ita, mex, neu, oeu, pak, usa)

# Significance
all$sign <- ifelse(all$p.value <0.05, 1, 0)

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

# Check
global <- global[complete.cases(global$ln_ely_q), ]
global <- global[complete.cases(global$ac), ]
global <- global[complete.cases(global$ln_total_exp_usd_2011), ]
global <- global[complete.cases(global$mean_CDD18_db), ]
global <- global[complete.cases(global$ownership_d), ]
global <- global[complete.cases(global$n_members), ]
global <- global[complete.cases(global$age_head), ]
global <- global[complete.cases(global$country), ]
global <- global[complete.cases(global$weight), ]
global <- global[complete.cases(global$sex_head), ]
global <- global[complete.cases(global$urban_sh), ]
global <- global[complete.cases(global$ln_ely_p), ]
global <- global[complete.cases(global$curr_CDD18_db), ]
global <- global[complete.cases(global$curr_HDD18_db), ]
global <- global[complete.cases(global$adm1), ]
global <- global %>% filter(ln_ely_q > 0)
global <- global %>% filter(weight > 0)

# Create order by household expenditure per capita
exp <- global %>% 
  dplyr::group_by(country2) %>% 
  dplyr::summarise(mean = weighted.mean(total_exp_usd_2011/n_members, weight, na.rm = TRUE)) %>% # expenditure per capita
  dplyr::arrange(mean)

colnames(exp)[1] <- "country"
order <- exp$country

# Weight Germany inside OECD-EU
germany_pop <- 83000000
germany_ac <- 0.01

epic_eu_pop <- (67e6 + 17e6 + 47e6 + 10e6)
epic_eu_ac <- (0.1267 * 67e6 + 0.1302 * 17e6 + 0.5444 * 47e6  + 0.1688 * 10e6) / ( 67e6 + 17e6 + 47e6  +  10e6)

weight1 <- (germany_pop*germany_ac) 
weight2 <- (epic_eu_pop*epic_eu_ac) 
weight1_sh <- (weight1) / (weight1+weight2)
weight2_sh <- (weight2) / (weight1+weight2)

all[all$country=="OECD-EU",c(3:7)] <- all[all$country=="OECD-EU",c(3:7)]*weight2_sh + all[all$country=="Germany",c(3:7)]*weight1_sh
all <- filter(all, country!="Germany")
all[all$country=="OECD-EU",6] <- pnorm((all[all$country=="OECD-EU",3] / all[all$country=="OECD-EU",4]), lower.tail=FALSE)
all[all$country=="OECD-EU",10] <- ifelse(all[all$country=="OECD-EU",5] <0.05, 1, 0)

# Factor order
all <- all %>% filter(country != "Global")
all$country <- factor(all$country, levels = c(order))

# Test LM
reg <- lm(estimate*100 ~ as.numeric(country), data = all); summary(reg) # p-value 0.00335
reg <- summary(reg)
pval <- round(reg$coefficients[2,4], 3)

## Coefficient plot
ac_country <- all %>% 
  ggplot(aes(x= country, y=estimate*100)) + 
  geom_point(position = position_dodge(0.4)) + 
  geom_smooth(aes(as.numeric(country), estimate*100), method='lm', se=T) + 
  geom_errorbar(aes(ymin=conf.low*100, ymax=conf.high*100), # colour=as.factor(sign)
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
  annotate("text", x = "OECD-NonEU", family = "mono", y=90, label = paste0("P-value = ", pval, sep = ""), fontface = "bold")

ac_country

# Save
ggsave(paste(output, 'Figure4.png', sep = ''), last_plot(), scale=2.5, width = 3.5, height = 2)
