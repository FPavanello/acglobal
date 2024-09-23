
##########################################

#                 Figure 3

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
library(ggpubr)

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


# Load input
load(paste(interm, 'global_csquintiles.RData', sep=''))
ac_eff_csqnt <- ac_eff
ac_eff_cscddqnt <- ac_cdd


## Coefficients Plot - Average Effect
# Graph
ac_csqnt <- ac_eff_csqnt %>% 
  ggplot(aes(x= factor(qnt_inc, levels = c("1st", "2nd", "3rd", "4th", "5th")), y=estimate*100)) + geom_point() + 
  geom_errorbar(aes(ymin=conf.low*100, ymax=conf.high*100), alpha=0.6, size=0.6, width=0.4) +
  geom_hline(aes(yintercept=36.2), linetype="dashed", color = "red") + 
  ylab("Additional Electricity Demand due to AC \n(Marginal Effect, %)") +
  xlab("Expenditure Quintile") +
  coord_cartesian(ylim = c(20, 60)) +
  theme_classic() +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color="lightgray"),
        panel.grid.minor=element_blank(),
        plot.title = element_text(size = 12, family = "Tahoma", hjust = 0.5), 
        text = element_text(size = 12, family = "Tahoma"), 
        #                  axis.title.x = element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_blank())

ac_csqnt

# Save
ggsave(paste(output, 'Figure_3A.png', sep = ''), last_plot(), scale=2.5, width = 3.5, height = 2)


## Coefficients Plot - Effect varies with CDD
# Adjust
ac_eff_cscddqnt$cdd <- as.character(ac_eff_cscddqnt$cdd)
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd0"] <- "0"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd1"] <- "1"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd2"] <- "2"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd3"] <- "3"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd4"] <- "4"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd5"] <- "5"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd6"] <- "6"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd7"] <- "7"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd8"] <- "8"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd9"] <- "9"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd10"] <- "10"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd11"] <- "11"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd12"] <- "12"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd13"] <- "13"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd14"] <- "14"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd15"] <- "15"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd16"] <- "16"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd17"] <- "17"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd18"] <- "18"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd19"] <- "19"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd20"] <- "20"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd21"] <- "21"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd22"] <- "22"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd23"] <- "23"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd24"] <- "24"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd25"] <- "25"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd26"] <- "26"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd27"] <- "27"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd28"] <- "28"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd29"] <- "29"
ac_eff_cscddqnt$cdd[ac_eff_cscddqnt$cdd == "cdd30"] <- "30"
ac_eff_cscddqnt$cdd <- as.numeric(ac_eff_cscddqnt$cdd)

# Graph
ac_cddqnt  <- ac_eff_cscddqnt %>% 
  ggplot(aes(x= cdd*100, 
             y=AME*100)) + geom_point() + 
  geom_errorbar(aes(ymin=CI_lb*100, ymax=CI_ub*100), alpha=0.6, size=0.6, width=0.4) +
  geom_hline(aes(yintercept=36.2), linetype="dashed", color = "red") + 
  ylab("Additional Electricity Demand due to AC \n(Marginal Effect, %)") +
  xlab("Cooling Degree Days") +
  coord_cartesian(ylim = c(-30, 80)) +
  theme_classic() +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color="lightgray"),
        panel.grid.minor=element_blank(),
        plot.title = element_text(size = 12, family = "Tahoma", hjust = 0.5), 
        text = element_text(size = 12, family = "Tahoma"), 
        #                  axis.title.x = element_blank(),
        axis.ticks=element_blank(),
        legend.title = element_blank()) +
  facet_grid(~ qnt_inc)

ac_cddqnt

# Save
ggsave(paste(output, 'Figure_3B.png', sep = ''), last_plot(), scale=2.5, width = 3.5, height = 2)

# Combine plots
Figure3 <- ggarrange(ac_csqnt, ac_cddqnt, 
                     labels = c("A", "B"),
                     ncol = 1, nrow = 2)
Figure3
ggsave(paste(output, 'Figure3.png', sep = ''), last_plot(), scale=2.5, width = 4.5, height = 4)
