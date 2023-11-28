
## This R-script:
##      1) plots AC elasticities by income country-specific quintile - coefficient plot (Figure 2)

# Free memory
.rs.restartR()
rm(list=ls(all=TRUE)) # Removes all previously created variables
gc()                  # frees up memory resources

# Packages
library(tidyverse)
library(data.table)
library(patchwork)
library(ggsci)

# Set users
user <- 'fp'
#user <- 'gf'

if (user=='fp') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

if (user=='gf') {
  stub <- 'H:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

results <- paste(stub,'6-Projections/results/regressions/', sep='')
graphs <- paste(stub,'6-Projections/results/graphs/Paper1/', sep='')


########################################################

#    PLOT OF MFXs BY COUNTRY-SPECIFIC INCOME DECILES   #

########################################################

# Global - quintiles
load(paste(results, 'for_graphs/subsamples/global_csquintiles.RData', sep=''))
ac_eff_csqnt <- ac_eff
ac_eff_cscddqnt <- ac_cdd

# Significance
ac_eff_csqnt <- ac_eff_csqnt %>% mutate(CI_lb = AME - 1.96*SE,
                                        CI_ub = AME + 1.96*SE)

## Coefficients Plot - Average Effect
# Graph
ac_csqnt <- ac_eff_csqnt %>% 
            ggplot(aes(x= factor(qnt_inc, levels = c("1st", "2nd", "3rd", "4th", "5th")), y=AME*100)) + geom_point() + 
            geom_errorbar(aes(ymin=CI_lb*100, ymax=CI_ub*100), alpha=0.6, size=0.6, width=0.4) +
            geom_hline(aes(yintercept=33.6), linetype="dashed", color = "red") + 
            ylab("Additional Electricity Demand due to AC \n(Marginal Effect, %)") +
            xlab("Expenditure Quintile") +
            coord_cartesian(ylim = c(20, 50)) +
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
ggsave(paste(graphs, 'Figure_3A.png', sep = ''), last_plot(), scale=2.5, width = 3.5, height = 2)


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
              geom_hline(aes(yintercept=33.6), linetype="dashed", color = "red") + 
              ylab("Additional Electricity Demand due to AC \n(Marginal Effect, %)") +
              xlab("Cooling Degree Days") +
              coord_cartesian(ylim = c(-15, 75)) +
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
ggsave(paste(graphs, 'Figure_3B.png', sep = ''), last_plot(), scale=2.5, width = 3.5, height = 2)

# Combine plots
Figure2 <- ac_csqnt + ac_cddqnt + plot_layout(guides = "collect", ncol=1) + plot_annotation(tag_levels = 'A')
ggsave(paste(graphs, 'Figure3.png', sep = ''), last_plot(), scale=2.5, width = 4.5, height = 4)
