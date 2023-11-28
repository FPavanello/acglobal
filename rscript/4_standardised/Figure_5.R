
## This R-script:
##      1) plots cooling energy poverty share by quintile (Figure 5)

# Free memory
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

# Directory ------------------------------------------------------------------
if (user=='fp') {
  stub <- 'G:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

if (user=='gf') {
  stub <- 'H:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/'
}

house <- paste(stub,'6-Projections/data/household/', sep='')
results <- paste(stub,'6-Projections/results/regressions/', sep='')
graphs <- paste(stub,'6-Projections/results/graphs/Paper1/', sep='')

## Figure 5
# Load global data
global <- readRDS(paste(house,'global.rds', sep=''))

# Load main model results
load(paste(results, 'for_projections/global_wgt_dmcf.RData', sep=''))

# Save coefficients
global <- global %>% mutate(beta_ac = reg_ely$coefficients[1],
                            beta_acxcdd = reg_ely$coefficients[17],
                            beta_acxcddsq = reg_ely$coefficients[18]) 

# Compute electricity for cooling
global <- global %>% mutate(ely_q = exp(ln_ely_q),
                            pct = beta_ac + beta_acxcdd*curr_CDD18_db + 
                                  beta_acxcddsq*(curr_CDD18_db^2),
                            ely_q_ac = ely_q*pct)
global$pct[global$ac == 0] <- 0
global$pct[global$ely_q_ac == 0] <- 0
global$ely_q_ac[global$ac == 0] <- 0
global$ely_q_ac[global$ely_q_ac < 0] <- 0

# Quintile of total expenditure
setDT(global)[,qnt_inc := cut(total_exp_usd_2011, breaks = quantile(total_exp_usd_2011, probs = seq(0, 1, 0.2)),
                              labels = c(1, 2, 3, 4, 5), include.lowest = TRUE), by = country]

# Expenditure for cooling
global <- global %>% mutate(ely_ac_exp = ely_q_ac*ely_p_usd_2011) # in $2011

# Shares
global <- global %>% mutate(share_ac = ely_ac_exp/total_exp_usd_2011*100,
                            sh_ely = sh_ely*100)
global <- global %>% filter(sh_ely <= 100) %>% filter(share_ac <= 100)
summary(global$sh_ely)
summary(global$sh_ely[global$ac == 1])
summary(global$share_ac[global$ac == 1])

# Set geographical areas
global$country = as.character(global$country)

global$country[global$country=="Ghana" | global$country=="Niger" | global$country=="Nigeria" | 
                 global$country=="Tanzania"| global$country=="Malawi"| global$country=="Kenya" | global$country=="Burkina Faso"] = "AFRICA"
global$country[global$country=="Argentina"] = "ARG"
global$country[global$country=="Brazil"] = "BRA"
global$country[global$country=="China"] = "CHN"
global$country[global$country=="Indonesia"] = "IDN"
global$country[global$country=="India"] = "IND"
global$country[global$country=="Italy"] = "ITA"
global$country[global$country=="Mexico"] = "MEX"
global$country[global$country=="Pakistan"] = "PAK"
global$country[global$country=="United States"] = "USA"
global$country[global$country=="Australia" | global$country=="Japan" | global$country=="Canada"] = "EPIC-NONEU"
global$country[global$country=="France" | global$country=="Germany"| global$country=="Netherlands"| global$country=="Spain"| global$country=="Sweden" |
                  global$country=="Switzerland"] = "EPIC-EU"

# Split samples to create groups
global_tot <- global %>% dplyr::select(c(hhid, country, sh_ely, qnt_inc)) %>% 
  mutate(sh = sh_ely, cat = "Total Electricity Expenditure (All Households)") %>% dplyr::select(-c(sh_ely))
global_ac <- global %>% filter(ac == 1) %>% dplyr::select(c(hhid, country, sh_ely, qnt_inc)) %>% 
  mutate(sh = sh_ely, cat = "Total Electricity Expenditure (AC = 1)") %>% dplyr::select(-c(sh_ely))
global_ace <- global %>% filter(ac == 1) %>% dplyr::select(c(hhid, country, share_ac, qnt_inc)) %>% 
  mutate(sh = share_ac, cat = "Cooling Electricity") %>% dplyr::select(-c(share_ac))

# Append
data_graph <- rbind(global_tot, global_ac, global_ace)
data_graph$cat <- as.factor(data_graph$cat)

# filtering function - turns outliers into NAs to be removed
filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}
# Graph
Figure5 <-  data_graph %>% group_by(country, qnt_inc) %>% mutate(sh_clean = filter_lims(sh)) %>%
            ggplot(aes(x= factor(qnt_inc, levels = c(1, 2, 3, 4, 5)), y=sh_clean, fill = cat)) + 
            geom_boxplot(na.rm = TRUE, coef = 5) +
            facet_wrap(vars(country), scales = "free", ncol=3)+
            ylab("Proportion of household total expenditure for (AC) electricity consumption")+
            xlab("Total household expenditure quintile in the survey year")+
            scale_fill_npg(name="")+
                      theme_classic() +
                      theme(panel.background=element_blank(),
                            panel.border=element_blank(),
                            panel.grid.major=element_line(color="lightgray"),
                            panel.grid.minor=element_blank(),
                            plot.title = element_text(size = 12, family = "Tahoma", hjust = 0.5), 
                            text = element_text(size = 12, family = "Tahoma"), 
                            legend.position = "bottom", 
                            legend.direction = "horizontal", 
                            strip.background = element_blank(),
                            axis.ticks=element_blank(),
                            legend.title = element_blank())

ac_csqnt

# Save
ggsave(paste(graphs, 'Figure_5.png', sep = ''), last_plot(), scale=2.5, width = 3.5, height = 2)
