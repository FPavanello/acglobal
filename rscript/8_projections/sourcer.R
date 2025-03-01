
### 
# Master script to replicate projections and figures/tables that are based on the projections
# Note: this script calls individual script files and it should be run exactly in this order to replicate the paper results

# national projections, national models, without social drivers
lapply(list.files("C:/Users/Utente/Downloads/acglobal/rscript/8_projections/national_model/only_macro_drivers/", pattern = "script", recursive=F, full.names = T), source, echo=T)

# national projections, national models
lapply(list.files("C:/Users/Utente/Downloads/acglobal/rscript/8_projections/national_model/", pattern = "_dmf", recursive=F, full.names = T), source, echo=T)

# replicate figures
source("C:/Users/Utente/Downloads/acglobal/rscript/8_projections/figures/Figure_A5.R")
source("C:/Users/Utente/Downloads/acglobal/rscript/7_implications/Figure_7.R")

#####

# # national projections, global model, without social drivers
lapply(list.files("C:/Users/Utente/Downloads/acglobal/rscript/8_projections/global_model/only_macro_drivers", pattern = "script", recursive=F, full.names = T), source, echo=T)

# # global projections, global model
lapply(list.files("C:/Users/Utente/Downloads/acglobal/rscript/8_projections/global_model/", pattern = "script", recursive=F, full.names = T)[6], source, echo=T)

# national projections, global model
lapply(list.files("C:/Users/Utente/Downloads/acglobal/rscript/8_projections/global_model/", pattern = "script", recursive=F, full.names = T)[-6], source, echo=T)

#####

# replicate figures and tables

source("C:/Users/Utente/Downloads/acglobal/rscript/7_implications/Table_4_A15.R")
source("C:/Users/Utente/Downloads/acglobal/rscript/8_projections/figures/Figure_A6.R")

source("C:/Users/Utente/Downloads/acglobal/rscript/9_supplementary/Figure_S2.R")
source("C:/Users/Utente/Downloads/acglobal/rscript/9_supplementary/Figure_S3.R")

source("C:/Users/Utente/Downloads/acglobal/rscript/8_projections/figures/Figure_A4.R")

