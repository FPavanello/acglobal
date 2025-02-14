#lapply(list.files("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/dmcf/regressions/country/with_continuous_urbanisation", pattern = "dmcf", recursive=F, full.names = T), source)

library(fixest)

# national projections, national models, without social drivers
lapply(list.files("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/projections/ac_ely/cont_urb/wsd/", pattern = "script", recursive=F, full.names = T), source, echo=T)

# national projections, national models
lapply(list.files("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/projections/ac_ely/cont_urb/", pattern = "_dmf", recursive=F, full.names = T)[-c(6)], source, echo=T)

source("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/figures_paper/fig2.R")

source("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/figures_paper/compare_all_drivers_only_income_climate.R")

# global projections, global model
lapply(list.files("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/projections/ac_ely/cont_urb/", pattern = "script", recursive=F, full.names = T)[6], source, echo=T)

# # national projections, global model, without social drivers
lapply(list.files("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/projections/ac_ely/cont_urb/glob_mod/wsd", pattern = "script", recursive=F, full.names = T), source, echo=T)

# national projections, global model
lapply(list.files("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/projections/ac_ely/cont_urb/glob_mod/", pattern = "script", recursive=F, full.names = T), source, echo=T)

source("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/figures_paper/compare_all_drivers_only_income_climate_glomod.R")


###

source("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/figures_paper/table_summary.R")

source("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/figures_paper/fig3_expenditure.R")

source("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/figures_paper/energy_poverty_figure.R")

# check bias size figures

source("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/projections/ac_ely/check_bias_size.R")

source("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/projections/ac_ely/check_bias_size_ely.R")

# comparison

source("F:/.shortcut-targets-by-id/1JhN0qxmpnYQDoWQdBhnYKzbRCVGH_WXE/6-Projections/rscripts/figures_paper/compare_nat_glob_mods.R")

