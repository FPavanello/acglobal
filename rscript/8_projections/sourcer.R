
# national projections, national models, without social drivers
lapply(list.files("C:/Users/falchetta/Downloads/acglobal/rscript/8_projections/national_model/only_macro_drivers/", pattern = "script", recursive=F, full.names = T), source, echo=T)

# national projections, national models
lapply(list.files("C:/Users/falchetta/Downloads/acglobal/rscript/8_projections/national_model/", pattern = "_dmf", recursive=F, full.names = T), source, echo=T)

# # global projections, global model, without social drivers
lapply(list.files("C:/Users/falchetta/Downloads/acglobal/rscript/8_projections/global_model/only_macro_drivers", pattern = "script", recursive=F, full.names = T), source, echo=T)

# global projections, global model
lapply(list.files("C:/Users/falchetta/Downloads/acglobal/rscript/8_projections/global_model/", pattern = "script", recursive=F, full.names = T), source, echo=T)

