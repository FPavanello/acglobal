
# national projections, national models, without social drivers
lapply(list.files("./rscript/8_projections/national_model/cont_urb/only_macro_drivers/", pattern = "script", recursive=F, full.names = T), source, echo=T)

# national projections, national models
lapply(list.files("./rscript/8_projections/national_model/cont_urb/", pattern = "_dmf", recursive=F, full.names = T)[-c(6)], source, echo=T)

# global projections, global model
lapply(list.files("./rscript/8_projections/global_model/cont_urb/", pattern = "script", recursive=F, full.names = T)[6], source, echo=T)

# # national projections, global model, without social drivers
lapply(list.files("./rscript/8_projections/global_model/cont_urb/only_macro_drivers", pattern = "script", recursive=F, full.names = T), source, echo=T)

# national projections, global model
lapply(list.files("./rscript/8_projections/global_model/cont_urb/", pattern = "script", recursive=F, full.names = T), source, echo=T)
