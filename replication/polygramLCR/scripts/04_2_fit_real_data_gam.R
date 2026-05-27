library(polygramLCR)
library(mgcv)

args = list(
  data_file = file.path("data", "ATUS_employed_vs_not_employed.csv"),
  out_dir = file.path("..", "outputs", "real_data", "fits", "ilr_gam"),
  years = "2022,2023,2024",
  age_min = "18",
  scenario_id = "ATUS",
  gam_k_grid = "50,80,100,150,200"
)

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

########## Prepare ATUS data ##########

years = as.integer(strsplit(args$years, ",", fixed = TRUE)[[1]])
atus = read.csv(args$data_file, check.names = FALSE)
# Keep the same adult ATUS analysis sample across all real-data methods.
atus = atus[
  atus$TUYEAR %in% years &
    atus$TEAGE >= as.integer(args$age_min) &
    atus$total_3part == 1440 &
    atus$x1 > 0 &
    atus$x2 > 0 &
    atus$x3 > 0 &
    !is.na(atus$y_emp),
  ,
  drop = FALSE
]
ilr = comp_to_ilr(as.matrix(atus[, c("x1", "x2", "x3")]))
analysis_data = data.frame(
  TUCASEID = atus$TUCASEID,
  TUYEAR = atus$TUYEAR,
  TEAGE = atus$TEAGE,
  x1 = atus$x1,
  x2 = atus$x2,
  x3 = atus$x3,
  z1 = ilr[, 1],
  z2 = ilr[, 2],
  v1 = atus$v1,
  v2 = atus$v2,
  y = as.integer(atus$y_emp),
  stringsAsFactors = FALSE
)
gam_k_grid = as.numeric(strsplit(args$gam_k_grid, ",", fixed = TRUE)[[1]])

scenario_id = args$scenario_id
analysis = analysis_data

########## Fit real-data GAM model ##########

start_time = proc.time()[["elapsed"]]
k_grid = as.integer(gam_k_grid)

# REML gives stable smoothing-parameter fitting for the ATUS binary response.
fits = lapply(k_grid, function(k) {
  gam(y ~ s(z1, z2, k = k), data = analysis, family = binomial(), method = "REML")
})

aic = vapply(fits, AIC, numeric(1))

# Store the AIC-selected basis size together with its fitted probabilities.
best_index = which.min(aic)
fit = fits[[best_index]]
attr(fit, "selected_k") = k_grid[[best_index]]
attr(fit, "k_grid") = paste(k_grid, collapse = ",")
attr(fit, "selection") = "aic"
prob = as.numeric(predict(fit, newdata = analysis, type = "response"))
elapsed = proc.time()[["elapsed"]] - start_time

predictions = data.frame(
  scenario_id = scenario_id,
  method = "ilr_gam",
  x1 = analysis$x1,
  x2 = analysis$x2,
  x3 = analysis$x3,
  z1 = analysis$z1,
  z2 = analysis$z2,
  v1 = analysis$v1,
  v2 = analysis$v2,
  y = analysis$y,
  prob = prob,
  elapsed_seconds = elapsed,
  stringsAsFactors = FALSE
)
diagnostics = data.frame(
  scenario_id = scenario_id,
  method = "ilr_gam",
  selected_k = attr(fit, "selected_k"),
  k_grid = attr(fit, "k_grid"),
  selection = attr(fit, "selection"),
  stringsAsFactors = FALSE
)

########## Save fitted model outputs ##########

saveRDS(fit, file.path(args$out_dir, paste0(scenario_id, "_fit.rds")))
write.csv(predictions, file.path(args$out_dir, "predictions.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(args$out_dir, "diagnostics.csv"), row.names = FALSE)
message("Wrote real-data GAM fits to: ", args$out_dir)
