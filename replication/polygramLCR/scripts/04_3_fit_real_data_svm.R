library(polygramLCR)
library(e1071)

args = list(
  data_file = file.path("data", "ATUS_employed_vs_not_employed.csv"),
  out_dir = file.path("..", "outputs", "real_data", "fits", "ilr_svm_rbf"),
  years = "2022,2023,2024",
  age_min = "18",
  scenario_id = "ATUS",
  svm_cost_grid = "0.25,1,4,16",
  svm_gamma_grid = "0.5,2,4,8",
  svm_cv_folds = "3",
  seed = "20260509"
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

svm_cost_grid = as.numeric(strsplit(args$svm_cost_grid, ",", fixed = TRUE)[[1]])
svm_gamma_grid = as.numeric(strsplit(args$svm_gamma_grid, ",", fixed = TRUE)[[1]])

scenario_id = args$scenario_id
analysis = analysis_data

########## Fit real-data SVM model ##########

start_time = proc.time()[["elapsed"]]
set.seed(as.integer(args$seed))
# Cross-validation selects the radial-kernel cost and gamma for the ATUS fit.
tuned = tune.svm(
  y ~ .,
  data = data.frame(y = factor(analysis$y), z1 = analysis$z1, z2 = analysis$z2),
  kernel = "radial",
  cost = as.numeric(svm_cost_grid),
  gamma = as.numeric(svm_gamma_grid),
  probability = TRUE,
  tunecontrol = tune.control(cross = as.integer(args$svm_cv_folds))
)

fit = tuned$best.model
attr(fit, "selected_cost") = tuned$best.parameters$cost
attr(fit, "selected_gamma") = tuned$best.parameters$gamma
attr(fit, "cost_grid") = paste(svm_cost_grid, collapse = ",")
attr(fit, "gamma_grid") = paste(svm_gamma_grid, collapse = ",")
attr(fit, "selection") = paste0(args$svm_cv_folds, "-fold-cv")
pred = predict(
  fit,
  newdata = data.frame(z1 = analysis$z1, z2 = analysis$z2),
  probability = TRUE
)

prob = as.numeric(attr(pred, "probabilities")[, "1"])
elapsed = proc.time()[["elapsed"]] - start_time

predictions = data.frame(
  scenario_id = scenario_id,
  method = "ilr_svm_rbf",
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
  method = "ilr_svm_rbf",
  selected_cost = attr(fit, "selected_cost"),
  selected_gamma = attr(fit, "selected_gamma"),
  cost_grid = attr(fit, "cost_grid"),
  gamma_grid = attr(fit, "gamma_grid"),
  selection = attr(fit, "selection"),
  stringsAsFactors = FALSE
)

########## Save fitted model outputs ##########

saveRDS(fit, file.path(args$out_dir, paste0(scenario_id, "_fit.rds")))
write.csv(predictions, file.path(args$out_dir, "predictions.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(args$out_dir, "diagnostics.csv"), row.names = FALSE)
message("Wrote real-data SVM fits to: ", args$out_dir)
