library(polygramLCR)

args = list(
  data_file = file.path("data", "ATUS_employed_vs_not_employed.csv"),
  out_dir = file.path("..", "outputs", "real_data", "fits", "polygram"),
  years = "2022,2023,2024",
  age_min = "18",
  scenario_id = "ATUS",
  polygram_centers = "100",
  polygram_admm_lambdas = "100",
  polygram_admm_lambda_max = "5",
  polygram_admm_lambda_max_min_ratio = "0.01",
  polygram_admm_selection = "aic",
  polygram_admm_rho = "auto",
  polygram_admm_max_iterations = "500",
  polygram_admm_max_newton_iterations = "50",
  polygram_admm_hessian_epsilon = "0",
  polygram_admm_require_newton_convergence = "true",
  polygram_admm_convergence_criterion = "objective",
  polygram_admm_objective_tol = "1e-5",
  polygram_admm_require_admm_convergence = "true",
  polygram_admm_normalize_constraint_rows = "true",
  polygram_admm_warm_start = "true",
  polygram_verbose = "false",
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

scenario_id = args$scenario_id
analysis = analysis_data

########## Fit real-data polygram model ##########

start_time = proc.time()[["elapsed"]]
# The real-data polygram fit is the source for both probabilities and active edges.
fit = polygram(
  analysis,
  centers = as.integer(args$polygram_centers),
  number_lambdas = as.integer(args$polygram_admm_lambdas),
  lambda_max = as.numeric(args$polygram_admm_lambda_max),
  lambda_max_min_ratio = as.numeric(args$polygram_admm_lambda_max_min_ratio),
  selection = args$polygram_admm_selection,
  rho = if (identical(args$polygram_admm_rho, "auto")) NULL else as.numeric(args$polygram_admm_rho),
  max_admm_iter = as.integer(args$polygram_admm_max_iterations),
  max_newton_iter = as.integer(args$polygram_admm_max_newton_iterations),
  hessian_epsilon = as.numeric(args$polygram_admm_hessian_epsilon),
  require_newton_convergence = tolower(args$polygram_admm_require_newton_convergence) == "true",
  convergence_criterion = args$polygram_admm_convergence_criterion,
  objective_tol = as.numeric(args$polygram_admm_objective_tol),
  require_admm_convergence = tolower(args$polygram_admm_require_admm_convergence) == "true",
  normalize_constraint_rows = tolower(args$polygram_admm_normalize_constraint_rows) == "true",
  warm_start = tolower(args$polygram_admm_warm_start) == "true",
  seed = as.integer(args$seed),
  verbose = tolower(args$polygram_verbose) == "true"
)
prob = predict_polygram(fit, analysis)
elapsed = proc.time()[["elapsed"]] - start_time

predictions = data.frame(
  scenario_id = scenario_id,
  method = "polygram",
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
  method = "polygram",
  summarize_polygram(fit),
  stringsAsFactors = FALSE
)

edge_summary = polygram_active_edge_summary(fit, scenario_id = scenario_id)
# Save the fitted probability grid once; script 05 uses it for the decision-boundary plot.
grid = make_simplex_grid(grid_n = 260, eps = 1e-8)
grid$prob = predict_polygram(fit, grid)

########## Save fitted model outputs ##########

saveRDS(fit, file.path(args$out_dir, paste0(scenario_id, "_fit.rds")))
write.csv(edge_summary, file.path(args$out_dir, paste0(scenario_id, "_active_edge_summary.csv")), row.names = FALSE)
write.csv(grid, file.path(args$out_dir, paste0(scenario_id, "_fitted_probability_grid.csv")), row.names = FALSE)
write.csv(predictions, file.path(args$out_dir, "predictions.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(args$out_dir, "diagnostics.csv"), row.names = FALSE)
message("Wrote real-data polygram fits to: ", args$out_dir)
