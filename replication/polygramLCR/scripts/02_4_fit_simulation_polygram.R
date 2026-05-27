library(polygramLCR)

args = list(
  data_dir = file.path("..", "outputs", "simulation", "data"),
  out_dir = file.path("..", "outputs", "simulation", "fits", "polygram"),
  scenarios = "sim0,sim1,sim2,sim3",
  seed_base = "75",
  polygram_centers = "100",
  polygram_high_centers = "200",
  polygram_high_centers_scenarios = "SIM-03-complex-radial-peaks",
  polygram_admm_lambdas = "100",
  polygram_admm_lambda_max = "5",
  polygram_admm_lambda_max_min_ratio = "0.01",
  polygram_admm_selection = "aic",
  polygram_admm_rho = "auto",
  polygram_admm_max_iterations = "500",
  polygram_admm_max_newton_iterations = "50",
  polygram_admm_hessian_epsilon = "1e-8",
  polygram_admm_require_newton_convergence = "true",
  polygram_admm_convergence_criterion = "objective",
  polygram_admm_objective_tol = "1e-5",
  polygram_admm_require_admm_convergence = "true",
  polygram_admm_normalize_constraint_rows = "true",
  polygram_admm_warm_start = "true",
  polygram_verbose = "false"
)

polygram_high_centers_scenarios = trimws(strsplit(args$polygram_high_centers_scenarios, ",", fixed = TRUE)[[1]])

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
prediction_rows = list()
diagnostic_rows = list()
path_rows = list()
fit_index = 1L

scenarios = trimws(strsplit(args$scenarios, ",", fixed = TRUE)[[1]])

########## Fit simulation polygram models ##########

for (slug in scenarios) {
  files = list.files(file.path(args$data_dir, slug), pattern = "_train[.]csv$", full.names = TRUE)
  for (train_path in files) {
    # Each replicate has a paired train/test file generated in script 01.
    replicate = as.integer(sub("^rep0*([0-9]+)_train[.]csv$", "\\1", basename(train_path)))
    test_path = sub("_train[.]csv$", "_test.csv", train_path)
    train = read.csv(train_path)
    test = read.csv(test_path)
    scenario_index = as.integer(sub("^SIM-0*([0-9]+).*$", "\\1", train$scenario_id[[1]])) + 1L
    seed = as.integer(args$seed_base) + scenario_index * 10000L + replicate * 100L

    start_time = proc.time()[["elapsed"]]
    # SIM03 uses more polygram vertices because its response surface has more local structure.
    centers = if (train$scenario_id[[1]] %in% polygram_high_centers_scenarios) {
      as.integer(args$polygram_high_centers)
    } else {
      as.integer(args$polygram_centers)
    }
    fit = polygram(
      train,
      centers = centers,
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
      seed = seed,
      verbose = tolower(args$polygram_verbose) == "true"
    )
    prob = predict_polygram(fit, test)
    elapsed = proc.time()[["elapsed"]] - start_time
    saveRDS(fit, file.path(args$out_dir, sprintf("%s_rep%03d_fit.rds", slug, replicate)))

    prediction_rows[[fit_index]] = data.frame(
      scenario_id = test$scenario_id,
      replicate = replicate,
      method = "polygram",
      n_train = nrow(train),
      n_test = nrow(test),
      x1 = test$x1,
      x2 = test$x2,
      x3 = test$x3,
      z1 = test$z1,
      z2 = test$z2,
      v1 = test$v1,
      v2 = test$v2,
      true_prob = test$true_prob,
      y = test$y,
      prob = prob,
      elapsed_seconds = elapsed,
      stringsAsFactors = FALSE
    )
    diagnostic_rows[[fit_index]] = data.frame(
      scenario_id = train$scenario_id[[1]],
      replicate = replicate,
      method = "polygram",
      summarize_polygram(fit),
      stringsAsFactors = FALSE
    )
    # The full ADMM path is kept so convergence and lambda selection can be inspected later.
    path_rows[[fit_index]] = data.frame(
      scenario_id = train$scenario_id[[1]],
      replicate = replicate,
      method = "polygram",
      fit$path_summary,
      stringsAsFactors = FALSE
    )
    fit_index = fit_index + 1L
  }
}

########## Save fitted model outputs ##########

# Script 03 reads these three files to compute metrics and draw active-edge figures.
write.csv(do.call(rbind, prediction_rows), file.path(args$out_dir, "predictions.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostic_rows), file.path(args$out_dir, "diagnostics.csv"), row.names = FALSE)
write.csv(do.call(rbind, path_rows), file.path(args$out_dir, "admm_path.csv"), row.names = FALSE)
message("Wrote simulation polygram fits to: ", args$out_dir)
