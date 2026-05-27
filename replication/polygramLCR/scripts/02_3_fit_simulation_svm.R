library(e1071)

args = list(
  data_dir = file.path("..", "outputs", "simulation", "data"),
  out_dir = file.path("..", "outputs", "simulation", "fits", "ilr_svm_rbf"),
  scenarios = "sim0,sim1,sim2,sim3",
  seed_base = "75",
  svm_cost_grid = "0.25,1,4",
  svm_gamma_grid = "0.1,0.5,2",
  svm_high_cost_grid = "0.25,1,4,16",
  svm_high_gamma_grid = "0.5,2,4,8",
  svm_high_grid_scenarios = "SIM-03-complex-radial-peaks",
  svm_cv_folds = "3"
)

svm_cost_grid = as.numeric(strsplit(args$svm_cost_grid, ",", fixed = TRUE)[[1]])
svm_gamma_grid = as.numeric(strsplit(args$svm_gamma_grid, ",", fixed = TRUE)[[1]])
svm_high_cost_grid = as.numeric(strsplit(args$svm_high_cost_grid, ",", fixed = TRUE)[[1]])
svm_high_gamma_grid = as.numeric(strsplit(args$svm_high_gamma_grid, ",", fixed = TRUE)[[1]])
svm_high_grid_scenarios = trimws(strsplit(args$svm_high_grid_scenarios, ",", fixed = TRUE)[[1]])

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
prediction_rows = list()
diagnostic_rows = list()
fit_index = 1L

scenarios = trimws(strsplit(args$scenarios, ",", fixed = TRUE)[[1]])

########## Fit simulation SVM models ##########

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
    # The heterogeneous local-peak scenario uses a wider SVM tuning grid.
    cost_grid = if (train$scenario_id[[1]] %in% svm_high_grid_scenarios) svm_high_cost_grid else svm_cost_grid
    gamma_grid = if (train$scenario_id[[1]] %in% svm_high_grid_scenarios) svm_high_gamma_grid else svm_gamma_grid
    set.seed(seed)
    # Cross-validation selects the radial-kernel cost and gamma for this replicate.
    tuned = tune.svm(
      y ~ .,
      data = data.frame(y = factor(train$y), z1 = train$z1, z2 = train$z2),
      kernel = "radial",
      cost = as.numeric(cost_grid),
      gamma = as.numeric(gamma_grid),
      probability = TRUE,
      tunecontrol = tune.control(cross = as.integer(args$svm_cv_folds))
    )
    fit = tuned$best.model
    attr(fit, "selected_cost") = tuned$best.parameters$cost
    attr(fit, "selected_gamma") = tuned$best.parameters$gamma
    attr(fit, "cost_grid") = paste(cost_grid, collapse = ",")
    attr(fit, "gamma_grid") = paste(gamma_grid, collapse = ",")
    attr(fit, "selection") = paste0(args$svm_cv_folds, "-fold-cv")
    pred = predict(
      fit,
      newdata = data.frame(z1 = test$z1, z2 = test$z2),
      probability = TRUE
    )
    prob = as.numeric(attr(pred, "probabilities")[, "1"])
    elapsed = proc.time()[["elapsed"]] - start_time
    saveRDS(fit, file.path(args$out_dir, sprintf("%s_rep%03d_fit.rds", slug, replicate)))

    prediction_rows[[fit_index]] = data.frame(
      scenario_id = test$scenario_id,
      replicate = replicate,
      method = "ilr_svm_rbf",
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
      method = "ilr_svm_rbf",
      selected_cost = attr(fit, "selected_cost"),
      selected_gamma = attr(fit, "selected_gamma"),
      cost_grid = attr(fit, "cost_grid"),
      gamma_grid = attr(fit, "gamma_grid"),
      selection = attr(fit, "selection"),
      stringsAsFactors = FALSE
    )
    fit_index = fit_index + 1L
  }
}

########## Save fitted model outputs ##########

# Save both fitted probabilities and selected SVM tuning parameters.
write.csv(do.call(rbind, prediction_rows), file.path(args$out_dir, "predictions.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostic_rows), file.path(args$out_dir, "diagnostics.csv"), row.names = FALSE)
message("Wrote simulation SVM fits to: ", args$out_dir)
