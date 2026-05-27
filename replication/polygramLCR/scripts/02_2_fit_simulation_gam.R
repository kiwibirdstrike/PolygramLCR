library(mgcv)

args = list(
  data_dir = file.path("..", "outputs", "simulation", "data"),
  out_dir = file.path("..", "outputs", "simulation", "fits", "ilr_gam"),
  scenarios = "sim0,sim1,sim2,sim3",
  gam_k_grid = "10,20,30,40,60,80,100",
  gam_high_k_grid = "200",
  gam_high_k_scenarios = "SIM-02-simplex-region,SIM-03-complex-radial-peaks"
)

gam_k_grid = as.numeric(strsplit(args$gam_k_grid, ",", fixed = TRUE)[[1]])
gam_high_k_grid = as.numeric(strsplit(args$gam_high_k_grid, ",", fixed = TRUE)[[1]])
gam_high_k_scenarios = trimws(strsplit(args$gam_high_k_scenarios, ",", fixed = TRUE)[[1]])

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
prediction_rows = list()
diagnostic_rows = list()
fit_index = 1L

scenarios = trimws(strsplit(args$scenarios, ",", fixed = TRUE)[[1]])

########## Fit simulation GAM models ##########

for (slug in scenarios) {
  files = list.files(file.path(args$data_dir, slug), pattern = "_train[.]csv$", full.names = TRUE)
  for (train_path in files) {
    # Each replicate has a paired train/test file generated in script 01.
    replicate = as.integer(sub("^rep0*([0-9]+)_train[.]csv$", "\\1", basename(train_path)))
    test_path = sub("_train[.]csv$", "_test.csv", train_path)
    train = read.csv(train_path)
    test = read.csv(test_path)

    start_time = proc.time()[["elapsed"]]
    # The local-peak scenarios use a larger basis grid to give the GAM enough flexibility.
    k_grid = if (train$scenario_id[[1]] %in% gam_high_k_scenarios) gam_high_k_grid else gam_k_grid
    k_grid = as.integer(k_grid)
    fits = lapply(k_grid, function(k) {
      gam(y ~ s(z1, z2, k = k), data = train, family = binomial())
    })
    aic = vapply(fits, AIC, numeric(1))
    # Store the AIC-selected GAM and its tuning choice for later diagnostics.
    best_index = which.min(aic)
    fit = fits[[best_index]]
    attr(fit, "selected_k") = k_grid[[best_index]]
    attr(fit, "k_grid") = paste(k_grid, collapse = ",")
    attr(fit, "selection") = "aic"
    prob = as.numeric(predict(fit, newdata = test, type = "response"))
    elapsed = proc.time()[["elapsed"]] - start_time
    saveRDS(fit, file.path(args$out_dir, sprintf("%s_rep%03d_fit.rds", slug, replicate)))

    prediction_rows[[fit_index]] = data.frame(
      scenario_id = test$scenario_id,
      replicate = replicate,
      method = "ilr_gam",
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
      method = "ilr_gam",
      selected_k = attr(fit, "selected_k"),
      k_grid = attr(fit, "k_grid"),
      selection = attr(fit, "selection"),
      stringsAsFactors = FALSE
    )
    fit_index = fit_index + 1L
  }
}

########## Save fitted model outputs ##########

# Predictions and selected k values are consumed by script 03.
write.csv(do.call(rbind, prediction_rows), file.path(args$out_dir, "predictions.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostic_rows), file.path(args$out_dir, "diagnostics.csv"), row.names = FALSE)
message("Wrote simulation GAM fits to: ", args$out_dir)
