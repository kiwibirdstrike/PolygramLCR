args = list(
  data_dir = file.path("..", "outputs", "simulation", "data"),
  out_dir = file.path("..", "outputs", "simulation", "fits", "ilr_logistic"),
  scenarios = "sim0,sim1,sim2,sim3"
)

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
prediction_rows = list()
fit_index = 1L

scenarios = trimws(strsplit(args$scenarios, ",", fixed = TRUE)[[1]])

########## Fit simulation logistic models ##########

for (slug in scenarios) {
  files = list.files(file.path(args$data_dir, slug), pattern = "_train[.]csv$", full.names = TRUE)
  for (train_path in files) {
    # Each replicate has a paired train/test file generated in script 01.
    replicate = as.integer(sub("^rep0*([0-9]+)_train[.]csv$", "\\1", basename(train_path)))
    test_path = sub("_train[.]csv$", "_test.csv", train_path)
    train = read.csv(train_path)
    test = read.csv(test_path)

    start_time = proc.time()[["elapsed"]]
    # Baseline global logistic regression in ILR coordinates.
    fit = glm(y ~ z1 + z2, data = train, family = binomial())
    prob = as.numeric(predict(fit, newdata = test, type = "response"))
    elapsed = proc.time()[["elapsed"]] - start_time
    saveRDS(fit, file.path(args$out_dir, sprintf("%s_rep%03d_fit.rds", slug, replicate)))

    prediction_rows[[fit_index]] = data.frame(
      scenario_id = test$scenario_id,
      replicate = replicate,
      method = "ilr_logistic",
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
    fit_index = fit_index + 1L
  }
}

########## Save fitted model outputs ##########

# Predictions are saved in one long table because script 03 computes all metrics from this file.
write.csv(do.call(rbind, prediction_rows), file.path(args$out_dir, "predictions.csv"), row.names = FALSE)
message("Wrote simulation logistic fits to: ", args$out_dir)
