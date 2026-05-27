library(polygramLCR)

args = list(
  data_file = file.path("data", "ATUS_employed_vs_not_employed.csv"),
  out_dir = file.path("..", "outputs", "real_data", "fits", "ilr_logistic"),
  years = "2022,2023,2024",
  age_min = "18",
  scenario_id = "ATUS"
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
# Save this once because the later real-data scripts and figures use the same analysis rows.
composition_summary = aggregate(
  analysis[, c("x1", "x2", "x3")],
  by = list(response = analysis$y),
  FUN = mean
)
composition_summary = data.frame(
  scenario_id = scenario_id,
  composition_summary,
  stringsAsFactors = FALSE
)

write.csv(analysis, file.path(args$out_dir, paste0(scenario_id, "_analysis_data.csv")), row.names = FALSE)
write.csv(composition_summary, file.path(args$out_dir, paste0(scenario_id, "_composition_summary.csv")), row.names = FALSE)

########## Fit real-data logistic model ##########

start_time = proc.time()[["elapsed"]]
# Global ILR logistic regression is the simplest real-data benchmark.
fit = glm(y ~ z1 + z2, data = analysis, family = binomial())
prob = as.numeric(predict(fit, newdata = analysis, type = "response"))
elapsed = proc.time()[["elapsed"]] - start_time

predictions = data.frame(
  scenario_id = scenario_id,
  method = "ilr_logistic",
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
scenario_summary = data.frame(
  scenario_id = scenario_id,
  n = nrow(analysis),
  n_y0 = sum(analysis$y == 0),
  n_y1 = sum(analysis$y == 1),
  x1_mean = mean(analysis$x1),
  x2_mean = mean(analysis$x2),
  x3_mean = mean(analysis$x3),
  stringsAsFactors = FALSE
)

########## Save fitted model outputs ##########

saveRDS(fit, file.path(args$out_dir, paste0(scenario_id, "_fit.rds")))
write.csv(predictions, file.path(args$out_dir, "predictions.csv"), row.names = FALSE)
write.csv(scenario_summary, file.path(args$out_dir, "scenario_summary.csv"), row.names = FALSE)
message("Wrote real-data logistic fits to: ", args$out_dir)
