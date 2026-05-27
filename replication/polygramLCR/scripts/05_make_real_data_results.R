library(polygramLCR)

args = list(
  fit_dir = file.path("..", "outputs", "real_data", "fits"),
  results_dir = file.path("..", "outputs", "real_data", "results"),
  figures_dir = file.path("..", "outputs", "figures"),
  methods = "ilr_logistic,ilr_gam,ilr_svm_rbf,polygram",
  scenario_id = "ATUS"
)

add_simplex_component_labels = function(labels = c("Maintenance", "Obligation", "Leisure")) {
  text(0.00, -0.030, labels[1], adj = c(0, 1), cex = 0.72, xpd = NA)
  text(1.00, -0.030, labels[2], adj = c(1, 1), cex = 0.72, xpd = NA)
  text(0.50, sqrt(3) / 2 + 0.020, labels[3], adj = c(0.5, 0), cex = 0.72, xpd = NA)
}

draw_atus_points = function(data, alpha, cex) {
  cols = ifelse(
    data$y == 1,
    rgb(178 / 255, 24 / 255, 43 / 255, alpha),
    rgb(33 / 255, 102 / 255, 172 / 255, alpha)
  )
  set.seed(1)
  order_index = sample(seq_len(nrow(data)))
  points(data$v1[order_index], data$v2[order_index], pch = 16, cex = cex, col = cols[order_index])
}

dir.create(args$results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$figures_dir, recursive = TRUE, showWarnings = FALSE)

########## Calculate real-data metrics ##########

methods = trimws(strsplit(args$methods, ",", fixed = TRUE)[[1]])
prediction_files = file.path(args$fit_dir, methods, "predictions.csv")
# Real-data metrics are computed from the saved fitted probabilities for each method.
predictions = do.call(rbind, lapply(prediction_files, read.csv))

model_groups = unique(predictions[c("scenario_id", "method")])
model_rows = vector("list", nrow(model_groups))
for (i in seq_len(nrow(model_groups))) {
  # One metric row per fitted method.
  keep = predictions$scenario_id == model_groups$scenario_id[i] &
    predictions$method == model_groups$method[i]
  dat = predictions[keep, , drop = FALSE]
  model_rows[[i]] = data.frame(
    scenario_id = model_groups$scenario_id[i],
    method = model_groups$method[i],
    n = nrow(dat),
    n_y0 = sum(dat$y == 0),
    n_y1 = sum(dat$y == 1),
    elapsed_seconds = dat$elapsed_seconds[[1]],
    classification_metrics(dat$y, dat$prob),
    stringsAsFactors = FALSE
  )
}

model_summary = do.call(rbind, model_rows)

diagnostic_files = file.path(args$fit_dir, methods, "diagnostics.csv")
diagnostic_files = diagnostic_files[file.exists(diagnostic_files)]
diagnostic_rows = lapply(diagnostic_files, read.csv)
# Diagnostics have method-specific columns, so align them before stacking.
diagnostic_cols = unique(unlist(lapply(diagnostic_rows, names), use.names = FALSE))
diagnostic_rows = lapply(diagnostic_rows, function(row) {
  for (col in setdiff(diagnostic_cols, names(row))) {
    row[[col]] = NA
  }
  row[, diagnostic_cols, drop = FALSE]
})
diagnostics = do.call(rbind, diagnostic_rows)

scenario_summary = read.csv(file.path(args$fit_dir, "ilr_logistic", "scenario_summary.csv"))
composition_summary = read.csv(file.path(args$fit_dir, "ilr_logistic", paste0(args$scenario_id, "_composition_summary.csv")))

########## Save real-data result tables ##########

write.csv(model_summary, file.path(args$results_dir, "atus_model_summary.csv"), row.names = FALSE)
write.csv(scenario_summary, file.path(args$results_dir, "atus_scenario_summary.csv"), row.names = FALSE)
write.csv(composition_summary, file.path(args$results_dir, "atus_summary.csv"), row.names = FALSE)
write.csv(predictions, file.path(args$results_dir, "atus_predictions.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(args$results_dir, "atus_model_diagnostics.csv"), row.names = FALSE)

########## Draw ATUS simplex figures ##########

analysis = read.csv(file.path(args$fit_dir, "ilr_logistic", paste0(args$scenario_id, "_analysis_data.csv")))
fit = readRDS(file.path(args$fit_dir, "polygram", paste0(args$scenario_id, "_fit.rds")))
edge_summary = read.csv(file.path(args$fit_dir, "polygram", paste0(args$scenario_id, "_active_edge_summary.csv")))
# The probability grid was saved by script 04_4 and is reused for the decision boundary.
grid = read.csv(file.path(args$fit_dir, "polygram", paste0(args$scenario_id, "_fitted_probability_grid.csv")))

pdf(file.path(args$figures_dir, "fig03_atus_sample.pdf"), width = 5.8, height = 5.3)
old_par = par(mar = c(1.05, 0.45, 2.00, 0.45), oma = c(0, 0, 0, 0))
draw_simplex_frame("Observed ATUS compositions")
add_simplex_component_labels(c("Maintenance", "Obligation", "Leisure"))
draw_atus_points(analysis, alpha = 0.36, cex = 0.78)
legend(
  "topright",
  legend = c("Not employed", "Employed"),
  col = c(
    rgb(33 / 255, 102 / 255, 172 / 255, 0.78),
    rgb(178 / 255, 24 / 255, 43 / 255, 0.78)
  ),
  pch = 16,
  bty = "n",
  cex = 0.68
)
n0 = sum(analysis$y == 0, na.rm = TRUE)
n1 = sum(analysis$y == 1, na.rm = TRUE)
mtext(
  sprintf("n = %s; not employed = %s, employed = %s", format(nrow(analysis), big.mark = ","), format(n0, big.mark = ","), format(n1, big.mark = ",")),
  side = 1,
  line = -0.16,
  cex = 0.56
)
par(old_par)
dev.off()

pdf(file.path(args$figures_dir, "fig04_atus_active_edges.pdf"), width = 9.8, height = 4.95)
old_par = par(mfrow = c(1, 2), mar = c(1.05, 0.35, 2.00, 0.35), oma = c(0, 0, 1.10, 0))

draw_simplex_frame(expression(hat(p) == 0.5 ~ "decision boundary"))
add_simplex_component_labels(c("Maintenance", "Obligation", "Leisure"))
draw_atus_points(analysis, alpha = 0.055, cex = 1.05)

# Convert the saved probability grid into a matrix for the p-hat = 0.5 contour.
v1_seq = sort(unique(grid$v1))
v2_seq = sort(unique(grid$v2))
boundary_prob = matrix(NA_real_, nrow = length(v1_seq), ncol = length(v2_seq))
boundary_prob[cbind(match(grid$v1, v1_seq), match(grid$v2, v2_seq))] = grid$prob
contour(
  v1_seq,
  v2_seq,
  boundary_prob,
  levels = 0.5,
  add = TRUE,
  drawlabels = FALSE,
  col = "#000000",
  lwd = 1.85
)

legend(
  "topright",
  legend = c("Not employed", "Employed", expression(hat(p) == 0.5)),
  col = c(
    rgb(33 / 255, 102 / 255, 172 / 255, 0.58),
    rgb(178 / 255, 24 / 255, 43 / 255, 0.58),
    "#000000"
  ),
  pch = c(16, 16, NA),
  lty = c(NA, NA, 1),
  lwd = c(NA, NA, 2.25),
  bty = "n",
  cex = 0.66
)
mtext("classification boundary from the fitted probability surface", side = 1, line = -0.16, cex = 0.56)

top_percent = 50
draw_simplex_frame("Fitted polygram active edges")
add_simplex_component_labels(c("Maintenance", "Obligation", "Leisure"))
draw_atus_points(analysis, alpha = 0.055, cex = 1.05)

edge_summary_plot = edge_summary
edge_summary_plot$active_interior = edge_summary_plot$active_interior &
  polygram_active_edge_mask(fit, edge_trim = 0.06)
# Plot the upper half of active edges by gradient-jump strength.
active_index = which(edge_summary_plot$active_interior)
n_keep = max(1L, ceiling(length(active_index) * top_percent / 100))
selected_edges = active_index[order(edge_summary_plot$gradient_jump[active_index], decreasing = TRUE)[seq_len(n_keep)]]
active_count = length(selected_edges)

vertices_bary = comp_to_bary(ilr_to_comp(as.matrix(fit$vertices)))
edges = fit$setup$common_edge[selected_edges, , drop = FALSE]
strength = edge_summary_plot$gradient_jump[selected_edges]
scaled_strength = pmin(strength / quantile(strength, 0.95), 1)
edge_cols = colorRampPalette(c("#D7D7D7", "#9E9E9E", "#5A5A5A", "#000000"))(101)
edge_col = edge_cols[pmax(1, pmin(101, floor(scaled_strength * 100) + 1))]
edge_lwd = 0.45 + 2.35 * sqrt(scaled_strength)
segments(
  vertices_bary[edges[, 1], "v1"],
  vertices_bary[edges[, 1], "v2"],
  vertices_bary[edges[, 2], "v1"],
  vertices_bary[edges[, 2], "v2"],
  col = edge_col,
  lwd = edge_lwd
)

legend(
  "topright",
  legend = c("Not employed", "Employed", "active edges"),
  col = c(
    rgb(33 / 255, 102 / 255, 172 / 255, 0.58),
    rgb(178 / 255, 24 / 255, 43 / 255, 0.58),
    "#000000"
  ),
  pch = c(16, 16, NA),
  lty = c(NA, NA, 1),
  lwd = c(NA, NA, 2.25),
  bty = "n",
  cex = 0.66
)
mtext(
  sprintf("top %d%% by gradient-jump strength (%d active edges)", top_percent, active_count),
  side = 1,
  line = -0.16,
  cex = 0.56
)

mtext("ATUS decision boundary and fitted active-edge structure", outer = TRUE, side = 3, line = 0.05, cex = 1.05, font = 2)
par(old_par)
dev.off()

write.csv(edge_summary, file.path(args$results_dir, "atus_active_edge_summary.csv"), row.names = FALSE)
write.csv(grid, file.path(args$results_dir, "atus_fitted_probability_grid.csv"), row.names = FALSE)

message("Wrote real-data results to: ", args$results_dir)
