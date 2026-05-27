library(polygramLCR)

args = list(
  fit_dir = file.path("..", "outputs", "simulation", "fits"),
  results_dir = file.path("..", "outputs", "simulation", "results"),
  figures_dir = file.path("..", "outputs", "figures"),
  representative_replicate = "1",
  methods = "ilr_logistic,ilr_gam,ilr_svm_rbf,polygram"
)

add_simplex_component_labels = function(labels = c(expression(x[1]), expression(x[2]), expression(x[3]))) {
  text(0.00, -0.030, labels[1], adj = c(0, 1), cex = 0.72, xpd = NA)
  text(1.00, -0.030, labels[2], adj = c(1, 1), cex = 0.72, xpd = NA)
  text(0.50, sqrt(3) / 2 + 0.020, labels[3], adj = c(0.5, 0), cex = 0.72, xpd = NA)
}

scenario_truth = function(scenario_id, X) {
  # Same truth definitions as script 01, used here only for reference displays.
  Z = comp_to_ilr(X)
  if (scenario_id == "SIM-00-linear-ilr") {
    eta = -0.60 + 1.50 * Z[, "z2"]
  } else if (scenario_id == "SIM-01-ilr-piecewise-boundary") {
    eta = -2.50 + 5.00 * X[, "x2"]
  } else if (scenario_id == "SIM-02-simplex-region") {
    peaks = data.frame(
      center1 = c(-0.9187318, 0.9187318, 0.0000000),
      center2 = c(-0.5304301, -0.5304301, 1.0608601),
      width1 = c(0.35, 0.35, 0.35),
      width2 = c(0.35, 0.35, 0.35),
      coefficient = c(4.80, 4.80, 4.80)
    )
    scores = sapply(seq_len(nrow(peaks)), function(i) {
      exp(-0.5 * (((Z[, "z1"] - peaks$center1[[i]]) / peaks$width1[[i]])^2 +
        ((Z[, "z2"] - peaks$center2[[i]]) / peaks$width2[[i]])^2))
    })
    eta = -2.10 + as.numeric(scores %*% peaks$coefficient)
  } else if (scenario_id == "SIM-03-complex-radial-peaks") {
    peaks = data.frame(
      center1 = c(-1.45, -0.15, 1.20),
      center2 = c(-0.95, 1.10, -0.45),
      width1 = c(0.95, 0.90, 0.80),
      width2 = c(0.70, 0.65, 0.75),
      coefficient = c(5.80, 4.45, 3.55)
    )
    scores = sapply(seq_len(nrow(peaks)), function(i) {
      pmax(
        1 -
          abs((Z[, "z1"] - peaks$center1[[i]]) / peaks$width1[[i]]) -
          abs((Z[, "z2"] - peaks$center2[[i]]) / peaks$width2[[i]]),
        0
      )
    })
    eta = -2.10 + as.numeric(scores %*% peaks$coefficient)
  }
  list(eta = eta, prob = plogis(eta))
}

add_truth_probability_contours = function(scenario_id, levels = 0.5, col = "#222222", lty = 3, lwd = 2.2) {
  # Draw contours of the known data-generating probability surface.
  grid = make_simplex_grid(grid_n = 220, eps = 1e-8)
  truth = scenario_truth(scenario_id, as.matrix(grid[, c("x1", "x2", "x3")]))
  v1_seq = sort(unique(grid$v1))
  v2_seq = sort(unique(grid$v2))
  zmat = matrix(NA_real_, nrow = length(v1_seq), ncol = length(v2_seq))
  zmat[cbind(match(grid$v1, v1_seq), match(grid$v2, v2_seq))] = truth$prob
  keep_levels = levels[levels >= min(truth$prob, na.rm = TRUE) & levels <= max(truth$prob, na.rm = TRUE)]
  if (length(keep_levels) == 0) {
    return(invisible(FALSE))
  }
  contour(v1_seq, v2_seq, zmat, levels = keep_levels, add = TRUE, drawlabels = FALSE, col = col, lty = lty, lwd = lwd)
  invisible(TRUE)
}

add_true_simulation_structure = function(scenario_id, col = "#222222", lty = 3, lwd = 2.2) {
  if (scenario_id == "SIM-00-linear-ilr") {
    z = cbind(z1 = seq(-3.2, 3.2, length.out = 500), z2 = 0.40)
    v = comp_to_bary(ilr_to_comp(z))
    lines(v[, "v1"], v[, "v2"], col = col, lty = lty, lwd = lwd)
  } else if (scenario_id == "SIM-01-ilr-piecewise-boundary") {
    add_truth_probability_contours(scenario_id, levels = 0.5, col = col, lty = lty, lwd = lwd)
  } else if (scenario_id == "SIM-02-simplex-region") {
    add_truth_probability_contours(scenario_id, levels = c(0.25, 0.5, 0.75), col = col, lty = lty, lwd = lwd)
  } else if (scenario_id == "SIM-03-complex-radial-peaks") {
    add_truth_probability_contours(scenario_id, levels = c(0.25, 0.5, 0.75), col = col, lty = lty, lwd = lwd)
  }
}

draw_true_probability_background = function(scenario_id, n_grid = 85) {
  # Faint background shading helps compare fitted active edges to the truth.
  grid = make_simplex_grid(grid_n = n_grid, eps = 1e-8)
  truth = scenario_truth(scenario_id, as.matrix(grid[, c("x1", "x2", "x3")]))
  prob = pmin(pmax(truth$prob, 0), 1)
  cols = colorRampPalette(
    c("#FFF7EC", "#FDD49E", "#FC8D59", "#D7301F", "#7F0000")
  )(101)
  col = adjustcolor(cols[pmax(1, pmin(101, floor(prob * 100) + 1))], alpha.f = 0.34)
  points(grid$v1, grid$v2, pch = 15, cex = 0.48, col = col)
}

scenario_title = function(scenario_id) {
  titles = c(
    `SIM-00-linear-ilr` = "SIM00 Linear ILR",
    `SIM-01-ilr-piecewise-boundary` = "SIM01 High x2",
    `SIM-02-simplex-region` = "SIM02 Matched local peaks",
    `SIM-03-complex-radial-peaks` = "SIM03 Heterogeneous local peaks"
  )
  titles[[scenario_id]]
}

draw_polygram_active_edges = function(
  fit,
  edge_summary,
  active_epsilon = 1e-7,
  edge_trim = 0.015,
  strength_col = "gradient_jump",
  top_percent = NULL,
  scale_by_strength = FALSE,
  edge_color = "#D7301F",
  edge_lwd = 1.45,
  alpha = 0.9
) {
  edge_summary$active_interior = edge_summary$active_interior &
    polygram_active_edge_mask(fit, active_epsilon = active_epsilon, edge_trim = edge_trim)
  keep = edge_summary$active_interior
  if (!is.null(top_percent)) {
    top_percent = as.numeric(top_percent)
    if (top_percent <= 1) {
      top_percent = top_percent * 100
    }
    top_percent = min(top_percent, 100)

    active_index = which(keep)
    active_strength = edge_summary[[strength_col]][active_index]
    n_keep = max(1L, ceiling(length(active_index) * top_percent / 100))
    selected = active_index[order(active_strength, decreasing = TRUE)[seq_len(n_keep)]]
    keep = rep(FALSE, nrow(edge_summary))
    keep[selected] = TRUE
  }

  vertices_bary = comp_to_bary(ilr_to_comp(as.matrix(fit$vertices)))
  edges = fit$setup$common_edge[keep, , drop = FALSE]
  if (scale_by_strength) {
    strength = edge_summary[[strength_col]][keep]
    strength_scale = quantile(strength, 0.95)
    scaled = pmin(strength / strength_scale, 1)
    edge_cols = if (tolower(edge_color) %in% c("black", "#000000")) {
      colorRampPalette(c("#D7D7D7", "#9E9E9E", "#5A5A5A", "#000000"))(101)
    } else {
      colorRampPalette(c("#CFCFCF", "#FEE0D2", "#FC9272", edge_color))(101)
    }
    edge_col = adjustcolor(
      edge_cols[pmax(1, pmin(101, floor(scaled * 100) + 1))],
      alpha.f = alpha
    )
    edge_lwd = 0.45 + 2.35 * sqrt(scaled)
  } else {
    edge_col = adjustcolor(edge_color, alpha.f = alpha)
  }
  segments(
    vertices_bary[edges[, 1], "v1"],
    vertices_bary[edges[, 1], "v2"],
    vertices_bary[edges[, 2], "v1"],
    vertices_bary[edges[, 2], "v2"],
    col = edge_col,
    lwd = edge_lwd
  )
  invisible(TRUE)
}

dir.create(args$results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$figures_dir, recursive = TRUE, showWarnings = FALSE)

methods = trimws(strsplit(args$methods, ",", fixed = TRUE)[[1]])
prediction_files = file.path(args$fit_dir, methods, "predictions.csv")
# All methods write the same prediction schema, so metrics start from one combined table.
predictions = do.call(rbind, lapply(prediction_files, read.csv))

########## Calculate simulation metrics ##########

metric_groups = unique(predictions[c("scenario_id", "replicate", "method")])
metric_rows = vector("list", nrow(metric_groups))
for (i in seq_len(nrow(metric_groups))) {
  # Compute metrics separately for every scenario, replicate, and method.
  keep = predictions$scenario_id == metric_groups$scenario_id[i] &
    predictions$replicate == metric_groups$replicate[i] &
    predictions$method == metric_groups$method[i]
  dat = predictions[keep, , drop = FALSE]
  metric_rows[[i]] = data.frame(
    scenario_id = metric_groups$scenario_id[i],
    replicate = metric_groups$replicate[i],
    method = metric_groups$method[i],
    n_train = dat$n_train[[1]],
    n_test = dat$n_test[[1]],
    status = "ok",
    error_message = "",
    elapsed_seconds = dat$elapsed_seconds[[1]],
    classification_metrics(dat$y, dat$prob),
    stringsAsFactors = FALSE
  )
}

replicate_metrics = do.call(rbind, metric_rows)

summary_groups = unique(replicate_metrics[c("scenario_id", "method")])
numeric_cols = names(replicate_metrics)[vapply(replicate_metrics, is.numeric, logical(1))]
numeric_cols = setdiff(numeric_cols, "replicate")
summary_rows = vector("list", nrow(summary_groups))
for (i in seq_len(nrow(summary_groups))) {
  # Summaries average replicate-level metrics; fixed sample sizes keep only means.
  keep = replicate_metrics$scenario_id == summary_groups$scenario_id[[i]] &
    replicate_metrics$method == summary_groups$method[[i]]
  dat = replicate_metrics[keep, , drop = FALSE]
  out = data.frame(
    scenario_id = summary_groups$scenario_id[[i]],
    method = summary_groups$method[[i]],
    n_replicates = nrow(dat),
    stringsAsFactors = FALSE
  )
  for (col in numeric_cols) {
    out[[paste0(col, "_mean")]] = mean(dat[[col]], na.rm = TRUE)
    if (!col %in% c("n_train", "n_test")) {
      out[[paste0(col, "_sd")]] = sd(dat[[col]], na.rm = TRUE)
    }
  }
  summary_rows[[i]] = out
}
summary_metrics = do.call(rbind, summary_rows)

representative_predictions = predictions[
  predictions$replicate == as.integer(args$representative_replicate),
  ,
  drop = FALSE
]

diagnostic_files = file.path(args$fit_dir, methods, "diagnostics.csv")
diagnostic_files = diagnostic_files[file.exists(diagnostic_files)]
diagnostic_tables = lapply(diagnostic_files, read.csv)
# Diagnostics differ by method, so align columns before stacking.
diagnostic_cols = unique(unlist(lapply(diagnostic_tables, names), use.names = FALSE))
for (i in seq_along(diagnostic_tables)) {
  missing_cols = setdiff(diagnostic_cols, names(diagnostic_tables[[i]]))
  for (col in missing_cols) {
    diagnostic_tables[[i]][[col]] = NA
  }
  diagnostic_tables[[i]] = diagnostic_tables[[i]][, diagnostic_cols, drop = FALSE]
}
polygram_diagnostics = do.call(rbind, diagnostic_tables)
admm_path = file.path(args$fit_dir, "polygram", "admm_path.csv")
polygram_admm_path = read.csv(admm_path)

########## Save simulation result tables ##########

write.csv(replicate_metrics, file.path(args$results_dir, "replicate_metrics.csv"), row.names = FALSE)
write.csv(summary_metrics, file.path(args$results_dir, "summary_metrics.csv"), row.names = FALSE)
write.csv(representative_predictions, file.path(args$results_dir, "representative_predictions.csv"), row.names = FALSE)
write.csv(polygram_diagnostics, file.path(args$results_dir, "polygram_diagnostics.csv"), row.names = FALSE)
write.csv(polygram_admm_path, file.path(args$results_dir, "polygram_admm_path.csv"), row.names = FALSE)

scenario_order = unique(representative_predictions$scenario_id)
scenario_order = scenario_order[order(as.integer(sub("^SIM-0*([0-9]+).*$", "\\1", scenario_order)))]
scenario_slugs = setNames(sub("^SIM-0*([0-9]+).*$", "sim\\1", scenario_order), scenario_order)

########## Draw simulation active-edge figures ##########

active_rows = list()
sim03_fit = NULL
sim03_edge_summary = NULL
# Appendix figure: one representative active-edge display for each simulation setting.
pdf(file.path(args$figures_dir, "figA1_simulation_active_edges.pdf"), width = 7.6, height = 8.3)
old_par = par(mfrow = c(2, 2), mar = c(0.92, 0.32, 1.75, 0.32), oma = c(1.55, 0.1, 0.2, 0.1))
for (scenario_id in scenario_order) {
  slug = scenario_slugs[[scenario_id]]
  fit_path = file.path(args$fit_dir, "polygram", sprintf("%s_rep%03d_fit.rds", slug, as.integer(args$representative_replicate)))
  fit = readRDS(fit_path)
  edge_summary = polygram_active_edge_summary(fit, scenario_id = scenario_id)
  active_rows[[length(active_rows) + 1L]] = edge_summary
  if (scenario_id == "SIM-03-complex-radial-peaks") {
    sim03_fit = fit
    sim03_edge_summary = edge_summary
  }
  active_keep = edge_summary$active_interior &
    polygram_active_edge_mask(fit, edge_trim = 0.015)
  draw_simplex_frame(paste0(scenario_title(scenario_id), ": active edges"))
  add_true_simulation_structure(scenario_id, col = rgb(0, 0.25, 0.65, 0.42), lty = 3, lwd = 0.90)
  draw_polygram_active_edges(fit, edge_summary = edge_summary, edge_trim = 0.015, scale_by_strength = TRUE)
  add_simplex_component_labels(c(expression(x[1]), expression(x[2]), expression(x[3])))
  mtext(sprintf("%d active interior edges", sum(active_keep)), side = 1, line = -0.16, cex = 0.58)
}
par(old_par)
dev.off()

pdf(file.path(args$figures_dir, "fig02_sim03_active_edges.pdf"), width = 9.80, height = 4.95)
old_par = par(mfrow = c(1, 2), mar = c(1.05, 0.35, 2.00, 0.35), oma = c(0, 0, 1.10, 0))

# Main simulation figure: true SIM03 facet changes beside fitted active edges.
draw_simplex_frame("True SIM03 facet-change structure")
draw_true_probability_background("SIM-03-complex-radial-peaks")
centers = data.frame(
  center1 = c(-1.45, -0.15, 1.20),
  center2 = c(-0.95, 1.10, -0.45),
  width1 = c(0.95, 0.90, 0.80),
  width2 = c(0.70, 0.65, 0.75)
)
for (i in seq_len(nrow(centers))) {
  c1 = centers$center1[[i]]
  c2 = centers$center2[[i]]
  w1 = centers$width1[[i]]
  w2 = centers$width2[[i]]
  diamond = rbind(
    c(c1, c2 + w2),
    c(c1 + w1, c2),
    c(c1, c2 - w2),
    c(c1 - w1, c2),
    c(c1, c2 + w2)
  )
  colnames(diamond) = c("z1", "z2")
  diamond_bary = comp_to_bary(ilr_to_comp(diamond))
  lines(diamond_bary[, "v1"], diamond_bary[, "v2"], col = adjustcolor("#2166AC", alpha.f = 0.80), lwd = 1.30, lty = 2)

  vertical = rbind(c(c1, c2 - w2), c(c1, c2 + w2))
  horizontal = rbind(c(c1 - w1, c2), c(c1 + w1, c2))
  colnames(vertical) = colnames(horizontal) = c("z1", "z2")
  vertical_bary = comp_to_bary(ilr_to_comp(vertical))
  horizontal_bary = comp_to_bary(ilr_to_comp(horizontal))
  lines(vertical_bary[, "v1"], vertical_bary[, "v2"], col = adjustcolor("#2166AC", alpha.f = 0.95), lwd = 1.75, lty = 1)
  lines(horizontal_bary[, "v1"], horizontal_bary[, "v2"], col = adjustcolor("#2166AC", alpha.f = 0.95), lwd = 1.75, lty = 1)
}
add_simplex_component_labels(c(expression(x[1]), expression(x[2]), expression(x[3])))
legend("topright", legend = c("true L1 facet-change lines"), col = c("#2166AC"), lty = c(1), lwd = c(1.75), bty = "n", cex = 0.68)
mtext("background = true success probability", side = 1, line = -0.16, cex = 0.56)

draw_simplex_frame("Fitted polygram active edges")
draw_true_probability_background("SIM-03-complex-radial-peaks")
draw_polygram_active_edges(
  sim03_fit,
  edge_summary = sim03_edge_summary,
  edge_trim = 0.015,
  top_percent = 50,
  scale_by_strength = TRUE,
  alpha = 0.92
)
add_simplex_component_labels(c(expression(x[1]), expression(x[2]), expression(x[3])))
legend("topright", legend = c("fitted active edges"), col = c("#DE2D26"), lty = c(1), lwd = c(2.25), bty = "n", cex = 0.68)
mtext("red edges = upper 50% active edges by gradient-jump strength", side = 1, line = -0.16, cex = 0.56)
mtext("SIM03 true facet changes and fitted active edges", outer = TRUE, side = 3, line = 0.05, cex = 1.05, font = 2)

par(old_par)
dev.off()

write.csv(do.call(rbind, active_rows), file.path(args$results_dir, "simulation_active_edge_summary.csv"), row.names = FALSE)

message("Wrote simulation results to: ", args$results_dir)
