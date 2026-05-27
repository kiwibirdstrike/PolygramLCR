library(polygramLCR)

args = list(
  replicates = "30",
  replicate_start = "1",
  n_train = "5000",
  n_test = "5000",
  high_n_train = "10000",
  high_n_scenarios = "SIM-03-complex-radial-peaks",
  seed_base = "75",
  scenarios = "SIM-00-linear-ilr,SIM-01-ilr-piecewise-boundary,SIM-02-simplex-region,SIM-03-complex-radial-peaks",
  out_dir = file.path("..", "outputs", "simulation", "data"),
  figures_dir = file.path("..", "outputs", "figures")
)

scenario_truth = function(scenario_id, X) {
  # True logit surfaces used in the four simulation settings.
  Z = comp_to_ilr(X)
  if (scenario_id == "SIM-00-linear-ilr") {
    eta = -0.60 + 1.50 * Z[, "z2"]
    boundary = rep(0.4, nrow(X))
  } else if (scenario_id == "SIM-01-ilr-piecewise-boundary") {
    eta = -2.50 + 5.00 * X[, "x2"]
    boundary = rep(NA_real_, nrow(X))
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
    boundary = rep(NA_real_, nrow(X))
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
    boundary = rep(NA_real_, nrow(X))
  }
  list(eta = eta, prob = plogis(eta), boundary = boundary)
}

generate_simulation_data = function(scenario_id, n, seed, split) {
  # Generate compositions and attach both ILR and barycentric coordinates.
  set.seed(seed)
  X = matrix(rgamma(n * 3L, shape = 1, rate = 1), nrow = n, ncol = 3)
  X = X / rowSums(X)
  colnames(X) = c("x1", "x2", "x3")
  Z = comp_to_ilr(X)
  V = comp_to_bary(X)
  truth = scenario_truth(scenario_id, X)
  y = rbinom(n, 1, truth$prob)
  data.frame(
    scenario_id = scenario_id,
    split = split,
    x1 = X[, "x1"],
    x2 = X[, "x2"],
    x3 = X[, "x3"],
    z1 = Z[, "z1"],
    z2 = Z[, "z2"],
    v1 = V[, "v1"],
    v2 = V[, "v2"],
    true_eta = truth$eta,
    true_prob = truth$prob,
    true_boundary_z2 = truth$boundary,
    y = y
  )
}

scenarios = trimws(strsplit(args$scenarios, ",", fixed = TRUE)[[1]])
replicates = seq.int(
  as.integer(args$replicate_start),
  as.integer(args$replicate_start) + as.integer(args$replicates) - 1L
)
simulation_jobs = expand.grid(
  scenario_id = scenarios,
  replicate = replicates,
  stringsAsFactors = FALSE
)

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

########## Generate simulation data ##########

high_scenarios = trimws(strsplit(args$high_n_scenarios, ",", fixed = TRUE)[[1]])
for (job_id in seq_len(nrow(simulation_jobs))) {
  scenario_id = simulation_jobs$scenario_id[[job_id]]
  replicate = simulation_jobs$replicate[[job_id]]
  scenario_index = match(scenario_id, scenarios)
  n_train = if (scenario_id %in% high_scenarios) as.integer(args$high_n_train) else as.integer(args$n_train)
  slug = sub("^SIM-0*([0-9]+).*$", "sim\\1", scenario_id)
  scenario_dir = file.path(args$out_dir, slug)
  dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)

  # Train and test sets use separate seeds for an independent evaluation sample.
  train_seed = as.integer(args$seed_base) + scenario_index * 10000L + replicate
  test_seed = as.integer(args$seed_base) + scenario_index * 10000L + replicate + 500000L
  train = generate_simulation_data(scenario_id, n_train, train_seed, split = "train")
  test = generate_simulation_data(scenario_id, as.integer(args$n_test), test_seed, split = "test")

  train$replicate = replicate
  test$replicate = replicate
  write.csv(train, file.path(scenario_dir, sprintf("rep%03d_train.csv", replicate)), row.names = FALSE)
  write.csv(test, file.path(scenario_dir, sprintf("rep%03d_test.csv", replicate)), row.names = FALSE)
}

########## Draw simulation truth figure ##########

dir.create(args$figures_dir, recursive = TRUE, showWarnings = FALSE)
scenario_titles = c(
  "SIM00: Linear ILR",
  "SIM01: High x2",
  "SIM02: Matched local peaks",
  "SIM03: Heterogeneous local peaks"
)
cols = colorRampPalette(c(
  "#3B0F70", "#2878B5", "#2DB27D", "#FDE725", "#F46D43", "#A50026"
))(101)
z_scale = 0.42

truth_grid = do.call(
  rbind,
  lapply(scenarios, function(scenario_id) {
    # This grid is for Figure 1 only; model fitting uses the simulated samples above.
    grid = make_simplex_grid(grid_n = 170, eps = 0.012)
    truth = scenario_truth(scenario_id, as.matrix(grid[, c("x1", "x2", "x3")]))
    data.frame(
      scenario_id = scenario_id,
      grid,
      true_prob = truth$prob,
      true_eta = truth$eta,
      true_boundary_z2 = truth$boundary
    )
  })
)

draw_surface_3d = function(panel, scenario_id, title) {
  # Lift the same true probability surface from the simplex into a 3D display.
  v1_seq = sort(unique(panel$v1))
  v2_seq = sort(unique(panel$v2))
  prob_mat = matrix(NA_real_, nrow = length(v1_seq), ncol = length(v2_seq))
  prob_mat[cbind(match(panel$v1, v1_seq), match(panel$v2, v2_seq))] = pmin(pmax(panel$true_prob, 0), 1)
  zmat = prob_mat * z_scale

  transform = persp(
    v1_seq,
    v2_seq,
    matrix(0, nrow = length(v1_seq), ncol = length(v2_seq)),
    zlim = c(0, z_scale),
    theta = -38,
    phi = 27,
    expand = 0.64,
    col = NA,
    border = "#00000000",
    ticktype = "simple",
    xlab = "",
    ylab = "",
    zlab = "p(X)",
    main = title,
    cex.axis = 0.52,
    cex.lab = 0.62
  )

  facets = list()
  facet_index = 1L
  # Draw facets by hand so missing points outside the simplex are skipped cleanly.
  for (row in seq_len(nrow(zmat) - 1L)) {
    for (col_index in seq_len(ncol(zmat) - 1L)) {
      corners = rbind(
        c(v1_seq[row], v2_seq[col_index], zmat[row, col_index]),
        c(v1_seq[row + 1L], v2_seq[col_index], zmat[row + 1L, col_index]),
        c(v1_seq[row + 1L], v2_seq[col_index + 1L], zmat[row + 1L, col_index + 1L]),
        c(v1_seq[row], v2_seq[col_index + 1L], zmat[row, col_index + 1L])
      )
      if (anyNA(corners[, 3])) {
        next
      }
      triangles = list(corners[c(1L, 2L, 3L), , drop = FALSE], corners[c(1L, 3L, 4L), , drop = FALSE])
      for (tri in triangles) {
        facets[[facet_index]] = list(
          coords = tri,
          mean_prob = mean(tri[, 3]) / z_scale,
          mean_v2 = mean(tri[, 2])
        )
        facet_index = facet_index + 1L
      }
    }
  }
  if (length(facets) > 0) {
    draw_order = order(vapply(facets, function(facet) facet$mean_v2, numeric(1)), decreasing = TRUE)
    for (facet in facets[draw_order]) {
      projected = trans3d(facet$coords[, 1], facet$coords[, 2], facet$coords[, 3], transform)
      color_value = pmin(pmax(facet$mean_prob, 0), 1)^0.72
      fill = cols[pmax(1, pmin(101, floor(color_value * 100) + 1))]
      polygon(projected, col = adjustcolor(fill, alpha.f = 0.82), border = NA)
    }
  }

  bary_xy = function(x) c(v1 = x[2] + 0.5 * x[3], v2 = sqrt(3) / 2 * x[3])
  grid_families = list(
    function(a) cbind(x1 = a, x2 = seq(0, 1 - a, length.out = 150L), x3 = rev(seq(0, 1 - a, length.out = 150L))),
    function(a) cbind(x1 = seq(0, 1 - a, length.out = 150L), x2 = a, x3 = rev(seq(0, 1 - a, length.out = 150L))),
    function(a) cbind(x1 = seq(0, 1 - a, length.out = 150L), x2 = rev(seq(0, 1 - a, length.out = 150L)), x3 = a)
  )
  # Add ternary grid lines on the lifted surface.
  for (family in grid_families) {
    for (a in seq(0.1, 0.9, by = 0.1)) {
      comp = family(a)
      comp = pmax(comp, 1e-8)
      comp = comp / rowSums(comp)
      xy = t(apply(comp, 1L, bary_xy))
      colnames(xy) = c("v1", "v2")
      z = scenario_truth(scenario_id, comp)$prob * z_scale
      keep = is.finite(z)
      xy = xy[keep, , drop = FALSE]
      z = z[keep]
      if (nrow(xy) < 2L) {
        next
      }
      for (k in seq_len(nrow(xy) - 1L)) {
        segment_prob = mean(z[k:(k + 1L)]) / z_scale
        segment_col = cols[pmax(1, pmin(101, floor(segment_prob * 100) + 1))]
        projected = trans3d(
          xy[k:(k + 1L), "v1"],
          xy[k:(k + 1L), "v2"],
          z[k:(k + 1L)] + 0.004,
          transform
        )
        lines(projected, col = adjustcolor(segment_col, alpha.f = 0.88), lwd = 0.42)
      }
    }
  }

  base = trans3d(
    c(0, 1, 0.5, 0),
    c(0, 0, sqrt(3) / 2, 0),
    c(0, 0, 0, 0),
    transform
  )
  lines(base, col = "#222222", lwd = 1.45)
}

pdf(file.path(args$figures_dir, "fig01_simulation_truth.pdf"), width = 5.95, height = 9.60)
old_par = par(no.readonly = TRUE)
layout(
  rbind(
    c(1, 2),
    c(3, 4),
    c(5, 6),
    c(7, 8),
    c(9, 9)
  ),
  widths = c(1.02, 1.35),
  heights = c(1, 1, 1, 1, 0.18)
)
par(mar = c(0.3, 0.35, 1.35, 0.35), oma = c(0, 0, 0, 0))

for (i in seq_along(scenarios)) {
  panel = truth_grid[truth_grid$scenario_id == scenarios[i], , drop = FALSE]
  prob = pmin(pmax(panel$true_prob, 0), 1)
  point_cols = cols[pmax(1, pmin(101, floor(prob * 100) + 1))]
  draw_simplex_frame(paste0(scenario_titles[i], " (2D)"))
  points(panel$v1, panel$v2, pch = 15, cex = 0.24, col = point_cols)
  add_simplex_grid_lines(col = adjustcolor("white", alpha.f = 0.58), lwd = 0.42)
  lines(c(0, 1), c(0, 0), lwd = 1.45, col = "#1F1F1F")
  lines(c(1, 0.5), c(0, sqrt(3) / 2), lwd = 1.45, col = "#1F1F1F")
  lines(c(0.5, 0), c(sqrt(3) / 2, 0), lwd = 1.45, col = "#1F1F1F")

  draw_surface_3d(panel, scenarios[i], paste0(scenario_titles[i], " (3D)"))
}

par(mar = c(0.15, 0, 0.05, 0))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
xleft = 0.28
xright = 0.72
ybottom = 0.34
ytop = 0.58
xs = seq(xleft, xright, length.out = length(cols) + 1)
for (i in seq_along(cols)) {
  rect(xs[i], ybottom, xs[i + 1], ytop, col = cols[i], border = NA, xpd = NA)
}
text((xleft + xright) / 2, 0.80, "true p(X)", cex = 0.68)
text(xleft, 0.12, "0", cex = 0.62)
text(xright, 0.12, "1", cex = 0.62)

par(old_par)
dev.off()
message("Wrote simulation data to: ", args$out_dir)
