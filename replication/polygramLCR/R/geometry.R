check_composition = function(X, eps = 1e-12) {
  X = as.matrix(X)
  colnames(X) = c("x1", "x2", "x3")
  X
}

comp_to_ilr = function(X) {
  X = check_composition(X)
  x1 = X[, 1]
  x2 = X[, 2]
  x3 = X[, 3]
  out = cbind(
    z1 = sqrt(1 / 2) * log(x1 / x2),
    z2 = sqrt(2 / 3) * log(sqrt(x1 * x2) / x3)
  )
  out
}

ilr_to_comp = function(Z) {
  Z = as.matrix(Z)
  z1 = Z[, 1]
  z2 = Z[, 2]
  g1 = exp(z1 / sqrt(2) + z2 / sqrt(6))
  g2 = exp(-z1 / sqrt(2) + z2 / sqrt(6))
  g3 = exp(-2 * z2 / sqrt(6))
  s = g1 + g2 + g3
  out = cbind(x1 = g1 / s, x2 = g2 / s, x3 = g3 / s)
  out
}

comp_to_bary = function(X) {
  X = check_composition(X)
  out = cbind(
    v1 = X[, 2] + 0.5 * X[, 3],
    v2 = (sqrt(3) / 2) * X[, 3]
  )
  out
}

bary_to_comp = function(V, eps = 1e-12) {
  V = as.matrix(V)
  x3 = 2 * V[, 2] / sqrt(3)
  x2 = V[, 1] - 0.5 * x3
  x1 = 1 - x2 - x3
  out = cbind(x1 = x1, x2 = x2, x3 = x3)
  out[out < eps] = eps
  out / rowSums(out)
}

make_simplex_grid = function(grid_n = 120, eps = 1e-10) {
  v1_seq = seq(0, 1, length.out = grid_n)
  v2_seq = seq(0, sqrt(3) / 2, length.out = grid_n)
  grid = expand.grid(v1 = v1_seq, v2 = v2_seq)
  ymax = pmin(sqrt(3) * grid$v1, sqrt(3) * (1 - grid$v1))
  keep = grid$v2 >= 0 & grid$v2 <= ymax
  bary = as.matrix(grid[keep, , drop = FALSE])
  comp = bary_to_comp(bary, eps = eps)
  ilr = comp_to_ilr(comp)
  data.frame(
    v1 = bary[, 1],
    v2 = bary[, 2],
    x1 = comp[, 1],
    x2 = comp[, 2],
    x3 = comp[, 3],
    z1 = ilr[, 1],
    z2 = ilr[, 2]
  )
}

draw_simplex_frame = function(title = NULL) {
  plot(
    NA,
    xlim = c(-0.05, 1.05),
    ylim = c(-0.05, sqrt(3) / 2 + 0.05),
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = if (is.null(title)) "" else title
  )
  lines(c(0, 1), c(0, 0), lwd = 2)
  lines(c(1, 0.5), c(0, sqrt(3) / 2), lwd = 2)
  lines(c(0.5, 0), c(sqrt(3) / 2, 0), lwd = 2)
}

add_simplex_grid_lines = function(at = seq(0.1, 0.9, by = 0.1), col = "#D6D6D6", lwd = 0.35) {
  bary_xy = function(x) c(v1 = x[2] + 0.5 * x[3], v2 = sqrt(3) / 2 * x[3])
  for (a in at) {
    p1 = bary_xy(c(a, 1 - a, 0))
    p2 = bary_xy(c(a, 0, 1 - a))
    lines(c(p1["v1"], p2["v1"]), c(p1["v2"], p2["v2"]), col = col, lwd = lwd)

    p1 = bary_xy(c(1 - a, a, 0))
    p2 = bary_xy(c(0, a, 1 - a))
    lines(c(p1["v1"], p2["v1"]), c(p1["v2"], p2["v2"]), col = col, lwd = lwd)

    p1 = bary_xy(c(1 - a, 0, a))
    p2 = bary_xy(c(0, 1 - a, a))
    lines(c(p1["v1"], p2["v1"]), c(p1["v2"], p2["v2"]), col = col, lwd = lwd)
  }
}

safe_log_loss = function(y, prob, eps = 1e-12) {
  prob = pmin(pmax(prob, eps), 1 - eps)
  -mean(y * log(prob) + (1 - y) * log(1 - prob))
}
