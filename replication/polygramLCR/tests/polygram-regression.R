library(polygramLCR)

X = matrix(
  c(
    0.62, 0.25, 0.13,
    0.48, 0.37, 0.15,
    0.34, 0.48, 0.18,
    0.24, 0.56, 0.20,
    0.18, 0.32, 0.50,
    0.36, 0.22, 0.42,
    0.52, 0.18, 0.30,
    0.28, 0.18, 0.54,
    0.16, 0.62, 0.22,
    0.44, 0.12, 0.44,
    0.58, 0.32, 0.10,
    0.22, 0.42, 0.36
  ),
  ncol = 3,
  byrow = TRUE
)

Z = comp_to_ilr(X)
eta = -0.35 + 1.25 * Z[, 1] - 0.85 * Z[, 2] + 0.9 * (Z[, 1] > 0)
train = data.frame(
  z1 = Z[, 1],
  z2 = Z[, 2],
  y = as.integer(eta > stats::median(eta))
)

vertices = unique(round(rbind(Z, Z[grDevices::chull(Z), , drop = FALSE]), digits = 12))
colnames(vertices) = c("z1", "z2")

fit = polygram(
  train,
  fixed_vertices = vertices,
  lambdas = c(0.4, 0.12, 0.04),
  selection = "aic",
  rho = 0.2,
  max_admm_iter = 80,
  max_newton_iter = 40,
  objective_tol = 1e-7,
  require_admm_convergence = TRUE,
  hessian_epsilon = 1e-8,
  warm_start = TRUE
)

expected_beta = c(
  9.239705116048512,
  2.146505507807661,
  -4.213228137101902,
  -9.138664098203730,
  -4.163766950815061,
  7.170349089094360,
  12.789285742412010,
  7.413257194491569,
  -12.927046308168750,
  16.582144592129200,
  -2.535356676431193,
  -5.908211042119404
)

expected_prob = c(
  0.9999029032076311,
  0.8953417780057881,
  0.01458271687885732,
  0.0001074192104979026,
  0.01531080973703815,
  0.9992315367739387,
  0.9999972094965570,
  0.9993971608915776,
  0.000002431386859001622,
  0.9999999371267256,
  0.0734164172273972,
  0.002709680897939644
)

expected_path_aic = c(
  6.550908740413956,
  3.405223149370209,
  2.442430104874161
)

stopifnot(identical(fit$selected_index, 3L))
stopifnot(isTRUE(all.equal(fit$theta, fit$beta, tolerance = 0)))
stopifnot(isTRUE(all.equal(fit$w, fit$z, tolerance = 0)))
stopifnot(isTRUE(all.equal(fit$beta, expected_beta, tolerance = 1e-8)))
stopifnot(isTRUE(all.equal(fit$fitted_prob, expected_prob, tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(fit$path_summary$aic, expected_path_aic, tolerance = 1e-10)))
stopifnot(isTRUE(all(fit$path_summary$converged)))
