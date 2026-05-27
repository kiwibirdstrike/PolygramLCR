# Build the triangulation, hat-basis matrix, and edge-contrast matrix used by
# the polygram fit. This is the geometry setup shared by every lambda value.
polygram_setup = function(
  train,
  centers,
  seed,
  normalize_constraint_rows = TRUE,
  extra_vertices = NULL,
  fixed_vertices = NULL
) {
  predictors = as.matrix(train[, c("z1", "z2")])
  y = as.numeric(train$y)

  # Use supplied vertices for reproducible checks; otherwise build them from
  # k-means interior centers plus the convex hull of the observed ILR points.
  vertices = if (!is.null(fixed_vertices)) {
    fixed_vertices = as.matrix(fixed_vertices)
    colnames(fixed_vertices) = colnames(predictors)
    unique(round(fixed_vertices, digits = 12))
  } else {
    make_polygram_vertices(predictors, centers = centers, seed = seed, extra_vertices = extra_vertices)
  }

  setup = list(
    response = y,
    predictors = predictors,
    vertex = vertices,
    sample_size = length(y),
    number_vertex = nrow(vertices)
  )

  # RTriangle returns the triangular cells and their edges. Interior shared
  # edges define the local plane contrasts penalized by the polygram.
  setup$triangulate = RTriangle::triangulate(RTriangle::pslg(setup$vertex))
  setup$triangle = setup$triangulate$T
  setup$number_triangle = nrow(setup$triangle)
  setup$triangle[, 2:3] = setup$triangle[, 3:2]

  setup$common_edge = setup$triangulate$E[setup$triangulate$EB == 0, ]
  setup$boundary_edge = setup$triangulate$E[setup$triangulate$EB == 1, ]
  setup$common_edge = matrix(setup$common_edge, ncol = 2)
  setup$number_common_edge = nrow(setup$common_edge)
  setup$number_boundary_edge = nrow(setup$boundary_edge)

  # For each vertex, store the triangles that contain it. The vertex is rotated
  # into the first position so the hat-basis C++ routine has a stable layout.
  setup$star_vertex = list()
  for (j in seq_len(setup$number_vertex)) {
    index = which(rowSums(setup$triangle == j) == 1)
    setup$star_vertex[[j]] = matrix(setup$triangle[index, ], ncol = 3)
    for (k in seq_len(nrow(setup$star_vertex[[j]]))) {
      if (setup$star_vertex[[j]][k, 2] == j) {
        setup$star_vertex[[j]][k, ] = setup$star_vertex[[j]][k, c(2, 3, 1)]
      }
      if (setup$star_vertex[[j]][k, 3] == j) {
        setup$star_vertex[[j]][k, ] = setup$star_vertex[[j]][k, c(3, 1, 2)]
      }
    }
  }

  # For each interior edge, store the two neighboring triangles. These pairs are
  # used below to build one contrast row per shared edge.
  setup$star_common_edge = list()
  setup$ce_index = NULL
  for (j in seq_len(setup$number_common_edge)) {
    ce = setup$common_edge[j, ]
    index = which(rowSums(matrix(setup$triangle %in% ce, ncol = 3)) == 2)
    setup$star_common_edge[[j]] = setup$triangle[index, ]
    setup$ce_index = rbind(setup$ce_index, index)
  }

  # Boundary edges are retained for plotting and active-edge filtering, but they
  # do not define plane-continuity contrasts.
  setup$star_boundary_edge = list()
  for (j in seq_len(setup$number_boundary_edge)) {
    be = setup$boundary_edge[j, ]
    index = which(rowSums(matrix(setup$triangle %in% be, ncol = 3)) == 2)
    setup$star_boundary_edge[[j]] = setup$triangle[index, ]
  }

  setup$dimension = length(setup$star_vertex)
  setup$constraint = matrix(0, setup$dimension, setup$number_common_edge)

  # Edge contrast c_e'theta compares the slave triangle's fitted vertex value
  # with the value predicted there by the affine plane on the master triangle.
  for (j in seq_len(setup$number_common_edge)) {
    master_triangle = setup$star_common_edge[[j]][1, ]
    slave_triangle = setup$star_common_edge[[j]][2, ]
    target_index = slave_triangle[!slave_triangle %in% setup$common_edge[j, ]]
    b = rect2bary(
      setup$vertex[target_index, ],
      setup$vertex[master_triangle, ]
    )
    setup$constraint[master_triangle, j] = b
    setup$constraint[target_index, j] = -1
  }

  # B is the fitted hat-basis design matrix. basis_vertex is kept for geometry
  # diagnostics and downstream plotting.
  setup$basis = hat_basis_linear(predictors, vertices, setup$star_vertex)
  setup$basis_vertex = hat_basis_linear(vertices, vertices, setup$star_vertex)

  # D has one row per interior edge and maps theta to normalized plane contrasts.
  D_raw = t(as.matrix(setup$constraint))
  constraint_row_scale = rep(1, nrow(D_raw))
  if (normalize_constraint_rows && nrow(D_raw) > 0) {
    constraint_row_scale = sqrt(rowSums(D_raw^2))
  }

  list(
    B = as.matrix(setup$basis),
    D = D_raw / constraint_row_scale,
    y = y,
    setup = setup,
    vertices = vertices,
    constraint_row_scale = constraint_row_scale,
    normalize_constraint_rows = normalize_constraint_rows
  )
}


# Choose polygram vertices from observed ILR points: interior k-means centers
# plus hull vertices, with optional extra vertices for fixed grids or boundaries.
make_polygram_vertices = function(predictors, centers, seed = NULL, extra_vertices = NULL) {
  predictors = as.matrix(predictors)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  centers = as.integer(centers)
  center_vertices = if (centers >= nrow(predictors)) {
    predictors
  } else {
    stats::kmeans(predictors, centers = centers, nstart = 10, iter.max = 100)$centers
  }
  hull_vertices = predictors[grDevices::chull(predictors), , drop = FALSE]
  if (!is.null(extra_vertices)) {
    extra_vertices = as.matrix(extra_vertices)
    colnames(extra_vertices) = colnames(predictors)
    vertices = rbind(center_vertices, hull_vertices, extra_vertices)
  } else {
    vertices = rbind(center_vertices, hull_vertices)
  }
  vertices = unique(round(vertices, digits = 12))
  colnames(vertices) = colnames(predictors)
  vertices
}


# Numerically stable binomial negative log-likelihood for eta = B theta.
polygram_logistic_nll = function(B, y, theta) {
  eta = as.numeric(B %*% theta)
  log_partition = ifelse(eta > 0, eta + log1p(exp(-eta)), log1p(exp(eta)))
  sum(log_partition - y * eta)
}


# Inner ADMM update for theta at fixed w and u. The logistic part is smooth but
# non-quadratic, so this block solves the theta subproblem by Newton iterations.
polygram_admm_theta_update = function(
  B,
  y,
  D,
  theta,
  w,
  u,
  rho,
  DtD,
  max_newton_iter = 50,
  newton_tol = 1e-6,
  hessian_epsilon = 0,
  max_line_search = 30
) {
  converged = FALSE
  error_message = ""

  # Smooth ADMM theta subproblem:
  #   L(theta) + (rho / 2) ||D theta - w + u||_2^2
  # This is the objective checked during the Newton line search.
  augmented_objective = function(candidate_theta) {
    penalty_residual = as.numeric(D %*% candidate_theta) - w + u
    polygram_logistic_nll(B, y, candidate_theta) + 0.5 * rho * sum(penalty_residual^2)
  }

  for (iter in seq_len(max_newton_iter)) {
    # Newton system for the local quadratic approximation:
    #   H d = -g,
    # where
    #   g = B' (p - y) + rho D' (D theta - w + u)
    #   H = B' diag(p(1-p)) B + rho D'D.
    eta = as.numeric(B %*% theta)
    prob = stats::plogis(eta)
    logistic_weight = prob * (1 - prob)
    penalty_residual = as.numeric(D %*% theta) - w + u
    grad = as.numeric(crossprod(B, prob - y) + rho * crossprod(D, penalty_residual))
    hessian = crossprod(B, B * logistic_weight) + rho * DtD

    if (hessian_epsilon > 0) {
      diag(hessian) = diag(hessian) + hessian_epsilon
    }
    direction = tryCatch(
      as.numeric(solve(hessian, -grad)),
      error = function(e) {
        error_message <= conditionMessage(e)
        NULL
      }
    )

    if (is.null(direction) || any(!is.finite(direction))) {
      error_message = if (nzchar(error_message)) {
        paste("Newton linear solve failed:", error_message)
      } else {
        "Newton linear solve returned non-finite direction."
      }
      break
    }

    grad_direction = sum(grad * direction)
    if (!is.finite(grad_direction)) {
      error_message = "Newton direction has non-finite directional derivative."
      break
    }

    if (grad_direction > 1e-10) {
      error_message = "Newton direction is not a descent direction for the theta subproblem."
      break
    }

    # Newton-decrement stopping rule for the inner theta update.
    newton_decrement = -grad_direction
    if (newton_decrement / 2 <= newton_tol) {
      converged = TRUE
      break
    }

    # Backtracking line search. Start with the full Newton step and halve
    # the step size until the Armijo decrease condition is satisfied.
    current_obj = augmented_objective(theta)
    step = 1
    accepted = FALSE
    for (line_iter in seq_len(max_line_search)) {
      candidate = theta + step * direction
      candidate_obj = augmented_objective(candidate)
      if (is.finite(candidate_obj) && candidate_obj <= current_obj + 1e-4 * step * grad_direction) {
        theta = candidate
        accepted = TRUE
        break
      }
      step = step * 0.5
    }
    if (!accepted) {
      error_message = "Backtracking line search failed in the theta subproblem."
      break
    }
  }

  if (!converged && !nzchar(error_message)) {
    error_message = "Newton theta subproblem did not converge within max_newton_iter."
  }

  list(
    theta = theta,
    converged = converged,
    iterations = iter,
    error_message = error_message
  )
}


# Fit one point on the lambda path. This is the outer ADMM loop:
# theta-update, w soft-threshold update, u dual update, then convergence check.
polygram_admm_fit_one_lambda = function(
  prepared,
  lambda = 0.02,
  rho = NULL,
  rho_scale = 0.1,
  rho_min = 1e-3,
  max_admm_iter = 500,
  max_newton_iter = 50,
  abstol = 1e-4,
  reltol = 1e-3,
  newton_tol = 1e-6,
  hessian_epsilon = 0,
  require_newton_convergence = TRUE,
  convergence_criterion = c("objective", "residual"),
  objective_tol = 1e-5,
  theta_init = NULL,
  w_init = NULL,
  u_init = NULL,
  verbose = FALSE,
  lambda_index = NA_integer_
) {
  B = prepared$B
  D = prepared$D
  y = prepared$y
  n = nrow(B)
  p = ncol(B)
  k = nrow(D)
  convergence_criterion = match.arg(convergence_criterion)
  if (is.null(rho)) {
    rho = max(lambda * rho_scale, rho_min)
  }
  if (!is.finite(rho) || rho <= 0) {
    stop("rho must be positive.", call. = FALSE)
  }

  # Warm starts pass the previous lambda's theta, w, and u into the next fit.
  theta = if (is.null(theta_init)) rep(0, p) else theta_init
  w = if (is.null(w_init)) as.numeric(D %*% theta) else w_init
  u = if (is.null(u_init)) rep(0, k) else u_init
  DtD = crossprod(D)

  history = vector("list", max_admm_iter)
  converged = FALSE
  objval_previous = Inf

  for (iter in seq_len(max_admm_iter)) {
    # 1. theta-update: smooth logistic subproblem solved by Newton steps.
    theta_update = polygram_admm_theta_update(
      B = B,
      y = y,
      D = D,
      theta = theta,
      w = w,
      u = u,
      rho = rho,
      DtD = DtD,
      max_newton_iter = max_newton_iter,
      newton_tol = newton_tol,
      hessian_epsilon = hessian_epsilon
    )
    if (require_newton_convergence && !isTRUE(theta_update$converged)) {
      stop(
        "ADMM theta subproblem failed at lambda_index=",
        lambda_index,
        ", iteration=",
        iter,
        ": ",
        theta_update$error_message,
        call. = FALSE
      )
    }

    theta = theta_update$theta

    # 2. w-update: soft-threshold the current edge contrasts D theta + u.
    w_old = w
    Dtheta = as.numeric(D %*% theta)
    w_input = Dtheta + u
    w = sign(w_input) * pmax(abs(w_input) - lambda / rho, 0)

    # 3. u-update: scaled ADMM dual variable.
    u = u + Dtheta - w

    # Diagnostics used for convergence checks, path summaries, and AIC.
    r_norm = sqrt(sum((Dtheta - w)^2))
    s_norm = sqrt(sum((rho * as.numeric(crossprod(D, w - w_old)))^2))
    eps_pri = sqrt(k) * abstol + reltol * max(sqrt(sum(Dtheta^2)), sqrt(sum(w^2)))
    eps_dual = sqrt(p) * abstol + reltol * sqrt(sum((rho * as.numeric(crossprod(D, u)))^2))
    nll = polygram_logistic_nll(B, y, theta)
    penalty = lambda * sum(abs(Dtheta))
    objval = nll + penalty
    active_jumps = sum(abs(w) > 1e-7)
    common_edges = sum(abs(w) <= 1e-7)
    effective_dimension = max(1, p - common_edges)
    history[[iter]] = data.frame(
      iteration = iter,
      objval = objval,
      nll = nll,
      penalty = penalty,
      r_norm = r_norm,
      s_norm = s_norm,
      eps_pri = eps_pri,
      eps_dual = eps_dual,
      active_jumps = active_jumps,
      common_edges = common_edges,
      effective_dimension = effective_dimension,
      objective_change = abs(objval - objval_previous),
      theta_newton_iterations = theta_update$iterations,
      theta_newton_converged = theta_update$converged,
      theta_newton_error = theta_update$error_message,
      stringsAsFactors = FALSE
    )

    if (verbose && (iter == 1 || iter %% 10 == 0)) {
      message(
        "lambda_index=", lambda_index,
        " lambda=", signif(lambda, 6),
        "ADMM iter=", iter,
        " obj=", signif(objval, 6),
        " r=", signif(r_norm, 4),
        " s=", signif(s_norm, 4)
      )
    }
    residual_converged = r_norm < eps_pri && s_norm < eps_dual
    objective_converged = abs(objval - objval_previous) < objective_tol
    if (
      (identical(convergence_criterion, "residual") && residual_converged) ||
        (identical(convergence_criterion, "objective") && objective_converged)
    ) {
      converged = TRUE
      break
    }
    objval_previous = objval
  }

  history = do.call(rbind, history[seq_len(iter)])
  eta = as.numeric(B %*% theta)
  loglik = -polygram_logistic_nll(B, y, theta)
  df_beta = sum(abs(theta) > 1e-7)
  df_jumps = sum(abs(w) > 1e-7)
  common_edges = sum(abs(w) <= 1e-7)
  df_effective = max(1, p - common_edges)

  # beta/z are retained as aliases for older scripts; theta/w are the names that
  # match the paper notation.
  list(
    theta = theta,
    w = w,
    beta = theta,
    z = w,
    u = u,
    lambda = lambda,
    rho = rho,
    history = history,
    converged = converged,
    iterations = iter,
    max_admm_iter = max_admm_iter,
    fitted_prob = stats::plogis(eta),
    loglik = loglik,
    aic = -2 * loglik + 2 * df_effective,
    bic = -2 * loglik + log(n) * df_effective,
    df_beta = df_beta,
    df_jumps = df_jumps,
    common_edges = common_edges,
    df_effective = df_effective,
    convergence_criterion = convergence_criterion,
    objective_tol = objective_tol,
    method = "polygram"
  )
}


# Main user-facing fit. It prepares the geometry once, fits every lambda value,
# and returns the AIC/BIC-selected polygram model with the full path attached.
polygram = function(
  train,
  centers = 70,
  number_lambdas = 20,
  lambda_max = 0.2,
  lambda_max_min_ratio = 1e-2,
  lambdas = NULL,
  selection = c("bic", "aic"),
  rho = NULL,
  rho_scale = 0.1,
  rho_min = 1e-3,
  max_admm_iter = 500,
  max_newton_iter = 50,
  abstol = 1e-4,
  reltol = 1e-3,
  newton_tol = 1e-6,
  hessian_epsilon = 0,
  require_newton_convergence = TRUE,
  convergence_criterion = c("objective", "residual"),
  objective_tol = 1e-5,
  require_admm_convergence = TRUE,
  normalize_constraint_rows = TRUE,
  extra_vertices = NULL,
  fixed_vertices = NULL,
  warm_start = TRUE,
  seed = 1,
  verbose = FALSE
) {
  selection = match.arg(selection)
  convergence_criterion = match.arg(convergence_criterion)
  prepared = polygram_setup(
    train = train,
    centers = centers,
    seed = seed,
    normalize_constraint_rows = normalize_constraint_rows,
    extra_vertices = extra_vertices,
    fixed_vertices = fixed_vertices
  )

  # Default path decreases geometrically from lambda_max.
  if (is.null(lambdas)) {
    lambdas = if (number_lambdas <= 1) {
      lambda_max
    } else {
      exp(seq(log(lambda_max), log(lambda_max * lambda_max_min_ratio), length.out = number_lambdas))
    }
  }
  lambdas = as.numeric(lambdas)
  rho_path = rho

  if (is.null(rho_path)) {
    rho_path = max(max(lambdas) * rho_scale, rho_min)
  }

  fit_list = vector("list", length(lambdas))
  path_summary = vector("list", length(lambdas))
  theta_init = w_init = u_init = NULL

  # Fit the path sequentially. Warm starts make each lambda use the previous
  # solution as its initial ADMM state.
  for (i in seq_along(lambdas)) {
    fit_i = polygram_admm_fit_one_lambda(
      prepared = prepared,
      lambda = lambdas[[i]],
      rho = rho_path,
      rho_scale = rho_scale,
      rho_min = rho_min,
      max_admm_iter = max_admm_iter,
      max_newton_iter = max_newton_iter,
      abstol = abstol,
      reltol = reltol,
      newton_tol = newton_tol,
      hessian_epsilon = hessian_epsilon,
      require_newton_convergence = require_newton_convergence,
      convergence_criterion = convergence_criterion,
      objective_tol = objective_tol,
      theta_init = if (warm_start) theta_init else NULL,
      w_init = if (warm_start) w_init else NULL,
      u_init = if (warm_start) u_init else NULL,
      verbose = verbose,
      lambda_index = i
    )
    fit_list[[i]] = fit_i
    last = fit_i$history[nrow(fit_i$history), , drop = FALSE]
    path_summary[[i]] = data.frame(
      lambda_index = i,
      lambda = fit_i$lambda,
      rho = fit_i$rho,
      converged = fit_i$converged,
      iterations = fit_i$iterations,
      loglik = fit_i$loglik,
      aic = fit_i$aic,
      bic = fit_i$bic,
      df_beta = fit_i$df_beta,
      df_jumps = fit_i$df_jumps,
      common_edges = fit_i$common_edges,
      df_effective = fit_i$df_effective,
      final_objval = last$objval,
      final_r_norm = last$r_norm,
      final_s_norm = last$s_norm,
      final_eps_pri = last$eps_pri,
      final_eps_dual = last$eps_dual,
      final_objective_change = last$objective_change,
      stringsAsFactors = FALSE
    )
    theta_init = fit_i$theta
    w_init = fit_i$w
    u_init = fit_i$u
  }

  path_summary = do.call(rbind, path_summary)
  selection_values = path_summary[[selection]]

  # AIC/BIC selection can optionally ignore lambda values that did not meet the
  # ADMM convergence criterion.
  if (require_admm_convergence) {
    selection_values[!path_summary$converged] = Inf
  }
  if (!any(is.finite(selection_values))) {
    stop(
      "No ADMM lambda fit satisfied the convergence criterion.",
      call. = FALSE
    )
  }

  selected_index = which.min(ifelse(is.finite(selection_values), selection_values, Inf))
  selected = fit_list[[selected_index]]
  selected$setup = prepared$setup
  selected$vertices = prepared$vertices
  selected$constraint_row_scale = prepared$constraint_row_scale
  selected$normalize_constraint_rows = prepared$normalize_constraint_rows
  selected$selected_index = selected_index
  selected$selection = selection
  selected$require_admm_convergence = require_admm_convergence
  selected$convergence_criterion = convergence_criterion
  selected$objective_tol = objective_tol
  selected$lambdas = lambdas
  selected$path = fit_list
  selected$path_summary = path_summary
  selected$method = "polygram"
  selected
}

# Compact one-row summary used by scripts and result tables.
summarize_polygram = function(model) {
  last = model$history[nrow(model$history), , drop = FALSE]
  data.frame(
    selected_index = if (is.null(model$selected_index)) NA_integer_ else model$selected_index,
    selection = if (is.null(model$selection)) NA_character_ else model$selection,
    lambda = model$lambda,
    rho = model$rho,
    converged = model$converged,
    iterations = model$iterations,
    final_objval = last$objval,
    final_r_norm = last$r_norm,
    final_s_norm = last$s_norm,
    final_eps_pri = last$eps_pri,
    final_eps_dual = last$eps_dual,
    final_objective_change = last$objective_change,
    convergence_criterion = if (is.null(model$convergence_criterion)) NA_character_ else model$convergence_criterion,
    objective_tol = if (is.null(model$objective_tol)) NA_real_ else model$objective_tol,
    active_jumps = last$active_jumps,
    common_edges = if (is.null(model$common_edges)) last$common_edges else model$common_edges,
    effective_dimension = if (is.null(model$df_effective)) last$effective_dimension else model$df_effective,
    df_beta = if (is.null(model$df_beta)) sum(abs(model$theta) > 1e-7) else model$df_beta,
    df_jumps = if (is.null(model$df_jumps)) last$active_jumps else model$df_jumps,
    aic = model$aic,
    bic = model$bic,
    min_path_aic = min(model$path_summary$aic),
    min_path_bic = min(model$path_summary$bic),
    n_lambdas = length(model$lambdas),
    n_lambdas_converged = sum(model$path_summary$converged),
    n_lambdas_at_max_iterations = sum(
      model$path_summary$iterations >= if (is.null(model$max_admm_iter)) model$iterations else model$max_admm_iter
    ),
    number_vertices = nrow(model$vertices),
    number_triangles = model$setup$number_triangle,
    number_common_edges = model$setup$number_common_edge,
    normalized_constraint_rows = model$normalize_constraint_rows,
    stringsAsFactors = FALSE
  )
}

# Predict probabilities at new ILR coordinates by evaluating the fitted hat
# basis and applying the logistic link.
predict_polygram = function(model, newdata) {
  Z = as.matrix(newdata[, c("z1", "z2")])
  basis = hat_basis_linear(Z, model$vertices, model$setup$star_vertex)
  eta = as.numeric(basis %*% model$theta)
  stats::plogis(eta)
}

# Gradient of one fitted affine plane over a triangle in ILR coordinates.
polygram_triangle_gradient = function(vertices, theta, triangle) {
  X = cbind(1, vertices[triangle, , drop = FALSE])
  coef = qr.coef(qr(X), theta[triangle])
  c(z1 = coef[[2]], z2 = coef[[3]])
}


# Edge strength used for active-edge displays: Euclidean jump between the two
# neighboring triangle gradients that share an interior edge.
polygram_edge_gradient_jumps = function(fit) {
  vertices = as.matrix(fit$vertices)
  edges = fit$setup$common_edge
  theta = if (is.null(fit$theta)) fit$beta else fit$theta
  vapply(seq_len(nrow(edges)), function(edge_index) {
    triangles = fit$setup$star_common_edge[[edge_index]]
    g1 = polygram_triangle_gradient(vertices, theta, as.integer(triangles[1, ]))
    g2 = polygram_triangle_gradient(vertices, theta, as.integer(triangles[2, ]))
    sqrt(sum((g1 - g2)^2))
  }, numeric(1))
}


# Active edges are nonzero fitted edge contrasts after optional trimming of
# simplex-boundary-adjacent edges for display.
polygram_active_edge_mask = function(fit, active_epsilon = 1e-7, edge_trim = 0.015) {
  vertices_comp = ilr_to_comp(as.matrix(fit$vertices))
  edges = fit$setup$common_edge
  edge_contrast = if (is.null(fit$w)) fit$z else fit$w
  interior = rowSums(vertices_comp[edges[, 1], , drop = FALSE] > edge_trim) == 3 &
    rowSums(vertices_comp[edges[, 2], , drop = FALSE] > edge_trim) == 3
  interior & abs(edge_contrast) > active_epsilon
}


# Table form of all shared edges, their contrast magnitudes, gradient jumps, and
# whether they pass the display-oriented active/interior filter.
polygram_active_edge_summary = function(fit, scenario_id = NA_character_, active_epsilon = 1e-7, edge_trim = 0.015) {
  edges = fit$setup$common_edge
  edge_contrast = if (is.null(fit$w)) fit$z else fit$w
  data.frame(
    scenario_id = scenario_id,
    edge_index = seq_len(nrow(edges)),
    v_start = edges[, 1],
    v_end = edges[, 2],
    active_contrast = abs(edge_contrast),
    gradient_jump = polygram_edge_gradient_jumps(fit),
    active_interior = polygram_active_edge_mask(
      fit,
      active_epsilon = active_epsilon,
      edge_trim = edge_trim
    ),
    stringsAsFactors = FALSE
  )
}
