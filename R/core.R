# Generated from the manuscript's tested method kernel.

psd_parts <- function(G, tol = 1e-10) {
  G <- as.matrix(G)
  stopifnot(nrow(G) == ncol(G), all(is.finite(G)))
  if (max(abs(G - t(G))) > tol * max(1, max(abs(G)))) stop("G is not symmetric")
  e <- eigen((G + t(G)) / 2, symmetric = TRUE)
  scale <- max(1, max(e$values))
  if (min(e$values) < -tol * scale) stop("G is not positive semidefinite")
  keep <- e$values > tol * scale
  U <- e$vectors[, keep, drop = FALSE]
  d <- e$values[keep]
  pinv <- if (length(d)) tcrossprod(sweep(U, 2, sqrt(d), "/")) else G * 0
  P <- tcrossprod(U)
  # The retained eigenspace is only a computational projection. Keep the full
  # spectrum separately: projection must not remove determinant penalties.
  list(U = U, d = d, pinv = pinv, P = P, rank = sum(keep), tol = tol,
       full_U = e$vectors, full_d = e$values)
}

joint_geometry <- function(G, sigma = 1, tau = 0.1, alpha = 0.05,
                           tol = 1e-10) {
  stopifnot(sigma > 0, tau > 0, alpha > 0, alpha < 1)
  z <- psd_parts(G, tol)
  lambda <- sigma^2 / tau^2
  V <- if (z$rank) {
    tcrossprod(sweep(z$U, 2, sigma * sqrt(z$d + lambda) / z$d, "*"))
  } else G * 0
  if (any(z$full_d + lambda <= 0)) {
    stop("Regularized information is not positive definite")
  }
  D <- sum(log1p(z$full_d / lambda))
  list(G = G, sigma = sigma, tau = tau, alpha = alpha, V = V,
       D = D, B = 2 * log(1 / alpha) + D, parts = z)
}

joint_fit <- function(G, s, sigma = 1, tau = 0.1, alpha = 0.05,
                      tol = 1e-10) {
  z <- joint_geometry(G, sigma, tau, alpha, tol)
  stopifnot(length(s) == nrow(G), all(is.finite(s)))
  # Do not reject a score merely because it has a discarded target component.
  # The full score is preserved; only the lower quadratic bound is projected.
  z$center <- drop(z$parts$pinv %*% s)
  z$s <- s
  z
}

mixture_precision <- function(fit) {
  # Precision of the projected quadratic lower bound, not a replacement for
  # the full-information mixture. Its determinant penalty remains unprojected.
  z <- fit$parts
  if (!z$rank) return(fit$V * 0)
  precision <- z$d^2 / (fit$sigma^2 * (z$d + fit$sigma^2 / fit$tau^2))
  tcrossprod(sweep(z$U, 2, sqrt(precision), "*"))
}

guard_dual <- function(raw, primal, scale, label) {
  # A dual objective is a lower bound in exact arithmetic. Subtract a
  # roundoff allowance before using it. If the guarded value still exceeds a
  # primal upper bound by more than numerical tolerance, stop. Zero is always
  # dual feasible and is the conservative fallback for a residual tiny breach.
  floating_guard <- 256 * .Machine$double.eps * pmax(1, scale)
  guarded <- raw - floating_guard
  bad <- guarded > primal + 1e-6 * pmax(1, primal)
  if (any(bad)) stop(label, " numerical duality check failed")
  guarded[guarded > primal] <- 0
  list(value = pmax(0, guarded), floating_guard = floating_guard)
}

joint_log_e <- function(G, s, beta, sigma = 1, tau = 0.1) {
  z <- joint_geometry(G, sigma, tau)
  stopifnot(length(s) == nrow(G), length(beta) == nrow(G),
            all(is.finite(s)), all(is.finite(beta)))
  r <- drop(s - G %*% beta)
  a <- drop(crossprod(z$parts$full_U, r))
  -z$D / 2 + sum(a^2 / (z$parts$full_d + sigma^2 / tau^2)) / (2 * sigma^2)
}

contrast_estimable <- function(fit, contrast) {
  # No residual tolerance can certify an original contrast with an unrestricted
  # kernel coefficient. Use full retained rank or exact Gram-column identities.
  p <- length(fit$center)
  stopifnot(length(contrast) == p, all(is.finite(contrast)))
  if (fit$parts$rank == p || all(contrast == 0)) return(TRUE)
  G <- as.matrix(fit$G)
  active <- which(colSums(abs(G)) > 0)
  if (any(contrast[setdiff(seq_len(p), active)] != 0)) return(FALSE)
  if (!length(active)) return(FALSE)
  signs <- vapply(active, function(j) sign(G[which(G[, j] != 0)[1], j]),
                  numeric(1))
  canonical <- sweep(G[, active, drop = FALSE], 2, signs, "*")
  representatives <- which(!duplicated(t(canonical)))
  # Additional dependencies or discarded positive directions have no certified
  # contrast range in this interface. Return an unbounded interval in that case.
  if (length(representatives) != fit$parts$rank) return(FALSE)
  for (j in seq_along(active)) {
    k <- representatives[vapply(representatives, function(k) {
      all(canonical[, j] == canonical[, k])
    }, logical(1))][1]
    if (contrast[active[j]] * signs[j] != contrast[active[k]] * signs[k]) {
      return(FALSE)
    }
  }
  TRUE
}

contrast_variance <- function(fit, contrast) {
  stopifnot(length(contrast) == length(fit$center), all(is.finite(contrast)))
  coordinates <- drop(crossprod(fit$parts$U, contrast))
  sum(coordinates^2 * fit$sigma^2 *
        (fit$parts$d + fit$sigma^2 / fit$tau^2) / fit$parts$d^2)
}

contrast_interval <- function(fit, contrast) {
  stopifnot(length(contrast) == length(fit$center), all(is.finite(contrast)))
  if (!contrast_estimable(fit, contrast)) {
    return(c(lower = -Inf, upper = Inf, center = NA, radius = Inf))
  }
  m <- sum(contrast * fit$center)
  w <- sqrt(max(0, fit$B * contrast_variance(fit, contrast)))
  c(lower = m - w, upper = m + w, center = m, radius = w)
}

profile_log_e <- function(fit, contrast, delta) {
  stopifnot(length(delta) == 1, delta >= 0, sum(contrast^2) > 0)
  ci <- contrast_interval(fit, contrast)
  if (!is.finite(ci["radius"])) return(c(non_negligible = -fit$D / 2,
                                         negligible = -fit$D / 2))
  v <- contrast_variance(fit, contrast)
  m <- abs(ci["center"])
  c(non_negligible = unname(-fit$D / 2 + pmax(m - delta, 0)^2 / (2 * v)),
    negligible = unname(-fit$D / 2 + pmax(delta - m, 0)^2 / (2 * v)))
}

sgpv <- function(lower, upper, delta) {
  stopifnot(delta > 0, all(lower <= upper))
  width <- upper - lower
  ans <- rep(NA_real_, length(width))
  ok <- is.finite(width) & width > 0
  overlap <- pmax(0, pmin(upper[ok], delta) - pmax(lower[ok], -delta))
  ans[ok] <- overlap / width[ok] * pmax(width[ok] / (4 * delta), 1)
  # The numerical SGPV is not assigned to an unbounded or degenerate interval.
  ans
}

decisions <- function(center, radius, delta) {
  list(non_negligible = abs(center) - radius > delta,
       negligible = abs(center) + radius < delta)
}

nuisance_inverse <- function(G, outside, tol = 1e-10) {
  # Only exact zero, duplicate, or sign-duplicate full Gram columns are removed.
  # A small positive nuisance eigenvalue must never be silently set to zero:
  # an unrestricted outside coefficient can amplify the resulting mean bias.
  G <- as.matrix(G)
  q <- length(outside)
  answer <- matrix(0, q, q)
  if (!q) return(list(inverse = answer, representatives = integer()))
  columns <- G[, outside, drop = FALSE]
  nonzero <- which(colSums(abs(columns)) > 0)
  if (!length(nonzero)) {
    return(list(inverse = answer, representatives = integer()))
  }
  canonical <- columns[, nonzero, drop = FALSE]
  signs <- vapply(seq_len(ncol(canonical)), function(j) {
    sign(canonical[which(canonical[, j] != 0)[1], j])
  }, numeric(1))
  canonical <- sweep(canonical, 2, signs, "*")
  representatives <- nonzero[!duplicated(t(canonical))]
  D <- G[outside[representatives], outside[representatives], drop = FALSE]
  parts <- psd_parts(D, tol)
  if (parts$rank != nrow(D)) {
    stop("Unresolved nuisance rank: no e-value is returned. Supply a design-justified full-rank nuisance representation; do not truncate small positive nuisance eigenvalues.")
  }
  answer[representatives, representatives] <- chol2inv(chol(D))
  list(inverse = answer, representatives = outside[representatives])
}

local_group_increment <- function(G, s, group, tol = 1e-10) {
  # Efficient score for beta_group after eliminating every coefficient outside
  # the group. The Schur complement is formed within each independent cohort.
  G <- as.matrix(G); p <- nrow(G)
  stopifnot(ncol(G) == p, length(s) == p, length(group) > 0,
            !anyDuplicated(group), all(group %in% seq_len(p)),
            all(is.finite(G)), all(is.finite(s)))
  if (max(abs(G - t(G))) > tol * max(1, max(abs(G)))) {
    stop("G is not symmetric")
  }
  outside <- setdiff(seq_len(p), group)
  if (!length(outside)) {
    H <- G[group, group, drop = FALSE]
    u <- s[group]
  } else {
    z <- nuisance_inverse(G, outside, tol)
    cross <- G[group, outside, drop = FALSE]
    H <- G[group, group, drop = FALSE] - cross %*% z$inverse %*%
      G[outside, group, drop = FALSE]
    u <- s[group] - drop(cross %*% z$inverse %*% s[outside])
  }
  H <- (H + t(H)) / 2
  psd_parts(H, tol) # Validate without projecting either returned quantity.
  list(H = H, u = u, group = group)
}

local_group_update <- function(state, G, s, group, tol = 1e-10) {
  inc <- local_group_increment(G, s, group, tol)
  if (is.null(state)) {
    state <- list(H = inc$H * 0, u = inc$u * 0, group = group)
  }
  stopifnot(identical(state$group, group))
  list(H = state$H + inc$H, u = state$u + inc$u, group = group)
}

heterogeneous_information <- function(G_list, group, tol = 1e-10) {
  # Singular version of the pooled-minus-local Schur identity. Each nuisance
  # block may be rank deficient. The identity is checked before return.
  stopifnot(length(G_list) > 0)
  G_list <- lapply(G_list, as.matrix)
  p <- nrow(G_list[[1]])
  stopifnot(all(vapply(G_list, function(G) {
    nrow(G) == p && ncol(G) == p
  }, logical(1))), length(group) > 0, !anyDuplicated(group),
  all(group %in% seq_len(p)))
  outside <- setdiff(seq_len(p), group)
  k <- length(group)
  if (!length(outside)) {
    zero <- matrix(0, k, k)
    summed <- Reduce(`+`, G_list)[group, group, drop = FALSE]
    return(list(local = summed, pooled = summed,
                loss = zero,
                K = replicate(length(G_list), matrix(0, 0, k), simplify = FALSE),
                K_pooled = matrix(0, 0, k)))
  }
  pieces <- lapply(G_list, function(G) {
    psd_parts(G, tol)
    D <- G[outside, outside, drop = FALSE]
    C <- G[outside, group, drop = FALSE]
    K <- nuisance_inverse(G, outside, tol)$inverse %*% C
    H <- G[group, group, drop = FALSE] - t(C) %*% K
    list(D = D, C = C, K = K, H = (H + t(H)) / 2)
  })
  D <- Reduce(`+`, lapply(pieces, `[[`, "D"))
  C <- Reduce(`+`, lapply(pieces, `[[`, "C"))
  K_pooled <- nuisance_inverse(Reduce(`+`, G_list), outside, tol)$inverse %*% C
  local <- Reduce(`+`, lapply(pieces, `[[`, "H"))
  G_pooled <- Reduce(`+`, G_list)
  pooled <- local_group_increment(G_pooled, rep(0, p), group, tol)$H
  loss <- Reduce(`+`, lapply(pieces, function(x) {
    dK <- x$K - K_pooled
    t(dK) %*% x$D %*% dK
  }))
  scale <- max(1, max(abs(pooled)), max(abs(local)), max(abs(loss)))
  if (max(abs(pooled - local - loss)) > 1e-7 * scale) {
    stop("Singular heterogeneous-information identity failed")
  }
  list(local = local, pooled = pooled, loss = loss,
       K = lapply(pieces, `[[`, "K"), K_pooled = K_pooled)
}

local_group_fit <- function(state, sigma = 1, tau = 0.1, alpha = 0.05,
                            weight = 1, tol = 1e-10) {
  stopifnot(weight > 0, weight <= 1)
  fit <- joint_fit(state$H, state$u, sigma, tau, alpha * weight, tol)
  fit$group <- state$group
  fit$weight <- weight
  fit
}

box_certificate <- function(fit, delta, tolerance = 1e-8) {
  # Computes the composite practical-null e-value. The feasible dual value
  # gives a conservative e-value, so optimization error cannot create a false
  # certificate. Singular information, including exact collinearity, is permitted.
  stopifnot(delta > 0, length(fit$center) > 0)
  m <- fit$center
  W <- mixture_precision(fit)
  precision <- fit$parts$d^2 /
    (fit$sigma^2 * (fit$parts$d + fit$sigma^2 / fit$tau^2))
  objective <- function(z) {
    coordinates <- drop(crossprod(fit$parts$U, z - m))
    sum(coordinates^2 * precision)
  }
  gradient <- function(z) drop(2 * W %*% (z - m))
  opt <- optim(pmax(-delta, pmin(delta, m)), objective, gradient,
               method = "L-BFGS-B", lower = rep(-delta, length(m)),
               upper = rep(delta, length(m)),
               control = list(factr = 1e5, pgtol = 1e-11, maxit = 2000))
  displacement_coordinates <- drop(crossprod(fit$parts$U, m - opt$par))
  u_coordinates <- precision * displacement_coordinates
  u <- drop(fit$parts$U %*% u_coordinates)
  primal <- objective(opt$par)
  complementarity_slack <- opt$par - delta * sign(u)
  raw_dual <- primal + 2 * sum(u * complementarity_slack)
  guarded <- guard_dual(raw_dual, primal,
                        1 + abs(primal) +
                          2 * sum(abs(u * complementarity_slack)),
                        "Box")
  dual <- guarded$value
  conservative_log_evalue <- (dual - fit$D) / 2
  primal_log_evalue <- (primal - fit$D) / 2
  list(conservative_evalue = exp(conservative_log_evalue),
       conservative_log_evalue = conservative_log_evalue,
       primal_evalue = exp(primal_log_evalue),
       primal_log_evalue = primal_log_evalue,
       certified = dual > fit$B + tolerance * max(1, fit$B),
       dual = dual, primal = primal, gap = max(0, primal - dual),
       floating_guard = guarded$floating_guard,
       convergence = opt$convergence,
       function_evaluations = unname(opt$counts["function"]),
       gradient_evaluations = unname(opt$counts["gradient"]))
}

pair_box_distance <- function(center, W, delta, V = NULL) {
  # Exact active-set solution for many two-dimensional centers with one common
  # positive-semidefinite precision matrix. The returned dual is a lower bound.
  center <- as.matrix(center); W <- as.matrix(W)
  stopifnot(ncol(center) == 2, all(dim(W) == c(2, 2)), delta > 0)
  psd_parts(W)
  B <- nrow(center)
  clamp <- function(x) pmax(-delta, pmin(delta, x))
  candidates <- list(cbind(clamp(center[, 1]), clamp(center[, 2])))
  scale <- max(1, max(abs(W)))
  for (a in c(-delta, delta)) {
    z2 <- if (W[2, 2] > 1e-12 * scale) {
      clamp(center[, 2] - W[1, 2] / W[2, 2] * (a - center[, 1]))
    } else clamp(center[, 2])
    candidates[[length(candidates) + 1]] <- cbind(rep(a, B), z2)
  }
  for (a in c(-delta, delta)) {
    z1 <- if (W[1, 1] > 1e-12 * scale) {
      clamp(center[, 1] - W[1, 2] / W[1, 1] * (a - center[, 2]))
    } else clamp(center[, 1])
    candidates[[length(candidates) + 1]] <- cbind(z1, rep(a, B))
  }
  values <- vapply(candidates, function(z) {
    d <- center - z
    rowSums((d %*% W) * d)
  }, numeric(B))
  if (is.null(dim(values))) values <- matrix(values, nrow = B)
  best <- max.col(-values, ties.method = "first")
  z <- matrix(0, B, 2)
  for (j in seq_along(candidates)) z[best == j, ] <- candidates[[j]][best == j, ]
  d <- center - z
  u <- d %*% W
  if (!is.null(V)) {
    V <- as.matrix(V)
    stopifnot(all(dim(V) == c(2, 2)), all(is.finite(V)),
              max(abs(V - t(V))) <= 1e-8 * max(1, max(abs(V))))
    psd_parts(V)
  }
  primal <- rowSums((d %*% W) * d)
  complementarity_slack <- z - delta * sign(u)
  raw_dual <- primal + 2 * rowSums(u * complementarity_slack)
  guarded <- guard_dual(raw_dual, primal,
                        1 + abs(primal) +
                          2 * rowSums(abs(u * complementarity_slack)),
                        "Pair-box")
  dual <- guarded$value
  list(dual = dual, primal = primal, gap = pmax(0, primal - dual),
       floating_guard = guarded$floating_guard, minimizer = z)
}
