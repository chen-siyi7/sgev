.scalar <- function(x, name, lower = 0, upper = Inf,
                    lower_closed = FALSE, upper_closed = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (ok) ok <- if (lower_closed) x >= lower else x > lower
  if (ok) ok <- if (upper_closed) x <= upper else x < upper
  if (!ok) stop("Invalid ", name, call. = FALSE)
  invisible(x)
}

.target <- function(target, p) {
  if (!is.numeric(target) || !length(target) || anyNA(target) ||
      any(!is.finite(target)) || any(target != as.integer(target)) ||
      any(target < 1 | target > p) || anyDuplicated(target)) {
    stop("target must contain distinct predictor indices", call. = FALSE)
  }
  as.integer(target)
}

.vector <- function(x, n, name) {
  if (!is.numeric(x) || length(x) != n || any(!is.finite(x))) {
    stop(name, " must be a finite numeric vector of length ", n, call. = FALSE)
  }
  as.numeric(x)
}

.information <- function(x, p, tolerance) {
  if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(p, p)) ||
      any(!is.finite(x))) stop("Invalid information matrix", call. = FALSE)
  psd_parts(x, tolerance)
  (x + t(x)) / 2
}

.fit <- function(fit) {
  if (!inherits(fit, "sgev_fit")) stop("Use a fit returned by sgev_fit()", call. = FALSE)
  invisible(fit)
}

sgev_start <- function(p, target, sigma = 1, tau = 0.1, alpha = 0.05,
                       weight = 1, tolerance = 1e-10) {
  .scalar(p, "p")
  if (p != as.integer(p)) stop("p must be a positive integer", call. = FALSE)
  p <- as.integer(p)
  target <- .target(target, p)
  .scalar(sigma, "sigma"); .scalar(tau, "tau")
  .scalar(alpha, "alpha", upper = 1)
  .scalar(weight, "weight", upper = 1, upper_closed = TRUE)
  .scalar(tolerance, "tolerance", upper = 1)
  k <- length(target)
  structure(list(p = p, target = target, H = matrix(0, k, k), u = numeric(k),
                 cohorts = 0L, sigma = sigma, tau = tau, alpha = alpha,
                 weight = weight, tolerance = tolerance), class = "sgev_state")
}

sgev_update <- function(object, information, score) {
  if (!inherits(object, "sgev_state")) stop("Use sgev_start() first", call. = FALSE)
  G <- .information(information, object$p, object$tolerance)
  s <- .vector(score, object$p, "score")
  zero <- colSums(abs(G)) == 0
  if (any(s[zero] != 0)) stop("A zero-information coordinate has nonzero score", call. = FALSE)
  increment <- local_group_increment(G, s, object$target, object$tolerance)
  object$H <- object$H + increment$H
  object$u <- object$u + increment$u
  object$cohorts <- object$cohorts + 1L
  object
}

sgev_fit <- function(object) {
  if (!inherits(object, "sgev_state")) stop("Use sgev_start() first", call. = FALSE)
  fit <- local_group_fit(list(H = object$H, u = object$u, group = object$target),
                         object$sigma, object$tau, object$alpha,
                         object$weight, object$tolerance)
  fit$basis <- fit$parts$U
  fit$rank <- fit$parts$rank
  fit$cohorts <- object$cohorts
  fit$target <- object$target
  class(fit) <- "sgev_fit"
  fit
}

sgev_box <- function(fit, margin = 0.05) {
  .fit(fit); .scalar(margin, "margin")
  result <- box_certificate(fit, margin)
  list(log_evalue = result$conservative_log_evalue,
       evalue = result$conservative_evalue, reject = result$certified,
       threshold = 1 / fit$alpha, primal = result$primal,
       dual = result$dual, gap = result$gap, convergence = result$convergence)
}

sgev_contrast <- function(fit, contrast, margin = 0.05) {
  .fit(fit); .scalar(margin, "margin")
  contrast <- .vector(contrast, length(fit$center), "contrast")
  if (all(contrast == 0)) stop("contrast must be nonzero", call. = FALSE)
  interval <- contrast_interval(fit, contrast)
  profiles <- profile_log_e(fit, contrast, margin)
  decision <- "unresolved"
  if (is.finite(interval["radius"])) {
    flags <- decisions(unname(interval["center"]), unname(interval["radius"]), margin)
    if (flags$non_negligible) decision <- "non-negligible"
    if (flags$negligible) decision <- "negligible"
  }
  list(interval = interval, log_evalues = profiles, decision = decision,
       sgpv = unname(sgpv(interval["lower"], interval["upper"], margin)))
}

sgev_convex <- function(fit, coordinates, support) {
  .fit(fit)
  coordinates <- .vector(coordinates, fit$rank, "coordinates")
  if (!is.function(support)) stop("support must be a function", call. = FALSE)
  direction <- drop(fit$basis %*% coordinates)
  h <- support(direction)
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h)) {
    stop("support must return one finite value", call. = FALSE)
  }
  if (all(direction == 0) && h < 0) stop("Invalid support value at zero", call. = FALSE)
  terms <- c(2 * sum(direction * fit$center),
             -sum(direction * drop(fit$V %*% direction)), -2 * h)
  allowance <- 256 * .Machine$double.eps * max(1, sum(abs(terms)))
  dual <- max(0, sum(terms) - allowance)
  log_evalue <- (dual - fit$D) / 2
  list(log_evalue = log_evalue, evalue = exp(log_evalue), dual = dual,
       direction = direction, reject = dual > fit$B + 1e-8 * max(1, fit$B))
}

sgev_sgpv <- function(lower, upper, margin = 0.05) {
  .scalar(margin, "margin")
  if (!is.numeric(lower) || !is.numeric(upper) || !length(lower) ||
      length(lower) != length(upper) || anyNA(lower) || anyNA(upper) ||
      any(lower > upper)) stop("Invalid interval endpoints", call. = FALSE)
  sgpv(lower, upper, margin)
}

sgev_information <- function(information, alternative, margin = 0.05,
                             sigma = 1, tolerance = 1e-10) {
  .scalar(margin, "margin"); .scalar(sigma, "sigma")
  .scalar(tolerance, "tolerance", upper = 1)
  p <- length(alternative)
  if (!p) stop("alternative must be nonempty", call. = FALSE)
  H <- .information(information, as.integer(p), tolerance)
  alternative <- .vector(alternative, p, "alternative")
  W <- H / sigma^2
  objective <- function(z) sum((z - alternative) * drop(W %*% (z - alternative)))
  gradient <- function(z) drop(2 * W %*% (z - alternative))
  opt <- stats::optim(pmax(-margin, pmin(margin, alternative)), objective, gradient,
                      method = "L-BFGS-B", lower = rep(-margin, p), upper = rep(margin, p),
                      control = list(factr = 1e5, pgtol = 1e-11, maxit = 2000))
  direction <- drop(W %*% (alternative - opt$par))
  primal <- objective(opt$par)
  slack <- opt$par - margin * sign(direction)
  raw <- primal + 2 * sum(direction * slack)
  dual <- guard_dual(raw, primal, 1 + abs(primal) + 2 * sum(abs(direction * slack)),
                     "Information")$value
  list(lower = dual / 2, upper = primal / 2, gap = max(0, primal - dual) / 2,
       minimizer = opt$par, convergence = opt$convergence)
}

sgev_partition <- function(information_list, target, tolerance = 1e-10) {
  .scalar(tolerance, "tolerance", upper = 1)
  if (!is.list(information_list) || !length(information_list)) {
    stop("information_list must be a nonempty list", call. = FALSE)
  }
  p <- nrow(information_list[[1L]])
  .scalar(p, "matrix dimension")
  p <- as.integer(p)
  target <- .target(target, p)
  matrices <- lapply(information_list, .information, p = p, tolerance = tolerance)
  heterogeneous_information(matrices, target, tolerance)
}
