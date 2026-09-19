# Vectorized scalar profiles for reduced Gaussian scores.
sgev_scalar <- function(score, information, margin = 0.05, sigma = 1,
                        tau = 0.1, alpha = 0.05) {
  .scalar(margin, "margin")
  .scalar(sigma, "sigma")
  .scalar(tau, "tau")
  .scalar(alpha, "alpha", upper = 1)
  if (!is.numeric(score) || !length(score) || any(!is.finite(score))) {
    stop("score must be a nonempty finite numeric vector", call. = FALSE)
  }
  if (!is.numeric(information) ||
      !(length(information) %in% c(1L, length(score))) ||
      any(!is.finite(information)) || any(information < 0)) {
    stop("information must be nonnegative and scalar or match score", call. = FALSE)
  }
  H <- rep(information, length.out = length(score))
  positive <- H > 0
  if (any(score[!positive] != 0)) {
    stop("Zero information requires zero score", call. = FALSE)
  }
  estimate <- rep(NA_real_, length(score))
  radius <- rep(Inf, length(score))
  log_out <- log_in <- numeric(length(score))
  lower <- rep(-Inf, length(score))
  upper <- rep(Inf, length(score))
  lambda <- sigma^2 / tau^2
  if (!is.finite(lambda) || lambda <= 0) {
    stop("sigma and tau produce an unrepresentable mixture scale", call. = FALSE)
  }
  if (any(positive)) {
    h <- H[positive]
    m <- score[positive] / h
    D <- log1p(h / lambda)
    v <- sigma^2 * (1 + lambda / h) / h
    r <- sqrt(v * (2 * log(1 / alpha) + D))
    estimate[positive] <- m
    radius[positive] <- r
    lower[positive] <- m - r
    upper[positive] <- m + r
    log_out[positive] <- -D / 2 + pmax(abs(m) - margin, 0)^2 / (2 * v)
    log_in[positive] <- -D / 2 + pmax(margin - abs(m), 0)^2 / (2 * v)
  }
  # Decisions use strict interval separation. Equality at a margin is unresolved.
  meaningful <- lower > margin | upper < -margin
  negligible <- lower > -margin & upper < margin
  decision <- rep("unresolved", length(score))
  decision[meaningful] <- "non-negligible"
  decision[negligible] <- "negligible"
  data.frame(estimate = estimate, radius = radius, lower = lower, upper = upper,
             log_out = log_out, log_in = log_in,
             meaningful = meaningful, negligible = negligible,
             decision = decision, sgpv = sgev_sgpv(lower, upper, margin),
             stringsAsFactors = FALSE)
}
