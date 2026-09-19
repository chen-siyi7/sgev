library(sgev)

near <- function(a, b, tolerance = 1e-9) {
  stopifnot(isTRUE(all.equal(unname(a), unname(b), tolerance = tolerance)))
}
fails <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))

fails(sgev_scalar(numeric(), 1))
fails(sgev_scalar(NA_real_, 1))
fails(sgev_scalar(1, 0))
fails(sgev_scalar(0, -1))
fails(sgev_scalar(c(0, 1), c(1, 2, 3)))
fails(sgev_scalar(0, 1, margin = 0))
fails(sgev_scalar(0, 1, sigma = 0))
fails(sgev_scalar(0, 1, tau = Inf))
fails(sgev_scalar(0, 1, alpha = 1))
zero <- sgev_scalar(0, 0)
stopifnot(zero$log_out == 0, zero$log_in == 0,
          zero$lower == -Inf, zero$upper == Inf,
          is.na(zero$sgpv), zero$decision == "unresolved")

# Independently evaluate the scalar formulas and compare the regression API.
cases <- 0L
for (sigma in c(.5, 1, 2)) for (tau in c(.01, .1, 1)) {
  for (alpha in c(.001, .05)) for (H in c(1, 100, 10000)) {
    m <- c(-2, -.2, -.05, 0, .05, .2, 2)
    result <- sgev_scalar(H * m, H, sigma = sigma, tau = tau, alpha = alpha)
    lambda <- sigma^2 / tau^2
    D <- log1p(H / lambda)
    v <- sigma^2 * (H + lambda) / H^2
    r <- sqrt(v * (2 * log(1 / alpha) + D))
    outward <- -D / 2 + pmax(abs(m) - .05, 0)^2 / (2 * v)
    inward <- -D / 2 + pmax(.05 - abs(m), 0)^2 / (2 * v)
    near(result$log_out, outward)
    near(result$log_in, inward)
    near(result$lower, m - r)
    near(result$upper, m + r)
    stopifnot(identical(result$meaningful, outward > log(1 / alpha)),
              identical(result$negligible, inward > log(1 / alpha)),
              !any(result$meaningful & result$negligible),
              all(result$sgpv >= 0 & result$sgpv <= 1))
    for (j in seq_along(m)) {
      state <- sgev_start(1, 1, sigma, tau, alpha)
      state <- sgev_update(state, matrix(H), H * m[j])
      contrast <- sgev_contrast(sgev_fit(state), 1)
      near(c(result$lower[j], result$upper[j]), contrast$interval[1:2])
      near(c(result$log_out[j], result$log_in[j]), contrast$log_evalues)
      stopifnot(result$decision[j] == contrast$decision)
      cases <- cases + 1L
    }
  }
}
heterogeneous <- sgev_scalar(c(0, 20, -600), c(0, 100, 5000))
stopifnot(nrow(heterogeneous) == 3, heterogeneous$decision[1] == "unresolved")
cat("PASS:", cases, "scalar formula and regression API comparisons\n")
