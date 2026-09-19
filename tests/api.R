library(sgev)
near <- function(a, b, tolerance = 1e-8) {
  stopifnot(length(a) == length(b), max(abs(a - b)) <= tolerance * max(1, abs(a), abs(b)))
}
fails <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
fails(sgev_start(2, c(1, 1)))
fails(sgev_start(2, 1, sigma = NA_real_))
fails(sgev_start(2, 1, alpha = 1))
fails(sgev_start(2, 1, weight = c(0.5, 0.5)))
state <- sgev_start(2, 1)
fails(sgev_update(state, matrix(c(1, 2, 2, 1), 2), c(0, 0)))
fails(sgev_update(state, diag(3), rep(0, 3)))
fails(sgev_update(state, matrix(0, 2, 2), c(1, 0)))
fit0 <- sgev_fit(state)
near(sgev_box(fit0)$log_evalue, 0)
stopifnot(!sgev_box(fit0)$reject, sgev_contrast(fit0, 1)$decision == "unresolved")

# Nuisance means cancel before accumulation, even with alternating designs.
matrices <- list()
for (b in 1:8) {
  rho <- (-1)^b * 0.8
  G <- 300 * matrix(c(1, rho, rho, 1), 2)
  score <- drop(G %*% c(0.12, 100 * sin(b)))
  previous <- state
  state <- sgev_update(state, G, score)
  stopifnot(previous$cohorts == b - 1L)
  near(state$H, b * 300 * (1 - 0.8^2))
  near(state$u, drop(state$H) * 0.12)
  matrices[[b]] <- G
}
fit <- sgev_fit(state)
H <- drop(state$H); u <- drop(state$u); lambda <- 100
expected <- -log1p(H / lambda) / 2 + max(abs(u) - H * 0.05, 0)^2 / (2 * (H + lambda))
near(sgev_box(fit)$log_evalue, expected)
info <- sgev_information(state$H, 0.12)
near(info$lower, H * (0.12 - 0.05)^2 / 2)
near(info$upper, H * (0.12 - 0.05)^2 / 2)
partition <- sgev_partition(matrices, 1)
near(partition$pooled - partition$local, partition$loss)

# Output order is preserved when all predictors are targets.
G <- matrix(c(3, 0.2, 0.2, 1), 2)
permuted <- sgev_partition(list(G), c(2, 1))
near(permuted$local, G[c(2, 1), c(2, 1)])

# The same local score history must produce the same mixture decisions.
direct <- sgev_update(sgev_start(1, 1), state$H, state$u)
near(sgev_box(sgev_fit(direct))$log_evalue, sgev_box(fit)$log_evalue)

# Singular targets: a mean of aliases is estimable, an individual coefficient is not.
alias <- sgev_update(sgev_start(2, 1:2), 1000 * matrix(1, 2, 2), c(500, 500))
alias_fit <- sgev_fit(alias)
stopifnot(sgev_contrast(alias_fit, c(1, 0))$decision == "unresolved")
stopifnot(sgev_contrast(alias_fit, c(0.5, 0.5))$decision == "non-negligible")
tiny <- sgev_start(2, 1)
fails(sgev_update(tiny, diag(c(1, 1e-14)), c(0, 0)))

# Boundary contact remains unresolved on the strict decision scale.
near(sgev_sgpv(0.05, 0.1), 0)
near(sgev_sgpv(-0.05, 0.05), 1)
stopifnot(is.na(sgev_sgpv(-Inf, Inf)), is.na(sgev_sgpv(0, 0)))
fails(sgev_sgpv(c(0, 1), 2))
fails(sgev_contrast(fit, 0))

# A supplied support function gives a lower convex e-value.
coordinates <- drop(crossprod(fit$basis, fit$center))
convex <- sgev_convex(fit, coordinates, function(v) 0.05 * sum(abs(v)))
stopifnot(convex$log_evalue <= sgev_box(fit)$log_evalue + 1e-8)
zero <- sgev_convex(fit0, numeric(), function(v) 0)
near(zero$log_evalue, 0)
fails(sgev_convex(fit, coordinates, function(v) NA_real_))

# Fixed-time mean under a boundary null, using reproducible Gaussian draws.
set.seed(129)
values <- replicate(2000, {
  H <- matrix(80, 1, 1)
  s <- 80 * 0.05 + sqrt(80) * rnorm(1)
  f <- sgev_fit(sgev_update(sgev_start(1, 1), H, s))
  exp(sgev_box(f)$log_evalue)
})
stopifnot(mean(values) < 1 + 5 * sd(values) / sqrt(length(values)))
cat("PASS: package API, scalar formulas, nuisance drift, target order, singular information, profiles, dual bounds, and null mean\n")
