# Shared simulation summaries and scalar comparator formulas.
binomial_interval <- function(rate, reps, level = .95) {
  count <- round(rate * reps)
  tail <- (1 - level) / 2
  cbind(lower = ifelse(count == 0, 0, qbeta(tail, count, reps - count + 1)),
        upper = ifelse(count == reps, 1, qbeta(1 - tail, count + 1, reps - count)))
}

paired_summary <- function(x, reference) {
  d <- as.numeric(x) - as.numeric(reference)
  c(difference = mean(d), mcse = sd(d) / sqrt(length(d)))
}

scalar_box_log_e <- function(score, information, delta, tau = .1) {
  lambda <- 1 / tau^2
  .5 * (pmax(abs(score / information) - delta, 0)^2 *
          information^2 / (information + lambda) - log1p(information / lambda))
}

stitched_boundary <- function(information, alpha = .05) {
  stopifnot(all(information >= 1), alpha > 0, alpha < 1)
  epoch <- floor(log2(information))
  allocation <- alpha / 2 * 6 / (pi^2 * (epoch + 1)^2)
  tilt <- sqrt(2 * log(1 / allocation) / 2^(epoch + 1))
  log(1 / allocation) / tilt + tilt * information / 2
}

mean_path_crossing <- function(gap, alpha = .05, tau = .1) {
  stopifnot(gap > 0)
  lambda <- 1 / tau^2
  f <- function(log_n) {
    n <- exp(log_n)
    n * gap^2 * n / (n + lambda) - log1p(n / lambda) - 2 * log(1 / alpha)
  }
  upper <- max(10, -4 * log(gap))
  while (f(upper) < 0) upper <- upper + 5
  exp(uniroot(f, c(-20, upper), tol = 1e-11)$root)
}
