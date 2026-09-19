# Fixed-alpha operating characteristics. No tuning uses performance paths.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(26090981)
reps <- if (getOption("sgev.simulation.mode") == "smoke") 200L else 5000L
alpha <- .05; delta <- .05; horizon <- 4000000L
times <- sort(unique(pmin(horizon, ceiling(100 * 1.02^(0:536)))))
stopifnot(tail(times, 1) == horizon)
gaps <- c(.0025, .005, .01, .02, .04, .08)
effects <- c(zero = 0, boundary_p = delta, boundary_m = -delta,
             setNames(delta + gaps, paste0("gap_", gaps)),
             negative_01 = -delta - .01, negative_04 = -delta - .04)
scales <- .1 * 2^(-6:3)
methods <- c("Fixed scale", "Scale mixture", "Stitched boundary")
first <- array(NA_integer_, c(reps, length(effects), length(methods)))
noise <- numeric(reps); previous <- 0L
for (n in times) {
  noise <- noise + rnorm(reps, sd = sqrt(n - previous))
  for (ei in seq_along(effects)) {
    score <- noise + n * effects[ei]
    fixed <- scalar_box_log_e(score, n, delta)
    logs <- vapply(scales, function(tau) scalar_box_log_e(score, n, delta, tau), numeric(reps))
    largest <- apply(logs, 1, max)
    mixed <- largest + log(rowMeans(exp(logs - largest)))
    decisions <- cbind(fixed > log(1 / alpha), mixed > log(1 / alpha),
      abs(score) - n * delta > stitched_boundary(n, alpha))
    for (mi in seq_along(methods)) {
      take <- is.na(first[, ei, mi]) & decisions[, mi]
      first[take, ei, mi] <- n
    }
  }
  previous <- n
}
results <- paired <- list()
for (ei in seq_along(effects)) for (mi in seq_along(methods)) {
  crossed <- !is.na(first[, ei, mi]); count <- sum(crossed)
  tt <- ifelse(crossed, first[, ei, mi], horizon)
  results[[length(results) + 1L]] <- data.frame(
    scenario = names(effects)[ei], beta = effects[ei],
    gap = abs(effects[ei]) - delta, null = abs(effects[ei]) <= delta,
    method = methods[mi], crossing_rate = mean(crossed),
    crossing_mcse = sd(crossed) / sqrt(reps),
    crossing_upper95 = if (count == reps) 1 else qbeta(.95, count + 1, reps - count),
    mean_truncated_n = mean(tt), mean_truncated_mcse = sd(tt) / sqrt(reps),
    censored_fraction = mean(!crossed), reps = reps, horizon = horizon,
    looks = length(times), alpha = alpha, delta = delta)
  if (mi > 1) {
    ref <- first[, ei, 1]
    dp <- paired_summary(crossed, !is.na(ref))
    dn <- paired_summary(tt, ifelse(is.na(ref), horizon, ref))
    paired[[length(paired) + 1L]] <- data.frame(scenario = names(effects)[ei],
      method = methods[mi], power_difference = dp[1], power_mcse = dp[2],
      mean_n_difference = dn[1], mean_n_mcse = dn[2], reps = reps, row.names = NULL)
  }
}
write.csv(do.call(rbind, results), "results/near_margin.csv", row.names = FALSE)
write.csv(do.call(rbind, paired), "results/near_margin_paired.csv", row.names = FALSE)
saveRDS(list(first = first, methods = methods, effects = effects, scales = scales,
  times = times, seed = 26090981, horizon = horizon),
  "results/near_margin_replicates.rds", version = 3)
gap_grid <- 10^(-seq(1, 12, length.out = 100))
nstar <- vapply(gap_grid, mean_path_crossing, numeric(1))
write.csv(data.frame(gap = gap_grid, mean_path_n = nstar,
  normalized = nstar * gap_grid^2 / (2 * log(1 / gap_grid))),
  "results/near_margin_mean_path.csv", row.names = FALSE)
cat("Small-gap study complete: 11 targets, 3 methods,", reps, "paths each.\n")
