# Stopping-time comparison under predictable nuisance effects.
reps <- if (getOption("sgev.simulation.mode") == "smoke") 200L else 3000L
horizon <- 6000L
levels <- c(.05, .01, 1e-4, 1e-8, 1e-16, 1e-32)
tau <- .1
seed <- 26091129L
specs <- stopping_specs()
records <- summaries <- diagnostics <- list()
for (g in seq_along(specs)) {
  spec <- specs[[g]]
  set.seed(seed + g)
  parts <- eigen(spec$H, symmetric = TRUE)
  root <- parts$vectors %*% diag(sqrt(pmax(parts$values, 0)))
  zstar <- stopping_projection(spec)
  direction <- spec$a - zstar
  information <- drop(crossprod(direction, spec$H %*% direction)) / 2
  intercept <- drop(crossprod(spec$a, spec$H %*% spec$a) -
                    crossprod(zstar, spec$H %*% zstar)) / 2
  U <- recovered_U <- matrix(0, reps, 2)
  V <- numeric(reps)
  times <- array(NA_integer_, c(reps, length(levels), 2),
                 dimnames = list(NULL, NULL, c("Local mixture", "Supporting oracle")))
  max_error <- max_log_error <- 0
  for (b in seq_len(horizon)) {
    # Each policy uses its own path's past, before the current innovations.
    D <- 1 + .5 * tanh(V / (1 + b))
    c1 <- .7 * tanh(V / (1 + b))
    c2 <- ifelse(U[, 1] >= U[, 2], .6, -.6)
    loading <- cbind(c1, c2)
    nuisance <- .4 * sin(b / 9) + .3 * tanh(V / (1 + b))
    u <- sweep(matrix(rnorm(reps * 2), reps, 2) %*% t(root),
               2, drop(spec$H %*% spec$a), "+")
    outside <- D * (nuisance + drop(loading %*% spec$a)) + sqrt(D) * rnorm(reps)
    target <- u + loading * outside
    recovered <- target - loading * outside
    U <- U + u
    recovered_U <- recovered_U + recovered
    V <- V + outside
    max_error <- max(max_error, abs(recovered_U - U))
    log_e <- stopping_log_e(recovered_U, b, spec, tau)
    oracle <- drop(recovered_U %*% direction) - b * intercept
    # Coupling to the homogeneous local experiment must preserve all decisions.
    direct <- stopping_log_e(U, b, spec, tau)
    max_log_error <- max(max_log_error, abs(log_e - direct))
    for (k in seq_along(levels)) {
      boundary <- log(1 / levels[k])
      stopifnot(identical(log_e >= boundary, direct >= boundary))
      hit <- is.na(times[, k, 1]) & log_e >= boundary
      times[hit, k, 1] <- b
      hit <- is.na(times[, k, 2]) & oracle >= boundary
      times[hit, k, 2] <- b
    }
    if (!anyNA(times)) break
  }
  records[[spec$id]] <- list(spec = spec, times = times, last_cohort = b,
                             information = information, null_projection = zstar)
  diagnostics[[g]] <- data.frame(geometry = spec$id, seed = seed + g,
    max_score_error = max_error, max_log_e_error = max_log_error,
    simulated_cohorts = b, coupled_decision_disagreements = 0L)
  for (method in seq_len(2)) for (k in seq_along(levels)) {
    t <- times[, k, method]
    censored <- is.na(t)
    t[censored] <- horizon
    factor <- information / log(1 / levels[k])
    summaries[[length(summaries) + 1L]] <- data.frame(
      geometry = spec$id, method = dimnames(times)[[3]][method], level = levels[k],
      information = information, reps = reps, horizon = horizon,
      censored = sum(censored), mean_truncated_time = mean(t),
      time_mcse = sd(t) / sqrt(reps), median_time = median(t),
      ratio = mean(t) * factor, ratio_mcse = sd(t) / sqrt(reps) * factor)
  }
  message(spec$id, ": ", b, " cohorts; I = ", signif(information, 6))
}
dir.create("results", showWarnings = FALSE)
write.csv(do.call(rbind, summaries), "results/stopping_efficiency.csv", row.names = FALSE)
write.csv(do.call(rbind, diagnostics), "results/stopping_efficiency_diagnostics.csv", row.names = FALSE)
saveRDS(list(reps = reps, horizon = horizon, levels = levels, tau = tau,
             seed = seed, records = records), "results/stopping_efficiency_paths.rds", version = 3)
writeLines(capture.output(sessionInfo()), "results/stopping_efficiency_session.txt")
