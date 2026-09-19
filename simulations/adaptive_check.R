# Bounded-stopping likelihood identity under adaptive designs and nuisance.
# This diagnostic is separate from the rejection rules.
set.seed(26090921)
check_reps <- if (getOption("sgev.simulation.mode") == "smoke") 1000L else 30000L
max_cohorts <- 30L
nb <- 40
a <- .12
z <- .05
x <- a - z
active <- rep(TRUE, check_reps)
sum_B <- local_U <- local_H <- ll_full <- ll_local <- compensator <- numeric(check_reps)
stopped_at <- rep(max_cohorts, check_reps)
for (b in seq_len(max_cohorts)) {
  ids <- which(active)
  if (!length(ids)) break
  # Both quantities depend only on the common observed history.
  rho <- ifelse(sum_B[ids] >= 0, .85, -.2)
  gamma <- .25 * tanh(sum_B[ids] / 40)
  gamma0 <- gamma + rho * x
  Z <- matrix(rnorm(2 * length(ids)), length(ids), 2)
  noise_A <- sqrt(nb) * Z[, 1]
  noise_B <- sqrt(nb) * (rho * Z[, 1] + sqrt(1 - rho^2) * Z[, 2])
  sA <- nb * (a + rho * gamma) + noise_A
  sB <- nb * (rho * a + gamma) + noise_B
  null_A <- nb * (z + rho * gamma0)
  null_B <- nb * (rho * z + gamma0)
  hb <- nb * (1 - rho^2)
  ub <- sA - rho * sB
  ll_full[ids] <- ll_full[ids] +
    x * (sA - null_A) - rho * x * (sB - null_B) - .5 * x^2 * hb
  ll_local[ids] <- ll_local[ids] + x * ub - .5 * (a^2 - z^2) * hb
  compensator[ids] <- compensator[ids] + .5 * x^2 * hb
  local_U[ids] <- local_U[ids] + ub
  local_H[ids] <- local_H[ids] + hb
  sum_B[ids] <- sum_B[ids] + sB
  now_stop <- abs(sum_B[ids]) > 100 |
    abs(local_U[ids] - z * local_H[ids]) > 2.5 * sqrt(local_H[ids] + 1)
  stopped_at[ids[now_stop]] <- b
  active[ids[now_stop]] <- FALSE
}
martingale_part <- ll_full - compensator
out <- data.frame(reps = check_reps, horizon_cohorts = max_cohorts,
  mean_stopped_cohort = mean(stopped_at),
  mean_log_likelihood_ratio = mean(ll_full), mean_KL_compensator = mean(compensator),
  mean_difference = mean(martingale_part),
  difference_mcse = sd(martingale_part) / sqrt(check_reps),
  maximum_path_identity_error = max(abs(ll_full - ll_local)))
stopifnot(out$maximum_path_identity_error < 1e-9)
write.csv(out, "results/adaptive_drift_kl.csv", row.names = FALSE)
cat("Adaptive stopped-likelihood identity checked on", check_reps, "paths.\n")
