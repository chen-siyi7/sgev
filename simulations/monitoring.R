# Same-target comparisons. All settings are fixed before Monte Carlo sampling.
# The invariant comparator implements Corollary 4.3 of Lindon et al. (2025),
# arXiv:2210.08589v5 (revised 7 July 2025), using fixed Phi.
# Its pointwise likelihood formula is Theorem 4.1, equation (9), in that version.
# https://arxiv.org/abs/2210.08589v5
# The related JASA article appeared in 2026; result numbers here refer to v5.
# No asymptotic sample g-prior is used.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(26090819)
reps <- if (getOption("sgev.simulation.mode") == "smoke") 200L else 3000L
alpha <- .05; delta <- .05; tau <- .1
p <- 20L; horizon <- 8000L; increment <- 40L
times <- seq(increment, horizon, by = increment)
schedules <- c(5L, 20L, 200L)
effects <- rbind(c(0, 0), c(delta, delta), c(.08, 0), c(.12, 0), c(.20, 0))
effect_names <- c("zero", "boundary", "near", "moderate", "strong")
method_names <- c("Local mixture", "Scheduled joint SGPV",
                  "Invariant t-mixture", "Spending joint SGPV",
                  "Calibrated joint region")

# Finite-schedule Gaussian maximum-radius calibration. This is an independent
# Monte Carlo tolerance bound, not tuning against the performance paths.
# With probability >=1-eta the fixed boundary's noncoverage is <=alpha-eta.
# Accounting for calibration failure gives unconditional error <=alpha.
set.seed(26090951)
calibration_reps <- if (getOption("sgev.simulation.mode") == "smoke") 10000L else 200000L
eta <- .001
order_index <- qbinom(1 - eta, calibration_reps, 1 - alpha + eta) + 1L
stopifnot(order_index <= calibration_reps,
  pbinom(order_index - 1, calibration_reps, 1 - alpha + eta,
         lower.tail = FALSE) <= eta)
W <- matrix(0, calibration_reps, 2)
maximum <- matrix(0, calibration_reps, length(schedules))
for (k in seq_along(times)) {
  W <- W + matrix(rnorm(2 * calibration_reps), calibration_reps, 2)
  q <- rowSums(W^2) / k
  for (j in seq_along(schedules)) {
    if (k %% (length(times) / schedules[j]) == 0)
      maximum[, j] <- pmax(maximum[, j], q)
  }
}
critical <- apply(maximum, 2, function(x) sort(x, partial = order_index)[order_index])
write.csv(data.frame(looks = schedules, critical = critical,
  calibration_reps = calibration_reps, order_index = order_index,
  conditional_error = alpha - eta, failure_budget = eta, total_budget = alpha,
  seed = 26090951), "results/sequential_calibration.csv", row.names = FALSE)
rm(W, maximum)
set.seed(26090819)

# Return the common confidence-region geometry. A nonpositive direction makes
# the ellipsoid shortcut unavailable; returning no certificate is conservative.
invariant_geometry <- function(H, n, parameters, level, phi = 100) {
  d <- nrow(H); nu <- n - parameters
  stopifnot(nu > 0, all(eigen(H, symmetric = TRUE)$values > 0))
  det_penalty <- as.numeric(determinant(diag(d) + H / phi,
                                        logarithm = TRUE)$modulus)
  log_a <- (2 * log(level) - det_penalty) / (nu + d)
  a <- exp(log_a)
  precision <- a * H - phi * H %*% solve(phi * diag(d) + H)
  precision <- (precision + t(precision)) / 2
  list(precision = precision, radius_factor = -expm1(log_a),
       usable = min(eigen(precision, symmetric = TRUE)$values) > 0,
       nu = nu)
}

# Verify the confidence-boundary inversion against the authors' pointwise
# formula, independently of the simulation outcomes.
Hcheck <- matrix(c(1000, 250, 250, 1000), 2)
acheck <- invariant_geometry(Hcheck, 1200, 20, alpha)
direction <- c(1, -.4)
rss_check <- 1180 * 1.1
edge <- direction * sqrt(rss_check * acheck$radius_factor /
                          drop(crossprod(direction, acheck$precision %*% direction)))
log_point <- -.5 * determinant(diag(2) + Hcheck / 100,
                               logarithm = TRUE)$modulus +
  (1180 + 2) / 2 * (log1p(drop(crossprod(edge, Hcheck %*% edge)) / rss_check) -
                    log1p(drop(crossprod(edge,
                      (Hcheck - Hcheck %*% solve(Hcheck + 100 * diag(2), Hcheck)) %*% edge)) /
                      rss_check))
stopifnot(abs(log_point - log(1 / alpha)) < 1e-10)

results <- paired <- paths <- list(); ri <- 0L
for (rho in c(0, .9)) {
  R <- diag(p)
  R[1:2, 1:2] <- matrix(c(1, rho, rho, 1), 2)
  Ginc <- increment * R
  inv_inc <- solve(Ginc)
  chol_inc <- chol(Ginc)
  noise_score <- matrix(0, reps, p)
  noise_yty <- numeric(reps)
  first <- array(NA_integer_, c(reps, nrow(effects), length(schedules),
                                length(method_names)))
  for (index in seq_along(times)) {
    n <- times[index]
    new_score <- matrix(rnorm(reps * p), reps, p) %*% chol_inc
    # Exact joint Gaussian sufficient statistics from a fixed full-rank design.
    # The residual sum of squares in a fresh cohort is independent chi-square.
    noise_yty <- noise_yty + rowSums((new_score %*% inv_inc) * new_score) +
      rchisq(reps, increment - p)
    noise_score <- noise_score + new_score
    G <- n * R
    rss <- noise_yty - rowSums((noise_score %*% solve(G)) * noise_score)
    stopifnot(all(rss > 0))
    H <- G[1:2, 1:2]
    centers_noise <- noise_score[, 1:2] %*% solve(H)
    local_geometry <- joint_geometry(H, tau = tau, alpha = alpha)
    local_precision <- psd_parts(local_geometry$V)$pinv
    inv_geometry <- invariant_geometry(H, n, p, alpha, phi = 1 / tau^2)
    for (ei in seq_len(nrow(effects))) {
      centers <- sweep(centers_noise, 2, effects[ei, ], "+")
      qlocal <- pair_box_distance(centers, local_precision, delta)$dual
      qfixed <- pair_box_distance(centers, H, delta)$dual
      cinv <- if (inv_geometry$usable) {
        pair_box_distance(centers, inv_geometry$precision, delta)$dual >
          rss * inv_geometry$radius_factor * (1 + 1e-8)
      } else rep(FALSE, reps)
      for (si in seq_along(schedules)) {
        stride <- length(times) / schedules[si]
        if (index %% stride != 0) next
        look <- index / stride
        decisions <- cbind(
          qlocal > local_geometry$B * (1 + 1e-8),
          qfixed > qchisq(1 - alpha / schedules[si], 2) * (1 + 1e-8),
          cinv,
          qfixed > qchisq(1 - alpha / (look * (look + 1)), 2) * (1 + 1e-8),
          qfixed > critical[si] * (1 + 1e-8))
        for (mi in seq_along(method_names)) {
          newly <- is.na(first[, ei, si, mi]) & decisions[, mi]
          first[newly, ei, si, mi] <- n
        }
      }
    }
  }
  paths[[as.character(rho)]] <- first
  for (ei in seq_len(nrow(effects))) for (si in seq_along(schedules)) {
    for (mi in seq_along(method_names)) {
      stop_time <- first[, ei, si, mi]
      crossed <- !is.na(stop_time)
      truncated <- ifelse(crossed, stop_time, horizon)
      count <- sum(crossed)
      ri <- ri + 1L
      results[[ri]] <- data.frame(
        rho = rho, scenario = effect_names[ei], beta1 = effects[ei, 1],
        beta2 = effects[ei, 2], looks = schedules[si], method = method_names[mi],
        crossing_rate = mean(crossed), crossing_mcse = sd(crossed) / sqrt(reps),
        crossing_upper95 = if (count == reps) 1 else qbeta(.95, count + 1, reps - count),
        censored_fraction = mean(!crossed),
        mean_truncated_n = mean(truncated), mean_truncated_mcse = sd(truncated) / sqrt(reps),
        median_n_among_crossings = if (count) median(stop_time[crossed]) else NA_real_,
        reps = reps, horizon = horizon, alpha = alpha, delta = delta, tau = tau)
      if (mi > 1L) {
        ref <- first[, ei, si, 1L]
        a <- paired_summary(crossed, !is.na(ref))
        b <- paired_summary(truncated, ifelse(is.na(ref), horizon, ref))
        paired[[length(paired) + 1L]] <- data.frame(
          rho = rho, scenario = effect_names[ei], looks = schedules[si],
          method = method_names[mi], reference = method_names[1], reps = reps,
          crossing_difference = a[1], crossing_difference_mcse = a[2],
          mean_n_difference = b[1], mean_n_difference_mcse = b[2], row.names = NULL)
      }
    }
  }
  cat("Same-target sequential comparison complete: rho =", rho, "\n")
}
write.csv(do.call(rbind, results), "results/sequential_comparisons.csv", row.names = FALSE)
write.csv(do.call(rbind, paired), "results/sequential_paired.csv", row.names = FALSE)
saveRDS(list(first = paths, methods = method_names, schedules = schedules,
  effects = effects, effect_names = effect_names, horizon = horizon),
  "results/sequential_replicates.rds", version = 3)

# Fixed-total-sample batching study. Every partition uses the same design and
# Gaussian outcome paths. Zero information at saturation is a structural fact,
# not a small-eigenvalue approximation.
set.seed(26090820)
n_batch_total <- 640L
q <- 8L
XB <- matrix(rnorm(n_batch_total * q), n_batch_total, q)
XA <- .8 * XB[, 1] + .6 * rnorm(n_batch_total)
X <- cbind(XA, XB)
errors <- matrix(rnorm(reps * n_batch_total), reps, n_batch_total)
beta_target <- .35
batch_sizes <- c(4L, 8L, 16L, 32L, 64L, 128L, 640L)
pool_residual <- qr.resid(qr(XB), XA)
pool_information <- sum(pool_residual^2)
pool_score <- drop(errors %*% pool_residual) + pool_information * beta_target
pool_fixed <- pmax(abs(pool_score / pool_information) - delta, 0)^2 *
  pool_information > qchisq(1 - alpha, 1)
batch_results <- batch_paths <- list()
for (bi in seq_along(batch_sizes)) {
  size <- batch_sizes[bi]
  chunks <- split(seq_len(n_batch_total), rep(seq_len(n_batch_total / size), each = size))
  information <- 0; score <- numeric(reps)
  first_local <- first_scheduled <- rep(NA_integer_, reps)
  for (ci in seq_along(chunks)) {
    rows <- chunks[[ci]]
    if (size <= q) {
      stopifnot(qr(t(XB[rows, , drop = FALSE]))$rank == size)
      residual <- rep(0, size)
    } else {
      stopifnot(qr(XB[rows, , drop = FALSE])$rank == q)
      residual <- qr.resid(qr(XB[rows, , drop = FALSE]), XA[rows])
      # Cross-check the score-summary implementation on each supported cohort.
      hscore <- local_group_increment(crossprod(X[rows, , drop = FALSE]),
                                      rep(0, q + 1), 1L)$H
      stopifnot(abs(hscore - sum(residual^2)) < 1e-8)
    }
    information <- information + sum(residual^2)
    score <- score + drop(errors[, rows, drop = FALSE] %*% residual) +
      sum(residual^2) * beta_target
    if (information <= 0) next
    dist <- pmax(abs(score / information) - delta, 0)^2
    qlocal <- dist * information^2 / (information + 1 / tau^2)
    threshold <- 2 * log(1 / alpha) + log1p(information * tau^2)
    local <- qlocal > threshold * (1 + 1e-8)
    scheduled <- dist * information >
      qchisq(1 - alpha / length(chunks), 1) * (1 + 1e-8)
    first_local[is.na(first_local) & local] <- max(rows)
    first_scheduled[is.na(first_scheduled) & scheduled] <- max(rows)
  }
  batch_paths[[as.character(size)]] <- list(local = first_local,
    scheduled = first_scheduled, pooled_final = pool_fixed)
  for (method in c("Local mixture", "Scheduled joint SGPV")) {
    ft <- if (method == "Local mixture") first_local else first_scheduled
    crossed <- !is.na(ft)
    batch_results[[length(batch_results) + 1L]] <- data.frame(
      batch_size = size, looks = length(chunks), outside_dimension = q,
      local_rank = as.integer(information > 0), information = information,
      pooled_information = pool_information, retained_fraction = information / pool_information,
      method = method, crossing_rate = mean(crossed),
      crossing_mcse = sd(crossed) / sqrt(reps),
      mean_truncated_n = mean(ifelse(crossed, ft, n_batch_total)),
      mean_truncated_mcse = sd(ifelse(crossed, ft, n_batch_total)) / sqrt(reps),
      pooled_fixed_power = mean(pool_fixed), reps = reps,
      horizon = n_batch_total, beta1 = beta_target, alpha = alpha, delta = delta, tau = tau)
  }
}
batch_results <- do.call(rbind, batch_results)
stopifnot(all(diff(subset(batch_results, method == "Local mixture")$information) >= -1e-8),
          all(subset(batch_results, batch_size <= q)$crossing_rate == 0))
write.csv(batch_results, "results/batching_comparisons.csv", row.names = FALSE)
saveRDS(batch_paths, "results/batching_replicates.rds", version = 3)
writeLines(c("Seed: 26090819 (monitoring), 26090820 (batching)",
             "Invariant comparator: Lindon et al. (2025), arXiv:2210.08589v5, revised 7 July 2025.",
             "Inversion: Corollary 4.3. Pointwise likelihood: Theorem 4.1, equation (9). Fixed Phi = 100 I.",
             "https://arxiv.org/abs/2210.08589v5",
             "Result numbers refer to this preprint, not the related 2026 JASA article.",
             "Also checked against authors' R/e_variables.R at the commit below.",
             "Authors' avlm commit: e4d98c5e7ffcf804ab1105baa65dc79109a1653a.",
             "https://github.com/michaellindon/avlm/tree/e4d98c5e7ffcf804ab1105baa65dc79109a1653a",
             "No sample-dependent g-prior substitution is used.",
             capture.output(sessionInfo())), "results/sequential_session_info.txt")
cat("Batching comparison complete.\n")
