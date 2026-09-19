# Nuisance-drift comparison using base R.

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(26090920)

reps <- if (getOption("sgev.simulation.mode") == "smoke") 500L else 5000L

alpha <- 0.05
delta <- 0.05
alternative <- 0.12
sigma <- 1
tau <- 0.1
lambda <- sigma^2 / tau^2
cohort_n <- 300L
cohorts <- 20L
horizon <- cohort_n * cohorts
rho_grid <- c(0, 0.5, 0.8, 0.9, 0.95, 0.99, 1)
scenario_names <- c("Common boundary null", "Matched drift null",
                    "Common alternative")
method_names <- c("Local mixture", "Pooled common-model mixture",
                  "Oracle fixed-horizon test")

summarize_decision <- function(first, rho, scenario, method, information) {
  crossed <- !is.na(first)
  truncated <- ifelse(crossed, first, horizon)
  rate <- mean(crossed)
  count <- sum(crossed)
  upper <- if (count == reps) 1 else qbeta(0.95, count + 1, reps - count)
  data.frame(
    rho = rho,
    scenario = scenario,
    method = method,
    crossing_rate = rate,
    crossing_mcse = sqrt(rate * (1 - rate) / reps),
    crossing_upper95 = upper,
    mean_truncated_n = mean(truncated),
    censored_fraction = mean(!crossed),
    reps = reps,
    horizon = horizon,
    alpha = alpha,
    delta = delta,
    alternative = alternative,
    local_information = information,
    pooled_information = horizon,
    drift_KL = (alternative - delta)^2 * information / (2 * sigma^2),
    common_KL = (alternative - delta)^2 * horizon / (2 * sigma^2),
    robust_power_envelope = pnorm(qnorm(alpha) +
      (alternative - delta) * sqrt(information) / sigma),
    guarantee_under_drift = method != "Pooled common-model mixture",
    monitoring = if (method == "Oracle fixed-horizon test") "Final only" else "All 20 looks"
  )
}

rows <- list()
row_index <- 0L
for (rho in rho_grid) {
  first <- array(NA_integer_, c(reps, length(scenario_names),
                                length(method_names)))
  score_noise <- matrix(0, reps, 2)
  score_means <- matrix(0, length(scenario_names), 2)
  local_noise <- numeric(reps)
  H <- 0
  G <- matrix(0, 2, 2)

  for (b in seq_len(cohorts)) {
    signed_rho <- if (b %% 2L) rho else -rho
    Gb <- cohort_n * matrix(c(1, signed_rho, signed_rho, 1), 2)
    Z <- matrix(rnorm(2 * reps), reps, 2)
    noise <- sqrt(cohort_n) * cbind(
      Z[, 1],
      signed_rho * Z[, 1] + sqrt(1 - rho^2) * Z[, 2]
    )
    score_noise <- score_noise + noise
    hb <- cohort_n * (1 - rho^2)
    local_noise <- local_noise + noise[, 1] - signed_rho * noise[, 2]
    H <- H + hb
    G <- G + Gb

    gamma_null <- signed_rho * (alternative - delta)
    beta_by_scenario <- rbind(c(delta, 0), c(delta, gamma_null),
                              c(alternative, 0))
    score_means <- score_means + beta_by_scenario %*% Gb
    geom <- joint_geometry(G, sigma = sigma, tau = tau, alpha = alpha)
    pooled_estimable <- max(abs(c(1, 0) -
      geom$parts$P %*% c(1, 0))) < 1e-8

    for (scenario_index in seq_along(scenario_names)) {
      target <- beta_by_scenario[scenario_index, 1]
      U <- local_noise + H * target
      local_log_e <- 0.5 * (
        pmax(abs(U) - delta * H, 0)^2 /
          (sigma^2 * (H + lambda)) - log1p(H / lambda)
      )

      S <- sweep(score_noise, 2, score_means[scenario_index, ], "+")
      pooled_log_e <- rep(-Inf, reps)
      if (pooled_estimable) {
        centers <- S %*% geom$parts$pinv
        pooled_log_e <- 0.5 * (
          pmax(abs(centers[, 1]) - delta, 0)^2 / geom$V[1, 1] - geom$D
        )
      }

      crossed <- cbind(local_log_e, pooled_log_e) > log(1 / alpha) + 1e-8
      for (method_index in 1:2) {
        new <- crossed[, method_index] &
          is.na(first[, scenario_index, method_index])
        first[new, scenario_index, method_index] <- b * cohort_n
      }
      # Compare the vectorized scalar reduction with full-score elimination.
      mean_b <- drop(Gb %*% beta_by_scenario[scenario_index, ])
      direct <- local_group_increment(Gb, noise[1, ] + mean_b, 1L)
      expected <- noise[1, 1] - signed_rho * noise[1, 2] + hb * target
      stopifnot(abs(direct$u - expected) < 1e-9, abs(direct$H - hb) < 1e-9)
    }
  }

  stopifnot(abs(H - horizon * (1 - rho^2)) < 1e-8,
            max(abs(G - horizon * diag(2))) < 1e-8)

  random_at_zero <- runif(reps) < alpha

  for (scenario_index in seq_along(scenario_names)) {
    target <- if (scenario_index == 3L) alternative else delta
    oracle <- if (H > 0) {
      (local_noise + H * (target - delta)) / (sigma * sqrt(H)) > qnorm(1 - alpha)
    } else random_at_zero
    first[oracle, scenario_index, 3] <- horizon
    for (method_index in seq_along(method_names)) {
      row_index <- row_index + 1L
      rows[[row_index]] <- summarize_decision(
        first[, scenario_index, method_index], rho,
        scenario_names[scenario_index], method_names[method_index], H
      )
    }
  }

  stopifnot(identical(first[, 1, 1], first[, 2, 1]))
  if (rho == 1) {
    stopifnot(max(abs(score_means[2, ] - score_means[3, ])) < 1e-9,
              identical(first[, 2, 2], first[, 3, 2]),
              all(is.na(first[, , 1])))
  }
}

drift <- do.call(rbind, rows)
stopifnot(nrow(drift) == 63L,
          max(abs(drift$local_information -
            drift$horizon * (1 - drift$rho^2))) < 1e-8,
          max(abs(drift$drift_KL -
            drift$common_KL * (1 - drift$rho^2))) < 1e-10,
          all(drift$crossing_rate >= 0 & drift$crossing_rate <= 1))

if (reps == 5000L) {
  reference_rates <- data.frame(
    rho = c(0.8, 0.8, 0.9, 0.9, 1, 1),
    scenario = rep(c("Common alternative", "Matched drift null"), 3),
    method = rep(c("Local mixture", "Pooled common-model mixture"), 3),
    rate = c(0.6508, 0.4516, 0.3218, 0.7722, 0, 0.9624)
  )
  observed <- vapply(seq_len(nrow(reference_rates)), function(i) {
    z <- subset(drift, rho == reference_rates$rho[i] &
      scenario == reference_rates$scenario[i] & method == reference_rates$method[i])
    stopifnot(nrow(z) == 1L)
    z$crossing_rate
  }, numeric(1))
  stopifnot(identical(observed, reference_rates$rate))
}

write.csv(drift, "results/drift_comparisons.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/drift_session.txt")

binomial_interval <- function(rate, n, level = 0.95) {
  count <- round(rate * n)
  tail <- (1 - level) / 2
  cbind(
    lower = ifelse(count == 0, 0, qbeta(tail, count, n - count + 1)),
    upper = ifelse(count == n, 1, qbeta(1 - tail, count + 1, n - count))
  )
}

blue <- "#2B5C8A"
ink <- "#111111"
reference <- "#777777"
pdf("figures/drift_robustness.pdf", width = 7.1, height = 5.50,
    family = "Helvetica", pointsize = 11.5, colormodel = "srgb",
    bg = "white", useDingbats = FALSE)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
par(mar = c(2.35, 4.55, 2.55, .85), oma = c(2.55, .1, .1, .1),
    mgp = c(2.35, .62, 0), las = 1, cex.axis = .94, cex.lab = 1,
    bty = "l", family = "sans", lend = "butt", ljoin = "mitre",
    tcl = -.20, fg = ink, col.axis = "#222222", col.lab = "#222222",
    col.main = ink, xaxs = "r", yaxs = "r")

panel <- function(letter, title) {
  mtext(paste0("(", letter, ") ", title), side = 3, line = .10,
        adj = 0, cex = 1)
}

for (scenario_name in c("Common alternative", "Matched drift null")) {
  is_power <- scenario_name == "Common alternative"
  plot(NA, xlim = c(-.035, 1.035), ylim = c(0, 1.02), xaxp = c(0, 1, 2),
       xlab = "",
       ylab = if (is_power) "Rejection probability" else "False certification")
  panel(if (is_power) "a" else "b",
        if (is_power) "Common alternative" else "Matched drift null")
  if (is_power) {
    fraction <- seq(0, 1, length.out = 301)
    upper <- pnorm(qnorm(alpha) +
      (alternative - delta) * sqrt(horizon * fraction))
    lines(fraction, upper, col = reference, lty = 3, lwd = .8)
  } else {
    abline(h = alpha, col = reference, lty = 3, lwd = .75)
  }

  for (method_name in method_names) {
    z <- drift[drift$scenario == scenario_name &
               drift$method == method_name, , drop = FALSE]
    x <- 1 - z$rho^2
    order_index <- order(x)
    x <- x[order_index]
    z <- z[order_index, ]
    interval <- binomial_interval(z$crossing_rate, z$reps)
    visible <- interval[, 2] > interval[, 1]
    arrows(x[visible], interval[visible, 1], x[visible], interval[visible, 2],
           angle = 90, code = 3, length = .025, lwd = .65,
           col = if (method_name == "Local mixture") blue else ink)
    local <- method_name == "Local mixture"
    lines(x, z$crossing_rate, type = "b",
          col = if (local) blue else ink,
          pch = if (local) 16 else 2,
          lty = if (local) 1 else 2,
          lwd = .90, cex = .70)
  }
}

z <- subset(drift, scenario == "Common alternative" &
            method == "Local mixture")
fraction <- 1 - z$rho^2
order_index <- order(fraction)
plot(NA, xlim = c(-.035, 1.035), ylim = c(0, 17.5), xaxp = c(0, 1, 2),
     xlab = "", ylab = "KL distance to the null")
lines(fraction[order_index], z$drift_KL[order_index], col = blue, lwd = .95)
abline(h = unique(z$common_KL), col = reference, lty = 3, lwd = .75)
text(.04, 15.6, "Common null", adj = 0, cex = .80, col = reference)
text(.70, 3.5, "Drift null", col = blue, cex = .80)
panel("c", "Information cost")

par(mar = rep(.5, 4))
plot.new()
legend("center",
       c("Local mixture", "Pooled mixture: common model", "Sharp power bound"),
       col = c(blue, ink, reference), pch = c(16, 2, NA),
       lty = c(1, 2, 3), lwd = .90, bty = "n", ncol = 1,
       cex = .92, pt.cex = .78, x.intersp = .8, y.intersp = 1.25,
       seg.len = 2.4, text.col = "#222222")
mtext("Information retained", side = 1, outer = TRUE, line = .45,
      cex = .94, col = "#222222")
invisible(dev.off())

cat("Reproduced Figure 2 with ", reps,
    " common paths per correlation.\n", sep = "")
