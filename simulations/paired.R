# First practical declaration from a common Gaussian-mixture confidence sequence.
# Run through simulations/run.R.
options(stringsAsFactors = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(250919)
reps <- if (getOption("sgev.simulation.mode") == "smoke") 200L else 5000L
cohort_size <- 40L
looks <- 500L
horizon <- cohort_size * looks
delta <- .05
alpha <- .05
sigma <- 1
tau <- .1
lambda <- sigma^2 / tau^2
ratios <- c(0, .5, 1, 1.5, 2, 3)
correlations <- c(0, .8)

# Standardized score innovations are shared by all effects and correlations.
innovations <- matrix(rnorm(reps * looks), nrow = reps)
interval <- function(count) {
  c(if (count == 0) 0 else qbeta(.025, count, reps - count + 1),
    if (count == reps) 1 else qbeta(.975, count + 1, reps - count))
}
summary_rows <- list()
checks <- 0L
for (rho in correlations) {
  h <- cohort_size * (1 - rho^2)
  for (ratio in ratios) {
    beta <- ratio * delta
    score <- numeric(reps)
    first <- integer(reps) # 0 unresolved, 1 meaningful, 2 negligible
    sample_size <- rep(horizon, reps)
    any_out <- any_in <- rep(FALSE, reps)
    for (k in seq_len(looks)) {
      score <- score + h * beta + sigma * sqrt(h) * innovations[, k]
      information <- k * h
      estimate <- score / information
      D <- log1p(information / lambda)
      v <- sigma^2 * (information + lambda) / information^2
      radius <- sqrt(v * (2 * log(1 / alpha) + D))
      out <- -D / 2 + pmax(abs(estimate) - delta, 0)^2 / (2 * v) > log(1 / alpha)
      inward <- -D / 2 + pmax(delta - abs(estimate), 0)^2 / (2 * v) > log(1 / alpha)
      implemented <- sgev_scalar(score, information, delta, sigma, tau, alpha)
      # Independent interval implementation, with the same strict conventions.
      lower <- estimate - radius
      upper <- estimate + radius
      stopifnot(identical(out, lower > delta | upper < -delta),
                identical(inward, lower > -delta & upper < delta),
                identical(out, implemented$meaningful),
                identical(inward, implemented$negligible),
                !any(out & inward))
      checks <- checks + reps
      any_out <- any_out | out
      any_in <- any_in | inward
      hit <- first == 0L & (out | inward)
      sample_size[hit] <- cohort_size * k
      first[hit & out] <- 1L
      first[hit & inward] <- 2L
    }
    # At the margin both strict declarations are false.
    false_ever <- if (beta < delta) any_out else if (beta > delta) any_in else any_out | any_in
    false_first <- if (beta < delta) first == 1L else if (beta > delta) first == 2L else first != 0L
    counts <- c(sum(first == 1L), sum(first == 2L), sum(first == 0L))
    ci <- t(vapply(counts, interval, numeric(2)))
    stopifnot(sum(counts) == reps, all(sample_size %% cohort_size == 0))
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      rho = rho, ratio = ratio, beta = beta, reps = reps,
      horizon = horizon, information_per_cohort = h,
      meaningful = counts[1] / reps, meaningful_lower = ci[1, 1], meaningful_upper = ci[1, 2],
      negligible = counts[2] / reps, negligible_lower = ci[2, 1], negligible_upper = ci[2, 2],
      unresolved = counts[3] / reps, unresolved_lower = ci[3, 1], unresolved_upper = ci[3, 2],
      mean_truncated_n = mean(sample_size), mcse_n = sd(sample_size) / sqrt(reps),
      false_first = mean(false_first), false_ever = mean(false_ever))
  }
}
results <- do.call(rbind, summary_rows)
write.csv(results, "results/paired_decisions_summary.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/paired_decisions_session.txt")

pdf("figures/paired_decisions.pdf", width = 7.1, height = 5.5,
    family = "Helvetica", pointsize = 11.5, useDingbats = FALSE)
par(mfrow = c(2, 2), mar = c(3.6, 4.1, 2, 1), oma = c(2.1, 0, 0, 0),
    family = "sans", bty = "l", las = 1, tcl = -.2, mgp = c(2.5, .6, 0),
    cex.axis = .9, cex.lab = .94, lwd = .7)
colors <- c("#2B5C8A", "#222222", "#888888")
shapes <- c(16, 2, 0)
types <- c(1, 2, 3)
outcomes <- c("meaningful", "negligible", "unresolved")
for (i in seq_along(correlations)) {
  d <- results[results$rho == correlations[i], ]
  plot(NA, xlim = range(ratios), ylim = c(-.025, 1.025), xlab = expression(beta / delta),
       ylab = "First-declaration probability", xaxt = "n", yaxt = "n", yaxs = "i")
  axis(1, at = ratios, labels = ratios)
  axis(2, at = seq(0, 1, .2), labels = sprintf("%.1f", seq(0, 1, .2)))
  abline(v = 1, col = "#BBBBBB", lty = 3, lwd = .6)
  mtext(bquote(.(letters[(i - 1) * 2 + 1]) * ")  " * abs(rho) == .(correlations[i])),
        side = 3, line = .25, adj = 0, cex = .96)
  for (j in seq_along(outcomes)) {
    outcome <- outcomes[j]
    lines(d$ratio, d[[outcome]], type = "b", pch = shapes[j], lty = types[j],
          col = colors[j], lwd = .85, cex = .65)
    arrows(d$ratio, d[[paste0(outcome, "_lower")]],
           d$ratio, d[[paste0(outcome, "_upper")]], angle = 90, code = 3,
           length = .025, lwd = .55, col = colors[j])
  }
  plot(d$ratio, d$mean_truncated_n / 1000, type = "b", ylim = c(0, 20.5),
       xlab = expression(beta / delta), ylab = "Mean observations (thousands)",
       pch = 16, col = "#222222", lwd = .85, cex = .65, xaxt = "n", yaxt = "n")
  axis(1, at = ratios, labels = ratios)
  axis(2, at = seq(0, 20, 5))
  abline(v = 1, col = "#BBBBBB", lty = 3, lwd = .6)
  arrows(d$ratio, (d$mean_truncated_n - 1.96 * d$mcse_n) / 1000,
         d$ratio, (d$mean_truncated_n + 1.96 * d$mcse_n) / 1000,
         angle = 90, code = 3, length = .025, lwd = .55)
  mtext(bquote(.(letters[(i - 1) * 2 + 2]) * ")  " * abs(rho) == .(correlations[i])),
        side = 3, line = .25, adj = 0, cex = .96)
}
par(fig = c(0, 1, 0, .085), mar = rep(0, 4), oma = rep(0, 4), new = TRUE)
plot.new()
legend("center", c("Meaningful", "Negligible", "Unresolved"),
       col = colors, pch = shapes, lty = types, lwd = .85, pt.cex = .7,
       cex = .93, bty = "n", horiz = TRUE, seg.len = 2.3)
invisible(dev.off())
print(results[, c("rho", "ratio", "meaningful", "negligible", "unresolved", "mean_truncated_n", "false_ever")])
cat("Profile/interval equality checked at", checks, "path-look combinations.\n")
