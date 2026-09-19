# Restyle the five manuscript figures from saved numerical summaries.
# No simulation is run and no numerical value is changed.

options(stringsAsFactors = FALSE)

data_dir <- "results"
figure_dir <- "figures"

seq_results <- read.csv(file.path(data_dir, "sequential_comparisons.csv"))
batch_results <- read.csv(file.path(data_dir, "batching_comparisons.csv"))
drift <- read.csv(file.path(data_dir, "drift_comparisons.csv"))
stopping <- read.csv(file.path(data_dir, "stopping_efficiency.csv"))
near <- read.csv(file.path(data_dir, "near_margin.csv"))

binomial_interval <- function(rate, reps, level = .95) {
  count <- round(rate * reps)
  tail <- (1 - level) / 2
  cbind(
    lower = ifelse(count == 0, 0, qbeta(tail, count, reps - count + 1)),
    upper = ifelse(count == reps, 1, qbeta(1 - tail, count + 1, reps - count))
  )
}

ink <- "#111111"
axis_ink <- "#222222"
reference_ink <- "#777777"
blue <- "#2B5C8A"
dark_gray <- "#3F3F3F"
mid_gray <- "#777777"
light_gray <- "#AAAAAA"

open_figure <- function(name, height) {
  pdf(file.path(figure_dir, name), width = 7.1, height = height,
      family = "Helvetica", pointsize = 11.5, colormodel = "srgb",
      bg = "white", useDingbats = FALSE)
}

plot_theme <- function(...) {
  par(family = "sans", fg = ink, col = ink, col.axis = axis_ink,
      col.lab = axis_ink, col.main = ink, bty = "l", las = 1,
      tcl = -.20, mgp = c(2.35, .62, 0), cex.axis = .94,
      cex.lab = 1.00, lend = "butt", ljoin = "mitre", lwd = .70,
      xaxs = "r", yaxs = "r", ...)
}

panel_grid <- function() {
  invisible(NULL)
}

panel_title <- function(label) {
  mtext(label, side = 3, line = .10, adj = 0, cex = 1.00)
}

error_bars <- function(x, lower, upper, color) {
  visible <- is.finite(lower) & is.finite(upper) & upper > lower
  arrows(x[visible], lower[visible], x[visible], upper[visible],
         angle = 90, code = 3, length = .025,
         col = color, lwd = .65)
}

series <- function(x, y, color, pch, lty) {
  lines(x, y, type = "b", col = color, pch = pch, lty = lty,
        lwd = .90, cex = .70)
}

bottom_legend <- function(height, labels, colors, pch, lty, ncol) {
  par(fig = c(0, 1, 0, height), mar = rep(0, 4), oma = rep(0, 4),
      new = TRUE, xpd = NA)
  plot.new()
  legend("center", legend = labels, col = colors, pch = pch, lty = lty,
         lwd = .90, pt.cex = .78, cex = .92, bty = "n", ncol = ncol,
         x.intersp = .65, y.intersp = 1.0, seg.len = 2.0,
         text.col = axis_ink)
}

legend_panel <- function(labels, colors, pch, lty) {
  par(mar = rep(.5, 4))
  plot.new()
  legend("center", legend = labels, col = colors, pch = pch, lty = lty,
         lwd = .90, pt.cex = .78, cex = .92, bty = "n", ncol = 1,
         x.intersp = .8, y.intersp = 1.25, seg.len = 2.4,
         text.col = axis_ink)
}

common_x_label <- function(label) {
  mtext(label, side = 1, outer = TRUE, line = .45,
        cex = .94, col = axis_ink)
}

# Figure 1: monitoring schedules.
method_colors <- c(
  "Local mixture" = blue,
  "Scheduled joint SGPV" = ink,
  "Invariant t-mixture" = dark_gray,
  "Spending joint SGPV" = mid_gray,
  "Calibrated joint region" = light_gray
)
method_shapes <- c(
  "Local mixture" = 16,
  "Scheduled joint SGPV" = 2,
  "Invariant t-mixture" = 0,
  "Spending joint SGPV" = 5,
  "Calibrated joint region" = 4
)
method_lines <- c(
  "Local mixture" = 1,
  "Scheduled joint SGPV" = 2,
  "Invariant t-mixture" = 4,
  "Spending joint SGPV" = 3,
  "Calibrated joint region" = 5
)
method_offsets <- setNames(seq(-.034, .034, length.out = 5), names(method_colors))

schedule_panel <- function(ylim, ylab, title) {
  plot(NA, xlim = log10(c(4.2, 245)), ylim = ylim, xaxt = "n",
       xlab = "", ylab = ylab)
  panel_grid()
  axis(1, at = log10(c(5, 20, 200)), labels = c(5, 20, 200),
       lwd = .65, lwd.ticks = .65)
  panel_title(title)
}

open_figure("sequential_comparisons.pdf", 5.35)
plot_theme(mfrow = c(2, 2), mar = c(2.55, 4.35, 1.95, .75),
           oma = c(4.45, 0, .05, 0))
for (scenario_id in c("near", "moderate", "strong")) {
  is_rate <- scenario_id != "strong"
  lim <- switch(scenario_id, near = c(0, .13), moderate = c(0, .85),
                strong = c(700, 2100))
  title <- switch(
    scenario_id,
    near = expression("(a) Near:" ~ beta[1] == 0.08),
    moderate = expression("(b) Moderate:" ~ beta[1] == 0.12),
    strong = expression("(c) Strong:" ~ beta[1] == 0.20)
  )
  schedule_panel(lim, if (is_rate) "Crossing probability" else "Mean sample size", title)
  for (method_id in names(method_colors)) {
    z <- subset(seq_results, rho == .9 & scenario == scenario_id & method == method_id)
    z <- z[order(z$looks), ]
    x <- log10(z$looks) + method_offsets[method_id]
    y <- if (is_rate) z$crossing_rate else z$mean_truncated_n
    ci <- if (is_rate) binomial_interval(y, z$reps) else
      cbind(pmax(0, y - 1.96 * z$mean_truncated_mcse),
            y + 1.96 * z$mean_truncated_mcse)
    error_bars(x, ci[, 1], ci[, 2], method_colors[method_id])
    series(x, y, method_colors[method_id], method_shapes[method_id],
           method_lines[method_id])
  }
}
schedule_panel(c(0, .06), "Largest null rate",
               "(d) Null calibration")
abline(h = .05, lty = 3, col = reference_ink, lwd = .75)
text(log10(230), .052, "0.05", adj = 1, cex = .78, col = reference_ink)
for (method_id in names(method_colors)) {
  z <- subset(seq_results,
              scenario %in% c("zero", "boundary") & method == method_id)
  means <- aggregate(crossing_rate ~ looks, z, max)
  uppers <- aggregate(crossing_upper95 ~ looks, z, max)
  x <- log10(means$looks) + method_offsets[method_id]
  error_bars(x, means$crossing_rate, uppers$crossing_upper95,
             method_colors[method_id])
  series(x, means$crossing_rate, method_colors[method_id],
         method_shapes[method_id], method_lines[method_id])
}
common_x_label("Planned looks")
bottom_legend(.12, names(method_colors), method_colors, method_shapes,
              method_lines, 3)
dev.off()

# Figure 2: nuisance heterogeneity.
drift_colors <- c("Local mixture" = blue,
                  "Pooled common-model mixture" = ink)
drift_shapes <- c("Local mixture" = 16,
                  "Pooled common-model mixture" = 2)
drift_lines <- c("Local mixture" = 1,
                 "Pooled common-model mixture" = 2)

open_figure("drift_robustness.pdf", 5.50)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot_theme(mar = c(2.35, 4.55, 2.55, .85), oma = c(2.55, .1, .1, .1))
for (scenario_id in c("Common alternative", "Matched drift null")) {
  is_power <- scenario_id == "Common alternative"
  plot(NA, xlim = c(-.035, 1.035), ylim = c(0, 1.02),
       xaxp = c(0, 1, 2), xlab = "",
       ylab = if (is_power) "Rejection probability" else "False certification")
  panel_grid()
  panel_title(if (is_power)
    "(a) Common alternative" else
    "(b) Matched drift null")
  if (is_power) {
    f <- seq(0, 1, length.out = 301)
    bound <- pnorm(qnorm(.05) + (.12 - .05) * sqrt(6000 * f))
    lines(f, bound, col = reference_ink, lty = 3, lwd = .8)
  } else {
    abline(h = .05, col = reference_ink, lty = 3, lwd = .75)
  }
  for (method_id in names(drift_colors)) {
    z <- subset(drift, scenario == scenario_id & method == method_id)
    x <- 1 - z$rho^2
    o <- order(x)
    x <- x[o]
    z <- z[o, ]
    ci <- binomial_interval(z$crossing_rate, z$reps)
    error_bars(x, ci[, 1], ci[, 2], drift_colors[method_id])
    series(x, z$crossing_rate, drift_colors[method_id],
           drift_shapes[method_id], drift_lines[method_id])
  }
}
z <- subset(drift,
            scenario == "Common alternative" & method == "Local mixture")
f <- 1 - z$rho^2
o <- order(f)
plot(NA, xlim = c(-.035, 1.035), ylim = c(0, 17.5), xaxp = c(0, 1, 2),
     xlab = "", ylab = "KL distance to the null")
panel_grid()
lines(f[o], z$drift_KL[o], col = blue, lwd = .95)
abline(h = unique(z$common_KL), col = reference_ink, lty = 3, lwd = .75)
text(.04, 15.6, "Common null", adj = 0, cex = .80, col = reference_ink)
text(.70, 3.5, "Drift null", cex = .80, col = blue)
panel_title("(c) Information cost")
legend_panel(
  c("Local mixture", "Pooled mixture: common model", "Sharp power bound"),
  c(drift_colors, reference_ink), c(drift_shapes, NA),
  c(drift_lines, 3))
common_x_label("Information retained")
dev.off()

# Figure S1: cohort size.
z_batch <- subset(batch_results, method == "Local mixture")
batch_ticks <- c(4, 16, 64, 640)
batch_axis <- function() {
  axis(1, at = log2(batch_ticks), labels = batch_ticks,
       lwd = .65, lwd.ticks = .65, gap.axis = -1)
}

open_figure("batching_comparisons.pdf", 5.50)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot_theme(mar = c(2.35, 4.55, 2.55, .85), oma = c(2.55, .1, .1, .1))
plot(NA, xlim = range(log2(z_batch$batch_size)) + c(-.12, .12),
     ylim = c(0, 1), xaxt = "n", xlab = "",
     ylab = "Retained information")
panel_grid()
series(log2(z_batch$batch_size), z_batch$retained_fraction, blue, 16, 1)
batch_axis()
panel_title("(a) Information")

plot(NA, xlim = range(log2(z_batch$batch_size)) + c(-.12, .12),
     ylim = c(0, 1.02), xaxt = "n", xlab = "",
     ylab = "Crossing probability")
panel_grid()
abline(h = unique(z_batch$pooled_fixed_power), col = reference_ink,
       lty = 3, lwd = .75)
for (method_id in c("Local mixture", "Scheduled joint SGPV")) {
  z <- subset(batch_results, method == method_id)
  ci <- binomial_interval(z$crossing_rate, z$reps)
  color <- method_colors[method_id]
  error_bars(log2(z$batch_size), ci[, 1], ci[, 2], color)
  series(log2(z$batch_size), z$crossing_rate, color,
         method_shapes[method_id], method_lines[method_id])
}
batch_axis()
panel_title("(b) Power")

plot(NA, xlim = range(log2(z_batch$batch_size)) + c(-.12, .12),
     ylim = c(395, 655), xaxt = "n", xlab = "",
     ylab = "Mean sample size")
panel_grid()
error_bars(log2(z_batch$batch_size),
           z_batch$mean_truncated_n - 1.96 * z_batch$mean_truncated_mcse,
           z_batch$mean_truncated_n + 1.96 * z_batch$mean_truncated_mcse, blue)
series(log2(z_batch$batch_size), z_batch$mean_truncated_n, blue, 16, 1)
batch_axis()
panel_title("(c) Monitoring delay")
legend_panel(
  c("Local mixture", "Scheduled joint SGPV", "Pooled final"),
  c(blue, ink, reference_ink), c(16, 2, NA), c(1, 2, 3))
common_x_label("Cohort size")
dev.off()

# Figure S2: stopping efficiency.
geometry_labels <- c(
  diagonal_box = "Full-rank box",
  singular_box = "Rank-one box",
  ellipsoid = "Rotated ellipse"
)
stop_colors <- c("Local mixture" = blue, "Supporting oracle" = ink)
stop_shapes <- c("Local mixture" = 16, "Supporting oracle" = 2)
stop_lines <- c("Local mixture" = 1, "Supporting oracle" = 2)

open_figure("stopping_efficiency.pdf", 5.50)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot_theme(mar = c(2.35, 4.55, 2.55, .85), oma = c(2.55, .1, .1, .1))
for (g in seq_along(geometry_labels)) {
  geometry_id <- names(geometry_labels)[g]
  plot(NA, xlim = c(-1, 33), ylim = c(.92, 2.55),
       xlab = "",
       ylab = if (g == 1) "Normalized mean time" else "")
  panel_grid()
  abline(h = 1, col = reference_ink, lty = 3, lwd = .75)
  panel_title(paste0("(", letters[g], ") ", geometry_labels[g]))
  for (method_id in names(stop_colors)) {
    z <- stopping[stopping$geometry == geometry_id &
                  stopping$method == method_id, ]
    x <- -log10(z$level)
    error_bars(x, z$ratio - 1.96 * z$ratio_mcse,
               z$ratio + 1.96 * z$ratio_mcse, stop_colors[method_id])
    series(x, z$ratio, stop_colors[method_id], stop_shapes[method_id],
           stop_lines[method_id])
  }
}
legend_panel(
  c("Local mixture", "Supporting oracle", "First-order bound"),
  c(stop_colors, reference_ink), c(stop_shapes, NA), c(stop_lines, 3))
common_x_label(expression(-log[10](epsilon)))
dev.off()

# Figure S3: behavior near the practical margin.
near_colors <- c("Fixed scale" = blue, "Scale mixture" = ink,
                 "Stitched boundary" = reference_ink)
near_shapes <- c("Fixed scale" = 16, "Scale mixture" = 2,
                 "Stitched boundary" = 0)
near_lines <- c("Fixed scale" = 1, "Scale mixture" = 2,
                "Stitched boundary" = 4)
positive <- subset(near, grepl("^gap_", scenario))

open_figure("near_margin.pdf", 5.50)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot_theme(mar = c(2.35, 4.55, 2.55, .85), oma = c(2.55, .1, .1, .1))
for (panel in 1:3) {
  plot(NA, xlim = c(.0022, .09),
       ylim = switch(panel, c(0, 1.02), c(1000, 4000000), c(.8, 4)),
       log = if (panel == 2) "xy" else "x", xaxt = "n",
       yaxt = if (panel == 2) "n" else "s",
       xlab = "",
       ylab = c("Crossing probability", "Mean information",
                "Normalized information")[panel])
  panel_grid()
  axis(1, at = c(.0025, .01, .04), labels = c(".0025", ".01", ".04"),
       lwd = .65, lwd.ticks = .65)
  if (panel == 2) {
    axis(2, at = 10^(3:6), labels = expression(10^3, 10^4, 10^5, 10^6),
         lwd = .65, lwd.ticks = .65)
  }
  if (panel == 3) {
    abline(h = 1, col = reference_ink, lty = 3, lwd = .75)
  }
  for (method_id in names(near_colors)) {
    z <- subset(positive, method == method_id)
    z <- z[order(z$gap), ]
    if (panel == 1) {
      value <- z$crossing_rate
      ci <- binomial_interval(value, z$reps)
    } else {
      divisor <- if (panel == 2) 1 else 2 * log(1 / .05) / z$gap^2
      value <- z$mean_truncated_n / divisor
      ci <- cbind(value - 1.96 * z$mean_truncated_mcse / divisor,
                  value + 1.96 * z$mean_truncated_mcse / divisor)
    }
    error_bars(z$gap, ci[, 1], ci[, 2], near_colors[method_id])
    series(z$gap, value, near_colors[method_id], near_shapes[method_id],
           near_lines[method_id])
  }
  title <- c("Power", "Mean information", "Normalized cost")[panel]
  panel_title(paste0("(", letters[panel], ") ", title))
}
legend_panel(names(near_colors), near_colors, near_shapes, near_lines)
common_x_label("Gap above margin")
dev.off()

message("Restyled five figures from saved numerical summaries.")
