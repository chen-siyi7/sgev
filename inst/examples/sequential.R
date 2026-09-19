# A self-contained simulation using the installed package.
library(sgev)
set.seed(20260912)
state <- sgev_start(2, target = 1, sigma = 1, tau = 0.1)
history <- data.frame(cohort = integer(), log_evalue = numeric(), reject = logical())
for (b in seq_len(50)) {
  rho <- if (b %% 2) 0.6 else -0.6
  G <- 100 * matrix(c(1, rho, rho, 1), 2)
  beta <- c(0.12, 0.15 * sin(b / 4))
  score <- drop(G %*% beta + t(chol(G)) %*% rnorm(2))
  state <- sgev_update(state, G, score)
  result <- sgev_box(sgev_fit(state))
  history[nrow(history) + 1L, ] <- list(b, result$log_evalue, result$reject)
  if (result$reject) break
}
print(history)
print(sgev_contrast(sgev_fit(state), 1))
