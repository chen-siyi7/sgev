# Geometry used only in the constant-local-information experiment.
# sigma = 1. All rows of m are target estimates.
stopping_specs <- function() {
  v <- c(1, -1) / sqrt(2)
  Q <- matrix(c(1, 1, -1, 1), 2) / sqrt(2)
  list(
    list(id = "diagonal_box", label = "Full-rank box", H = diag(c(1, .25)),
         a = c(.35, 0), kind = "box", delta = .05),
    list(id = "singular_box", label = "Rank-one box", H = tcrossprod(v),
         a = c(.30, -.30), kind = "rank_one", delta = .05, v = v),
    list(id = "ellipsoid", label = "Rotated ellipse", H = diag(2),
         a = c(.30, .15), kind = "ellipse", Q = Q, radii = c(.06, .03)))
}

ellipse_project <- function(m, radii) {
  m <- as.matrix(m)
  z <- m
  outside <- rowSums(sweep(m, 2, radii, "/")^2) > 1
  if (any(outside)) {
    x <- m[outside, , drop = FALSE]
    lo <- numeric(nrow(x))
    hi <- sqrt(rowSums(sweep(x, 2, radii, "*")^2))
    # The upper endpoint is feasible. Fixed iterations avoid data-dependent
    # accuracy settings. A feasible support-function dual is used for rejection.
    for (i in seq_len(45)) {
      mid <- (lo + hi) / 2
      norm2 <- rowSums(sweep(x, 2, radii, "*")^2 /
                        outer(mid, radii^2, "+")^2)
      larger <- norm2 > 1
      lo[larger] <- mid[larger]
      hi[!larger] <- mid[!larger]
    }
    z[outside, ] <- sweep(x, 2, radii^2, "*") / outer(hi, radii^2, "+")
  }
  z
}

stopping_projection <- function(spec) {
  if (spec$kind == "box") return(pmax(-spec$delta, pmin(spec$delta, spec$a)))
  if (spec$kind == "rank_one") return(spec$delta * sign(spec$v))
  drop(ellipse_project(matrix(spec$a, 1) %*% spec$Q, spec$radii) %*% t(spec$Q))
}

stopping_log_e <- function(U, b, spec, tau = .1) {
  lambda <- 1 / tau^2
  if (spec$kind == "box") {
    h <- diag(spec$H)
    m <- sweep(U, 2, b * h, "/")
    w <- (b * h)^2 / (b * h + lambda)
    q <- drop(pmax(abs(m) - spec$delta, 0)^2 %*% w)
    determinant <- sum(log1p(b * h / lambda))
  } else if (spec$kind == "rank_one") {
    m <- drop(U %*% spec$v) / b
    q <- b^2 / (b + lambda) * pmax(abs(m) - spec$delta * sum(abs(spec$v)), 0)^2
    determinant <- log1p(b / lambda)
  } else {
    m <- U %*% spec$Q / b
    z <- ellipse_project(m, spec$radii)
    d <- m - z
    support <- sqrt(rowSums(sweep(d, 2, spec$radii, "*")^2))
    # Evaluating any feasible dual direction gives a lower bound on q.
    q <- b^2 / (b + lambda) * pmax(0, 2 * rowSums(d * m) - rowSums(d^2) - 2 * support)
    determinant <- 2 * log1p(b / lambda)
  }
  (q - determinant) / 2
}
