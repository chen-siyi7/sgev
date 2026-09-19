\name{sgev}
\alias{sgev}
\alias{sgev-package}
\alias{sgev_start}
\alias{sgev_update}
\alias{sgev_fit}
\alias{sgev_box}
\alias{sgev_contrast}
\alias{sgev_convex}
\alias{sgev_sgpv}
\alias{sgev_information}
\alias{sgev_partition}
\title{Second-Generation E-Values for Sequential Regression}
\description{
Sequential normal-mixture evidence for practical-null hypotheses after
cohortwise elimination of nuisance coefficients. The target coefficients are
common across cohorts; outside coefficients may change predictably.
}
\usage{
sgev_start(p, target, sigma = 1, tau = 0.1, alpha = 0.05,
           weight = 1, tolerance = 1e-10)
sgev_update(object, information, score)
sgev_fit(object)
sgev_box(fit, margin = 0.05)
sgev_contrast(fit, contrast, margin = 0.05)
sgev_convex(fit, coordinates, support)
sgev_sgpv(lower, upper, margin = 0.05)
sgev_information(information, alternative, margin = 0.05,
                 sigma = 1, tolerance = 1e-10)
sgev_partition(information_list, target, tolerance = 1e-10)
}
\arguments{
\item{p}{Positive integer number of full-model predictors.}
\item{target}{Distinct indices of the stable target coefficients, in output order.}
\item{sigma}{Known positive error standard deviation, not variance.}
\item{tau}{Positive normal-mixture coefficient scale, fixed before monitoring.}
\item{alpha}{Error level strictly between zero and one.}
\item{weight}{Prespecified envelope weight in (0, 1]. Weights across envelopes
must sum to at most one. Each fit uses local level \code{alpha * weight}.}
\item{tolerance}{Positive relative spectral threshold, less than one.}
\item{object}{An object returned by \code{sgev_start} or \code{sgev_update}.}
\item{information}{For an update, the new cohort Gram matrix, with dimensions
\code{p} by \code{p}. For \code{sgev_information}, the target's reduced
information matrix. Matrices must be symmetric positive semidefinite.}
\item{score}{The full score vector for the new cohort. In a regression this is
\code{crossprod(X, y)} and information is \code{crossprod(X)}.}
\item{fit}{An object returned by \code{sgev_fit}.}
\item{margin}{Positive common practical-null margin.}
\item{contrast}{Nonzero numeric vector in target-coordinate order.}
\item{coordinates}{Dual direction coordinates in \code{fit$basis}, a vector of
length \code{fit$rank}. They may depend on the observed data.}
\item{support}{Function evaluating the support function of a fixed nonempty
compact convex null at a target-coordinate vector. A certified upper bound
on the support function also yields a lower e-value.}
\item{lower, upper}{Numeric interval endpoints of equal positive length.
Infinite endpoints are allowed. Missing values are not allowed.}
\item{alternative}{Finite target vector for the information-distance calculation.}
\item{information_list}{Nonempty list of full cohort information matrices.}
}
\details{
Supply each disjoint cohort exactly once. Supplying cumulative statistics as
successive increments counts observations repeatedly and is invalid. Cohort
designs and nuisance coefficients must be predictable. Under the Gaussian
model, the conditional score covariance is \eqn{\sigma^2 G}. The conditional
sub-Gaussian moment bound in the manuscript also suffices for e-process validity.

The target, scale, variance, level, and weight are stored at initialization.
Do not modify state or fit fields during monitoring. For data-selected
contrasts, every coefficient in the containing envelope must remain stable.
Joint confidence coverage protects selected contrasts and margins. Individual
profile e-process guarantees use fixed contrasts and margins.

\code{sgev_box} minimizes over the common-margin target box and uses a feasible
dual lower bound for rejection. \code{sgev_convex} evaluates a feasible dual
bound at a supplied direction; it does not optimize that direction. The support
function must describe the intended null. Both retain the full determinant penalty.

Small positive target directions may be omitted from the quadratic term.
Positive nuisance directions are not silently truncated. Exact zero, duplicate,
and sign-duplicate nuisance columns are reduced; unresolved nuisance rank stops
evaluation. Finite contrast intervals require full retained rank or exact column
identities. Other contrasts return an unbounded interval, even if a mathematical
range calculation could establish estimability by another argument.

Dual lower bounds use numerical roundoff allowances. These are floating-point
implementations, not interval-arithmetic certificates. Inspect the primal-dual
gap and optimizer convergence. Use log e-values to avoid overflow at large evidence.
The primal objective is a diagnostic upper bound and is not used for rejection.

The information distance is half the minimum squared distance to the box in
the metric \eqn{H/\sigma^2}. Its lower and upper bounds describe numerical
optimization error. The manuscript's sharp power result uses deterministic
designs; its stopping-time attainment result uses constant local information.
}
\value{
\code{sgev_start} and \code{sgev_update} return an \code{sgev_state} containing
accumulated reduced score \code{u}, information \code{H}, cohort count, and
the fixed analysis settings. Updates return a new state.

\code{sgev_fit} returns an \code{sgev_fit} with \code{center}, the retained
\code{basis} and \code{rank}, and the covariance geometry used by the decisions.

\code{sgev_box} returns \code{log_evalue}, \code{evalue}, \code{reject},
\code{threshold}, and optimization diagnostics \code{primal}, \code{dual},
\code{gap}, and \code{convergence}. The rejection rule includes a numerical
tolerance above the theoretical threshold.

\code{sgev_contrast} returns the interval, two log profile e-values, an SGPV,
and a \code{decision}: non-negligible, negligible, or unresolved. Boundary
contacts are unresolved. SGPVs for unbounded or degenerate intervals are
\code{NA}; unbounded intervals have unresolved decisions.

\code{sgev_convex} returns the dual direction, lower objective, log e-value,
e-value, and rejection indicator. \code{sgev_sgpv} returns a numeric vector.

\code{sgev_information} returns lower and upper information-distance bounds,
their gap, a primal minimizer, and the optimizer convergence code.
\code{sgev_partition} returns local and pooled information, their loss, and
the corresponding nuisance-loading matrices in the requested target order.
}
\references{
Chen, S. Second-generation e-values for sequential assessment of
practically meaningful effects. Manuscript, 2026.

Blume, J. D., D'Agostino McGowan, L., Dupont, W. D. and Greevy, R. A. Jr.
(2018). Second-generation p-values. PLOS ONE 13(3), e0188299.
\doi{10.1371/journal.pone.0188299}.
}
\examples{
state <- sgev_start(p = 2, target = 1, sigma = 1)
for (b in 1:20) {
  rho <- if (b \%\% 2) 0.6 else -0.6
  G <- 300 * matrix(c(1, rho, rho, 1), 2)
  score <- drop(G \%*\% c(0.12, 0.1 * sin(b)))
  state <- sgev_update(state, G, score)
}
fit <- sgev_fit(state)
sgev_box(fit, margin = 0.05)
sgev_contrast(fit, contrast = 1, margin = 0.05)
sgev_information(state$H, alternative = 0.12)

# A feasible dual bound for the scalar interval null.
v <- drop(crossprod(fit$basis, fit$center))
sgev_convex(fit, coordinates = v, support = function(x) 0.05 * sum(abs(x)))
sgev_sgpv(c(0.05, -0.05), c(0.1, 0.05))
}
\keyword{htest}
