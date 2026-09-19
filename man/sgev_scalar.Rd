\name{sgev_scalar}
\alias{sgev_scalar}
\title{Scalar Second-Generation E-Values and Practical Decisions}
\description{Vectorized outward and inward normal-mixture profiles from
accumulated reduced Gaussian scores.}
\usage{
sgev_scalar(score, information, margin = 0.05, sigma = 1,
            tau = 0.1, alpha = 0.05)
}
\arguments{
\item{score}{Nonempty finite numeric vector of accumulated reduced scores.}
\item{information}{Accumulated nonnegative scalar Gram information, or a vector
of the same length as \code{score}. It is not divided by the error variance.}
\item{margin}{Positive practical-null half-width, fixed before monitoring.}
\item{sigma}{Known positive error standard deviation.}
\item{tau}{Positive mixture scale, fixed before monitoring.}
\item{alpha}{Error level strictly between zero and one.}
}
\details{
At each time, the score has conditional increments with mean
\eqn{h\beta} and variance \eqn{\sigma^2 h}. For sequential validity, inputs
must arise from one cumulative score process with predictable information.
The same confidence sequence gives strict exclusion of the interval null
and strict containment within it. Both decisions share the coverage error
bound. Unrelated rows are not automatically protected against multiplicity.

Zero information requires zero score. It returns an unbounded interval,
zero log e-values, an unresolved decision, and an undefined (\code{NA}) SGPV.
All calculations use the known variance. Log e-values avoid exponentiation
at large evidence. Values at a floating-point threshold should be treated
with numerical caution; strict interval comparisons determine the decision.
}
\value{A data frame with one row per input score. Columns are \code{estimate},
\code{radius}, \code{lower}, \code{upper}, \code{log_out}, \code{log_in},
\code{meaningful}, \code{negligible}, \code{decision}, and \code{sgpv}.
The decision labels are \code{non-negligible}, \code{negligible}, and
\code{unresolved}. The outward profile tests absolute effects at most the
margin. The inward profile tests absolute effects at least the margin.}
\seealso{\code{\link{sgev_contrast}}, \code{\link{sgev_start}}}
\examples{
sgev_scalar(score = c(0, 600), information = 5000)
sgev_scalar(score = 0, information = 0)
}
