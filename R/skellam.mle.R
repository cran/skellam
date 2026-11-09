## --------------------------------------------------------------------------
## skellam.mle.R -- Maximum likelihood estimation
##
## Public functions:
##   skellam.mle()   MLE of (lambda1, lambda2) from a sample of differences
##
## Original author: Michail Tsagris
## Corrections:     Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' Maximum Likelihood Estimation for the Skellam Distribution
#'
#' Estimates the parameters of a Skellam distribution using maximum likelihood.
#'
#' @param x A vector of integers (positive or negative).
#'
#' @details
#' At the maximum likelihood estimate the fitted mean equals the sample mean,
#' \eqn{\hat\lambda_1 - \hat\lambda_2 = \bar{x}}{lambda1 - lambda2 = mean(x)}, so the
#' function maximizes the profile log-likelihood over \eqn{\lambda_2}{lambda2} alone,
#' with \eqn{\lambda_1 = \lambda_2 + \bar{x}}{lambda1 = lambda2 + mean(x)}, using
#' \code{\link[stats]{optimize}} on the admissible range
#' \eqn{\lambda_2 \ge \max(0, -\bar{x})}{lambda2 >= max(0, -mean(x))}, which keeps both
#' rates non-negative.
#'
#' The estimate can lie on the boundary, with one rate equal to 0 (the data are
#' then best described by a Poisson or a reflected Poisson). This happens when
#' the sample variance does not exceed the absolute sample mean, which a Skellam
#' with two positive rates cannot produce (its variance is
#' \eqn{\lambda_1 + \lambda_2 \ge |\lambda_1 - \lambda_2|}{lambda1 + lambda2 >= |lambda1 - lambda2|});
#' for example, a sample of zeros gives \eqn{(0, 0)} and a constant sample
#' \eqn{c > 0} gives \eqn{(c, 0)}.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{iters}{Number of log-likelihood evaluations used by the optimizer.}
#'   \item{loglik}{Maximized log-likelihood value.}
#'   \item{param}{Estimated parameters (\eqn{\hat{\lambda}_1}, \eqn{\hat{\lambda}_2}).}
#' }
#'
#' @references
#' \itemize{
#'   \item Butler, R. (2007) \emph{Saddlepoint Approximations with Applications},
#'   Cambridge University Press.
#'
#'   \item Johnson, N. L. (1959) On an extension of the connection between Poisson
#'   and \eqn{\chi^2} distributions. \emph{Biometrika} \bold{46}, 352-362.
#'
#'   \item Johnson, N. L.; Kotz, S.; Kemp, A. W. (1993) \emph{Univariate Discrete
#'   Distributions}, 2nd ed., John Wiley and Sons.
#'
#'   \item Skellam, J. G. (1946) The frequency distribution of the difference between
#'   two Poisson variates belonging to different populations. \emph{Journal of the
#'   Royal Statistical Society, Series A} \bold{109}(3), 296.
#'
#'   \item Strackee, J.; van der Gon, J. J. D. (1962) The frequency distribution of
#'   the difference between two Poisson variates. \emph{Statistica Neerlandica}
#'   \bold{16}(1), 17-23.
#'
#'   \item Abdulhamid, A. A.; Maha, A. O. (2010) On The Poisson Difference
#'   Distribution Inference and Applications. \emph{Bulletin of the Malaysian
#'   Mathematical Sciences Society} \bold{33}(1), 17-45.
#'
#'   \item Wikipedia: Skellam distribution \url{https://en.wikipedia.org/wiki/Skellam_distribution}
#' }
#'
#' @examples
#' # Basic example
#' set.seed(1)
#' x1 <- rpois(1000, 10)
#' x2 <- rpois(1000, 6)
#' x <- x1 - x2
#' skellam.mle(x)
#'
#' # Larger sample size
#' x1 <- rpois(10000, 10)
#' x2 <- rpois(10000, 6)
#' x <- x1 - x2
#' skellam.mle(x)
#'
#' @keywords distribution models
#' @author Michail Tsagris
#' @export
skellam.mle <- function(x) {

  if (!is.numeric(x)) stop("'x' must be a numeric vector of integers")
  if (anyNA(x)) stop("'x' contains missing values; remove them before fitting")
  if (any(!is.finite(x)) || any(.skellam_nonint(x)))
    stop("'x' must contain finite integers (differences of counts)")
  n <- length(x)
  if (n < 2L) stop("at least two observations are needed")
  x <- round(x)
  mx <- mean(x)

  # Profile log-likelihood in lambda2, with lambda1 = lambda2 + mean(x)
  # (the MLE matches the sample mean). Both rates must be >= 0, so
  # lambda2 >= lo = max(0, -mean(x)). theta = lambda2 - lo >= 0.
  lo <- max(0, -mx)
  evals <- 0L
  prof <- function(theta) {
    evals <<- evals + 1L
    l2 <- lo + theta
    sum(dskellam(x, l2 + mx, l2, log = TRUE))
  }

  # Search interval: the moment estimate of lambda2 is (var - mean)/2; the
  # optimum lies well inside [0, 10 * (var + |mean|) + 1].
  upper <- 10 * (stats::var(x) + abs(mx)) + 1
  opt <- stats::optimize(prof, c(0, upper), maximum = TRUE, tol = 1e-10 * (1 + upper))
  # optimize() never evaluates the end points; the MLE can sit on theta = 0
  # (one rate zero), e.g. for samples with var(x) <= |mean(x)|.
  at0 <- prof(0)
  if (at0 >= opt$objective) {
    theta <- 0; loglik <- at0
  } else {
    theta <- opt$maximum; loglik <- opt$objective
  }

  l2 <- lo + theta
  param <- c(l2 + mx, l2)
  names(param) <- c("mu 1", "mu 2")
  list(iters = evals, loglik = loglik, param = param)

}
