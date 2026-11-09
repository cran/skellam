## --------------------------------------------------------------------------
## dskellam.r -- Skellam probability mass function
##
## Public functions:
##   dskellam()      P(X = x), on the natural or log scale
##
## Original author: Jerry W. Lewis
## Corrections:     Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' The Skellam Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the Skellam distribution, the distribution of the difference of two
#' independent Poisson variables, together with saddlepoint approximations to
#' the density and distribution function.
#'
#' @details
#' If \eqn{Y_1 \sim \mathrm{Poisson}(\lambda_1)}{Y1 ~ Poisson(lambda1)} and
#' \eqn{Y_2 \sim \mathrm{Poisson}(\lambda_2)}{Y2 ~ Poisson(lambda2)} are independent, then
#' \eqn{X = Y_1 - Y_2} has the Skellam distribution with probability mass
#' function
#' \deqn{P(X = x) = e^{-(\lambda_1+\lambda_2)} (\lambda_1/\lambda_2)^{x/2}
#'   I_{|x|}(2\sqrt{\lambda_1 \lambda_2}), \qquad x = \ldots, -1, 0, 1, \ldots}{
#'   P(X = x) = exp(-(lambda1+lambda2)) (lambda1/lambda2)^(x/2) I_|x|(2 sqrt(lambda1 lambda2)),  x = ..., -1, 0, 1, ...}
#' where \eqn{I_\nu} is the modified Bessel function of the first kind. The
#' mean is \eqn{\lambda_1 - \lambda_2} and the variance
#' \eqn{\lambda_1 + \lambda_2}. (If \eqn{Y_1} and \eqn{Y_2} are correlated
#' Poisson variables with means \eqn{\mu_i} and correlation \eqn{\rho}, the
#' difference is still Skellam, with
#' \eqn{\lambda_i = \mu_i - \rho\sqrt{\mu_1 \mu_2}}{lambda_i = mu_i - rho sqrt(mu1 mu2)}.)
#'
#' All functions follow the conventions of the \pkg{stats} distribution
#' functions: arguments are recycled to a common length (with a warning when
#' the lengths are not multiples), missing values propagate, invalid rates and
#' probabilities give \code{NaN} with a warning, a non-integer \code{x} has
#' density 0 (with a warning), and \code{qskellam} returns the smallest integer
#' \eqn{x} with \eqn{P(X \le x) \ge p}{P(X <= x) >= p} (for \code{lower.tail = TRUE}), so that
#' \code{qskellam(pskellam(x, ...), ...)} recovers \code{x}.
#'
#' \strong{Numerical methods.} \code{dskellam} evaluates the log of the Bessel
#' function directly, so that \code{log = TRUE} is finite and accurate for any
#' input: through \code{\link{besselI}} where that is reliable, through the
#' uniform asymptotic expansion for large orders, a stable downward recurrence
#' for large arguments, and the power series where \code{besselI} underflows.
#' The relative error is below \eqn{10^{-12}}{1e-12} for rates up to at least
#' \eqn{10^7}{1e7}, and below \eqn{10^{-9}}{1e-9} for extremely unequal large rates such
#' as \eqn{(10^6, 10^3)}{(1e6, 1e3)}. \code{pskellam} uses the non-central chi-square
#' identity of Johnson (1959) through \code{\link{pchisq}} when
#' \eqn{\lambda_1 + \lambda_2 \le 10^4}{lambda1 + lambda2 <= 1e4}, and otherwise, as well as in tails
#' where \code{pchisq} underflows, sums the log-pmf directly toward the nearer
#' tail, so that far tails are exact on the log scale. \code{qskellam} starts
#' from a Cornish-Fisher approximation and searches the integers with
#' \code{pskellam}. \code{rskellam} returns the difference of two
#' \code{\link{rpois}} draws.
#'
#' \code{dskellam.sp} and \code{pskellam.sp} are the second-order saddlepoint
#' approximation to the pmf and the Lugannani-Rice approximation to the cdf
#' with Butler's second continuity correction (Butler 2007). Their relative
#' error is of order \eqn{1/(\lambda_1+\lambda_2)^2}{1/(lambda1+lambda2)^2} and
#' \eqn{1/(\lambda_1+\lambda_2)}{1/(lambda1+lambda2)} respectively; they are provided for users who
#' want the approximations themselves and are not used by \code{dskellam} or
#' \code{pskellam}.
#'
#' @param x,q Numeric vector of quantiles (integers for the pmf).
#' @param p Numeric vector of probabilities.
#' @param n Number of observations; if \code{length(n) > 1}, the length is
#'   taken to be the number required.
#' @param lambda1,lambda2 Numeric vectors of non-negative rates;
#'   \code{lambda2} defaults to \code{lambda1}.
#' @param log,log.p Logical; if \code{TRUE}, probabilities are given as
#'   \code{log(p)}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le x)}{P(X <= x)}, otherwise \eqn{P(X > x)}.
#'
#' @return
#' \code{dskellam} gives the (log) probability mass function, \code{pskellam}
#' the (log) distribution function, \code{qskellam} the quantile function and
#' \code{rskellam} generates random deviates (an integer vector).
#' \code{dskellam.sp} and \code{pskellam.sp} give the saddlepoint
#' approximations. Names of the first argument are kept.
#'
#' Invalid \code{lambda} values return \code{NaN} with a warning; a zero rate
#' is allowed and reduces the distribution to a (reflected) Poisson.
#'
#' @references
#' \itemize{
#'   \item Butler, R. (2007) \emph{Saddlepoint Approximations with Applications},
#'         Cambridge University Press.
#'   \item Johnson, N. L. (1959) On an extension of the connection between Poisson and \eqn{\chi^2} distributions.
#'         \emph{Biometrika} 46, 352-362.
#'   \item Johnson, N. L., Kotz, S., & Kemp, A. W. (1993) \emph{Univariate Discrete Distributions},
#'         2nd ed., John Wiley and Sons.
#'   \item Olver, F. W. J. et al. (eds.) \emph{NIST Digital Library of Mathematical Functions},
#'         Section 10.41 (uniform asymptotic expansions of Bessel functions). \url{https://dlmf.nist.gov/10.41}
#'   \item Skellam, J. G. (1946) The frequency distribution of the difference between two Poisson variates.
#'         \emph{Journal of the Royal Statistical Society, Series A} 109(3), 296.
#'   \item Strackee, J. & van der Gon, J. J. D. (1962) The frequency distribution of the difference
#'         between two Poisson variates. \emph{Statistica Neerlandica} 16(1), 17-23.
#' }
#'
#' @seealso \code{\link{skellam.mle}} and \code{\link{skellam.reg}} for
#'   estimation; \code{\link{dpois}} for the Poisson distribution.
#'
#' @examples
#' # Compare with Poisson when one lambda = 0
#' dskellam(0:10, 5, 0)
#' dpois(0:10, 5)
#'
#' # Both lambdas non-zero
#' dskellam(c(-1,1), c(12,10), c(10,12))
#' pskellam(c(-1,0), c(12,10), c(10,12))
#'
#' # Quantile function inverts the distribution function
#' qskellam(c(0.05, 0.95), 3, 4)
#' qskellam(pskellam(-3:3, 3, 4), 3, 4)
#'
#' # Log-scale values stay finite far in the tails and at large rates
#' dskellam(500, 3, 2, log = TRUE)
#' pskellam(0, 1e6, 1e6)
#'
#' # Random generation
#' set.seed(1)
#' rskellam(10, 8.5, 10.25)
#'
#' @name skellam
#' @aliases dskellam pskellam qskellam rskellam dskellam.sp pskellam.sp
NULL

#' @export
#' @rdname skellam
dskellam <- function(x, lambda1, lambda2=lambda1, log=FALSE){
 # density (PMF) of Skellam distriubition (difference of Poissons)
    log <- isTRUE(as.logical(log[1]))
    nms <- names(x)
    nx <- length(x)
   # recycle to a common length; warn when lengths are not multiples (N13)
    a <- .skellam_recycle(x, lambda1, lambda2)
    x <- as.numeric(a[[1]]); lambda1 <- as.numeric(a[[2]]); lambda2 <- as.numeric(a[[3]])
    len <- length(x)
    ret <- rep(NaN, length.out = len)
   # invalid rates -> NaN with a warning, as dpois does
    bad <- .skellam_bad_rate(lambda1) | .skellam_bad_rate(lambda2)
    if (any(bad)) warning("NaNs produced")
   # missing values propagate (N1): NA stays NA, NaN stays NaN
    miss <- !bad & (is.na(x) | is.na(lambda1) | is.na(lambda2))
    ret[miss] <- (x + lambda1 + lambda2)[miss]
    ok <- !bad & !miss
   # outside the integer support the pmf is 0 (N11); infinite x likewise
    nonint <- ok & .skellam_nonint(x)
    if (any(nonint)) {
        xn <- x[nonint]
        warning("non-integer x = ", paste(format(xn[seq_len(min(5L, length(xn)))]), collapse = ", "),
                if (length(xn) > 5) ", ..." else "", call. = FALSE)
    }
    off <- nonint | (ok & is.infinite(x))
    ret[off] <- if (log) -Inf else 0
    ok <- ok & !off
    x[ok] <- round(x[ok])
   # handle a zero lambda separately (avoids division by zero & accuracy issues for large values of lambda or x)
    z1 <- ok & lambda1 == 0
    z2 <- ok & lambda2 == 0
    ret[z1] <- stats::dpois(-x[z1], lambda2[z1], log = log)
    ret[z2] <- stats::dpois( x[z2], lambda1[z2], log = log)    # corrects bug in VGAM 0.7-9
   # both rates positive: log-space core (see R/bessel.R), finite for every input
    nz <- ok & (lambda1 > 0) & (lambda2 > 0)
    if (any(nz)) {
        lp <- .skellam_logpmf(x[nz], lambda1[nz], lambda2[nz])
        ret[nz] <- if (log) lp else exp(lp)
    }
    if (length(nms) && nx == len) names(ret) <- nms
    ret
}
