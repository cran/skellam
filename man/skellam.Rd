\name{Skellam}
\alias{Skellam}
\alias{skellam}
\alias{dskellam}
\alias{pskellam}
\alias{qskellam}
\alias{rskellam}
\alias{dskellam.sp}
\alias{pskellam.sp}
\title{ The Skellam Distribution }
\description{
  Density, distribution function, quantile function and random number
  generation for the Skellam distribution.
}
\usage{
dskellam(x, lambda1, lambda2 = lambda1, log = FALSE)
pskellam(q, lambda1, lambda2 = lambda1, 
    lower.tail = TRUE, log.p = FALSE)
qskellam(p, lambda1, lambda2 = lambda1, 
    lower.tail = TRUE, log.p = FALSE)
rskellam(n, lambda1, lambda2 = lambda1)
dskellam.sp(x, lambda1, lambda2 = lambda1, log = FALSE)
pskellam.sp(q, lambda1, lambda2 = lambda1, 
    lower.tail = TRUE, log.p = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1},
    the length is taken to be the number required.}
  \item{lambda1, lambda2}{vectors of (non-negative) means.}
  \item{log, log.p}{ logical; if \code{TRUE}, probabilities p are given 
  as log(p). }
  \item{lower.tail}{ logical; if \code{TRUE} 
  (default), probabilities are \eqn{P[X \le x]},
    otherwise, \eqn{P[X > x]}. }
}
\details{
If \eqn{Y_1} and \eqn{Y_2} are Poisson variables with means 
  \eqn{\mu_1} and \eqn{\mu_2}
  and correlation \eqn{\rho}, then \eqn{X = Y_1 - Y_2} is Skellam with parameters
  \eqn{\lambda_1 = \mu_1 - \rho \sqrt{\mu_1 \mu_2}} and 
  \eqn{\lambda_2 = \mu_2 - \rho \sqrt{\mu_1 \mu_2}}.
  
  \code{dskellam} returns a value equivalent to
  \deqn{I(2 \sqrt{\lambda_1 \lambda_2}, |x|)  
    (\lambda_1/\lambda_2)^{x/2} \exp(-\lambda_1-\lambda_2)}
  where \eqn{I(y,\nu)} is the modified Bessel function of the first kind.
  The \eqn{|x|} differs from most Skellam expressions in the literature, but is correct since \eqn{x} is an integer, 
  resulting in improved portability and (in R versions prior to 2.9) improved accuracy for \eqn{x<0}.
  Exponential scaling is used in \code{\link[base]{besselI}} to postpone numerical problems. 
  When numerical problems do occur, a saddlepoint approximation is substituted,
  which typically gives at least 4-figure accuracy.
  An alternative representation is
  \eqn{dchisq(2 \lambda_1,2(x+1),2\lambda_2) 2}
  for \eqn{x \ge 0}, and
  \eqn{dchisq(2 \lambda_2,2(1-x),2\lambda_1) 2}
  for \eqn{x \le 0}; but in R \code{besselI} appears to be more accurately implemented
  (for very small probabilities) than \code{dchisq}.

  \code{pskellam(x,lambda1,lambda2)} returns
  \code{pchisq(2*lambda2, -2*x, 2*lambda1)}
  for \eqn{x<=0} and
  \code{1 - pchisq(2*lambda1, 2*(x+1), 2*lambda2)}
  for \eqn{x>=0}.
  When \code{\link[stats]{pchisq}} incorrectly returns 0, a saddlepoint approximation is substituted,
  which typically gives at least 2-figure accuracy.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{F(x) \ge p}, where \eqn{F} is the distribution function.
  For \code{lower.tail=FALSE}, the quantile is defined as the largest value \eqn{x} such that
  \eqn{F(x,}\code{lower.tail=FALSE}\eqn{) \le p}.

  \code{rskellam} is calculated as \code{rpois(n,lambda1)-rpois(n,lambda2)}

  \code{dskellam.sp} and \code{pskellam.sp} return saddlepoint approximations to the pmf and cdf.
  They are called by \code{dskellam} and \code{pskellam} when results from primary methods are in doubt.
}


\note{
  The \pkg{VGAM} package
   also contains \link[VGAM]{Skellam} functions, which are syntactically similar;
   independently developed versions are included here for completeness.
  Moreover, this \code{dskellam} function offers a broader working range,
   correct handling of cases where at least one rate parameter is zero,
   enhanced argument checking,
   and (in R versions prior to 2.9) improved accuracy for \eqn{x<0}.
  If both packages are loaded,
   \code{skellam::dskellam}
   or
   \code{VGAM::dskellam}
   can unambiguously specify which implementation to use.
}
\value{
  \code{dskellam} gives the (log) density, \code{pskellam} gives the (log) distribution function,
  \code{qskellam} gives the quantile function, and \code{rskellam} generates random deviates. 
  Invalid \code{lambda}s will result in return value \code{NaN}, with a warning.
}
\source{
  The relation of \code{dgamma} to the modified Bessel function of the first kind
   was given by Skellam (1946).
  The relation of \code{pgamma} to the noncentral chi-square was given by Johnson (1959).
  Tables are given by Strackee and van der Gon (1962), which can be used to verify this implementation 
   (cf. direct calculation in the examples below).

  \code{qskellam} uses the Cornish--Fisher expansion to include skewness and kurtosis
  corrections to a normal approximation, followed by a search.
  If \code{getOption("verbose")==TRUE}, then \code{qskellam} will not use \code{\link[stats]{qpois}}
  when one of the \code{lambda}s is zero,
  in order to verify that this search algorithm has been implemented properly.
}
\references{
  Butler, R. (2007) \emph{Saddlepoint Approximations with Applications},
   Cambridge University Press, Cambridge & New York, p.17.

  Johnson, N. L. (1959) On an extension of the connection between Poisson and \eqn{\chi^2} distributions.
   \emph{Biometrika} \bold{46}, 352--362.

  Johnson, N. L.; Kotz, S.; Kemp, A. W. (1993) \emph{Univariate Discrete Distributions}, 2nd ed.,
   John Wiley and Sons, New York, pp.190-192.

  Skellam, J. G. (1946) The frequency distribution of the difference between two Poisson variates
   belonging to different populations. \emph{Journal of the Royal Statistical Society, series A} \bold{109}/3, 26.

  Strackee, J.; van der Gon, J. J. D. (1962) The frequency distribution of the difference between two Poisson variates.
   \emph{Statistica Neerlandica} \bold{16}/1, 17-23.

  Wikipedia. \emph{Skellam distribution} \url{https://en.wikipedia.org/wiki/Skellam_distribution}
}
\examples{
require('skellam')

# one lambda = 0 ~ Poisson
c(dskellam(0:10,5,0), dpois(0:10,5))
c(dskellam(-(0:10),0,5), dpois(0:10,5))
c(pskellam(0:10,5,0,lower.tail=TRUE), 
    ppois(0:10,5,lower.tail=TRUE))
c(pskellam(0:10,5,0,lower.tail=FALSE),
    ppois(0:10,5,lower.tail=FALSE))
c(pskellam(-(0:10),0,5,lower.tail=FALSE),
    ppois(0:10-1,5,lower.tail=TRUE))
c(pskellam(-(0:10),0,5,lower.tail=TRUE),
    ppois(0:10-1,5,lower.tail=FALSE))

# both lambdas != 0 ~ convolution of Poissons
dskellam(1,0.5,0.75)  # sum(dpois(1+0:10,0.5)*dpois(0:10,0.75))
pskellam(1,0.5,0.75)  # sum(dskellam(-10:1,0.5,0.75))
dskellam(c(-1,1),c(12,10),c(10,12))  # c(0.0697968,0.0697968)
dskellam(c(-1,1),c(12,10),c(10,12),log=TRUE) 
 # log(dskellam(c(-1,1),c(12,10),c(10,12)))
dskellam(256,257,1)  
# 0.024829348733183769  
# exact result for comparison with saddlepoint
dskellam(-3724,2000,3000)  
# 3.1058145363400105e-308  
# exact result for comparison with saddlepoint 
# (still accurate in extreme tail)
pskellam(c(-1,0),c(12,10),c(10,12))  # c(0.2965079,0.7034921)
pskellam(c(-1,0),c(12,10),c(10,12),lower.tail=FALSE) 
 # 1-pskellam(c(-1,0),c(12,10),c(10,12))
pskellam(-2:2,8.5,10.25,log.p=TRUE)  # log(pskellam(-2:2,8.5,10.25))
qskellam(c(0.05,0.95),3,4)  # c(-5,3); pskellam(cbind(-6:-5,2:3),3,4)
qskellam(c(0.05,0.95),3,0)  # c(1,6); qpois(c(0.05,0.95),3)
rskellam(35,8.5,10.25)
}
\keyword{distribution}
\author{ Jerry W. Lewis, Patrick E. Brown}
